"""
GUI for displaying LSL logger messages in real-time.

This module provides a Tkinter-based GUI for monitoring log messages from the
LSL logger. It's designed to be integrated into applications that use LSL,
particularly those where console output might be hidden or inconvenient.

Usage Examples
--------------
Basic integration into an LSL script::
    
    from lsl.logger_gui import LoggerGUI
    from lsl.common import sdf
    
    # Create and start the GUI
    gui = LoggerGUI()
    
    # Your LSL processing code - log messages will appear in the GUI
    station = sdf.parse_sdf('myfile.sdf')
    # ... do work ...
    
    # Run the GUI event loop (this blocks until window is closed)
    gui.mainloop()

Running in a separate thread (non-blocking)::
    
    import threading
    from lsl.logger_gui import LoggerGUI
    from lsl.imaging import selfcal
    
    # Start GUI in background thread
    gui = LoggerGUI()
    gui_thread = threading.Thread(target=gui.mainloop, daemon=True)
    gui_thread.start()
    
    # Continue with LSL work while GUI runs
    # All log messages will appear in the GUI window
    gains, delays = selfcal.self_cal(...)
    
    # Keep script running
    input("Press Enter to exit...")

Using individual components::
    
    import tkinter as tk
    from lsl.logger_gui import LoggerFrame, FilterFrame
    
    # Create custom window with logger components
    root = tk.Tk()
    frame = tk.Frame(root)
    frame.pack()
    
    logger_display = LoggerFrame(frame)
    logger_display.start()
    
    root.mainloop()
"""

import queue
import signal
import logging
import threading

try:
    have_tk = True
    
    import tkinter as tk
    from tkinter import ttk, filedialog
    from tkinter.scrolledtext import ScrolledText
except ImportError:
    have_tk = False


from lsl import logger as lsl_logger


class LoggerFrame(object):
    """
    Simple frame to show log messages from the LSL logger/ThreadedHandler for
    easy inspection outside of a terminal.
    
    .. note:: Initializing this class will attach a new ThreadedHandler to the
              main LSL logger.
    """
    
    def __init__(self, frame, update_interval_ms=100, max_lines=10000):
        self._frame = frame
        self._update_interval_ms = update_interval_ms
        self._max_lines = max_lines
        self._line_count = 0
        
        # Setup the text area
        self._text = ScrolledText(frame, state='disabled', height=12)
        self._text.grid(row=0, column=0, sticky=(tk.N, tk.S, tk.W, tk.E))
        self._text.configure(font='TkFixedFont')
        self._text.tag_config('DEBUG', foreground='gray')
        self._text.tag_config('INFO', foreground='black')
        self._text.tag_config('WARNING', foreground='orange')
        self._text.tag_config('ERROR', foreground='red')
        self._text.tag_config('CRITICAL', foreground='red', underline=1)
        
        # Attach the ThreadedHandler so we that can access messages from other
        # threads
        self._handler = lsl_logger.ThreadedHandler()
        lsl_logger.add_handler(self._handler)
        
    def _display(self, record):
        """
        Method to take a LogRecord, format it, and display it in the text area.
        """
        
        # Format
        msg = self._handler.format(record)
        
        # Add
        self._text.configure(state='normal')
        self._text.insert(tk.END, msg + '\n', record.levelname)
        self._line_count += 1
        
        # Trim buffer if it exceeds max_lines
        if self._line_count > self._max_lines:
            # Delete the first line
            self._text.delete(1.0, 2.0)
            self._line_count -= 1
            
        self._text.configure(state='disabled')
        
        # Scroll
        self._text.yview(tk.END)
        
    def _poll_log_queue(self):
        """
        Method to poll the main LSL logger queue for new LogRecords and update
        the text area.  You only need to call this once since it reschedules itself.
        """
        
        while True:
            try:
                record = lsl_logger.LSL_LOG_QUEUE.get(block=False)
            except queue.Empty:
                break
            else:
                self._display(record)
                
        self._id = self._frame.after(self._update_interval_ms, self._poll_log_queue)
        
    def start(self):
        """
        Start up the background poller.
        """
        
        # Start up the poller
        self._id = self._frame.after(self._update_interval_ms, self._poll_log_queue)
        
    def stop(self):
        """
        Stop the background poller and remove the ThreadedHandler handler from
        the main LSL logger.
        """
        
        # Stop the polling
        try:
            self._frame.after_cancel(self._id)
        except AttributeError:
            pass
            
        # Remove the handler
        try:
            lsl_logger.remove_handler(self._handler)
        except AttributeError:
            pass
            
    def show_at_level(self, level):
        """
        Show/hide log messages at different levels by eliding text tags.
        
        Parameters
        ----------
        level : str
            The minimum level to display ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL').
        """
        
        found = False
        for l in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'):
            if level == l:
                found = True
                
            try:
                self._text.tag_config(l, elide=not found)
            except Exception as e:
                print(str(e))
                
    def clear(self):
        """
        Clear all text from the display buffer.
        """
        
        self._text.configure(state='normal')
        try:
            self._text.delete(1.0, tk.END)
            self._line_count = 0
        except Exception as e:
            print(str(e))
        self._text.configure(state='disabled')


class FilterFrame(object):
    """
    Simple frame that can be combined with LoggerFrame to control what kinds of
    log levels are displayed going forward, enable file logging, console output,
    and pattern-based filtering.
    """
    
    _values = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    _level_map = {logging.DEBUG: 0, logging.INFO: 1, logging.WARNING: 2,
                  logging.ERROR: 3, logging.CRITICAL: 4}
    _reverse_level_map = {0: logging.DEBUG, 1: logging.INFO, 2: logging.WARNING,
                          3: logging.ERROR, 4: logging.CRITICAL}
    
    def __init__(self, frame, logger_frame):
        self._frame = frame
        self._logger_frame = logger_frame
        self._logging_file = None
        self._active_patterns = []
        
        # Setup the combo box with the various options for logger level
        self._level = tk.StringVar()
        ttk.Label(self._frame, text='Logger Level:').grid(column=0, row=0, sticky=tk.W, padx=5)
        self._combobox = ttk.Combobox(self._frame,
                                      textvariable=self._level,
                                      width=10,
                                      state='readonly',
                                      values=self._values
                                     )
        current_level = lsl_logger.get_log_level()
        self._combobox.current(self._level_map.get(current_level, 1))
        self._combobox.grid(column=1, row=0, sticky=tk.W, padx=5)
        self._combobox.bind('<<ComboboxSelected>>', self._on_level_change)
        
        # Console output toggle
        self._console_enabled = tk.BooleanVar(value=False)
        self._console_check = ttk.Checkbutton(self._frame,
                                              text="Console Output",
                                              variable=self._console_enabled,
                                              command=self._on_console_toggle)
        self._console_check.grid(column=2, row=0, sticky=tk.W, padx=5)
        
        # File logging button
        self._file_button = ttk.Button(self._frame, text="Log to File...",
                                       command=self._on_file_logging)
        self._file_button.grid(column=3, row=0, sticky=tk.W, padx=5)
        
        # Clear button
        self._clear = ttk.Button(self._frame, text="Clear Buffer", command=self._on_clear)
        self._clear.grid(column=4, row=0, sticky=tk.E, padx=5)
        
        # Pattern filtering - second row
        ttk.Label(self._frame, text='Module Filter:').grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)
        self._pattern_var = tk.StringVar()
        self._pattern_entry = ttk.Entry(self._frame, textvariable=self._pattern_var, width=25)
        self._pattern_entry.grid(column=1, row=1, columnspan=2, sticky=(tk.W, tk.E), padx=5, pady=5)
        self._pattern_entry.bind('<Return>', lambda e: self._on_add_filter())
        
        self._add_filter_button = ttk.Button(self._frame, text="Add Filter",
                                             command=self._on_add_filter)
        self._add_filter_button.grid(column=3, row=1, sticky=tk.W, padx=5, pady=5)

        self._clear_filters_button = ttk.Button(self._frame, text="Clear Filters",
                                                command=self._on_clear_filters)
        self._clear_filters_button.grid(column=4, row=1, sticky=tk.E, padx=5, pady=5)
        
        # Active filters display
        ttk.Label(self._frame, text='Active Filters:').grid(column=0, row=2, sticky=tk.W, padx=5)
        self._filters_label = ttk.Label(self._frame, text='None', foreground='gray')
        self._filters_label.grid(column=1, row=2, columnspan=4, sticky=tk.W, padx=5)
        
    def _on_level_change(self, event):
        """
        Method to propagate the change in the log level to both the logger and
        the display filtering.
        """
        
        level_name = self._combobox.get()
        level_index = self._values.index(level_name)
        level = self._reverse_level_map[level_index]
        
        # Set the actual logger level
        lsl_logger.set_log_level(level)
        
        # Also update display filtering
        self._logger_frame.show_at_level(level_name)
        
    def _on_console_toggle(self):
        """
        Enable or disable console logging based on checkbox state.
        """
        
        if self._console_enabled.get():
            lsl_logger.enable_console_logging()
        else:
            lsl_logger.disable_console_logging()
            
    def _on_file_logging(self):
        """
        Enable or disable file logging with file dialog.
        """
        
        if self._logging_file is None:
            filename = filedialog.asksaveasfilename(
                defaultextension=".log",
                filetypes=[("Log files", "*.log"), ("All files", "*.*")],
                title="Save Log File As"
            )
            if filename:
                lsl_logger.enable_file_logging(filename)
                self._logging_file = filename
                self._file_button.config(text="Stop File Logging")
        else:
            lsl_logger.disable_file_logging()
            self._logging_file = None
            self._file_button.config(text="Log to File...")
            
    def _on_add_filter(self):
        """
        Add a module pattern filter to the logger.
        """
        
        pattern = self._pattern_var.get().strip()
        if pattern and pattern not in self._active_patterns:
            lsl_logger.add_filter(pattern)
            self._active_patterns.append(pattern)
            self._update_filters_display()
            self._pattern_var.set('')
            
    def _on_clear_filters(self):
        """
        Clear all active filters from the logger.
        """
        
        lsl_logger.clear_filters()
        self._active_patterns.clear()
        self._update_filters_display()
        
    def _update_filters_display(self):
        """
        Update the display of active filters.
        """
        
        if self._active_patterns:
            filters_text = ', '.join(self._active_patterns)
            self._filters_label.config(text=filters_text, foreground='black')
        else:
            self._filters_label.config(text='None', foreground='gray')
            
    def _on_clear(self):
        """
        Method to clear the LoggerFrame via its clear() method.
        """
        
        self._logger_frame.clear()
        
    def start(self):
        """
        Stub function to make FilterFrame look like LoggerFrame.
        """
        
        pass
        
    def stop(self):
        """
        Stub function to make FilterFrame look like LoggerFrame.
        """
        
        pass


class LoggerGUI(object):
    """
    Tk GUI to help display messages from the main LSL logger in a "not on the
    console" way.
    
    Features
    --------
    - Real-time display of log messages with color-coded levels
    - Adjustable logger level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    - Optional console output toggle
    - File logging with save dialog
    - Module-based pattern filtering (e.g., 'lsl.imaging.*')
    - Automatic buffer management (max 10,000 lines)
    - Clear buffer button
    
    Parameters
    ----------
    root : tk.Tk, optional
        Existing Tk root window. If None, a new root window is created.
    
    Examples
    --------
    Simple usage with blocking GUI::
        
        from lsl.logger_gui import LoggerGUI
        gui = LoggerGUI()
        # ... do LSL work that generates log messages ...
        gui.mainloop()
    
    Non-blocking background GUI::
        
        import threading
        from lsl.logger_gui import LoggerGUI
        
        gui = LoggerGUI()
        threading.Thread(target=gui.mainloop, daemon=True).start()
        # ... continue with other work ...
    
    Notes
    -----
    The GUI attaches a ThreadedHandler to the main LSL logger, allowing it to
    capture messages from any thread. When the GUI is closed, the handler is
    automatically removed.
    """
    
    def __init__(self, root=None):
        if not have_tk:
            raise RuntimeError("Cannot create GUI because tk module is not available")
            
        if root is None:
            root = tk.Tk()
            
        self._root = root
        self._root.title('LSL Logger')
        self._root.columnconfigure(0, weight=1)
        self._root.rowconfigure(0, weight=1)
        
        # Create the panels and frames needed for LoggerFrame and FilterFrame
        vertical_pane = ttk.PanedWindow(self._root, orient=tk.VERTICAL)
        vertical_pane.grid(row=0, column=0, sticky=(tk.N, tk.S, tk.W, tk.E))
        ## LoggerFrame holder
        display_frame = ttk.Labelframe(vertical_pane, text="")
        display_frame.columnconfigure(0, weight=1)
        display_frame.rowconfigure(0, weight=1)
        vertical_pane.add(display_frame, weight=1)
        ## FilterFrame holder
        filter_frame = ttk.Labelframe(vertical_pane, text="")
        filter_frame.columnconfigure(0, weight=1)
        filter_frame.columnconfigure(1, weight=1)
        filter_frame.columnconfigure(2, weight=1)
        filter_frame.columnconfigure(3, weight=1)
        filter_frame.columnconfigure(4, weight=1)
        filter_frame.rowconfigure(0, weight=1)
        filter_frame.rowconfigure(1, weight=1)
        filter_frame.rowconfigure(2, weight=1)
        vertical_pane.add(filter_frame, weight=1)
        
        # Create the frame classes themselves
        self._display = LoggerFrame(display_frame)
        self._filter = FilterFrame(filter_frame, self._display)
        
        # Bind
        self._root.protocol('WM_DELETE_WINDOW', self.quit)
        try:
            signal.signal(signal.SIGINT, self.quit)
        except (ValueError, OSError):
            # Signal handling may not work in all contexts (e.g., threads, Windows)
            pass
            
        # Start
        self._display.start()
        self._filter.start()
        
    def mainloop(self):
        """
        Call the mainloop of Tk.
        """
        
        self._root.mainloop()
        
    def quit(self, *args):
        """
        Quit out of the Tk window.
        """
        
        self._filter.stop()
        self._display.stop()
        self._root.destroy()


if __name__ == '__main__':
    import time
    import threading
    
    print("Starting logger GUI demo mode...")
    print("This will generate test log messages to demonstrate the GUI features.")
    print("Close the GUI window to exit.\n")
    
    # Create the GUI
    gui = LoggerGUI()
    
    # Set logger to DEBUG level to see all messages
    lsl_logger.set_log_level(logging.DEBUG)
    
    # Demo function that generates log messages
    def generate_demo_messages():
        time.sleep(2)  # Wait for GUI to initialize
        
        lsl_logger.LSL_LOGGER.debug("Debug message: Detailed diagnostic information")
        time.sleep(1)
        
        lsl_logger.LSL_LOGGER.info("Info message: Normal operational message")
        time.sleep(1)
        
        lsl_logger.LSL_LOGGER.warning("Warning message: Something unexpected happened")
        time.sleep(1)
        
        lsl_logger.LSL_LOGGER.error("Error message: Operation failed")
        time.sleep(1)
        
        lsl_logger.LSL_LOGGER.critical("Critical message: System in critical state")
        time.sleep(2)
        
        lsl_logger.LSL_LOGGER.info("Try adjusting the Logger Level dropdown")
        time.sleep(1)
        lsl_logger.LSL_LOGGER.info("Enable Console Output to see messages in terminal")
        time.sleep(1)
        lsl_logger.LSL_LOGGER.info("Click 'Log to File...' to save messages to disk")
        time.sleep(1)
        lsl_logger.LSL_LOGGER.info("Use Module Filter to filter by pattern (e.g., 'lsl_logger')")
        time.sleep(2)
        
        # Simulate some work with progress messages
        for i in range(5):
            lsl_logger.LSL_LOGGER.debug(f"Processing iteration {i+1} of 5...")
            time.sleep(1)
            
        lsl_logger.LSL_LOGGER.info("Demo complete! Close window to exit.")
        
    # Start demo message generator in background thread
    demo_thread = threading.Thread(target=generate_demo_messages, daemon=True)
    demo_thread.start()
    
    # Run GUI (blocks until window closed)
    gui.mainloop()
