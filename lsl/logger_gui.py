"""
Tk-based GUI widgets for displaying LSL logger messages in real-time.
Provides LoggerFrame for embedding in existing Tk applications,
FilterFrame for controlling log display, and LoggerGUI for standalone
use.
"""

import queue
import signal
import logging
import threading

from lsl import logger as lsl_logger
from lsl.misc import telemetry
telemetry.track_module()

try:
    have_tk = True

    import tkinter as tk
    from tkinter import ttk, filedialog
    from tkinter.scrolledtext import ScrolledText
except ImportError:
    have_tk = False


class LoggerFrame(object):
    """
    Tk frame to show log messages from the LSL logger.  Attaches a
    ThreadedHandler to the main LSL logger on initialization.
    """
    
    def __init__(self, frame, update_interval_ms=100, max_lines=10000):
        if not have_tk:
            raise RuntimeError("Cannot create LoggerFrame because tk module is not available")
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
        self._text.tag_config('NOTE', foreground='blue')
        
        # Attach the ThreadedHandler so we that can access messages from other
        # threads
        self._handler = lsl_logger.ThreadedHandler()
        lsl_logger.add_handler(self._handler)
        
    def _display(self, record):
        """Format and display a LogRecord in the text area."""

        # Format
        msg = self._handler.format(record)
        levelname = record.levelname

        # Add
        self._text.configure(state='normal')
        self._text.insert(tk.END, msg + '\n', levelname)
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
        """Poll the log queue and reschedule."""
        
        while True:
            try:
                record = lsl_logger.LSL_LOG_QUEUE.get(block=False)
            except queue.Empty:
                break
            else:
                self._display(record)
                
        self._id = self._frame.after(self._update_interval_ms, self._poll_log_queue)
        
    def start(self):
        """Start the background queue poller."""

        self._id = self._frame.after(self._update_interval_ms, self._poll_log_queue)

    def stop(self):
        """Stop the poller and remove the handler from the LSL logger."""
        
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
        """Show/hide log messages by eliding text tags below `level`."""
        
        found = False
        for l in ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'):
            if level == l:
                found = True
                
            try:
                self._text.tag_config(l, elide=not found)
            except Exception as e:
                lsl_logger.LSL_LOGGER.warning(f"Failed to update tag elide state: {e}")

    def clear(self):
        """Clear the display buffer."""

        self._text.configure(state='normal')
        try:
            self._text.delete(1.0, tk.END)
            self._line_count = 0
        except Exception as e:
            lsl_logger.LSL_LOGGER.warning(f"Failed to clear display buffer: {e}")
        self._text.configure(state='disabled')


class FilterFrame(object):
    """
    Tk frame for controlling log level, console/file output, and
    module pattern filtering.  Designed to pair with a LoggerFrame.
    """
    
    _values = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    _level_map = {logging.DEBUG: 0, logging.INFO: 1, logging.WARNING: 2,
                  logging.ERROR: 3, logging.CRITICAL: 4}
    _reverse_level_map = {0: logging.DEBUG, 1: logging.INFO, 2: logging.WARNING,
                          3: logging.ERROR, 4: logging.CRITICAL}
    
    def __init__(self, frame, logger_frame):
        if not have_tk:
            raise RuntimeError("Cannot create FilterFrame because tk module is not available")
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
        """Propagate the level change to the logger and display."""
        
        level_name = self._combobox.get()
        level_index = self._values.index(level_name)
        level = self._reverse_level_map[level_index]
        
        # Set the actual logger level
        lsl_logger.set_log_level(level)
        
        # Also update display filtering
        self._logger_frame.show_at_level(level_name)
        
    def _on_console_toggle(self):
        """Toggle console logging."""
        
        if self._console_enabled.get():
            lsl_logger.enable_console_logging()
        else:
            lsl_logger.disable_console_logging()
            
    def _on_file_logging(self):
        """Toggle file logging via a save dialog."""
        
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
        """Add a module pattern filter."""
        
        pattern = self._pattern_var.get().strip()
        if pattern and pattern not in self._active_patterns:
            lsl_logger.add_filter(pattern)
            self._active_patterns.append(pattern)
            self._update_filters_display()
            self._pattern_var.set('')
            
    def _on_clear_filters(self):
        """Clear all active filters."""
        
        lsl_logger.clear_filters()
        self._active_patterns.clear()
        self._update_filters_display()
        
    def _update_filters_display(self):
        """Update the active filters label."""
        
        if self._active_patterns:
            filters_text = ', '.join(self._active_patterns)
            self._filters_label.config(text=filters_text, foreground='black')
        else:
            self._filters_label.config(text='None', foreground='gray')
            
    def _on_clear(self):
        """Clear the LoggerFrame display."""
        
        self._logger_frame.clear()
        
    def start(self):
        """No-op; provided for interface parity with LoggerFrame."""
        pass

    def stop(self):
        """No-op; provided for interface parity with LoggerFrame."""
        pass


class LoggerGUI(object):
    """
    Standalone Tk window combining a LoggerFrame and FilterFrame for
    displaying LSL logger messages.  Creates a new Tk root or a
    Toplevel if one already exists.
    """

    def __init__(self, title='LSL Logger'):
        if not have_tk:
            raise RuntimeError("Cannot create GUI because tk module is not available")

        # Automatically detect if a root window exists
        if tk._default_root is None:
            self._root = tk.Tk()
        else:
            self._root = tk.Toplevel()

        self._root.title(title)
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
        self._root.mainloop()

    def quit(self, *args):
        """Close the window and clean up handlers."""
        
        self._filter.stop()
        self._display.stop()
        self._root.destroy()
