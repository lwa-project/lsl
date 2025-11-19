"""
Unit tests for the lsl.logger_gui module.
"""

import os
import sys
import time
import queue
import logging
import unittest
import threading

# Check if Tk is available
try:
    import tkinter as tk
    have_tk = True
    # Try to create a test root to see if display is available
    try:
        test_root = tk.Tk()
        test_root.destroy()
        have_display = True
    except tk.TclError:
        have_display = False
except ImportError:
    have_tk = False
    have_display = False

from lsl import logger as lsl_logger

# Only import GUI if Tk is available
if have_tk:
    from lsl.logger_gui import LoggerFrame, FilterFrame, LoggerGUI

__version__ = "0.1"
__author__ = "Jayce Dowell"


@unittest.skipIf(not have_tk, "Tkinter not available")
@unittest.skipIf(not have_display, "No display available for GUI tests")
class logger_gui_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.logger_gui module."""

    def setUp(self):
        """Set up test fixtures."""
        # Store original logger state
        self._original_level = lsl_logger.get_log_level()
        # Clear the log queue
        while not lsl_logger.LSL_LOG_QUEUE.empty():
            try:
                lsl_logger.LSL_LOG_QUEUE.get_nowait()
            except queue.Empty:
                break

    def tearDown(self):
        """Clean up after tests."""
        # Restore original log level
        lsl_logger.set_log_level(self._original_level)
        # Clear console/file logging
        lsl_logger.disable_console_logging()
        lsl_logger.disable_file_logging()
        # Clear filters
        lsl_logger.clear_filters()

    def test_logger_frame_creation(self):
        """Test creating a LoggerFrame."""
        root = tk.Tk()
        try:
            frame = tk.Frame(root)
            frame.pack()

            logger_frame = LoggerFrame(frame)
            self.assertIsNotNone(logger_frame)
            self.assertIsInstance(logger_frame._handler, lsl_logger.ThreadedHandler)

            # Verify handler was added to logger
            self.assertIn(logger_frame._handler, lsl_logger.LSL_LOGGER.handlers)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_logger_frame_display(self):
        """Test displaying messages in LoggerFrame."""
        root = tk.Tk()
        try:
            frame = tk.Frame(root)
            frame.pack()

            logger_frame = LoggerFrame(frame)
            logger_frame.start()

            # Generate a log message
            lsl_logger.LSL_LOGGER.info("Test GUI message")

            # Give GUI time to process (poll interval is 100ms by default)
            root.update()
            time.sleep(0.15)
            root.update()

            # Check if message appears in text widget
            text_content = logger_frame._text.get(1.0, tk.END)
            self.assertIn("Test GUI message", text_content)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_logger_frame_max_lines(self):
        """Test buffer management with max_lines."""
        root = tk.Tk()
        try:
            frame = tk.Frame(root)
            frame.pack()

            # Create frame with small buffer
            logger_frame = LoggerFrame(frame, max_lines=5)
            logger_frame.start()

            # Generate more messages than max_lines
            for i in range(10):
                lsl_logger.LSL_LOGGER.info(f"Message {i}")
                root.update()
                time.sleep(0.02)

            # Give GUI time to process all messages
            time.sleep(0.2)
            root.update()

            # Check line count doesn't exceed max
            self.assertLessEqual(logger_frame._line_count, 5)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_logger_frame_clear(self):
        """Test clearing the LoggerFrame buffer."""
        root = tk.Tk()
        try:
            frame = tk.Frame(root)
            frame.pack()

            logger_frame = LoggerFrame(frame)
            logger_frame.start()

            # Generate messages
            lsl_logger.LSL_LOGGER.info("Message before clear")
            root.update()
            time.sleep(0.15)
            root.update()

            # Clear the buffer
            logger_frame.clear()
            root.update()

            # Check buffer is empty
            text_content = logger_frame._text.get(1.0, tk.END).strip()
            self.assertEqual(text_content, "")
            self.assertEqual(logger_frame._line_count, 0)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_logger_frame_show_at_level(self):
        """Test filtering display by log level."""
        root = tk.Tk()
        try:
            frame = tk.Frame(root)
            frame.pack()

            logger_frame = LoggerFrame(frame)
            logger_frame.start()

            # Set logger to DEBUG to capture all messages
            lsl_logger.set_log_level(logging.DEBUG)

            # Generate messages at different levels
            lsl_logger.LSL_LOGGER.debug("Debug message")
            lsl_logger.LSL_LOGGER.info("Info message")
            lsl_logger.LSL_LOGGER.warning("Warning message")

            # Give GUI time to process
            root.update()
            time.sleep(0.2)
            root.update()

            # Show only WARNING and above
            logger_frame.show_at_level('WARNING')
            root.update()

            # Note: show_at_level uses tag eliding, so messages are still in the
            # text widget but DEBUG and INFO are hidden. We can't easily test
            # visibility, but we can verify the method doesn't crash.

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_filter_frame_creation(self):
        """Test creating a FilterFrame."""
        root = tk.Tk()
        try:
            display_frame = tk.Frame(root)
            display_frame.pack()
            logger_frame = LoggerFrame(display_frame)

            control_frame = tk.Frame(root)
            control_frame.pack()
            filter_frame = FilterFrame(control_frame, logger_frame)

            self.assertIsNotNone(filter_frame)
            self.assertEqual(filter_frame._logger_frame, logger_frame)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_filter_frame_level_change(self):
        """Test changing logger level via FilterFrame."""
        root = tk.Tk()
        try:
            display_frame = tk.Frame(root)
            display_frame.pack()
            logger_frame = LoggerFrame(display_frame)

            control_frame = tk.Frame(root)
            control_frame.pack()
            filter_frame = FilterFrame(control_frame, logger_frame)

            # Change to WARNING level
            filter_frame._combobox.set('WARNING')
            filter_frame._on_level_change(None)
            root.update()

            # Check logger level was changed
            self.assertEqual(lsl_logger.get_log_level(), logging.WARNING)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_filter_frame_console_toggle(self):
        """Test console output toggle."""
        root = tk.Tk()
        try:
            display_frame = tk.Frame(root)
            display_frame.pack()
            logger_frame = LoggerFrame(display_frame)

            control_frame = tk.Frame(root)
            control_frame.pack()
            filter_frame = FilterFrame(control_frame, logger_frame)

            # Initially console should be disabled
            self.assertIsNone(lsl_logger._console_handler)

            # Enable console
            filter_frame._console_enabled.set(True)
            filter_frame._on_console_toggle()
            root.update()

            # Check console handler was added
            self.assertIsNotNone(lsl_logger._console_handler)

            # Disable console
            filter_frame._console_enabled.set(False)
            filter_frame._on_console_toggle()
            root.update()

            # Check console handler was removed
            self.assertIsNone(lsl_logger._console_handler)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_filter_frame_add_clear_filters(self):
        """Test adding and clearing pattern filters."""
        root = tk.Tk()
        try:
            display_frame = tk.Frame(root)
            display_frame.pack()
            logger_frame = LoggerFrame(display_frame)

            control_frame = tk.Frame(root)
            control_frame.pack()
            filter_frame = FilterFrame(control_frame, logger_frame)

            # Add a filter
            filter_frame._pattern_var.set('lsl.imaging.*')
            filter_frame._on_add_filter()
            root.update()

            # Check filter was added
            self.assertIn('lsl.imaging.*', filter_frame._active_patterns)
            self.assertEqual(len(lsl_logger._active_filters), 1)

            # Add another filter
            filter_frame._pattern_var.set('lsl_logger')
            filter_frame._on_add_filter()
            root.update()

            self.assertIn('lsl_logger', filter_frame._active_patterns)
            self.assertEqual(len(lsl_logger._active_filters), 2)

            # Clear filters
            filter_frame._on_clear_filters()
            root.update()

            self.assertEqual(len(filter_frame._active_patterns), 0)
            self.assertEqual(len(lsl_logger._active_filters), 0)

            # Clean up
            logger_frame.stop()
        finally:
            root.destroy()

    def test_logger_gui_creation(self):
        """Test creating a complete LoggerGUI."""
        root = tk.Tk()
        try:
            gui = LoggerGUI(root=root)
            self.assertIsNotNone(gui)
            self.assertIsNotNone(gui._display)
            self.assertIsNotNone(gui._filter)

            # Clean up
            gui.quit()
        finally:
            # Ensure root is destroyed even if quit fails
            try:
                root.destroy()
            except tk.TclError:
                pass

    def test_logger_gui_with_messages(self):
        """Test LoggerGUI with actual log messages."""
        root = tk.Tk()
        try:
            gui = LoggerGUI(root=root)

            # Generate log messages
            lsl_logger.LSL_LOGGER.info("Test message 1")
            lsl_logger.LSL_LOGGER.warning("Test warning")

            # Give GUI time to process
            root.update()
            time.sleep(0.2)
            root.update()

            # Check messages appear in display
            text_content = gui._display._text.get(1.0, tk.END)
            self.assertIn("Test message 1", text_content)
            self.assertIn("Test warning", text_content)

            # Clean up
            gui.quit()
        finally:
            try:
                root.destroy()
            except tk.TclError:
                pass


class logger_gui_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.logger_gui unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(logger_gui_tests))


if __name__ == '__main__':
    unittest.main()
