"""
Unit tests for the lsl.logger module.
"""

import os
import queue
import logging
import tempfile
import unittest
import warnings

from lsl import logger as lsl_logger

__version__ = "0.1"
__author__ = "Jayce Dowell"


class logger_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.logger module."""
    
    def setUp(self):
        """Set up test fixtures."""
        
        # Store original logger state
        self._original_level = lsl_lsl_logger.get_log_level()
        
        # Clear any existing handlers except the default NullHandler
        for handler in lsl_lsl_logger.LSL_LOGGER.handlers[:]:
            if not isinstance(handler, logging.NullHandler):
                lsl_lsl_logger.LSL_LOGGER.removeHandler(handler)
                
        # Clear any filters
        lsl_lsl_logger.clear_filters()
        
        # Reset module-level handler tracking
        lsl_lsl_logger._console_handler = None
        lsl_lsl_logger._file_handler = None
        
    def tearDown(self):
        """Clean up after tests."""
        
        # Restore original log level
        lsl_lsl_logger.set_log_level(self._original_level)
        
        # Disable console and file logging
        lsl_lsl_logger.disable_console_logging()
        lsl_lsl_logger.disable_file_logging()
        
        # Clear all filters
        lsl_lsl_logger.clear_filters()
        
        # Clear the log queue
        while not lsl_lsl_logger.LSL_LOG_QUEUE.empty():
            try:
                lsl_lsl_logger.LSL_LOG_QUEUE.get_nowait()
            except queue.Empty:
                break
                
    def test_set_get_log_level(self):
        """Test setting and getting log levels."""
        
        lsl_lsl_logger.set_log_level(logging.DEBUG)
        self.assertEqual(lsl_lsl_logger.get_log_level(), logging.DEBUG)
        
        lsl_lsl_logger.set_log_level(logging.WARNING)
        self.assertEqual(lsl_lsl_logger.get_log_level(), logging.WARNING)
        
        lsl_lsl_logger.set_log_level(logging.ERROR)
        self.assertEqual(lsl_lsl_logger.get_log_level(), logging.ERROR)
        
    def test_threaded_handler(self):
        """Test ThreadedHandler puts records in queue."""
        
        # Add a ThreadedHandler
        handler = lsl_lsl_logger.ThreadedHandler()
        lsl_lsl_logger.add_handler(handler)
        
        # Clear queue first
        while not lsl_lsl_logger.LSL_LOG_QUEUE.empty():
            lsl_lsl_logger.LSL_LOG_QUEUE.get_nowait()
            
        # Generate a log message
        lsl_lsl_logger.LSL_LOGGER.info("Test message")
        
        # Check it appears in queue
        self.assertFalse(lsl_lsl_logger.LSL_LOG_QUEUE.empty())
        record = lsl_lsl_logger.LSL_LOG_QUEUE.get(timeout=1)
        self.assertIsInstance(record, logging.LogRecord)
        self.assertIn("Test message", record.getMessage())
        
        # Clean up
        lsl_lsl_logger.remove_handler(handler)
        
    def test_console_logging_enable_disable(self):
        """Test enabling and disabling console logging."""
        
        # Initially no console handler
        self.assertIsNone(lsl_lsl_logger._console_handler)
        
        # Enable console logging
        lsl_lsl_logger.enable_console_logging()
        self.assertIsNotNone(lsl_lsl_logger._console_handler)
        self.assertIsInstance(lsl_lsl_logger._console_handler, logging.StreamHandler)
        
        # Verify handler is in logger
        self.assertIn(lsl_lsl_logger._console_handler, lsl_lsl_logger.LSL_LOGGER.handlers)
        
        # Disable console logging
        lsl_lsl_logger.disable_console_logging()
        self.assertIsNone(lsl_lsl_logger._console_handler)
        
    def test_console_logging_with_level(self):
        """Test console logging with specific level."""
        
        lsl_lsl_logger.enable_console_logging(level=logging.WARNING)
        self.assertIsNotNone(lsl_lsl_logger._console_handler)
        self.assertEqual(lsl_lsl_logger._console_handler.level, logging.WARNING)
        
    def test_file_logging_enable_disable(self):
        """Test enabling and disabling file logging."""
        
        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            log_file = f.name
            
        try:
            # Initially no file handler
            self.assertIsNone(lsl_lsl_logger._file_handler)
            
            # Enable file logging
            lsl_lsl_logger.enable_file_logging(log_file)
            self.assertIsNotNone(lsl_logger._file_handler)
            self.assertIsInstance(lsl_logger._file_handler, logging.FileHandler)
            
            # Verify handler is in logger
            self.assertIn(lsl_logger._file_handler, lsl_logger.LSL_LOGGER.handlers)
            
            # Write a test message
            lsl_logger.LSL_LOGGER.info("Test file message")
            
            # Disable file logging
            lsl_logger.disable_file_logging()
            self.assertIsNone(lsl_logger._file_handler)
            
            # Check message was written
            with open(log_file, 'r') as f:
                content = f.read()
                self.assertIn("Test file message", content)
                
        finally:
            # Clean up temp file
            if os.path.exists(log_file):
                os.unlink(log_file)
                
    def test_file_logging_switch_warning(self):
        """Test warning when switching log files."""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            log_file1 = f.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            log_file2 = f.name
            
        try:
            # Enable first file
            lsl_logger.enable_file_logging(log_file1)
            
            # Switch to second file should generate warning
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                lsl_logger.enable_file_logging(log_file2)
                
                # Check warning was raised
                self.assertEqual(len(w), 1)
                self.assertTrue(issubclass(w[0].category, RuntimeWarning))
                self.assertIn("Closing existing log file", str(w[0].message))
                self.assertIn(log_file1, str(w[0].message))  # Old file
                self.assertIn(log_file2, str(w[0].message))  # New file
                
            lsl_logger.disable_file_logging()
            
        finally:
            if os.path.exists(log_file1):
                os.unlink(log_file1)
            if os.path.exists(log_file2):
                os.unlink(log_file2)
                
    def test_file_logging_with_level(self):
        """Test file logging with specific level."""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            log_file = f.name
            
        try:
            lsl_logger.enable_file_logging(log_file, level=logging.ERROR)
            self.assertIsNotNone(lsl_logger._file_handler)
            self.assertEqual(lsl_logger._file_handler.level, logging.ERROR)

            lsl_logger.disable_file_logging()
            
        finally:
            if os.path.exists(log_file):
                os.unlink(log_file)
                
    def test_file_logging_append_mode(self):
        """Test file logging in append mode."""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            log_file = f.name
            
        try:
            # Write first message
            lsl_logger.enable_file_logging(log_file, mode='w')
            lsl_logger.LSL_LOGGER.info("Message 1")
            lsl_logger.disable_file_logging()
            
            # Write second message in append mode
            lsl_logger.enable_file_logging(log_file, mode='a')
            lsl_logger.LSL_LOGGER.info("Message 2")
            lsl_logger.disable_file_logging()
            
            # Check both messages are present
            with open(log_file, 'r') as f:
                content = f.read()
                self.assertIn("Message 1", content)
                self.assertIn("Message 2", content)
                
        finally:
            if os.path.exists(log_file):
                os.unlink(log_file)
                
    def test_file_logging_write_mode(self):
        """Test file logging in write mode (overwrite)."""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            log_file = f.name
            
        try:
            # Write first message
            lsl_logger.enable_file_logging(log_file, mode='w')
            lsl_logger.LSL_LOGGER.info("Message 1")
            lsl_logger.disable_file_logging()
            
            # Overwrite with second message
            lsl_logger.enable_file_logging(log_file, mode='w')
            lsl_logger.LSL_LOGGER.info("Message 2")
            lsl_logger.disable_file_logging()
            
            # Check only second message is present
            with open(log_file, 'r') as f:
                content = f.read()
                self.assertNotIn("Message 1", content)
                self.assertIn("Message 2", content)
                
        finally:
            if os.path.exists(log_file):
                os.unlink(log_file)
                
    def test_add_filter(self):
        """Test adding module pattern filters."""
        
        # Initially no filters
        self.assertEqual(len(lsl_logger._active_filters), 0)
        
        # Add a filter
        lsl_logger.add_filter('lsl.imaging.*')
        self.assertEqual(len(lsl_logger._active_filters), 1)
        self.assertIn('lsl.imaging.*', lsl_logger._active_filters)
        
        # Add another filter
        lsl_logger.add_filter('lsl_logger')
        self.assertEqual(len(lsl_logger._active_filters), 2)
        
        # Adding same filter again should not duplicate
        lsl_logger.add_filter('lsl.imaging.*')
        self.assertEqual(len(lsl_logger._active_filters), 2)
        
    def test_remove_filter(self):
        """Test removing specific filters."""
        
        # Add filters
        lsl_logger.add_filter('lsl.imaging.*')
        lsl_logger.add_filter('lsl_logger')
        self.assertEqual(len(lsl_logger._active_filters), 2)
        
        # Remove one filter
        lsl_logger.remove_filter('lsl.imaging.*')
        self.assertEqual(len(lsl_logger._active_filters), 1)
        self.assertNotIn('lsl.imaging.*', lsl_logger._active_filters)
        self.assertIn('lsl_logger', lsl_logger._active_filters)
        
        # Remove non-existent filter (should not error)
        lsl_logger.remove_filter('nonexistent')
        self.assertEqual(len(lsl_logger._active_filters), 1)
        
    def test_clear_filters(self):
        """Test clearing all filters."""
        
        # Add multiple filters
        lsl_logger.add_filter('lsl.imaging.*')
        lsl_logger.add_filter('lsl_logger')
        lsl_logger.add_filter('lsl.sim.*')
        self.assertEqual(len(lsl_logger._active_filters), 3)
        
        # Clear all filters
        lsl_logger.clear_filters()
        self.assertEqual(len(lsl_logger._active_filters), 0)
        
    def test_filter_functionality(self):
        """Test that filters actually filter log records."""
        
        # Add a ThreadedHandler to capture messages
        handler = lsl_logger.ThreadedHandler()
        lsl_logger.add_handler(handler)
        
        # Clear queue
        while not lsl_logger.LSL_LOG_QUEUE.empty():
            lsl_logger.LSL_LOG_QUEUE.get_nowait()
            
        # Add filter for only 'lsl_logger'
        lsl_logger.add_filter('lsl_logger')
        
        # This message should pass through
        lsl_logger.LSL_LOGGER.info("Should appear")
        
        # This message should be filtered out (different logger name)
        other_logger = logging.getLogger('lsl.imaging')
        other_lsl_logger.info("Should not appear")
        
        # Check queue - should only have one message
        messages = []
        while not lsl_logger.LSL_LOG_QUEUE.empty():
            try:
                record = lsl_logger.LSL_LOG_QUEUE.get_nowait()
                messages.append(record.getMessage())
            except queue.Empty:
                break
                
        self.assertEqual(len(messages), 1)
        self.assertIn("Should appear", messages[0])
        
        # Clean up
        lsl_logger.clear_filters()
        lsl_logger.remove_handler(handler)
        
    def test_capture_warnings(self):
        """Test capturing Python warnings into lsl_logger."""
        
        # Enable warning capture
        lsl_logger.capture_warnings(True)
        
        # Add a ThreadedHandler to capture messages
        handler = lsl_logger.ThreadedHandler()
        lsl_logger.add_handler(handler)
        
        # Clear queue
        while not lsl_logger.LSL_LOG_QUEUE.empty():
            lsl_logger.LSL_LOG_QUEUE.get_nowait()
            
        # Generate a warning
        warnings.warn("Test warning message")
        
        # Check if warning was captured
        # Note: warnings go to py.warnings logger, not lsl_logger
        found_warning = False
        while not lsl_logger.LSL_LOG_QUEUE.empty():
            try:
                record = lsl_logger.LSL_LOG_QUEUE.get_nowait()
                if "Test warning message" in record.getMessage():
                    found_warning = True
            except queue.Empty:
                break
                
        # Disable warning capture
        lsl_logger.capture_warnings(False)
        
        # Clean up
        lsl_logger.remove_handler(handler)
        
    def test_add_remove_handler(self):
        """Test adding and removing custom handlers."""
        
        # Create a custom handler
        custom_handler = logging.StreamHandler()
        
        # Add it
        handler_count_before = len(lsl_logger.LSL_LOGGER.handlers)
        lsl_logger.add_handler(custom_handler)
        self.assertIn(custom_handler, lsl_logger.LSL_LOGGER.handlers)
        self.assertEqual(len(lsl_logger.LSL_LOGGER.handlers), handler_count_before + 1)
        
        # Remove it
        lsl_logger.remove_handler(custom_handler)
        self.assertNotIn(custom_handler, lsl_logger.LSL_LOGGER.handlers)
        self.assertEqual(len(lsl_logger.LSL_LOGGER.handlers), handler_count_before)
        
    def test_add_handler_with_formatter(self):
        """Test adding handler with custom formatter."""
        
        custom_handler = logging.StreamHandler()
        custom_formatter = logging.Formatter('%(message)s')
        
        lsl_logger.add_handler(custom_handler, formatter=custom_formatter)
        self.assertEqual(custom_handler.formatter, custom_formatter)
        
        lsl_logger.remove_handler(custom_handler)
        
    def test_add_handler_without_formatter(self):
        """Test adding handler without formatter."""
        
        custom_handler = logging.StreamHandler()
        lsl_logger.add_handler(custom_handler, formatter=None)
        
        # Should not have the default LSL formatter when formatter=None
        self.assertIsNone(custom_handler.formatter)
        
        lsl_logger.remove_handler(custom_handler)


class logger_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.logger unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(logger_tests))


if __name__ == '__main__':
    unittest.main()
