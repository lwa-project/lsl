import sys
import queue
import logging

from lsl.common.color import uncolorfy

__version__ = '0.1'
__all__ = ['LSL_LOGGER', 'LSL_LOG_FORMAT', 'LSL_LOG_QUEUE', 'set_log_level',
           'get_log_level', 'ThreadedHandler', 'add_handler', 'remove_handler',
           'capture_warnings']


#: The LSL logger instance that users should use
LSL_LOGGER = logging.getLogger('lsl_logger')


#: Default LSL_LOGGER entry formatter
LSL_LOG_FORMAT = logging.Formatter('%(asctime)s [%(levelname)-8s] %(message)s',
                                    datefmt='%Y-%m-%d %H:%M:%S')


#: Default queue for added log entries across threads with the ThreadedHandler
LSL_LOG_QUEUE = queue.Queue()


# Basic setup of the logger and default logging level
_LSL_LOG_HANDLER = logging.NullHandler()
_LSL_LOG_HANDLER.setFormatter(LSL_LOG_FORMAT)
LSL_LOGGER.addHandler(_LSL_LOG_HANDLER)
LSL_LOGGER.setLevel(logging.INFO)


def set_log_level(logging_level):
    """
    Set the logging level for the LSL logger.
    """
    
    global LSL_LOGGER
    
    LSL_LOGGER.setLevel(logging_level)


def get_log_level():
    """
    Get the current logging level for the LSL logger.
    """
    
    global LSL_LOGGER
    
    return LSL_LOGGER.level


class ThreadedHandler(logging.Handler):
    """
    logging.Handler class that allows messages to be captured in one thread and
    displayed in another.  Useful for powering a GUI a al the CASA logger.
    """
    
    def emit(self, record):
        global LSL_LOG_QUEUE
        
        LSL_LOG_QUEUE.put(record)


def add_handler(logging_handler, formatter=LSL_LOG_FORMAT):
    """
    Add a new handler to the LSL logger.  If the `formatter` keyword is not None
    then the handler's formatter will be set before it is added.
    """
    
    global LSL_LOGGER
    
    if formatter is not None:
        logging_handler.setFormatter(formatter)
        
    LSL_LOGGER.addHandler(logging_handler)


def remove_handler(logging_handler):
    """
    Remove the specified handler from the LSL logger.
    """
    
    global LSL_LOGGER
    
    LSL_LOGGER.removeHandler(logging_handler)


class _WarningHandler(logging.Handler):
    """
    Handler class that can be added to py.warnings to get those messages added
    to the main LSL_LOGGER.
    """
    
    def emit(self, record):
        global LSL_LOGGER
        
        record.msg = uncolorfy(record.msg)
        LSL_LOGGER.handle(record)


# Setup the py.warnings logger in case the user wants to redirect those as well
_warning_logger = logging.getLogger('py.warnings')
_warning_handler = _WarningHandler()
_warning_handler.setFormatter(LSL_LOG_FORMAT)
_warning_logger.addHandler(_warning_handler)


def capture_warnings(enable_capture):
    """
    Enable/disable capturing `warning.warn()` calls into the main LSL logger.
    """
    
    logging.captureWarnings(enable_capture)
