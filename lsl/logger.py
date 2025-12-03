import sys
import queue
import fnmatch
import logging
import warnings

__version__ = '0.1'
__all__ = ['LSL_LOGGER', 'LSL_LOG_FORMAT', 'LSL_LOG_QUEUE', 'set_log_level',
           'get_log_level', 'ThreadedHandler', 'add_handler', 'remove_handler',
           'capture_warnings', 'enable_console_logging', 'disable_console_logging',
           'enable_file_logging', 'disable_file_logging', 'add_filter',
           'remove_filter', 'clear_filters']


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

# Track console and file handlers
_console_handler = None
_file_handler = None
_active_filters = {}


def set_log_level(logging_level):
    """
    Set the logging level for the LSL logger.

    Parameters
    ----------
    logging_level : int
        The logging level (e.g., logging.DEBUG, logging.INFO, logging.WARNING).
    """
    
    global LSL_LOGGER
    
    LSL_LOGGER.setLevel(logging_level)


def get_log_level():
    """
    Get the current logging level for the LSL logger.

    Returns
    -------
    int
        The current logging level.
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

    Parameters
    ----------
    logging_handler : logging.Handler
        The handler to add to the logger.
    formatter : logging.Formatter, optional
        The formatter to apply to the handler. If None, no formatter is set.
    """
    
    global LSL_LOGGER
    
    if formatter is not None:
        logging_handler.setFormatter(formatter)
        
    LSL_LOGGER.addHandler(logging_handler)


def remove_handler(logging_handler):
    """
    Remove the specified handler from the LSL logger.

    Parameters
    ----------
    logging_handler : logging.Handler
        The handler to remove from the logger.
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
        
        from lsl.common.color import uncolorfy
        
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

    Parameters
    ----------
    enable_capture : bool
        If True, capture warnings into the logger. If False, disable capture.
    """

    logging.captureWarnings(enable_capture)


def enable_console_logging(level=None, stream=None):
    """
    Enable logging to the console (stdout/stderr).

    Parameters
    ----------
    level : int, optional
        The logging level for console output. If None, inherits from the logger's level.
    stream : file-like object, optional
        The stream to write to. Defaults to sys.stderr.
    """

    global _console_handler, LSL_LOGGER

    # Remove existing console handler if present
    if _console_handler is not None:
        LSL_LOGGER.removeHandler(_console_handler)

    # Create and configure new console handler
    if stream is None:
        stream = sys.stderr

    _console_handler = logging.StreamHandler(stream)
    _console_handler.setFormatter(LSL_LOG_FORMAT)

    if level is not None:
        _console_handler.setLevel(level)

    LSL_LOGGER.addHandler(_console_handler)


def disable_console_logging():
    """
    Disable logging to the console.
    """

    global _console_handler, LSL_LOGGER

    if _console_handler is not None:
        LSL_LOGGER.removeHandler(_console_handler)
        _console_handler = None


def enable_file_logging(filename, level=None, mode='a'):
    """
    Enable logging to a file.

    Parameters
    ----------
    filename : str
        The path to the log file.
    level : int, optional
        The logging level for file output. If None, inherits from the logger's level.
    mode : str, optional
        The file mode ('a' for append, 'w' for overwrite). Defaults to 'a'.
    """

    global _file_handler, LSL_LOGGER

    # Remove existing file handler if present
    if _file_handler is not None:
        warnings.warn(f"Closing existing log file '{_file_handler.baseFilename}' and switching to '{filename}'", RuntimeWarning)
        LSL_LOGGER.removeHandler(_file_handler)
        _file_handler.close()

    # Create and configure new file handler
    _file_handler = logging.FileHandler(filename, mode=mode)
    _file_handler.setFormatter(LSL_LOG_FORMAT)

    if level is not None:
        _file_handler.setLevel(level)

    LSL_LOGGER.addHandler(_file_handler)


def disable_file_logging():
    """
    Disable logging to a file and close the file handler.
    """

    global _file_handler, LSL_LOGGER

    if _file_handler is not None:
        LSL_LOGGER.removeHandler(_file_handler)
        _file_handler.close()
        _file_handler = None


class _ModuleFilter(logging.Filter):
    """
    Filter log records based on module name pattern matching.
    """

    def __init__(self, pattern):
        super().__init__()
        self.pattern = pattern

    def filter(self, record):
        """
        Return True if the record's module name matches the pattern.
        Supports wildcards: '*' matches any sequence, '?' matches one character.
        """
        
        return fnmatch.fnmatch(record.name, self.pattern)


def add_filter(pattern):
    """
    Add a filter to the logger based on a module name pattern.

    Only log records matching the pattern will be processed. Supports shell-style
    wildcards: '*' matches any sequence, '?' matches one character.

    Parameters
    ----------
    pattern : str
        The module name pattern to filter on (e.g., 'lsl.imaging.*', 'lsl.sim.vis').

    Examples
    --------
    >>> from lsl.logger import add_filter
    >>> add_filter('lsl.imaging.*')  # Only show logs from imaging modules
    >>> add_filter('lsl_logger')     # Only show logs from the main logger
    """

    global _active_filters, LSL_LOGGER

    if pattern not in _active_filters:
        filter_obj = _ModuleFilter(pattern)
        _active_filters[pattern] = filter_obj
        LSL_LOGGER.addFilter(filter_obj)


def remove_filter(pattern):
    """
    Remove a filter from the logger.

    Parameters
    ----------
    pattern : str
        The module name pattern to remove.
    """

    global _active_filters, LSL_LOGGER

    if pattern in _active_filters:
        filter_obj = _active_filters.pop(pattern)
        LSL_LOGGER.removeFilter(filter_obj)


def clear_filters():
    """
    Remove all filters from the logger.
    """

    global _active_filters, LSL_LOGGER

    for filter_obj in _active_filters.values():
        LSL_LOGGER.removeFilter(filter_obj)

    _active_filters.clear()
