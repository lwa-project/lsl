import sys
import logging

__all__ = ['LSL_LOGGER',]

LSL_LOGGER = logging.getLogger('lsl_logger')
LSL_WARNING_LOGGER = logging.getLogger("py.warnings")


LSL_LOG_FORMAT = logging.Formatter('%(asctime)s [%(levelname)-8s] %(message)s',
                                    datefmt='%Y-%m-%d %H:%M:%S')

_LSL_LOG_HANDLER = logging.NullHandler()
_LSL_LOG_HANDLER.setFormatter(LSL_LOG_FORMAT)

LSL_LOGGER.addHandler(_LSL_LOG_HANDLER)
LSL_WARNING_LOGGER.addHandler(_LSL_LOG_HANDLER)

LSL_LOGGER.setLevel(logging.INFO)


def set_log_level(logging_level):
    global LSL_LOGGER
    global LSL_WARNING_LOGGER
    
    LSL_LOGGER.setLevel(logging_level)
    LSL_WARNING_LOGGER.setLevel(logging_level)


def add_handler(logging_handler, formatter=LSL_LOG_FORMAT):
    global LSL_LOGGER
    global LSL_WARNING_LOGGER
    
    logging_handler.setFormatter(formatter)
    
    LSL_LOGGER.addHandler(logging_handler)
    LSL_WARNING_LOGGER.addHandler(logging_handler)


def capture_warnings(enable_capture):
    logging.captureWarnings(enable_capture)
