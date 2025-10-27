import sys
import logging

__all__ = ['LSL_LOGGER',]

LSL_LOGGER = logging.getLogger('lsl_logger')
logging.captureWarnings(True)

LSL_WARNING_LOGGER = logging.getLogger("py.warnings")


_LSL_LOG_FORMAT = logging.Formatter('%(asctime)s [%(levelname)-8s] %(message)s',
                                    datefmt='%Y-%m-%d %H:%M:%S')

_LSL_LOG_HANDLER = logging.StreamHandler(sys.stdout)
_LSL_LOG_HANDLER.setFormatter(_LSL_LOG_FORMAT)

LSL_LOGGER.addHandler(_LSL_LOG_HANDLER)
LSL_WARNING_LOGGER.addHandler(_LSL_LOG_HANDLER)

LSL_LOGGER.setLevel(logging.INFO)
