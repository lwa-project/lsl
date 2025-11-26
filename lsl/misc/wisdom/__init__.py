"""
Module for building and saving LSL-specific FFTW wisdom.

.. versionadded:: 1.0.1
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import

import os
import logging
import numpy
from datetime import datetime

from lsl.common.paths import WISDOM as wisdomPath
from lsl.common.busy import BusyIndicator
from lsl.misc import _wisdom

from lsl.logger import LSL_LOGGER

from lsl.misc import telemetry
telemetry.track_module()


__version__ = "0.4"
__all__ = ["make", "show"]


# Path to the LSL-specific FFTW wisdom file
_WISDOM_FFTW = os.path.join(wisdomPath, 'fftw_wisdom.txt')


def make():
    """
    Build a new set of LSL-specific FFTW wisdom.
    """
    
    bi = BusyIndicator(message="Building FFTW wisdom")
    bi.start()
    
    _wisdom.buildWisdom(_WISDOM_FFTW)
    
    bi.stop()
    
    return True


def show():
    """
    List information about the current LSL-specific FFTW wisdom.
    """
    
    if not os.path.exists(_WISDOM_FFTW):
        print("No LSL-specific FFTW wisdom file found, consider running 'python -m lsl.misc.wisdom'")
        LSL_LOGGER.warning("No LSL-specific FFTW wisdom file found, consider running 'python -m lsl.misc.wisdom'")
        return False
        
    fh = open(_WISDOM_FFTW, 'r')
    lines = fh.readlines()
    fh.close()
    
    print("LSL FFTW Wisdom:")
    print(f" Lines: {len(lines)}")
    print(f" Size: {os.path.getsize(_WISDOM_FFTW)} bytes")
    print(f" Last Modified: {datetime.utcfromtimestamp(os.stat(_WISDOM_FFTW)[8])}")
    LSL_LOGGER.info("LSL FFTW Wisdom:")
    LSL_LOGGER.info(f" Lines: {len(lines)}")
    LSL_LOGGER.info(f" Size: {os.path.getsize(_WISDOM_FFTW)} bytes")
    LSL_LOGGER.info(f" Last Modified: {datetime.utcfromtimestamp(os.stat(_WISDOM_FFTW)[8])}")
    
    return True
