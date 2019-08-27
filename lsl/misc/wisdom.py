# -*- coding: utf-8 -*-

"""
Module for building and saving LSL-specific FFTW and PyFFTW wisdom.

.. versionadded:: 1.0.1
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import numpy
import pickle
from datetime import datetime

from lsl.common.paths import data as dataPath
from lsl.common.busy import BusyIndicator
from lsl.misc import _wisdom


__version__ = "0.2"
__revision__ = "$Rev$"
__all__ = ["make", "show"]


# Path to the LSL-specific FFTW wisdom file
_WISDOM_FFTW = os.path.join(dataPath, 'fftw_wisdom.txt')

# Path to the LSL-specific PyFFTW wisdom file
_WISDOM_PYFFTW = os.path.join(dataPath, 'pyfftw_wisdom.pkl')


def _make_wisdom_fftw():
    """
    Build a new set of LSL-specific FFTW wisdom.
    """
    
    bi = BusyIndicator(message="Building FFTW wisdom")
    bi.start()
    
    _wisdom.buildWisdom(_WISDOM_FFTW)
    
    bi.stop()
    
    return True


def _make_wisdom_pyfftw():
    """
    Build a new set of LSL-specific PyFFTW wisdom.
    """
    
    MAXTRANSFORM = 262144
    
    try:
        import pyfftw
        
        bi = BusyIndicator(message="Building PyFFTW wisdom")
        bi.start()
        
        # Enable the PyFFTW cache
        if not pyfftw.interfaces.cache.is_enabled():
            pyfftw.interfaces.cache.enable()
            pyfftw.interfaces.cache.set_keepalive_time(60)
            
        fftFunction = lambda x: pyfftw.interfaces.numpy_fft.fft(x, planner_effort='FFTW_PATIENT')
        ifftFunction = lambda x: pyfftw.interfaces.numpy_fft.ifft(x, planner_effort='FFTW_PATIENT', overwrite_input=True)
        
        # Setup
        fftlen = 2
        while fftlen <= MAXTRANSFORM:
            data = numpy.ones(fftlen, dtype=numpy.complex64)
            fftFunction(data)
            ifftFunction(data)
            
            fftlen *= 2
            
        fftlen = 10
        while fftlen <= MAXTRANSFORM:
            data = numpy.ones(fftlen, dtype=numpy.complex64)
            fftFunction(data)
            ifftFunction(data)
            
            fftlen *= 10
            
        fh = open(_WISDOM_PYFFTW, 'wb')
        pickle.dump(pyfftw.export_wisdom(), fh)
        fh.close()
        
        bi.stop()
        
    except ImportError:
        print("PyFFTW is not installed, skipping")
        return False
        
    return True


def make(FFTW=True, PyFFTW=True):
    """
    Build a new set of LSL-specific FFTW and, optionally, PyFFTW wisdom.
    """
    
    # FFTW
    if FFTW:
        _make_wisdom_fftw()
        
    # PyFFTW
    if PyFFTW:
        _make_wisdom_pyfftw()
        
    return True


def _show_wisdom_fftw():
    """
    List information about the current LSL-specific FFTW wisdom.
    """
    
    if not os.path.exists(_WISDOM_FFTW):
        print("No LSL-specific FFTW wisdom file found, consider running make()")
        return False
        
    fh = open(_WISDOM_FFTW, 'r')
    lines = fh.readlines()
    fh.close()
    
    print("LSL FFTW Wisdom:")
    print(" Lines: %i" % len(lines))
    print(" Size: %i bytes" % os.path.getsize(_WISDOM_FFTW))
    print(" Last Modified: %s" % datetime.utcfromtimestamp(os.stat(_WISDOM_FFTW)[8]))
    
    return True


def _show_wisdom_pyfftw():
    """
    List information about the current LSL-specific FFTW wisdom.
    """
    
    if not os.path.exists(_WISDOM_PYFFTW):
        print("No LSL-specific PyFFTW wisdom file found, consider running make()")
        return False
        
    try:
        import pyfftw
        
        fh = open(_WISDOM_PYFFTW, 'r')
        d,s,ld = pickle.load(fh)
        fh.close()
        
        d = d.split('\n')
        s = s.split('\n')
        ld = ld.split('\n')
        
        print("LSL PyFFTW Wisdom:")
        print(" Lines: %i (single) %i (double) %i (long double)" % (len(s), len(d), len(ld)))
        print(" Size: %i bytes" % os.path.getsize(_WISDOM_PYFFTW))
        print(" Last Modified: %s" % datetime.utcfromtimestamp(os.stat(_WISDOM_PYFFTW)[8]))
        
    except ImportError:
        print("PyFFTW is not installed, skipping")
        return False
    
    return True


def show(fftw=True, pyfftw=True):
    """
    List information about the current LSL-specific FFTW and PyFFTW wisdom.
    """
    
    # FFTW
    if FFTW:
        _show_wisdom_fftw()
        
    # PyFFTW
    if PyFFTW:
        _show_wisdom_pyfftw()
        
    return True
