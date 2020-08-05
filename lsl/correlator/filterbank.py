"""
This module implements a uniform DFT filter bank for use in calculating 
spectra as an alternative to a simple FFT.  The implementation here is based 
on:  http://www.scribd.com/doc/20561850/6/Polyphase-Filter-Coef%EF%AC%81cients
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import numpy

from lsl.correlator.fx import null_window

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
__all__ = ['fft', 'fft2', 'fft4', 'fft8', 'fft16', 'fft32']

def __filterCoeff(N, P):
    """
    Private function to generate the filter bank coefficients for N 
    channels using P taps.
    """

    t = numpy.arange(N*P)
    return numpy.sinc((t - N*P/2.0 + 0.5)/N)


def fft(signal, N, P=1, window=null_window):
    """
    FFT-based poly-phase filter bank for creating N channels with P
    taps.  Optionally, a window function can be specified using the 
    'window' keyword.  See :mod:`lsl.correlator.fx.calcSpectra` for 
    details on using window functions.
    """
    
    filteredSignal = signal[0:N*P]*window(N*P)*__filterCoeff(N, P)
    
    for i in range(0, P):
        fbTemp = numpy.fft.fft(filteredSignal[i*N:(i+1)*N])
        try:
            fbOutput += fbTemp
        except NameError:
            fbOutput = fbTemp*1.0
            
    return fbOutput

def fft2(signal, N, window=null_window):
    """
    Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses two taps.
    """

    return fft(signal, N, P=2, window=window)

def fft4(signal, N, window=null_window):
    """
    Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses four taps.
    """

    return fft(signal, N, P=4, window=window)

def fft8(signal, N, window=null_window):
    """
    Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses eight taps.
    """

    return fft(signal, N, P=8, window=window)

def fft16(signal, N, window=null_window):
    """
    Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses 16 taps.
    """

    return fft(signal, N, P=16, window=window)

def fft32(signal, N, window=null_window):
    """
    Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses 32 taps.
    """

    return fft(signal, N, P=32, window=window)
