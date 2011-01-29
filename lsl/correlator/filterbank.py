# -*- coding: utf-8 -*-

"""This module implements a poly-phase filter bank for use in calculating spectra
as an alternative to a simple FFT.  The implementation here is based on:
http://www.scribd.com/doc/20561850/6/Polyphase-Filter-Coef%EF%AC%81cients
"""

import numpy
from fx import noWindow

__version__ = "0.1"
__revision__ = "$ Revision: 1 $"
__all__ = ['fft', 'fft2', 'fft4', 'fft8', 'fft16', 'fft32', '__version__', '__revision__', '__all__']

def __filterCoeff(N, P):
	"""Private function to generate the filter bank coefficients for N 
	channels using P taps."""

	t = numpy.arange(N*P)
	return numpy.sinc((t - N*P/2 + 0.5)/N)/N


def fft(signal, N, P=1, window=noWindow):
	"""FFT-based poly-phase filter bank for creating N channels with P
	taps.  Optionally, a window function can be specified using the 
	'window' keyword.  See :mod:`lsl.correlator.fx.calcSpectra` for 
	details on using window functions."""

	filteredSignal = signal[0:N*P]*__filterCoeff(N,P)
	
	fbOutput = numpy.fft.fft(window(N)*filteredSignal[0:N])
	for i in range(1,P):
		fbOutput += numpy.fft.fft(window(N)*filteredSignal[i*N:(i+1)*N])

	return fbOutput

def fft2(signal, N, window=noWindow):
	"""Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses two taps."""

	return fft(signal, N, P=2, window=window)

def fft4(signal, N, window=noWindow):
	"""Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses four taps."""

	return fft(signal, N, P=4, window=window)

def fft8(signal, N, window=noWindow):
	"""Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses eight taps."""

	return fft(signal, N, P=8, window=window)

def fft16(signal, N, window=noWindow):
	"""Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses 16 taps."""

	return fft(signal, N, P=16, window=window)

def fft32(signal, N, window=noWindow):
	"""Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses 32 taps."""

	return fft(signal, N, P=32, window=window)
