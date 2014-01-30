# -*- coding: utf-8 -*-

"""
This module implements a uniform DFT filter bank for use in calculating 
spectra as an alternative to a simple FFT.  The implementation here is based 
on:  http://www.scribd.com/doc/20561850/6/Polyphase-Filter-Coef%EF%AC%81cients

.. versionchanged:: 1.0.1
	Added support for using PyFFTW instead of NumPy for the FFTs
"""

import numpy
from fx import noWindow

try:
	import os
	import pickle
	import pyfftw
	
	from lsl.common.paths import data as dataPath
	
	# Enable the PyFFTW cache
	if not pyfftw.interfaces.cache.is_enabled():
		pyfftw.interfaces.cache.enable()
		pyfftw.interfaces.cache.set_keepalive_time(60)
		
	# Read in the wisdom (if it exists)
	wisdomFilename = os.path.join(dataPath, 'pyfftw-wisdom.pkl')
	if os.path.exists(wisdomFilename):
		fh = open(wisdomFilename, 'r')
		wisdom = pickle.load(fh)
		fh.close()
		
		pyfftw.import_wisdom(wisdom)
		useWisdom = True
	else:
		useWisdom = False
		
	usePyFFTW = True
	
except ImportError:
	usePyFFTW = False
	useWisdom = False


__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['fft', 'fft2', 'fft4', 'fft8', 'fft16', 'fft32', '__version__', '__revision__', '__all__']

def __filterCoeff(N, P):
	"""
	Private function to generate the filter bank coefficients for N 
	channels using P taps.
	"""

	t = numpy.arange(N*P)
	return numpy.sinc((t - N*P/2 + 0.5)/N)/N


def fft(signal, N, P=1, window=noWindow):
	"""
	FFT-based poly-phase filter bank for creating N channels with P
	taps.  Optionally, a window function can be specified using the 
	'window' keyword.  See :mod:`lsl.correlator.fx.calcSpectra` for 
	details on using window functions.
	"""
	
	filteredSignal = signal[0:N*P]*window(N*P)*__filterCoeff(N,P)
	
	if usePyFFTW and filteredSignal.dtype in (numpy.complex64, numpy.complex128):
		dd = filteredSignal.dtype
		di = numpy.empty(N, dtype=dd)
		do = numpy.empty(N, dtype=dd)
		
		forwardPlan = pyfftw.FFTW(di1, do1, direction='FFTW_FORWARD', flags=('FFTW_ESTIMATE', 'FFTW_UNALIGNED'))
		
		fbInput = numpy.empty(N, dtype=dd)
		fbTemp = numpy.empty(N, dtype=dd)
		for i in range(0,P):
			fbInput[:] = filteredSignal[i*N:(i+1)*N]
			forwardPlan(fbInput, fbTemp)
			try:
				fbOutput += fbTemp[:]
			except NameError:
				fbOutput = fbTemp[:]
				
	else:
		fbOutput = numpy.fft.fft(filteredSignal[0:N])
		for i in range(1,P):
			fbOutput += numpy.fft.fft(filteredSignal[i*N:(i+1)*N])
			
	return fbOutput

def fft2(signal, N, window=noWindow):
	"""
	Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses two taps.
	"""

	return fft(signal, N, P=2, window=window)

def fft4(signal, N, window=noWindow):
	"""
	Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses four taps.
	"""

	return fft(signal, N, P=4, window=window)

def fft8(signal, N, window=noWindow):
	"""
	Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses eight taps.
	"""

	return fft(signal, N, P=8, window=window)

def fft16(signal, N, window=noWindow):
	"""
	Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses 16 taps.
	"""

	return fft(signal, N, P=16, window=window)

def fft32(signal, N, window=noWindow):
	"""
	Sub-type of :mod:`lsl.correlator.filterbank.fft` that uses 32 taps.
	"""

	return fft(signal, N, P=32, window=window)

