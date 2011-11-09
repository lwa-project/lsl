# -*- coding: utf-8 -*-

"""
Module for calculating dispersion delay due to an ionized ISM and performing
incoherent dedispersion.
"""

import numpy

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['delay', 'incoherent', '__version__', '__revision__', '__all__']


def delay(freq, dm):
	"""
	Calculate the relative delay due to dispersion over a given frequnecy
	range in Hz for a particular dispersion measure in pc cm^-3.  Return 
	the dispersive delay in seconds.
	"""
	
	# Dispersion constant in MHz^2 s / pc cm^-3
	D = 4.148808e3
	
	# Delay in s
	tDelay = dm*D*((1e6/freq)**2 - (1e6/freq.max())**2)
	
	return tDelay


def incoherent(freq, waterfall, tInt, dm):
	"""
	Given a list of frequencies in Hz, a 2-D array of spectra as a function of
	time (time by frequency), and an integration time in seconds, perform 
	incoherent dedispersion on the data.
	"""
	
	# Compute the dispersive delay for the given frequency range
	tDelay = delay(freq, dm)
	
	# Convert the delays to integration periods
	tDelay = numpy.round(tDelay / tInt)
	tDelay = tDelay.astype(numpy.int32)
	
	# Roll the various frequency bins by the right amount
	ddWaterfall = waterfall*0.0
	for i,d in enumerate(tDelay):
		ddWaterfall[:,i] = numpy.roll(waterfall[:,i], -d)
		
	# Return
	return ddWaterfall


def __taperFunction(freq):
	"""
	Taper function based Equation (1) of "Pulsar Coherent De-dispersion 
	Experiment at Urumqi Observatory" CJA&A, 2006, S2, 53.
	"""

	freqMHz = freq / 1e6
	fMHz0 = freqMHz.mean()
	fMHz1 = freqMHz - fMHz0
	BW = fMHz1.max() - fMHz1.min()

	taper = 1.0 / numpy.sqrt( 1.0 + ( numpy.abs(fMHz1) / (0.47*BW) )**80 )

	return taper


def __chirpFunction(freq, dm, taper=False):
	"""
	Chip function for coherent dedispersion for a given set of frequencies (in Hz).  
	Based on Equation (6) of "Pulsar Observations II -- Coherent Dedispersion, 
	Polarimetry, and Timing" By Stairs, I. H.
	"""
	
	freqMHz = freq / 1e6
	fMHz0 = freqMHz.mean()
	fMHz1 = freqMHz - fMHz0
	BW = fMHz1.max() - fMHz1.min()
	
	chirp = numpy.exp(-2j*numpy.pi*dm/2.41033087e-10 * (fMHz1**2/ (fMHz0**2*(fMHz0 + fMHz1))))
	if taper:
		chirp *= __taperFunction(freq)
	
	return chirp


def coherent(timeseries, centralFreq, sampleRate, dm, taper=False):
	"""
	Simple coherent dedispersion of complex-valued time-series data at a given central
	frequency and sample rate.
	"""
	
	# Roughly estimate the number of points we need to look at to do the dedispersion 
	# correctrly.  Based on the GMRT coherent dedispersion pipeline
	N = 4*(202e6/centralFreq)**3*dm*(sampleRate/1e6)**2
	N = 2**int(numpy.ceil(numpy.log10(N)/numpy.log10(2.0)))
	if N < 2048:
		N = 2048
	
	# Compute the chirp function
	freq = numpy.fft.fftfreq(N, d=1/sampleRate) + centralFreq
	chirp = __chirpFunction(freq, dm, taper=taper)
	
	# Figure out the output array size
	nSets = len(timeseries) / N
	nDM = N / 8
	out = numpy.zeros(timeseries.size, dtype=timeseries.dtype)
	
	# Go!
	for i in xrange(nSets):
		start = i*N - nDM
		if start < 0:
			start = 0
		stop = start + N
		
		dataIn = timeseries[start:stop]
		dataOut = numpy.fft.ifft( numpy.fft.fft(dataIn) * chirp )
		
		out[(nDM/2 + i*(N-nDM)):(nDM/2 + (i+1)*(N-nDM))] = dataOut[nDM/2:N-nDM/2]
	
	return out
	
	