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
	
	K = 4.149e6
	
	# Delay in ms
	tDelay = dm*K*((1e6/freq)**2 - (1e6/freq.max())**2)
	
	# Conversion of s
	tDelay /= 1000.0
	
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


def __chirpFunction(freq, dm):
	"""
	Chip function for coherent dedispersion for a given set of frequencies (in Hz).  
	Based on Equation (6) of "Pulsar Observations II -- Coherent Dedispersion, 
	Polarimetry, and Timing" By Stairs, I. H.
	"""
	
	freqMHz = freq / 1e6
	fMHz0 = freqMHz.mean()
	fMHz1 = freqMHz - fMHz0
	BW = fMHz1.max() - fMHz1.min()
	
	chirp = numpy.exp(2j*numpy.pi*dm/2.41e-10 * (fMHz1**2 / (fMHz0**2*(fMHz0 + fMHz1)))) / len(freq)
	
	return chirp


def coherent(timeseries, centralFreq, sampleRate, dm):
	"""
	Simple coherent dedispersion of complex-valued time-series data at a given central
	frequency and sample rate.
	"""
	
	# Roughly estimate the number of points we need to look at to do the dedispersion 
	# correctrly.  Based on the GMRT coherent dedispersion pipeline
	N = 4*(centralFreq/202e6)**3*dm*(sampleRate/1e6)**2
	N = 2**int(numpy.ceil(numpy.log10(N)/numpy.log10(2.0)))
	if N < 2048:
		N = 2048
	
	# Compute the chirp function
	freq = numpy.fft.fftfreq(N, d=1/sampleRate) + centralFreq
	chirp = __chirpFunction(freq, dm)
	
	# Figure out the output array size
	nSets = len(timeseries) / N
	nDM = N / 4
	out = numpy.zeros(nSets*(N-nDM), dtype=timeseries.dtype)
	
	# Go!
	for i in xrange(nSets):
		start = i*N - nDM
		if start < 0:
			start = 0
		stop = start + N
		
		dataIn = timeseries[start:stop]
		dataOut = numpy.fft.fft(dataIn) #* chirp
		dataOut = numpy.fft.ifft(dataOut)
		
		out[i*(N-nDM):(i+1)*(N-nDM)] = dataOut[nDM/2:N-nDM/2]
	
	return out
	
	