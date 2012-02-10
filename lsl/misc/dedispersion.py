# -*- coding: utf-8 -*-

"""
Module for calculating dispersion delay due to an ionized ISM and performing
incoherent dedispersion.
"""

import numpy

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['delay', 'incoherent', 'getCoherentSampleSize', 'coherent', '__version__', '__revision__', '__all__']


# Dispersion constant in MHz^2 s / pc cm^-3
_D = 4.148808e3


def delay(freq, dm):
	"""
	Calculate the relative delay due to dispersion over a given frequnecy
	range in Hz for a particular dispersion measure in pc cm^-3.  Return 
	the dispersive delay in seconds.
	"""
	
	# Delay in s
	tDelay = dm*_D*((1e6/freq)**2 - (1e6/freq.max())**2)
	
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


def getCoherentSampleSize(centralFreq, sampleRate, dm):
	"""
	Estimate the number of samples needed to successfully apply coherent 
	dedispersion to a data stream.
	"""
	
	# Roughly estimate the number of points we need to look at to do the dedispersion 
	# correctrly.  Based on the the relative dispersion delay between the high and low
	# ends of an observational band.
	F0 = centralFreq
	BW = sampleRate

	delay = 4*dm*_D / (F0**3/BW - 2*BW*F0 + BW**3/F0) * (1e6)**2	# Dispersion delay across the band
	samples = delay*BW						# Conversion to samples
	samples = 2**(numpy.ceil(numpy.log(samples)/numpy.log(2)))	# Conversion to next largest power of 2
	samples *= 2							# Correction for the 'wings' of the convolution
	return int(samples)


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
	
	chirp = numpy.exp(-2j*numpy.pi*_D*1e6 / (fMHz0**2*(fMHz0 + fMHz1)) * dm*fMHz1**2)
	if taper:
		chirp *= __taperFunction(freq)
	
	return chirp


def coherent(timeseries, centralFreq, sampleRate, dm, taper=False):
	"""
	Simple coherent dedispersion of complex-valued time-series data at a given central
	frequency and sample rate.
	"""
	
	# Get an idea of how many samples we need to do the dedispersion correctly
	N = getCoherentSampleSize(centralFreq, sampleRate, dm)
	
	# Compute the chirp function
	freq = numpy.fft.fftfreq(N, d=1/sampleRate) + centralFreq
	chirp = __chirpFunction(freq, dm, taper=taper)
	
	# Figure out the output array size
	nSets = len(timeseries) / N
	out = numpy.zeros(timeseries.size, dtype=timeseries.dtype)
	
	# Go!
	last = False
	for i in xrange(2*nSets+1):
		start = i*N/2 - N/4
		stop = start + N

		if start < 0:
			dataIn = numpy.zeros(N, dtype=numpy.complex64)
			dataIn[-start:N] = timeseries[0:N+start]
		elif stop > timeseries.size:
			dataIn = numpy.zeros(N, dtype=numpy.complex64)
			dataIn[0:timeseries.size-start] = timeseries[start:]
			# If this is it, end after saving the data
			last = True
		else:
			dataIn = timeseries[start:stop]

		dataOut = numpy.fft.fft( dataIn )
		dataOut *= chirp
		dataOut = numpy.fft.ifft( dataOut )

		out[i*N/2:(i+1)*N/2] = dataOut[N/4:3*N/4]
		if last:
			break
	
	return out
	
	
