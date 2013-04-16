# -*- coding: utf-8 -*-

"""
Module for calculating dispersion delay due to an ionized ISM and performing
incoherent/coherent dedispersion.
"""

import numpy

__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['delay', 'incoherent', 'getCoherentSampleSize', 'coherent', '__version__', '__revision__', '__all__']


# Dispersion constant in MHz^2 s / pc cm^-3
_D = 4.148808e3


def delay(freq, dm):
	"""
	Calculate the relative delay due to dispersion over a given frequency
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
	# correctly.  Based on the the relative dispersion delay between the high and low
	# ends of an observational band.
	F0 = centralFreq
	BW = sampleRate

	delayBand = dm*_D *((1e6/(F0-BW/2.0))**2 - (1e6/(F0+BW/2.0))**2)	# Dispersion delay across the band
	samples = delayBand*BW									# Conversion to samples
	samples = 2**(numpy.ceil(numpy.log(samples)/numpy.log(2)))		# Conversion to next largest power of 2
	samples *= 2											# Correction for the 'wings' of the convolution
	
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


def coherent(t, timeseries, centralFreq, sampleRate, dm, taper=False, previousTime=None, previousData=None, nextTime=None, nextData=None):
	"""
	Simple coherent dedispersion of complex-valued time-series data at a given central
	frequency and sample rate.  A tapering function can also be applied to the chirp of 
	the form:
		
		:math:`\\sqrt{1 + \\left(\\frac{\\Delta f_{MHz}}{0.47 \\times \\mbox{BW}}\\right)^{80}}`, 
		
	where :math:`\\Delta f_{MHz}` is the frequency difference in MHz from the band 
	center and BW is the bandwidth in MHz.
	
	.. note::
		At the large fractional bandwidths of LWA, the window size needed for coherent 
		dedispersion can be prohibitive.  For example, at 74 MHz with 19.6 MS/s and a
		DM or 10 pc / cm^3 this function uses a window size of about 268 million points.
		
	.. versionchanged:: 0.6.4
		Added support for keeping track of time through the dedispersion process.
	"""
	
	# Get an idea of how many samples we need to do the dedispersion correctly
	N = getCoherentSampleSize(centralFreq, sampleRate, dm)
	
	# Compute the chirp function
	freq = numpy.fft.fftfreq(N, d=1/sampleRate) + centralFreq
	chirp = __chirpFunction(freq, dm, taper=taper)
	
	# Figure out the output array size
	nSets = len(timeseries) / N
	outT = numpy.zeros(timeseries.size, dtype=t.dtype)
	outD = numpy.zeros(timeseries.size, dtype=timeseries.dtype)
	
	if nSets == 0:
		RuntimeWarning("Too few data samples for proper dedispersion")
		
	# Go!
	endBuffered = False
	for i in xrange(2*nSets+1):
		start = i*N/2 - N/4
		stop = start + N
		print '@', i, start, stop, timeseries.size

		if start < 0:
			timeIn = numpy.zeros(N, dtype=t.dtype)
			dataIn = numpy.zeros(N, dtype=timeseries.dtype)
			
			if previousData is not None:
				try:
					timeIn[:-start] = previousTime[start:]
					dataIn[:-start] = previousData[start:]
				except ValueError:
					raise RuntimeError("Too few data samples for proper start buffering")
					
			timeIn[-start:N] = t[0:N+start]
			dataIn[-start:N] = timeseries[0:N+start]
			
		elif stop > timeseries.size:
			timeIn = numpy.zeros(N, dtype=t.dtype)
			dataIn = numpy.zeros(N, dtype=timeseries.dtype)
			
			if nextData is not None:
				ns = nextData.size
				df = dataIn.size - (timeseries.size-start)
				
				if endBuffered:
					raise RuntimeError("Too few data samples for proper end buffering")
					
				try:
					timeIn[timeseries.size-start:] = nextTime[:(dataIn.size-(timeseries.size-start))]
					dataIn[timeseries.size-start:] = nextData[:(dataIn.size-(timeseries.size-start))]
					
					endBuffered = True
				except ValueError:
					raise RuntimeError("Too few data samples for proper end buffering")
					
			if start < timeseries.size:
				timeIn[0:(timeseries.size-start)] = t[start:]
				dataIn[0:(timeseries.size-start)] = timeseries[start:]
				
		else:
			timeIn = t[start:stop]
			dataIn = timeseries[start:stop]
			
		timeOut = timeIn
		dataOut = numpy.fft.fft( dataIn )
		dataOut *= chirp
		dataOut = numpy.fft.ifft( dataOut )
		
		# Get the output data ranges
		outStart  = i*N/2
		outStop   = outStart + N/2
		dataStart = N/4
		dataStop  = dataStart + N/2
		
		# Make sure we don't fall off the end of the array
		if outStop >= outD.size:
			diff = outStop - outD.size
			outStop  -= diff
			dataStop -= diff
			
		if i == 0 and previousData is None:
			continue
		if i == (2*nSets) and nextData is None:
			continue
			
		outT[outStart:outStop] = timeOut[dataStart:dataStop]
		outD[outStart:outStop] = dataOut[dataStart:dataStop]
		
	return outT, outD
	
	
