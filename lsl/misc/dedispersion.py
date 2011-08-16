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
