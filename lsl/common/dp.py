# -*- coding: utf-8 -*-

"""
Module that contains common values found in the DP ICD, revision I.  The values 
are:
  * f_S - Sampleng rate in samples per second
  * T - Slot duration in seconds
  * T_2 - Sub-slot duration
  * N_MAX_UDP - Maximum UDP packet size
  
Also included are two functions to convert between frequencies and DP tuning 
words and functions for calculating the magnitude response of the TBN and DRX 
filters.
"""

import numpy
from scipy.signal import freqz
from scipy.interpolate import interp1d


__version__ = '0.4'
__revision__ = '$ Revision: 7 $'
__all__ = ['fS', 'T', 'T2', 'N_MAX', 'freq2word', 'word2freq', 'tbnFilter', 'drxFilter', '__version__', '__revision__', '__all__']

fS = 196.0e6	# Hz
T = 1.0		# seconds
T2 = 0.010	# seconds
N_MAX = 8192	# bytes


def freq2word(freq):
	"""
	Given a frequency in Hz, convert it to the closest DP tuning word.
	"""
	
	return int(round(freq/fS * 2**32))


def word2freq(word):
	"""
	Given a DP tuning word, convert it to a frequncy in Hz.
	"""
	
	return word*fS / 2**32


def tbnFilter(sampleRate=1e5):
	"""
	Return a function that will generate the shape of a TBN filter for a given sample
	rate.
	"""
	
	decimation = fS / sampleRate / 10
	
	# FIR coefficients
	tbnFIR = [-2737,    531,    516,    521,    543,    579,    625,    679,    738, 
			  798,    858,    915,    966,   1009,   1042,   1062,   1068,   1057,
			 1028,    980,    911,    820,    711,    574,    421,    247,     51,
			 -163,   -392,   -634,   -886,  -1145,  -1408,  -1671,  -1930,  -2179,
			-2414,  -2629,  -2820,  -2982,  -3110,  -3199,  -3244,  -3241,  -3186,
			-3073,  -2903,  -2671,  -2374,  -2013,  -1585,  -1092,   -532,     90,
			  775,   1517,   2314,   3160,   4051,   4980,   5941,   6927,   7929,
			 8940,   9952,  10955,  11943,  12904,  13831,  14715,  15549,  16323,
			17032,  17668,  18225,  18698,  19082,  19373,  19569,  19667,  19667,
			19569,  19373,  19082,  18698,  18225,  17668,  17032,  16323,  15549,
			14715,  13831,  12904,  11943,  10955,   9952,   8940,   7929,   6927,
			 5941,   4980,   4051,   3160,   2314,   1517,    775,     90,   -532,
			-1092,  -1585,  -2013,  -2374,  -2671,  -2903,  -3073,  -3186,  -3241,
			-3244,  -3199,  -3110,  -2982,  -2820,  -2629,  -2414,  -2179,  -1930,
			-1671,  -1408,  -1145,   -886,   -634,   -392,   -163,     51,    247,
			  421,    574,    711,    820,    911,    980,   1028,   1057,   1068,
			 1062,   1042,   1009,    966,    915,    858,    798,    738,    679,
			  625,    579,    543,    521,    516,    531,  -2737]
	
	# Part 1 - FIR filter
	h, wFIR = freqz(tbnFIR, 1, 200)
	
	# Mirror and filter magnitude
	h = numpy.concatenate([-h[::-1], h[1:]])
	h *= fS / decimation / numpy.pi
	w = numpy.concatenate([w[::-1], w[1:]])
	w = numpy.abs(w)**2
	
	# Return the interpolating function
	return interp1d(h, w/w.max(), kind='cubic')


def drxFilter(sampleRate=19.6e6):
	"""
	Return a function that will generate the shape of a DRX filter for a given sample
	rate.
	
	Based on memo DRX0001.
	"""
	
	decimation = fS / sampleRate
	decimationCIC = decimation / 2
	
	# CIC settings
	N = 3
	R = 5
	
	# FIR coefficients
	drxFIR = [ 7.6000000e+001, -5.6100000e+002, -1.3700000e+002,  5.1300000e+002,  1.4400000e+002, -7.6800000e+002, 
			-1.8700000e+002,  1.1230000e+003,  2.4300000e+002, -1.6520000e+003, -3.3200000e+002,  2.5350000e+003, 
			 5.3300000e+002, -4.3870000e+003, -1.3800000e+003,  1.1030000e+004,  1.8466000e+004,  1.1030000e+004, 
			-1.3800000e+003, -4.3870000e+003,  5.3300000e+002,  2.5350000e+003, -3.3200000e+002, -1.6520000e+003, 
			 2.4300000e+002,  1.1230000e+003, -1.8700000e+002, -7.6800000e+002,  1.4400000e+002,  5.1300000e+002, 
			-1.3700000e+002, -5.6100000e+002,  7.6000000e+001]
	     
	# Part 1 - CIC filter
	h = numpy.linspace(0, numpy.pi/decimationCIC, num=200, endpoint=False)
	wCIC = (numpy.sin(h*R)/numpy.sin(h/2))**N
	wCIC[0] = (2*R)**N
	
	# Part 2 - FIR filter
	h, wFIR = freqz(drxFIR, 1, 200)
	
	# Cascade
	w = numpy.abs(wCIC[:200]) * numpy.abs(wFIR)
	
	# Mirror and filter magnitude
	h = numpy.concatenate([-h[::-1], h[1:]])
	h *= fS / decimation / numpy.pi
	w = numpy.concatenate([w[::-1], w[1:]])
	w = numpy.abs(w)**2
	
	# Return the interpolating function
	return interp1d(h, w/w.max(), kind='cubic')
	