# -*- coding: utf-8 -*-

"""
Module that contains common values found in the ADP ICD.  The values 
are:
  * f_S - Sampling rate in samples per second
  * f_C - Width of each correlator frequency channel
  * T - Slot duration in seconds
  * T_2 - Sub-slot duration

Also included are two functions to convert between frequencies and ADP tuning 
words.
"""

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['fS', 'fC', 'T', 'T2', 'N_MAX', 'freq2word', 'word2freq', '__version__', '__revision__', '__all__']

fS = 196.0e6	# Hz
fC = 25e3		# Hz
T = 1.0		# seconds
T2 = 0.010	# seconds


def freq2word(freq):
	"""
	Given a frequency in Hz, convert it to the closest DP tuning word.
	"""
	
	return int(round(freq*2**32 / fS))


def word2freq(word):
	"""
	Given a DP tuning word, convert it to a frequncy in Hz.
	"""
	
	return word*fS / 2**32