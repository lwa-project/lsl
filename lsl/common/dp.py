# -*- coding: utf-8 -*-

"""
Module that contains common values found in the DP ICD, revision I.  The values 
are:
  * f_S - Sampleng rate in samples per second
  * T - Slot duration in seconds
  * T_2 - Sub-slot duration
  * N_MAX_UDP - Maximum UDP packet size
  
Also included are two functions to convert between frequencies and DP tuning 
words.
"""

__version__ = '0.3'
__revision__ = '$ Revision: 6 $'
__all__ = ['fS', 'T', 'T2', 'N_MAX', 'freq2word', 'word2freq', '__version__', '__revision__', '__all__']

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
