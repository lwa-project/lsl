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
__all__ = ['fS', 'fC', 'T', 'T2', 'N_MAX', 'freq2word', 'word2freq', 'delaytoDPD', 'DPDtodelay', 
		 'gaintoDPG', 'DPGtogain', '__version__', '__revision__', '__all__']

fS = 196.0e6	# Hz
fC = 25e3		# Hz
T = 1.0		# seconds
T2 = 0.010	# seconds
N_MAX = 8192	# bytes


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


def delaytoDPD(delay):
	"""Given a delay in ns, convert it to a course and fine portion and into the 
	final format expected by ADP (big endian 16.12 unsigned integer)."""
	
	# Convert the delay to a combination of FIFO delays (~5.1 ns) and 
	# FIR delays (~0.3 ns)
	sample = int(round(delay * fS * 16 / 1e9))
	course = sample // 16
	fine   = sample % 16
	
	# Combine into one value
	combined = (course << 4) | fine
	
	# Convert to big-endian
	combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
	
	return combined


def DPDtodelay(combined):
	"""Given a delay value in the final format expect by ADP, return the delay in ns."""
	
	# Convert to little-endian
	combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
	
	# Split
	fine = combined & 15
	course = (combined >> 4) & 4095
	
	# Convert to time
	delay = (course + fine/16.0) / fS
	delay *= 1e9
	
	return delay


def gaintoDPG(gain):
	"""Given a gain (between 0 and 1), convert it to a gain in the final form 
	expected by ADP (big endian 16.1 signed integer)."""
	
	# Convert
	combined = int(32767*gain)
	
	# Convert to big-endian
	combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
	
	return combined


def DPGtogain(combined):
	"""Given a gain value in the final format expected by ADP, return the gain
	as a decimal value (0 to 1)."""
	
	# Convert to little-endian
	combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
	
	# Convert back
	gain = combined / 32767.0
	
	return gain