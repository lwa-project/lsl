# -*- coding: utf-8 -*-

"""
Module that contains common values found in the DP ICD, revision I.  The values 
are:
  * f_S - Sampleng rate in samples per second
  * T - Slot duration in seconds
  * T_2 - Sub-slot duration
  * N_MAX_UDP - Maximum UDP packet size
"""

__version__ = '0.2'
__revision__ = '$ Revision: 6 $'
__all__ = ['fS', 'T', 'T2', 'N_MAX', '__version__', '__revision__', '__all__']

fS = 196.0e6	# Hz
T = 1.0		# seconds
T2 = 0.010	# seconds
N_MAX = 8192	# bytes



