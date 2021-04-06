"""
Module that contains common values found in the ADP ICD.  The values 
are:
 * f_S - Sampling rate in samples per second
 * f_C - Width of each correlator frequency channel
 * T - Slot duration in seconds
 * T_2 - Sub-slot duration

Also included are two functions to convert between frequencies and ADP tuning 
words and functions for calculating the magnitude response of the TBN and DRX 
filters.

.. versionchanged:: 1.2.1
    Added in functions to calculate the magnitude response of the TBN and DRX
    filters
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import numpy
from scipy.signal import freqz, lfilter
from scipy.interpolate import interp1d

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.4'
__all__ = ['fS', 'fC', 'T', 'T2', 'N_MAX', 'TBN_TUNING_WORD_MIN', 'TBN_TUNING_WORD_MAX',
           'DRX_TUNING_WORD_MIN', 'DRX_TUNING_WORD_MAX', 'DRX_BEAMS_MAX', 
           'freq_to_word', 'word_to_freq', 'delay_to_dpd', 'dpd_to_delay', 
           'gain_to_dpg', 'dpg_to_gain', 'tbn_filter', 'drx_filter']

fS = 196.0e6	# Hz
fC = 25e3		# Hz
T = 1.0		# seconds
T2 = 0.010	# seconds
N_MAX = 8192	# bytes

TBN_TUNING_WORD_MIN = 65739295         # Tuning word
TBN_TUNING_WORD_MAX = 2037918156       # Tuning word

DRX_TUNING_WORD_MIN = 222417950        # Tuning word
DRX_TUNING_WORD_MAX = 1928352663       # Tuning word

DRX_BEAMS_MAX = 3

# FIR Filters
## TBN
_TBN_FIR = [-3.22350e-05, -0.00021433,  0.0017756, -0.0044913,  0.0040327, 
             0.00735870,  -0.03218100,  0.0553980, -0.0398360, -0.0748920, 
             0.58308000,   0.58308000, -0.0748920, -0.0398360,  0.0553980,
            -0.03218100,   0.00735870,  0.0040327, -0.0044913,  0.0017756,
            -0.00021433,  -3.2235e-05]
## DRX
_DRX_FIR = [ 0.0111580, -0.0074330,  0.0085684, -0.0085984,  0.0070656, -0.0035905, 
            -0.0020837,  0.0099858, -0.0199800,  0.0316360, -0.0443470,  0.0573270, 
            -0.0696630,  0.0804420, -0.0888320,  0.0941650,  0.9040000,  0.0941650, 
            -0.0888320,  0.0804420, -0.0696630,  0.0573270, -0.0443470,  0.0316360, 
            -0.0199800,  0.0099858, -0.0020837, -0.0035905,  0.0070656, -0.0085984,  
             0.0085684, -0.0074330,  0.0111580]


_N_PTS = 1000 # Number of points to use in calculating the bandpasses


def freq_to_word(freq):
    """
    Given a frequency in Hz, convert it to the closest DP tuning word.
    """
    
    return int(round(freq*2**32 / fS))


def word_to_freq(word):
    """
    Given a DP tuning word, convert it to a frequncy in Hz.
    """
    
    return word*fS / 2**32


def delay_to_dpd(delay):
    """
    Given a delay in ns, convert it to a course and fine portion and into the 
    final format expected by ADP (big endian 16.12 unsigned integer)
    ."""
    
    # Convert the delay to a combination of FIFO delays (~5.1 ns) and 
    # FIR delays (~0.3 ns)
    sample = delay * (fS/1e9)
    course = int(sample)
    fine   = int(16*(sample - course))
    
    # Combine into one value
    combined = (course << 4) | fine
    
    # Convert to big-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    return combined


def dpd_to_delay(combined):
    """
    Given a delay value in the final format expect by ADP, return the delay in ns.
    """
    
    # Convert to little-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    # Split
    fine = combined & 15
    course = (combined >> 4) & 4095
    
    # Convert to time
    delay = (course + fine/16.0) * (1e9/fS)
    
    return delay


def gain_to_dpg(gain):
    """
    Given a gain (between 0 and 1), convert it to a gain in the final form 
    expected by ADP (big endian 16.1 signed integer).
    """
    
    # Convert
    combined = int(32767*gain)
    
    # Convert to big-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    return combined


def dpg_to_gain(combined):
    """
    Given a gain value in the final format expected by ADP, return the gain
    as a decimal value (0 to 1).
    """
    
    # Convert to little-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    # Convert back
    gain = combined / 32767.0
    
    return gain


def tbn_filter(sample_rate=1e5, npts=_N_PTS):
    """
    Return a function that will generate the shape of a TBN filter for a given sample
    rate.
    """
    
    # Part 0 - FIR filter
    h, wFIR = freqz(_TBN_FIR, 1, npts)
    w = numpy.abs(wFIR)
    
    # Convert to a "real" frequency and magnitude response
    h *= sample_rate / 2.0 / numpy.pi
    w = numpy.abs(w)**2
    
    # Mirror
    h = numpy.concatenate([-h[::-1], h[1:]])
    w = numpy.concatenate([ w[::-1], w[1:]])
    
    # Return the interpolating function
    return interp1d(h, w/w.max(), kind='cubic', bounds_error=False, fill_value=0.0)


def drx_filter(sample_rate=19.6e6, npts=_N_PTS):
    """
    Return a function that will generate the shape of a DRX filter for a given sample
    rate.
    
    Based on memo DRX0001.
    """
    # Part 0 - FIR filter
    h, wFIR = freqz(_DRX_FIR, 1, npts)
    w = numpy.abs(wFIR)
    
    # Convert to a "real" frequency and magnitude response
    h *= sample_rate / 2.0 / numpy.pi
    w = numpy.abs(w)**2
    
    # Mirror
    h = numpy.concatenate([-h[::-1], h[1:]])
    w = numpy.concatenate([w[::-1], w[1:]])
    
    # Return the interpolating function
    return interp1d(h, w/w.max(), kind='cubic', bounds_error=False, fill_value=0.0)
