"""
Module that contains common values found in the DP ICD, revision I.  The values 
are:
 * f_S - Sampling rate in samples per second
 * T - Slot duration in seconds
 * T_2 - Sub-slot duration
 * N_MAX_UDP - Maximum UDP packet size

Also included are two functions to convert between frequencies and DP tuning 
words and functions for calculating the magnitude response of the TBN and DRX 
filters and a software version of DP.
"""

import numpy as np
import concurrent.futures as cf
from scipy.signal import freqz, lfilter
from scipy.interpolate import interp1d

from lsl.common import _fir
from lsl.common.data_access import DataAccess

from lsl.misc import telemetry
telemetry.track_module()

from typing import Callable, List, Optional, Tuple


__version__ = '0.7'
__all__ = ['fS', 'T', 'T2', 'N_MAX', 'TBN_TUNING_WORD_MIN', 'TBN_TUNING_WORD_MAX',
           'DRX_TUNING_WORD_MIN', 'DRX_TUNING_WORD_MAX', 'DRX_BEAMS_MAX', 
           'freq_to_word', 'word_to_freq', 'delay_to_dpd', 'dpd_to_delay', 
           'gain_to_dpg', 'dpg_to_gain', 'tbn_filter', 'drx_filter', 'SoftwareDP']

#: Sample rate in Hz that is the basis for LWA time tags
fS = 196.0e6	# Hz

T = 1.0		# seconds
T2 = 0.010	# seconds
N_MAX = 8192	# bytes

#: Minimum TBN tuning word
TBN_TUNING_WORD_MIN = 109565492        # Tuning word

#: Maximum TBN tuning word
TBN_TUNING_WORD_MAX = 2037918156       # Tuning word

#: Minimum DRX tuning word
DRX_TUNING_WORD_MIN = 219130984        # Tuning word

#: Maximum DRX tuning word
DRX_TUNING_WORD_MAX = 1928352663       # Tuning word

#: Maximum number of beams
DRX_BEAMS_MAX = 4

with DataAccess.open('digital/dp_coeffs.npz', 'rb') as fh:
    dataDict = np.load(fh)
    
    # CIC Filters
    ## TBN CIC filter #7 with order 2, decimation by 98
    _TBN_CIC_7 = dataDict['TBN_CIC_7'][...]
    ## TBN CIC filter #6 with order 2, decimation by 196
    _TBN_CIC_6 = dataDict['TBN_CIC_6'][...]
    ## TBN CIC filter #5 with order 2, decimation by 392
    _TBN_CIC_5 = dataDict['TBN_CIC_5'][...]
    
    ## DRX CIC filter #7 with order 5, decimation by 5
    _DRX_CIC_7 = dataDict['DRX_CIC_7'][...]
    ## DRX CIC filter #6 with order 5, decimation by 10
    _DRX_CIC_6 = dataDict['DRX_CIC_6'][...]
    ## DRX CIC filter #5 with order 5, decimation by 20
    _DRX_CIC_5 = dataDict['DRX_CIC_5'][...]
    ## DRX CIC filter #4 with order 5, decimation by 49
    _DRX_CIC_4 = dataDict['DRX_CIC_4'][...]
    ## DRX CIC filter #3 with order 5, decimation by 98
    _DRX_CIC_3 = dataDict['DRX_CIC_3'][...]
    
    # FIR Filters
    ## Default beamformer delay FIR filters
    _DELAY_FIRS = dataDict['DELAY_FIRS'][...]
    
    ## TBN FIR filter with decimation of 20
    _TBN_FIR = dataDict['TBN_FIR'][...]
    
    ## DRX FIR filter with decimation of 2
    _DRX_FIR = dataDict['DRX_FIR'][...]
    
    try:
        dataDict.close()
    except AttributeError:
        pass
        
_N_PTS = 1000 # Number of points to use in calculating the bandpasses


def freq_to_word(freq: float) -> int:
    """
    Given a frequency in Hz, convert it to the closest DP tuning word.
    """
    
    return int(round(freq*2**32 / fS))


def word_to_freq(word: int) -> float:
    """
    Given a DP tuning word, convert it to a frequncy in Hz.
    """
    
    return word*fS / 2**32


def delay_to_dpd(delay: float) -> int:
    """
    Given a delay in ns, convert it to a course and fine portion and into the 
    final format expected by DP (big endian 16.12 unsigned integer).
    """
    
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


def dpd_to_delay(combined: int) -> float:
    """
    Given a delay value in the final format expect by DP, return the delay in ns.
    """
    
    # Convert to little-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    # Split
    fine = combined & 15
    course = (combined >> 4) & 4095
    
    # Convert to time
    delay = (course + fine/16.0) * (1e9/fS)
    
    return delay


def gain_to_dpg(gain: float) -> int:
    """
    Given a gain (between 0 and 1), convert it to a gain in the final form 
    expected by DP (big endian 16.1 signed integer).
    """
    
    # Convert
    combined = int(32767*gain)
    
    # Convert to big-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    return combined


def dpg_to_gain(combined: int) -> float:
    """
    Given a gain value in the final format expected by DP, return the gain
    as a decimal value (0 to 1).
    """
    
    # Convert to little-endian
    combined = ((combined & 0xFF) << 8) | ((combined >> 8) & 0xFF)
    
    # Convert back
    gain = combined / 32767.0
    
    return gain


def tbn_filter(sample_rate: float=1e5, npts: int=_N_PTS) -> Callable:
    """
    Return a function that will generate the shape of a TBN filter for a given sample
    rate.
    """
    
    decimation = fS / sample_rate / 10
    decimationCIC = decimation / 2
    
    # CIC settings
    N =  2
    R = decimationCIC
    
    # Part 1 - CIC filter
    h = np.linspace(0, np.pi/decimationCIC/2, num=npts, endpoint=True)
    wCIC = (np.sin(h*R)/np.sin(h/2))**N
    wCIC[0] = (2*R)**N
    
    # Part 2 - FIR filter
    h, wFIR = freqz(_TBN_FIR, 1, npts)
    
    # Cascade
    w = np.abs(wCIC) * np.abs(wFIR)
    
    # Convert to a "real" frequency and magnitude response
    h *= fS / decimation / np.pi
    w = np.abs(w)**2
    
    # Mirror
    h = np.concatenate([-h[::-1], h[1:]])
    w = np.concatenate([ w[::-1], w[1:]])
    
    # Return the interpolating function
    return interp1d(h, w/w.max(), kind='cubic', bounds_error=False, fill_value=0.0)


def drx_filter(sample_rate: float=19.6e6, npts: int=_N_PTS) -> Callable:
    """
    Return a function that will generate the shape of a DRX filter for a given sample
    rate.
    
    Based on memo DRX0001.
    """
    
    decimation = fS / sample_rate
    decimationCIC = decimation / 2
    
    # CIC settings
    N = 5
    R = decimationCIC
        
    # Part 1 - CIC filter
    h = np.linspace(0, np.pi/decimationCIC/2, num=npts, endpoint=True)
    wCIC = (np.sin(h*R)/np.sin(h/2))**N
    wCIC[0] = (2*R)**N
    
    # Part 2 - FIR filter
    h, wFIR = freqz(_DRX_FIR, 1, npts)
    
    # Cascade
    w = np.abs(wCIC) * np.abs(wFIR)
    
    # Convert to a "real" frequency and magnitude response
    h *= fS / decimation / np.pi
    w = np.abs(w)**2
    
    # Mirror
    h = np.concatenate([-h[::-1], h[1:]])
    w = np.concatenate([w[::-1], w[1:]])
    
    # Return the interpolating function
    return interp1d(h, w/w.max(), kind='cubic', bounds_error=False, fill_value=0.0)


def _process_stream_filter(time, data, filter_pack, central_freq):
    """
    Backend worker function for SoftwareDP for actually doing the DSP filtering.
    """
    
    # Mix with the NCO
    temp = data*np.exp(-2j*np.pi*central_freq*(time/fS))

    # CIC filter + decimation
    temp = lfilter(filter_pack['CIC'], 1, temp)[::filter_pack['cicD']] / filter_pack['cicD']
    scale = np.round( np.log10(np.array(filter_pack['CIC']).sum()) / np.log10(2.0) )
    temp /= 2**scale
    
    # FIR filter + decimation
    temp = lfilter(filter_pack['FIR'], 1, temp)[::filter_pack['firD']] / filter_pack['firD']
    scale = np.round( np.log10(np.array(filter_pack['FIR']).sum()) / np.log10(2.0) )
    temp /= 2**scale
    
    return temp


class SoftwareDP(object):
    """
    Class to deal with processing TBW data after the fact like DP would.  This 
    provides a means to recreate any other DP output from a TBW capture for a 
    variety of purposes.  For example, a TBW file could be processed with the
    DRX filter 4 to create a data stream that can be correlated and imaged.
    
    .. note::
        Not all DP filters are supported by this class.  Supported filters are:
         * TBN, filters 5, 6, and 7
         * DRX, filters 3, 4, 5, 6, and 7
    
    .. versionchanged:: 0.5.2
        Added support for beamforming using the DP FIR coefficients and renamed 
        SoftwareDP.apply() to SoftwareDP.applyFilter().
    """
    
    avaliableModes = {'TBN': {7: {'totalD': 1960, 'CIC': _TBN_CIC_7, 'cicD':  98, 'FIR': _TBN_FIR, 'firD': 20},
                              6: {'totalD': 3920, 'CIC': _TBN_CIC_6, 'cicD': 196, 'FIR': _TBN_FIR, 'firD': 20},
                              5: {'totalD': 7840, 'CIC': _TBN_CIC_5, 'cicD': 392, 'FIR': _TBN_FIR, 'firD': 20},
                             }, 
                     'DRX': {7: {'totalD':   10, 'CIC': _DRX_CIC_7, 'cicD':   5, 'FIR': _DRX_FIR, 'firD':  2}, 
                             6: {'totalD':   20, 'CIC': _DRX_CIC_6, 'cicD':  10, 'FIR': _DRX_FIR, 'firD':  2}, 
                             5: {'totalD':   40, 'CIC': _DRX_CIC_5, 'cicD':  20, 'FIR': _DRX_FIR, 'firD':  2}, 
                             4: {'totalD':   98, 'CIC': _DRX_CIC_4, 'cicD':  49, 'FIR': _DRX_FIR, 'firD':  2}, 
                             3: {'totalD':  196, 'CIC': _DRX_CIC_3, 'cicD':  98, 'FIR': _DRX_FIR, 'firD':  2},
                            },}
                        
    delayFIRs: List[List] = []
    for i in range(520):
        delayFIRs.append([])
        delayFIRs[-1].extend(_DELAY_FIRS)
    
    def __init__(self, mode: str='DRX', filter: int=7, central_freq: float=74e6):
        """
        Setup DP for processing an input TBW signal.  Keywords accepted are:
        * mode -> mode of operation (DRX or TBN)
        * filter -> filter code for the given mode
        * central_freq -> tuning frequency for the output
        """
        
        # Set the mode and make sure it is valid
        if mode not in self.avaliableModes:
            raise ValueError(f"Unknown mode '{mode}'")
        self.mode = mode
        
        # Set the filter and make sure it is valid
        filter = int(filter)
        if filter not in self.avaliableModes[self.mode]:
            raise ValueError(f"Unknown or unsupported filter for {self.mode}, '{filter}'")
        self.filter = filter
        
        # Set the tuning frequency and make sure it is valid
        central_freq = float(central_freq)
        central_word = freq_to_word(central_freq)
        
        if central_word < (DRX_TUNING_WORD_MIN if self.mode == 'DRX' else TBN_TUNING_WORD_MIN) \
           or central_word > (DRX_TUNING_WORD_MAX if self.mode == 'DRX' else TBN_TUNING_WORD_MAX):
            raise ValueError(f"Central frequency of {central_freq/1e6:.2f} MHz outside the DP tuning range.")
        self.central_freq = central_freq
        
    def __str__(self):
        return f"Sofware DP: {self.mode} with filter {self.filter} at {self.central_freq/1e6:.3f} MHz"
        
    def set_mode(self, mode: str):
        """
        Set the mode of operation for the software DP instance.
        """
        
        if mode not in self.avaliableModes:
            raise ValueError(f"Unknown mode '{mode}'")
        self.mode = mode
        
    def set_filter(self, filter: int):
        """
        Set the filter code for the current mode.
        """
        
        filter = int(filter)
        
        if filter not in self.avaliableModes[self.mode]:
            raise ValueError(f"Unknown or unsupported filter for {self.mode}, '{filter}'")
        self.filter = filter
        
    def set_tuning_freq(self, central_freq: float):
        """
        Set the tuning frequency for the current setup.
        """
        
        central_freq = float(central_freq)
        central_word = freq_to_word(central_freq)
        
        if central_word < (DRX_TUNING_WORD_MIN if self.mode == 'DRX' else TBN_TUNING_WORD_MIN) \
           or central_word > (DRX_TUNING_WORD_MAX if self.mode == 'DRX' else TBN_TUNING_WORD_MAX):
            raise ValueError(f"Central frequency of {central_freq/1e6:.2f} MHz outside the DP tuning range.")
        self.central_freq = central_freq
        
    def set_delay_firs(self, channel: int, coeffs: List[List]):
        """
        Set the delay FIR coefficients for a particular channel to the list of lists 
        provided (filter set by filter coefficients).  If channel is 0, the delay FIR 
        filters for all channels are set to the provided values.  If channel is -1, 
        the delay FIR filters for all channels are set to the DP default values.
        """
        
        # Make sure we have a list of lists
        try:
            len(coeffs[0])
        except TypeError:
            raise ValueError("Expected a list of lists for the coefficients.")
        
        if channel == -1:
            self.delayFIRs = []
            for i in range(520):
                self.delayFIRs.append([])
                self.delayFIRs[-1].extend(_DELAY_FIRS)
        
        if channel == 0:
            self.delayFIRs = []
            for i in range(520):
                self.delayFIRs.append([])
                self.delayFIRs[-1].extend(coeffs)
                
        else:
            self.delayFIRs[channel-1] = coeffs
            
    def form_beam(self, antennas: List, time: np.ndarray, data: np.ndarray, course_delays: Optional[np.ndarray]=None, fine_delays: Optional[np.ndarray]=None, gains: Optional[np.ndarray]=None) -> Tuple[np.ndarray,np.ndarray]:
        """
        Process a given batch of TBW data using the provided delay and gain information to
        form a beam.  Returns a two-element tuple, one for each beam polarization.
        """
        
        filters = np.array(self.delayFIRs, dtype=np.int16)
        course  = np.array(course_delays, dtype=np.int16)
        fine    = np.array(fine_delays, dtype=np.int16)
        gain    = (np.array(gains)*32767).astype(np.int16)		
        return _fir.integerBeamformer(data, filters, course, fine, gain)
        
    def apply_filter(self, time: np.ndarray, data: np.ndarray) -> np.ndarray:
        """
        Process a given batch of TBW data using the current mode of operation.  This 
        function requires both an array of times (int64 in fS since the UNIX epoch) 
        and data (1-D or 2-D).  If 2-D data are given, the first dimension should be 
        over inputs and the second over time.
        
        .. versionchanged:: 0.5.2
            Renamed SoftwareDP.apply() to SoftwareDP.applyFilter()
        """
        
        if len(data.shape) == 1:
            # Single input
            output = _process_stream_filter(time, data, self.avaliableModes[self.mode][self.filter], self.central_freq)
        else:
            # Multiple inputs - loop over time
            output = [None for i in range(data.shape[0])]
            with cf.ProcessPoolExecutor() as tpe:
                futures = {}
                for i in range(data.shape[0]):
                    task = tpe.submit(_process_stream_filter,
                                      time, data[i,:],
                                      self.avaliableModes[self.mode][self.filter],
                                      self.central_freq)
                    futures[task] = i
                    
                for task in cf.as_completed(futures):
                    i = futures[task]
                    output[i] = task.result()
                    
            output = np.array(output)
            
        return output
