"""
Module for calculating dispersion delay due to an ionized ISM and performing
incoherent/coherent dedispersion.
"""

import numpy as np

from lsl.misc import telemetry
telemetry.track_module()

from typing import Dict, Optional, Sequence, Tuple, Union


__version__ = '0.6'
__all__ = ['delay', 'incoherent', 'get_coherent_sample_size', 'coherent']


# Dispersion constant in MHz^2 s / pc cm^-3
_D = 4.148808e3


# Coherent dedispersion N and chirp cache
_coherentCache: Dict[Tuple[float,float,float,bool,np.dtype],Tuple[int,np.ndarray]] = {}


def delay(freq: Union[float,Sequence,np.ndarray], dm:float) -> Union[float,np.ndarray]:
    """
    Calculate the relative delay due to dispersion over a given frequency
    range in Hz for a particular dispersion measure in pc cm^-3.  Return 
    the dispersive delay in seconds.
    
    .. versionchanged:: 1.1.1
        If only a single frequency is provided, the returned delay is 
        relative to infinite frequency.
    """
    
    # Validate in input frequencies
    ## Right Type?
    if not isinstance(freq, np.ndarray):
        freq = np.array(freq, ndmin=1)
    ## Right size?
    singleFreq = False
    if freq.size == 1:
        singleFreq = True
        freq = np.append(freq, np.inf)
        
    # Delay in s
    tDelay = dm*_D*((1e6/freq)**2 - (1e6/freq.max())**2)
    
    # Cleanup
    if singleFreq:
        tDelay = tDelay[0]
        
    return tDelay


def incoherent(freq: np.ndarray, waterfall: np.ndarray, tInt: float, dm: float, boundary: str='wrap', fill_value: float=np.nan) -> np.ndarray:
    """
    Given a list of frequencies in Hz, a 2-D array of spectra as a function of
    time (time by frequency), and an integration time in seconds, perform 
    incoherent dedispersion on the data.
    
    The 'boundary' keyword is used to control how the boundary is treated.  The
    two options are:
      * wrap - Wrap the time boundary for each frequency (default)
      * fill - Fill the data not within the dedispersed time range with
               the value specified by the 'fill_value' keyword
    
    .. versionchanged:: 1.0.3
        Added in options for handling the boundary conditions
    """
    
    # Validate the boundary mode
    if boundary not in ('wrap', 'fill'):
        raise ValueError(f"Unknown boundary handling type '{boundary}'")
        
    # Compute the dispersive delay for the given frequency range
    tDelay = delay(freq, dm)
    
    # Convert the delays to integration periods
    tDelay = np.round(tDelay / tInt)
    tDelay = tDelay.astype(np.int32)
    
    # Roll the various frequency bins by the right amount
    ddWaterfall = waterfall*0.0
    for i,d in enumerate(tDelay):
        ddWaterfall[:,i] = np.roll(waterfall[:,i], -d)
        if boundary == 'fill' and d > 0:
            ddWaterfall[-d:,i] = fill_value
            
    # Return
    return ddWaterfall


def get_coherent_sample_size(central_freq: float, sample_rate: float, dm: float) -> int:
    """
    Estimate the number of samples needed to successfully apply coherent 
    dedispersion to a data stream.
    """
    
    # Roughly estimate the number of points we need to look at to do the dedispersion 
    # correctly.  Based on the the relative dispersion delay between the high and low
    # ends of an observational band.
    F0 = central_freq
    BW = sample_rate

    delayBand = dm*_D *((1e6/(F0-BW/2.0))**2 - (1e6/(F0+BW/2.0))**2)	# Dispersion delay across the band
    samples = delayBand*BW									# Conversion to samples
    samples = 2**(np.ceil(np.log(samples)/np.log(2)))		# Conversion to next largest power of 2
    samples *= 2											# Correction for the 'wings' of the convolution
    
    return int(samples)


def _taper(freq: np.ndarray) -> np.ndarray:
    """
    Taper function based Equation (1) of "Pulsar Coherent De-dispersion 
    Experiment at Urumqi Observatory" CJA&A, 2006, S2, 53.
    """

    freqMHz = freq / 1e6
    fMHz0 = freqMHz.mean()
    fMHz1 = freqMHz - fMHz0
    BW = fMHz1.max() - fMHz1.min()

    taper = 1.0 / np.sqrt( 1.0 + ( np.abs(fMHz1) / (0.47*BW) )**80 )

    return taper


def _chirp(freq: np.ndarray, dm: float, taper: bool=False) -> np.ndarray:
    """
    Chip function for coherent dedispersion for a given set of frequencies (in Hz).  
    Based on Equation (6) of "Pulsar Observations II -- Coherent Dedispersion, 
    Polarimetry, and Timing" By Stairs, I. H.
    """
    
    freqMHz = freq / 1e6
    fMHz0 = freqMHz.mean()
    fMHz1 = freqMHz - fMHz0
    
    chirp = np.exp(-2j*np.pi*_D*1e6 / (fMHz0**2*(fMHz0 + fMHz1)) * dm*fMHz1**2)
    if taper:
        chirp *= _taper(freq)
    
    return chirp


def coherent(t: np.ndarray, timeseries: np.ndarray, central_freq: float, sample_rate: float, dm: float, taper: bool=False, previous_time: Optional[np.ndarray]=None, previous_data: Optional[np.ndarray]=None, next_time: Optional[np.ndarray]=None, next_data: Optional[np.ndarray]=None, enable_caching: bool=True) -> Tuple[np.ndarray,np.ndarray]:
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
        
    .. note::
        Both previous_time and previous_data need to be provided for them to be
        used in the function call.  The same is true for next_time and next_data.
        
    .. versionchanged:: 1.0.1
        Added a cache for storing the chrip function between subsequent calls
        
    .. versionchanged:: 0.6.4
        Added support for keeping track of time through the dedispersion process.
    """
    
    # Caching for N and the chirp function
    try:
        pair = (central_freq*1.0, sample_rate*1.0, dm*1.0, taper, timeseries.dtype)
        
        # Get an idea of how many samples we need to do the dedispersion correctly
        # Compute the chirp function 
        N, chirp = _coherentCache[pair]
    except KeyError:
        # Get an idea of how many samples we need to do the dedispersion correctly
        N = get_coherent_sample_size(central_freq, sample_rate, dm)
        
        # Compute the chirp function 
        freq = np.fft.fftfreq(N, d=1.0/sample_rate) + central_freq 
        chirp = _chirp(freq, dm, taper=taper) 
        chirp = chirp.astype(timeseries.dtype)
        
        # Update the cache
        if enable_caching:
            _coherentCache[pair] = (N, chirp)
            
    # Figure out the output array size
    nSets = len(timeseries) // N
    outT = np.zeros(timeseries.size, dtype=t.dtype)
    outD = np.zeros(timeseries.size, dtype=timeseries.dtype)
    
    if nSets == 0:
        raise RuntimeWarning("Too few data samples for proper dedispersion")
        
    # Go!
    for i in range(2*nSets+1):
        start = N//2*i - N//4
        stop = start + N
        
        if start < 0:
            timeIn = np.zeros(N, dtype=t.dtype)
            dataIn = np.zeros(N, dtype=timeseries.dtype)
            
            if previous_time is not None and previous_data is not None:
                try:
                    timeIn[:-start] = previous_time[start:]
                    dataIn[:-start] = previous_data[start:]
                except ValueError:
                    raise RuntimeError("Too few data samples for proper start buffering")
                    
            timeIn[-start:N] = t[0:N+start]
            dataIn[-start:N] = timeseries[0:N+start]
            
        elif stop > timeseries.size:
            timeIn = np.zeros(N, dtype=t.dtype)
            dataIn = np.zeros(N, dtype=timeseries.dtype)
            
            if next_time is not None and next_data is not None:
                try:
                    timeIn[timeseries.size-start:] = next_time[:(dataIn.size-(timeseries.size-start))]
                    dataIn[timeseries.size-start:] = next_data[:(dataIn.size-(timeseries.size-start))]
                except ValueError:
                    raise RuntimeError("Too few data samples for proper end buffering")
                    
            if start < timeseries.size:
                timeIn[0:(timeseries.size-start)] = t[start:]
                dataIn[0:(timeseries.size-start)] = timeseries[start:]
                
        else:
            timeIn = t[start:stop]
            dataIn = timeseries[start:stop]
            
        timeOut = timeIn
        dataOut = np.fft.fft( dataIn )
        dataOut *= chirp
        dataOut = np.fft.ifft( dataOut )
        
        # Get the output data ranges
        outStart  = N//2*i
        outStop   = outStart + N//2
        dataStart = N//4
        dataStop  = dataStart + N//2
        
        # Make sure we don't fall off the end of the array
        if outStop >= outD.size:
            diff = outStop - outD.size
            outStop  -= diff
            dataStop -= diff
            
        if i == 0 and previous_data is None:
            continue
        if i == (2*nSets) and next_data is None:
            continue
            
        outT[outStart:outStop] = timeOut[dataStart:dataStop]
        outD[outStart:outStop] = dataOut[dataStart:dataStop]
        
    return outT, outD
