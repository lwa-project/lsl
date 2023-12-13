"""
Python module to handle the channelization and cross-correlation of TBW and
TBN data.  The main python functions in this module are:
  * calcSpectra - calculate power spectra for a collection of signals
  * FXCorrelator - calculate cross power spectra for a collection of signals
  * FXStokes - calculate Stokes cross power spectra for a collection of signals
               both of which have been deprecated in favor of the new C extension
               based  routines listed below.

The main python/C extension functions in this module are:
  * SpecMaster - similar to calcSpectra but uses the _spec module for all 
                 computations and does not support automatic sub-integration
  * StokesMaster - similar to SpecMaster but computes all four Stokes parameters
  * FXMaster - calculate cross power spectra for a collection of signals

Each function is set up to process the signals in parallel using the 
multiprocessing module and accepts a variety of options controlling the processing
of the data, including various window functions and time averaging.

.. versionchanged:: 1.0.1
    Removed SpecMasterP.

.. versionchanged:: 1.0.0
    All of the functions here now return all 'LFFT' channels.
"""

import ephem
import numpy
from astropy.constants import c as speedOfLight
from astropy.coordinates import AltAz as AstroAltAz

from lsl.reader.base import CI8
from lsl.common import dp as dp_common
from lsl.correlator import uvutils, _spec, _stokes, _core

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '1.1'
__all__ = ['pol_to_pols', 'null_window', 'SpecMaster', 'StokesMaster', 'FXMaster', 'FXStokes']


speedOfLight = speedOfLight.to('m/s').value


def pol_to_pols(pol):
    """
    Convert a polarization string, e.g., XX/XY or RR/LL, to a numeric :class:`lsl.common.stations.Antenna`
    instance polarization.
    """
    
    pol = pol.upper()
    out = []
    for p in pol:
        if p in ('X', 'R'):
            out.append(0)
        elif p in ('Y', 'L'):
            out.append(1)
        else:
            raise RuntimeError(f"Unknown polarization code '{pol}'")
        
    return out


def null_window(L):
    """
    Default "empty" windowing function for use with the various routines.  This
    function returned a numpy array of '1's of the specified length.
    """

    return numpy.ones(L)


def SpecMaster(signals, LFFT=64, window=null_window, pfb=False, verbose=False, sample_rate=None, central_freq=0.0, clip_level=0):
    """
    A more advanced version of calcSpectra that uses the _spec C extension 
    to handle all of the P.S.D. calculations in parallel.  Returns a two-
    element tuple of the frequencies (in Hz) and PSDs in linear power/RBW.
    
    .. note::
        SpecMaster currently average all data given and does not support the
        SampleAverage keyword that calcSpectra does.
        
    .. versionchanged:: 1.2.5
        Added the 'pfb' keyword to enable a 4-tap Hamming windowed PFB.  Enabling
        this overrides the 'window' keyword.
    """
    
    # Figure out if we are working with complex (I/Q) data or only real.  This
    # will determine how the FFTs are done since the real data mirrors the pos-
    # itive and negative Fourier frequencies.
    if signals.dtype.kind == 'c' or signals.dtype == CI8:
        lFactor = 1
        doFFTShift = True
        central_freq = float(central_freq)
        if signals.dtype == CI8:
            signals = signals.view(numpy.int8)
            signals = signals.reshape(signals.shape[:-1]+(-1, 2))
    else:
        lFactor = 2
        doFFTShift = False

    # Calculate the frequencies of the FFTs.  We do this for twice the FFT length
    # because the real-valued signal only occupies the positive part of the 
    # frequency space.
    if sample_rate is None:
        sample_rate = dp_common.fS
    freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/sample_rate)
    # Deal with TBW and TBN data in the correct way
    if doFFTShift:
        freq += central_freq
        freq = numpy.fft.fftshift(freq)
    freq = freq[:LFFT]
    
    if window is null_window:
        window = None
    if window is not None and pfb:
        raise RuntimeError("Cannot use a seperate window function with the PFB")
        
    if pfb:
        func = _spec.PFBPSD
    else:
        func = _spec.FPSD
    output = func(signals, LFFT=LFFT, overlap=1, clip_level=clip_level, window=window)
    
    return (freq, output)


def StokesMaster(signals, antennas, LFFT=64, window=null_window, pfb=False, verbose=False, sample_rate=None, central_freq=0.0, clip_level=0):
    """
    Similar to SpecMaster, but accepts an array of signals and a list of 
    antennas in order to compute the PSDs for the four Stokes parameters: 
    I, Q, U, and V.  Returns a two-element tuple of the frequencies (in Hz) 
    and PSDs in linear power/RBW.  The PSD are three dimensional with 
    dimensions Stokes parameter (0=I, 1=Q, 2=U, 3=V) by stand by channel).
    
    .. versionchanged:: 1.2.5
        Added the 'pfb' keyword to enable a 4-tap Hamming windowed PFB.  Enabling
        this overrides the 'window' keyword.
    """
    
    # Figure out if we are working with complex (I/Q) data or only real.  This
    # will determine how the FFTs are done since the real data mirrors the pos-
    # itive and negative Fourier frequencies.
    if signals.dtype.kind == 'c' or signals.dtype == CI8:
        lFactor = 1
        doFFTShift = True
        central_freq = float(central_freq)
        if signals.dtype == CI8:
            signals = signals.view(numpy.int8)
            signals = signals.reshape(signals.shape[:-1]+(-1, 2))
    else:
        lFactor = 2
        doFFTShift = False
        
    # Match the X and Y stand data
    signalsIndex1 = [i for (i, a) in enumerate(antennas) if a.pol == 0]
    signalsIndex2 = [i for (i, a) in enumerate(antennas) if a.pol == 1]
    if len(signalsIndex1) != len(signalsIndex2):
        raise RuntimeError("Supplied data does not contain an equal number of X and Y signals.")

    # Calculate the frequencies of the FFTs.  We do this for twice the FFT length
    # because the real-valued signal only occupies the positive part of the 
    # frequency space.
    if sample_rate is None:
        sample_rate = dp_common.fS
    freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/sample_rate)
    # Deal with TBW and TBN data in the correct way
    if doFFTShift:
        freq += central_freq
        freq = numpy.fft.fftshift(freq)
    freq = freq[:LFFT]
    
    if window is null_window:
        window = None
    if window is not None and pfb:
        raise RuntimeError("Cannot use a seperate window function with the PFB")
        
    if pfb:
        func = _stokes.PFBPSD
    else:
        func = _stokes.FPSD
    output = func(signals[signalsIndex1], signals[signalsIndex2], LFFT=LFFT, overlap=1, clip_level=clip_level, window=window)
    
    return (freq, output)


def FXMaster(signals, antennas, LFFT=64, overlap=1, include_auto=False, verbose=False, window=null_window, pfb=False, sample_rate=None, central_freq=0.0, pol='XX', gain_correct=False, return_baselines=False, clip_level=0, phase_center='z'):
    """
    A more advanced version of FXCorrelator for TBW and TBN data.  Given an 
    2-D array of signals (stands, time-series) and an array of stands, compute 
    the cross-correlation of the data for all baselines.  Return the frequencies 
    and visibilities as a two-elements tuple.
    
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.stations` module instead of a list of stand ID
        numbers.
        
    .. versionchanged:: 1.0.0
        Added a phase-center keyword that accept a two-element tuple of 
        azimuth and elelvation (in degrees) to change where the 
        correlations are phased to
        
    .. versionchanged:: 1.1.0
        Made the 'phase_center' keyword more flexible.  It can now be either:
         * 'z' to denote the zenith,
         * a ephem.Body instances which has been computed for the observer, or
         * a two-element tuple of azimuth, elevation in degrees.
         
    .. versionchanged:: 1.2.5
        Added the 'pfb' keyword to enable a 4-tap Hamming windowed PFB.  Enabling
        this overrides the 'window' keyword.
        
    .. versionchanged:: 2.0.1
        Added support for phase_center to be an astropy.coordinates.AltAz instance
    """
    
    # Decode the polarization product into something that we can use to figure 
    # out which antennas to use for the cross-correlation
    pol1, pol2 = pol_to_pols(pol)
    
    antennas1 = [a for a in antennas if a.pol == pol1]
    signalsIndex1 = [i for (i, a) in enumerate(antennas) if a.pol == pol1]
    antennas2 = [a for a in antennas if a.pol == pol2]
    signalsIndex2 = [i for (i, a) in enumerate(antennas) if a.pol == pol2]
    
    nStands = len(antennas1)
    baselines = uvutils.get_baselines(antennas1, antennas2=antennas2, include_auto=include_auto, indicies=True)
    
    # Figure out if we are working with complex (I/Q) data or only real.  This
    # will determine how the FFTs are done since the real data mirrors the pos-
    # itive and negative Fourier frequencies.
    if signals.dtype.kind == 'c' or signals.dtype == CI8:
        lFactor = 1
        doFFTShift = True
        central_freq = float(central_freq)
        if signals.dtype == CI8:
            signals = signals.view(numpy.int8)
            signals = signals.reshape(signals.shape[:-1]+(-1, 2))
    else:
        lFactor = 2
        doFFTShift = False
        
    if sample_rate is None:
        sample_rate = dp_common.fS
    freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/sample_rate)
    if doFFTShift:
        freq += central_freq
        freq = numpy.fft.fftshift(freq)
    freq = freq[:LFFT]
    
    # Get the location of the phase center in radians and create a 
    # pointing vector
    if phase_center == 'z':
        azPC = 0.0
        elPC = numpy.pi/2.0
    else:
        if isinstance(phase_center, ephem.Body):
            azPC = phase_center.az * 1.0
            elPC = phase_center.alt * 1.0
        elif isinstance(phase_center, AstroAltAz):
            azPC = phase_center.az.radian
            elPC = phase_center.alt.radian
        else:
            azPC = phase_center[0]*numpy.pi/180.0
            elPC = phase_center[1]*numpy.pi/180.0
            
    source = numpy.array([numpy.cos(elPC)*numpy.sin(azPC), 
                    numpy.cos(elPC)*numpy.cos(azPC), 
                    numpy.sin(elPC)])
                    
    # Define the cable/signal delay caches to help correlate along and compute 
    # the delays that we need to apply to align the signals
    dlyRef = len(freq)//2
    delays1 = numpy.zeros((nStands,LFFT))
    delays2 = numpy.zeros((nStands,LFFT))
    for i in list(range(nStands)):
        xyz1 = numpy.array([antennas1[i].stand.x, antennas1[i].stand.y, antennas1[i].stand.z])
        xyz2 = numpy.array([antennas2[i].stand.x, antennas2[i].stand.y, antennas2[i].stand.z])
        
        delays1[i,:] = antennas1[i].cable.delay(freq) - numpy.dot(source, xyz1) / speedOfLight
        delays2[i,:] = antennas2[i].cable.delay(freq) - numpy.dot(source, xyz2) / speedOfLight
    if not numpy.isfinite(delays1.max()):
        delays1[numpy.where( ~numpy.isfinite(delays1) )] = delays1[numpy.where( numpy.isfinite(delays1) )].max()
    if not numpy.isfinite(delays2.max()):
        delays2[numpy.where( ~numpy.isfinite(delays2) )] = delays2[numpy.where( numpy.isfinite(delays2) )].max()
    if delays1[:,dlyRef].min() < delays2[:,dlyRef].min():
        minDelay = delays1[:,dlyRef].min()
    else:
        minDelay = delays2[:,dlyRef].min()
    delays1 -= minDelay
    delays2 -= minDelay
    
    if window is null_window:
        window = None
    if window is not None and pfb:
        raise RuntimeError("Cannot use a seperate window function with the PFB")
        
    # F - defaults to running parallel in C via OpenMP
    if pfb:
        func = _core.PFBEngine
    else:
        func = _core.FEngine
    if signals.shape[0] != len(signalsIndex1):
        signalsF1, validF1 = func(signals[signalsIndex1,:], freq, delays1, LFFT=LFFT, overlap=overlap, sample_rate=sample_rate, clip_level=clip_level, window=window)
    else:
        signalsF1, validF1 = func(signals, freq, delays1, LFFT=LFFT, overlap=overlap, sample_rate=sample_rate, clip_level=clip_level, window=window)
        
    if pol2 == pol1:
        signalsF2 = signalsF1
        validF2 = validF1
    else:
        signalsF2, validF2 = func(signals[signalsIndex2,:], freq, delays2, LFFT=LFFT, overlap=overlap, sample_rate=sample_rate, clip_level=clip_level, window=window)
        
    # X
    output = _core.XEngine2(signalsF1, signalsF2, validF1, validF2)
    if not include_auto:
        # Remove auto-correlations from the output of the X engine if we don't 
        # need them.  To do this we need to first build the full list of baselines
        # (including auto-correlations) and then prune that.
        baselinesFull = uvutils.get_baselines(antennas1, antennas2=antennas2, include_auto=True, indicies=True)
        fom = numpy.array([a1-a2 for (a1,a2) in baselinesFull])
        nonAuto = numpy.where( fom != 0 )[0]
        output = output[nonAuto,:]
        
    # Apply cable gain corrections (if needed)
    if gain_correct:
        for bl in range(output.shape[0]):
            cableGain1 = antennas1[baselines[bl][0]].cable.gain(freq)
            cableGain2 = antennas2[baselines[bl][1]].cable.gain(freq)
            
            output[bl,:] /= numpy.sqrt(cableGain1*cableGain2)
            
    # Create antenna baseline list (if needed)
    if return_baselines:
        antennaBaselines = []
        for bl in range(output.shape[0]):
            antennaBaselines.append( (antennas1[baselines[bl][0]], antennas2[baselines[bl][1]]) )
        returnValues = (antennaBaselines, freq, output)
    else:
        returnValues = (freq, output)

    return returnValues


def FXStokes(signals, antennas, LFFT=64, overlap=1, include_auto=False, verbose=False, window=null_window, pfb=False, sample_rate=None, central_freq=0.0, gain_correct=False, return_baselines=False, clip_level=0, phase_center='z'):
    """
    A more advanced version of FXCorrelator for TBW and TBN data.  Given an 
    2-D array of signals (stands, time-series) and an array of stands, compute 
    the cross-correlation of the data for all baselines.  Return the frequencies 
    and visibilities as a two-elements tuple.
    
    .. versionchanged:: 2.0.1
        Added support for phase_center to be an astropy.coordinates.AltAz instance
    
    .. versionchanged:: 1.0.0
        Added a phase_center keyword that accept a two-element tuple of 
        azimuth and elelvation (in degrees) to change where the 
        correlations are phased to
        
    .. versionchanged:: 1.1.0
        Made the 'phase_center' keyword more flexible.  It can now be either:
         * 'z' to denote the zenith,
         * a ephem.Body instances which has been computed for the observer, or
         * a two-element tuple of azimuth, elevation in degrees.
         
    .. versionchanged:: 1.2.5
        Added the 'pfb' keyword to enable a 4-tap Hamming windowed PFB.  Enabling
        this overrides the 'window' keyword.
    """
    
    # Since we want to compute Stokes parameters, we need both pols
    pol1 = 0
    pol2 = 1
    
    antennas1 = [a for a in antennas if a.pol == pol1]
    signalsIndex1 = [i for (i, a) in enumerate(antennas) if a.pol == pol1]
    antennas2 = [a for a in antennas if a.pol == pol2]
    signalsIndex2 = [i for (i, a) in enumerate(antennas) if a.pol == pol2]
    
    nStands = len(antennas1)
    baselines = uvutils.get_baselines(antennas1, antennas2=antennas2, include_auto=include_auto, indicies=True)
    
    # Figure out if we are working with complex (I/Q) data or only real.  This
    # will determine how the FFTs are done since the real data mirrors the pos-
    # itive and negative Fourier frequencies.
    if signals.dtype.kind == 'c' or signals.dtype == CI8:
        lFactor = 1
        doFFTShift = True
        central_freq = float(central_freq)
        if signals.dtype == CI8:
            signals = signals.view(numpy.int8)
            signals = signals.reshape(signals.shape[:-1]+(-1, 2))
    else:
        lFactor = 2
        doFFTShift = False

    if sample_rate is None:
        sample_rate = dp_common.fS
    freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/sample_rate)
    if doFFTShift:
        freq += central_freq
        freq = numpy.fft.fftshift(freq)
    freq = freq[:LFFT]
    
    # Get the location of the phase center in radians and create a 
    # pointing vector
    if phase_center == 'z':
        azPC = 0.0
        elPC = numpy.pi/2.0
    else:
        if isinstance(phase_center, ephem.Body):
            azPC = phase_center.az * 1.0
            elPC = phase_center.alt * 1.0
        elif isinstance(phase_center, AstroAltAz):
            azPC = phase_center.az.radian
            elPC = phase_center.alt.radian
        else:
            azPC = phase_center[0]*numpy.pi/180.0
            elPC = phase_center[1]*numpy.pi/180.0
    source = numpy.array([numpy.cos(elPC)*numpy.sin(azPC), 
                    numpy.cos(elPC)*numpy.cos(azPC), 
                    numpy.sin(elPC)])
                    
    # Define the cable/signal delay caches to help correlate along and compute 
    # the delays that we need to apply to align the signals
    dlyRef = len(freq)//2
    delays1 = numpy.zeros((nStands,LFFT))
    delays2 = numpy.zeros((nStands,LFFT))
    for i in list(range(nStands)):
        xyz1 = numpy.array([antennas1[i].stand.x, antennas1[i].stand.y, antennas1[i].stand.z])
        xyz2 = numpy.array([antennas2[i].stand.x, antennas2[i].stand.y, antennas2[i].stand.z])
        
        delays1[i,:] = antennas1[i].cable.delay(freq) - numpy.dot(source, xyz1) / speedOfLight
        delays2[i,:] = antennas2[i].cable.delay(freq) - numpy.dot(source, xyz2) / speedOfLight
    if not numpy.isfinite(delays1.max()):
        delays1[numpy.where( ~numpy.isfinite(delays1) )] = delays1[numpy.where( numpy.isfinite(delays1) )].max()
    if not numpy.isfinite(delays2.max()):
        delays2[numpy.where( ~numpy.isfinite(delays2) )] = delays2[numpy.where( numpy.isfinite(delays2) )].max()
    if delays1[:,dlyRef].min() < delays2[:,dlyRef].min():
        minDelay = delays1[:,dlyRef].min()
    else:
        minDelay = delays2[:,dlyRef].min()
    delays1 -= minDelay
    delays2 -= minDelay
    
    if window is null_window:
        window = None
    if window is not None and pfb:
        raise RuntimeError("Cannot use a seperate window function with the PFB")
        
    # F - defaults to running parallel in C via OpenMP
    if pfb:
        func = _core.PFBEngine
    else:
        func = _core.FEngine
    signalsF1, validF1 = func(signals[signalsIndex1,:], freq, delays1, LFFT=LFFT, overlap=overlap, sample_rate=sample_rate, clip_level=clip_level, window=window)
    
    signalsF2, validF2 = func(signals[signalsIndex2,:], freq, delays2, LFFT=LFFT, overlap=overlap, sample_rate=sample_rate, clip_level=clip_level, window=window)
    
    # X
    output = _stokes.XEngine3(signalsF1, signalsF2, validF1, validF2)
    if not include_auto:
        # Remove auto-correlations from the output of the X engine if we don't 
        # need them.  To do this we need to first build the full list of baselines
        # (including auto-correlations) and then prune that.
        baselinesFull = uvutils.get_baselines(antennas1, antennas2=antennas2, include_auto=True, indicies=True)
        fom = numpy.array([a1-a2 for (a1,a2) in baselinesFull])
        nonAuto = numpy.where( fom != 0 )[0]
        output = output[:,nonAuto,:]
        
    # Apply cable gain corrections (if needed)
    if gain_correct:
        for bl in range(output.shape[0]):
            cableGain1 = antennas1[baselines[bl][0]].cable.gain(freq)
            cableGain2 = antennas2[baselines[bl][1]].cable.gain(freq)
            
            output[:,bl,:] /= numpy.sqrt(cableGain1*cableGain2)
            
    # Create antenna baseline list (if needed)
    if return_baselines:
        antennaBaselines = []
        for bl in range(output.shape[1]):
            antennaBaselines.append( (antennas1[baselines[bl][0]], antennas2[baselines[bl][1]]) )
        returnValues = (antennaBaselines, freq, output)
    else:
        returnValues = (freq, output)

    return returnValues	
