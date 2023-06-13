"""
Module to simulate observations made with the DP system.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import time
import numpy
from aipy import coord as aipycoord
from astropy.constants import c as speedOfLight

from lsl import astro
from lsl.common import dp as dp_common
from lsl.common import stations as lwa_common
from lsl.sim import tbn
from lsl.sim import drx
from lsl.sim import vis
from lsl.reader.tbn import FILTER_CODES as TBNFilters
from lsl.reader.drx import FILTER_CODES as DRXFilters

from lsl.misc import telemetry
telemetry.track_module()
speedOfLight = speedOfLight.to('m/s').value


__version__ = '0.5'
__all__ = ['basic_signal', 'point_source']


def _basic_tbn(fh, stands, nframes, **kwargs):
    """
    Private function for generating a basic TBN signal.
    """

    start_time = kwargs['start_time']
    filter = kwargs['filter']
    verbose = kwargs['verbose']
    noise_strength = kwargs['noise_strength']
    sample_rate = TBNFilters[filter]
    
    maxValue = 127
    samplesPerFrame = 512
    upperSpike = sample_rate / 4.0
    lowerSpike = -sample_rate / 4.0
    
    if verbose:
        print("Simulating %i frames of TBN Data @ %.2f kHz for %i stands:" % \
            (nframes, sample_rate/1e3, len(stands)))
    
    for i in range(nframes):
        if i % 1000 == 0 and verbose:
            print(" frame %i" % (i+1))
        t = int(start_time*dp_common.fS) + int(i*dp_common.fS*samplesPerFrame/sample_rate)
        tFrame = t/dp_common.fS - start_time + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sample_rate
        
        for j,stand in enumerate(stands):
            ## NB:  Stand/pol labels in the TBN data are based on the digitizer not
            ## the stand
            stand_id = (stand.digitizer - 1) // 2 + 1
            pol_id = (stand.digitizer - 1) % 2
            
            cFrame = tbn.SimFrame(stand=stand_id, pol=pol_id, central_freq=40e6, gain=20, frame_count=i+1, obs_time=t)
            cFrame.data = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
            cFrame.data += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
            cFrame.data *= maxValue*noise_strength
            cFrame.data += maxValue*numpy.exp(2j*numpy.pi*upperSpike*tFrame)
            cFrame.write_raw_frame(fh)


def _basic_drx(fh, stands, nframes, **kwargs):
    """
    Private function for generating a basic TBN signal.
    """

    start_time = kwargs['start_time']
    filter = kwargs['filter']
    ntuning = kwargs['ntuning']
    verbose = kwargs['verbose']
    noise_strength = kwargs['noise_strength']
    sample_rate = DRXFilters[filter]
    decimation = int(round(dp_common.fS / sample_rate))
    
    maxValue = 7
    samplesPerFrame = 4096
    upperSpike1 = sample_rate / 4.0
    lowerSpike1 = -sample_rate / 4.0
    upperSpike2 = sample_rate / 3.0
    lowerSpike2 = -sample_rate / 3.0

    if verbose:
        print("Simulating %i frames of DRX Data @ %.2f MHz for %i beams, %i tunings each:" % \
            (nframes, sample_rate/1e6, len(stands), ntuning))

    beams = stands
    for i in range(nframes):
        if i % 1000 == 0 and verbose:
            print(" frame %i" % i)
        t = int(start_time*dp_common.fS) + int(i*dp_common.fS*samplesPerFrame/sample_rate)
        tFrame = t/dp_common.fS - start_time + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sample_rate
        for beam in beams:
            for tune in range(1, ntuning+1):
                for pol in (0, 1):
                    if tune == 1:
                        if pol == 0:
                            spike = upperSpike1
                        else:
                            spike = lowerSpike1
                    else:
                        if pol == 0:
                            spike = lowerSpike2
                        else:
                            spike = upperSpike2
                            
                    cFrame = drx.SimFrame(beam=beam, tune=tune, pol=pol, frame_count=i+1, decimation=decimation, time_offset=0, obs_time=t, flags=0)
                    cFrame.data = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
                    cFrame.data += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
                    cFrame.data *= maxValue*noise_strength
                    cFrame.data += maxValue*numpy.exp(2j*numpy.pi*spike*tFrame)
                    cFrame.write_raw_frame(fh)


def basic_signal(fh, stands, nframes, station=lwa_common.lwa1, mode='DRX', filter=6, ntuning=2, start_time=0, noise_strength=0.1, verbose=False):
    """
    Generate a collection of frames with a basic test signal for TBN and 
    DRX.  The signals for the three modes are:
    
    TBN
     * noise + (sample_rate/4) kHz signal for x-pol. and noise + 
        (-sample_rate/4) for y-pol.

    DRX
     * noise + (sample_rate/4) kHz signal for x-pol. and noise + 
        (-sample_rate/4) for y-pol. -> tuning 1
     * noise + (-sample_rate/3) kHz signal for x-pol. and noise + 
        (sample_rate/3) for y-pol. -> tuning 2
        
    All modes need to have stands (a list of :class:`lsl.common.stations.Antenna`
    instances for TBN, a list of integer beams numbers for DRX) and number of
    frames to generate.  The TBN and DRX frames need the 'filter'
    keyword set to specify the filter width.  In addition, the 'stands' 
    argument is interpreted as beam numbers for DRX.
    
    .. versionchanged:: 0.4.4
        Added the `noise_strength` keyword to control how much noise is added to 
        the data.
        
    .. versionchanged:: 2.0.0
        Removed support for generating TBW data.
        
    .. versionchanged:: 2.1.8:
        Add the `station` keyword and documentation cleanup
        `stands` is now a list of :class:`lsl.common.stations.Antenna`
        instances for TBN
    """

    if start_time == 0:
        start_time = time.time()

    if mode == 'TBN':
        _basic_tbn(fh, stands, nframes, filter=filter, start_time=start_time, noise_strength=noise_strength, verbose=verbose)
    elif mode == 'DRX':
        _basic_drx(fh, stands, nframes, filter=filter, ntuning=ntuning, start_time=start_time, noise_strength=noise_strength, verbose=verbose)
    else:
        raise RuntimeError("Unknown observations mode: %s" % mode)


def _get_antennaarray(station, stands, utime, freqs):
    """
    Given a LWA station object, a list of stands, an observation time, and
    a list of frequencies in Hz, build an aipy AntennaArray object.
    """

    return vis.build_sim_array(station, stands, freqs/1e9, jd=astro.unix_to_utcjd(utime))


def _get_source_parameters(aa, timestamp, srcs):
    """
    Given an aipy AntennaArray object, an observation time, and aipy.src 
    object, return all of the parameters needed for a simulation.
    """
    
    # Set the time for the array
    aa.set_unixtime(timestamp)

    # Compute the source parameters
    srcs_tp = []
    srcs_mt = []
    srcs_jy = []
    srcs_fq = []
    for name in srcs:
        ## Update the source's coordinates
        src = srcs[name]
        src.compute(aa)

        ## Get parameters
        top = src.get_crds(crdsys='top', ncrd=3)	# topo. coords.
        mat = src.map							# equitorial -> topo. rotation matrix
        jys = src.get_jys()						# F_nu
        frq = aa.get_afreqs()					# nu

        ## Fix the lowest frequencies to avoid problems with the flux blowing up
        ## at nu = 0 Hz by replacing flux values below 1 MHz with the flux at 
        ## 1 MHz
        Jyat1MHz = jys[ numpy.where( numpy.abs(frq-0.001) == numpy.abs(frq-0.001).min() ) ]
        jys = numpy.where( frq >= 0.001, jys, Jyat1MHz )

        ## Filter out sources that are below the horizon or have no flux
        srcAzAlt = aipycoord.top2azalt(top)
        if srcAzAlt[1] <= 0 or jys.sum() <= 0:
            continue

        ## Save values into the source arrays
        srcs_tp.append( top )
        srcs_mt.append( mat )
        srcs_jy.append( jys )
        srcs_fq.append( frq )

    # Return the values as a dictionary
    return {'topo': srcs_tp, 'trans': srcs_mt, 'flux': srcs_jy, 'freq': srcs_fq}


def _build_signals(aa, stands, src_params, times):
    """
    Given an aipy AntennaArray, a list of stand numbers, a dictionary of source 
    parameters, and an array of times in ns, return a numpy array of the simulated 
    signals that is Nstands x Ntimes in shape
    ."""

    # Find out how many stands, srcs, and samples (times) we are working with
    Nstand = len(aa.ants)
    Ntime = len(times)
    
    # Setup a temporary array to hold the signals per stand and time.
    # This array is complex so that it can accommidate TBN data
    temp = numpy.zeros((Nstand, Ntime), dtype=numpy.complex128)

    # Loop over sources and stands to build up the signals
    for topo,trans,flux,freq in zip(src_params['topo'], src_params['trans'], src_params['flux'], src_params['freq']):
        # Random Guassian noise for seeding this source
        rv_samp = numpy.random.randn(2*freq.size).view(numpy.complex128)
        
        for j,(ant,std) in enumerate(zip(aa.ants, stands)):
            ## Zeroth, get the beam response in the direction of the current source for all frequencies
            antResponse = numpy.squeeze( ant.bm_response(topo, pol='x' if std.pol == 0 else 'y') )
            ## Create array of stand position for geometric delay calculations
            xyz = numpy.array([std.stand.x, std.stand.y, std.stand.z])

            ## First, do the geometric delay
            geoDelay = numpy.dot(topo, xyz) / speedOfLight * 1e9 # s -> ns 
            
            ## Second, do the cable delay
            delayAt1MHz = std.cable.delay(frequency=1.0e6) * 1e9 # s -> ns
            cblDelay = std.cable.delay(frequency=aa.get_afreqs()*1e9) * 1e9 # s -> ns
            ## NB: Replace the cable delays below 1 MHz with the 1 MHz value to keep the 
            ## delays from blowing up for small f
            cblDelay = numpy.where( freq >= 0.001, cblDelay, delayAt1MHz )
            
            ##Finally, load the cable gain
            cblGain = std.cable.gain(frequency=aa.get_afreqs()*1e9)
            
            ## Put it all together
            factor = numpy.sqrt(antResponse * cblGain * flux / (2*numpy.pi)) * rv_samp
            temp[j,:] += numpy.fft.ifft(factor * numpy.exp(-2j*numpy.pi*freq*(cblDelay - geoDelay)))
            
    # Scale temp to sqrt of half the FFT-Length
    temp /= numpy.sqrt(freq.size/2)
    # Done
    return temp


def _point_source_tbn(fh, stands, src, nframes, **kwargs):
    """
    Private function to build TBN point sources.
    """
    
    central_freq = kwargs['central_freq']
    filter = kwargs['filter']
    start_time = kwargs['start_time']
    verbose = kwargs['verbose']
    noise_strength = kwargs['noise_strength']
    
    sample_rate = TBNFilters[filter]
    maxValue = 127
    samplesPerFrame = 512
    freqs = (numpy.fft.fftfreq(samplesPerFrame, d=1.0/sample_rate)) + central_freq
    aa = _get_antennaarray(kwargs['station'], stands, start_time, freqs)
    
    if verbose:
        print("Simulating %i frames of TBN Data @ %.2f kHz for %i stands:" % \
            (nframes, sample_rate/1e3, len(stands)))
    
    for i in range(nframes):
        if i % 1000 == 0 and verbose:
            print(" frame %i" % (i+1))
        t = int(start_time*dp_common.fS) + int(i*dp_common.fS*samplesPerFrame/sample_rate)
        tFrame = t/dp_common.fS - start_time + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sample_rate
        
        # Get the source parameters
        src_params = _get_source_parameters(aa, t/dp_common.fS, src)
        
        # Generate the time series response of each signal at each frequency
        tdSignals = _build_signals(aa, stands, src_params, tFrame*1e9)
        
        for j,stand in enumerate(stands):
            ## NB:  Stand/pol labels in the TBN data are based on the digitizer not
            ## the stand
            stand_id = (stand.digitizer - 1) // 2 + 1
            pol_id = (stand.digitizer - 1) % 2
            
            cFrame = tbn.SimFrame(stand=stand_id, pol=pol_id, central_freq=central_freq, gain=19, frame_count=i+1, obs_time=t)
            cFrame.data = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
            cFrame.data += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
            cFrame.data *= maxValue*noise_strength
            cFrame.data += maxValue*tdSignals[j,:].astype(numpy.singlecomplex)
            
            cFrame.write_raw_frame(fh)


def point_source(fh, stands, src, nframes, station=lwa_common.lwa1, mode='TBN', central_freq=49.0e6, filter=7, start_time=0, noise_strength=0.1, verbose=False):
    """
    Generate a collection of frames with a point source signal for TBN.  
    The point source is specified as a aipy.src object.
        
    All modes need to have stands (a list of :class:`lsl.common.stations.Antenna`
    instances), a number of frames to generate, and the `filter' keyword set to
    specify the filter width.
    
    .. versionchanged:: 0.4.4
        Added the `noise_strength` keyword to control how much noise is added to 
        the data.
        
    .. versionchanged:: 2.0.0
        Removed support for generating TBW data.
        
    .. versionchanged:: 2.1.8
        Add the `station` keyword and documentation cleanup
        `stands` is now a list of :class:`lsl.common.stations.Antenna`
        instances
    """

    if start_time == 0:
        start_time = time.time()

    if mode == 'TBN':
        _point_source_tbn(fh, stands, src, nframes, station=station, central_freq=central_freq, filter=filter, start_time=start_time, noise_strength=noise_strength, verbose=verbose)
    else:
        raise RuntimeError("Unknown observations mode: %s" % mode)
