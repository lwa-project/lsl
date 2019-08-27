# -*- coding: utf-8 -*-

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
"""
Module to allow for post-acquisition delay-and-sum beamforming with integer 
sample delays for TBW time series data (int_delay_and_sum) and phase-and-sum 
beamforming for TBN time series data (delayAndSum).
"""

import os
import sys
import aipy
import math
import numpy
from astropy.constants import c as speedOfLight

from lsl.common.paths import DATA as dataPath
from lsl.common import dp as dp_common


__version__ = '0.6'
__revision__ = '$Rev$'
__all__ = ['calc_delay', 'int_delay_and_sum', 'int_beam_shape', 'phase_and_sum', 'phase_beam_shape', 'circularize']


speedOfLight = speedOfLight.to('m/s').value


def _load_stand_response(freq=49.0e6):
    """
    Create an aipy.amp.beam object that holds the response for a single 
    isolated stand.  The stand response is based on NEC4 models at a variety
    of frequencies within the LWA frequency range.
    """
    
    # Read in the spherical harmonic representation of the beam distributed with
    # LSL
    dd = numpy.load(os.path.join(dataPath, 'beam-shape.npz'))
    coeffs = dd['coeffs']
    try:
        dd.close()
    except AttributeError:
        pass
        
    # Calculate how many harmonics are stored in the data set and reorder the data
    # to AIPY's liking
    deg = coeffs.shape[0]-1
    lmax = int((math.sqrt(1+8*coeffs.shape[1])-3)/2)
    beam_shapeDict = {}
    for i in range(deg+1):
        beam_shapeDict[i] = numpy.squeeze(coeffs[-1-i,:])
        
    # Build the beam object and done
    return aipy.amp.BeamAlm(numpy.array([freq/1e9]), lmax=lmax, mmax=lmax, deg=deg, nside=128, coeffs=beam_shapeDict)


def calc_delay(antennas, freq=49.0e6, azimuth=0.0, elevation=90.0):
    """
    Calculate the time delays for delay-and-sum beam forming a collection of 
    stands looking in at a particular azimuth and elevation (both in degrees).  
    A numpy array of the geometric + cable delays in seconds is returned.
    
    .. versionchanged:: 0.5.0
        Changed the computed array center to exclude stands #257 through #260
    """
    
    # Make sure the pointing coordinates make sense
    if elevation < 0 or elevation > 90:
        raise ValueError("Pointing elevation (%.2f deg) is out of range [0, 90]" % elevation)
    if azimuth < 0 or azimuth > 360:
        raise ValueError("Pointing azimuth (%.2f deg) is out of range [0, 360]" % azimuth)
        
    # Get the positions of the stands and compute the mean center of the array
    xyz = numpy.zeros((len(antennas),3))
    i = 0
    good = []
    for ant in antennas:
        if ant.stand.id <= 256:
            good.append(i)
            
        xyz[i,0] = ant.stand.x
        xyz[i,1] = ant.stand.y
        xyz[i,2] = ant.stand.z
        i += 1
        
    arrayX = xyz[good,0].mean()
    arrayY = xyz[good,1].mean()
    arrayZ = xyz[good,2].mean()
    
    # Build up a unit vector that points in the direction azimuth,elevation
    rAz = azimuth*numpy.pi/180.0
    rEl = elevation*numpy.pi/180.0
    source = numpy.array([numpy.cos(rEl)*numpy.sin(rAz), 
                        numpy.cos(rEl)*numpy.cos(rAz), 
                        numpy.sin(rEl)])
                        
    # Compute the stand positions relative to the average and loop over stands
    # to compute the time delays in seconds
    arrayXYZ = xyz - numpy.array([arrayX, arrayY, arrayZ])
    delays = numpy.zeros((len(antennas),))
    for i in list(range(len(antennas))):
        delays[i] = numpy.dot(source, arrayXYZ[i,:]) / speedOfLight
        
    # Get the cable delays for each stand and add that in as well
    for i in list(range(len(antennas))):
        delays[i] = antennas[i].cable.delay(freq) - delays[i]
        
    # Done
    return delays


def int_delay_and_sum(antennas, data, sample_rate=dp_common.fS, freq=49e6, azimuth=0.0, elevation=90.0):
    """
    Given a list of antennas and a 2-D data stream with stands enumerated
    along the first axis and time series samples along the second axis, 
    delay and sum the data stream into one beam.  The delays applied are 
    integer sample delays.  Return a 1-D numpy array of the time series data 
    associated with the formed beam.
    
    .. note:
        This task is primarily intended for use with TBW data.  The time resolution of
        TBN data, even at the highest bandwidth, has an integer sample delay of zero
        samples for all pointings.
        
    .. note:
        "Bad" antennas (those with antenna.combined_status != 33) are automatically
        excluded from the beam.
        
    .. note:
        In order for this output to be properly processed by :func:`lsl.correlate.fx.calcSpectra`,
        The data need to be convert to a 2-D array via::
            
            >>> beamData = int_delay_and_sum(stands, data)
            >>> beamData.shape = (1,)+beamData.shape
            
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.stations` module instead of a list of stand ID
        numbers.
    """
    
    # Get the stand delays and convert the delay times from seconds to samples
    delays = calc_delay(antennas, freq=freq, azimuth=azimuth, elevation=elevation)
    delays = numpy.round(delays*sample_rate).astype(numpy.int16)
    
    # Figure out the polarizations
    pols = numpy.array([ant.pol for ant in antennas])
    pol0 = numpy.where( pols == 0 )[0]
    pol1 = numpy.where( pols == 1 )[0]
    
    # Make the delays into something meaningful for the shifting of the data 
    # streams
    delays -= delays.min()
    
    # Delay and sum by looping over stands inside of looping over times
    output = numpy.zeros((2, (data.shape[1]-delays.max())), dtype=data.dtype)
    for s,p in zip(range(len(antennas)), pols):
        if antennas[s].combined_status != 33:
            continue
            
        start = delays[s]
        stop = data.shape[1] - delays.max() + start
        output[p,:] += data[s,start:stop]
        
    # Check for empty polarization data.  Always return a 2-D array for the data
    if len(pol0) == 0:
        output = output[1,:]
        output.shape = (1,) + output.shape
    elif len(pol1) == 0:
        output = output[0,:]
        output.shape = (1,) + output.shape
    else:
        pass
        
    # Done
    return output


def _int_beep_and_sweep(antennas, arrayXYZ, t, freq, azimuth, elevation, beam_shape=1.0, sample_rate=dp_common.fS, direction=(0.0, 90.0)):
    """
    Worker function for int_beam_shape that 'beep's (makes a simulated signals) and
    'sweep's (delays it appropriately).
    """
    
    # Convert from degrees to radian
    rAz = azimuth*numpy.pi/180.0
    rEl = elevation*numpy.pi/180.0
    
    # Unit vector for the currect on-sky location
    currPos = numpy.array([numpy.cos(rEl)*numpy.sin(rAz), 
                        numpy.cos(rEl)*numpy.cos(rAz), 
                        numpy.sin(rEl)])
    # Stand response in this direction
    currResponse = beam_shape
    
    # Loop over stands to build the simulated singnals
    signals = currResponse * numpy.exp(2j*numpy.pi*freq*t).astype(numpy.complex64)
    signals.shape = (1,signals.size)
    signals = numpy.repeat(signals, len(antennas), axis=0)
    delays = numpy.zeros((len(antennas),1), dtype=numpy.complex64)
    for i in xrange(len(antennas)):
        currDelay = antennas[i].cable.delay(freq) - numpy.dot(currPos, arrayXYZ[i,:]) / speedOfLight
        delays[i,0] = numpy.exp(-2j*numpy.pi*freq*currDelay)
    signals *= delays
        
    # Beamform with delay-and-sum and store the RMS result
    beamHere = int_delay_and_sum(antennas, signals, sample_rate=sample_rate, freq=freq, azimuth=direction[0], elevation=direction[1])
    
    # Reduce the array dimensions
    beamHere = beamHere[0,:]
    
    # Return
    sigHere = (numpy.abs(beamHere)**2).mean()
    return sigHere


def int_beam_shape(antennas, sample_rate=dp_common.fS, freq=49e6, azimuth=0.0, elevation=90.0, progress=False, disable_pool=False):
    """
    Given a list of antennas, compute the on-sky response of the delay-and-sum
    scheme implemented in int_delay_and_sum.  A 360x90 numpy array spanning azimuth
    and elevation is returned.
    
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.stations` module instead of a list of stand ID
        numbers.
        
    .. versionchanged:: 0.4.2
        Allowed for multiple polarization data to be delayed-and-summed correctly 
        and insured that a 2-D array is always returned (pol(s) by samples)
    """
    
    # Get the stand delays and convert the delay times from seconds to samples
    delays = calc_delay(antennas, freq=freq, azimuth=azimuth, elevation=elevation)
    delays = numpy.round(delays*sample_rate).astype(numpy.int16)
    
    # Build up a base time array, load in the cable delays, and get the stand 
    # positions for geometric delay calculations.
    t = numpy.arange(0,1500)/sample_rate
    xyz = numpy.zeros((len(antennas),3))
    i = 0
    good = []
    for ant in antennas:
        if ant.stand.id <= 256:
            good.append(i)
            
        xyz[i,0] = ant.stand.x
        xyz[i,1] = ant.stand.y
        xyz[i,2] = ant.stand.z
        i += 1
        
    arrayX = xyz[good,0].mean()
    arrayY = xyz[good,1].mean()
    arrayZ = xyz[good,2].mean()
    arrayXYZ = xyz - numpy.array([arrayX, arrayY, arrayZ])
    
    # Load in the response of a single isolated stand
    standBeam = _load_stand_response(freq)
    
    # The multiprocessing module allows for the creation of worker pools to help speed
    # things along.  If the processing module is found, use it.  Otherwise, set
    # the 'usePool' variable to false and run single threaded.
    try:
        from multiprocessing import Pool, cpu_count
        
        # To get results pack from the pool, you need to keep up with the workers.  
        # In addition, we need to keep up with which workers goes with which 
        # baseline since the workers are called asynchronously.  Thus, we need a 
        # taskList array to hold tuples of baseline ('count') and workers.
        taskPool = Pool(processes=cpu_count())
        taskList = []
        
        usePool = True
        progress = False
    except ImportError:
        usePool = False
        
    # Turn off the thread pool if we are explicitly told not to use it.
    if disable_pool:
        usePool = False
        
    # Build up the beam shape over all azimuths and elevations
    beam_shape =  numpy.zeros((360,90))
    for az in list(range(360)):
        rAz = az*numpy.pi/180.0
        for el in list(range(90)):
            rEl = el*numpy.pi/180.0
            beam_shape[az,el] = standBeam.response(aipy.coord.azalt2top(numpy.concatenate([[rAz], [rEl]])))[0][0]
            
    # Build the output array and loop over all azimuths and elevations
    output = numpy.zeros((360,90))
    for az in list(range(360)):
        rAz = az*numpy.pi/180.0
        for el in list(range(90)):
            rEl = el*numpy.pi/180.0
            
            # Display the progress meter if the `progress' keyword is set to True.  The
            # progress meter displays a `.' every 2% complete and the percentages every
            # 10%.  At 100%, `Done' is displayed.
            if progress:
                fracDone = (az*90+el) / 32400.0 * 100
                if fracDone % 10 == 0 and round(fracDone, 2) != 100:
                    sys.stdout.write("%i%%" % fracDone)
                elif round(fracDone, 2) == 100:
                    sys.stdout.write("Done\n")
                elif round(fracDone, 3) % 2 == 0:
                    sys.stdout.write(".")
                else:
                    pass
                sys.stdout.flush()
                
            if usePool:
                task = taskPool.apply_async(_int_beep_and_sweep, args=(antennas, arrayXYZ, t, freq, az, el), kwds={'beam_shape': beam_shape[az,el], 'sample_rate': sample_rate, 'direction': (azimuth, elevation)})
                taskList.append((az,el,task))
            else:
                output[az,el] = _int_beep_and_sweep(antennas, arrayXYZ, t, freq, az, el, beam_shape=beam_shape[az,el], sample_rate=sample_rate, direction=(azimuth, elevation))
                
    # If pooling... Close the pool so that it knows that no ones else is joining.
    # Then, join the workers together and wait on the last one to finish before
    # saving the results.
    if usePool:
        taskPool.close()
        taskPool.join()
        
        # This is where he taskList list comes in handy.  We now know who did what
        # when we unpack the various results
        for az,el,task in taskList:
            output[az,el] = task.get()
            
    # Done
    return output


def phase_and_sum(antennas, data, sample_rate=dp_common.fS, central_freq=49.0e6, azimuth=0.0, elevation=90.0):
    """
    Given a list of antennas and a data stream of the form stands x times, 
    delay and sum the data stream into one beam.  Return a 1-D numpy array 
    of the time series data associated with the formed beam.
    
    .. note:
        This task is intended to be used with TBN data streams.
        
    .. note:
        "Bad" antennas (those with antenna.combined_status != 33) are automatically
        excluded from the beam.
    """
    
    # Get the stand delays in seconds
    delays = calc_delay(antennas, freq=central_freq, azimuth=azimuth, elevation=elevation)
    
    # Make the delays into something meaningful for the shifting of the data 
    # streams.  Then, get the beamforming coefficients (b^l_n (a la Steve's 
    # "Fun with TBN" memo).
    delays -= delays.min()
    bln = numpy.exp(2j*numpy.pi*central_freq*delays)
    
    # Figure out the polarizations
    pols = numpy.array([ant.pol for ant in antennas])
    stat = numpy.array([ant.combined_status for ant in antennas])
    pol0 = numpy.where( (pols == 0) & (stat == 33) )[0]
    pol1 = numpy.where( (pols == 1) & (stat == 33) )[0]
    
    # Loop over stands to compute the formed beam
    output = numpy.zeros((2, data.shape[1]), dtype=numpy.complex64)
    if len(pol0):
        output[0,:] = numpy.dot(bln[pol0], data[pol0,:]) / len(pol0)
    if len(pol1):
        output[1,:] = numpy.dot(bln[pol1], data[pol1,:]) / len(pol1)
        
    # Check for empty polarization data.  Always return a 3-D array for the data
    if len(pol0) == 0:
        output = output[1,:]
        output.shape = (1,) + output.shape
    elif len(pol1) == 0:
        output = output[0,:]
        output.shape = (1,) + output.shape
    else:
        pass
        
    # Done
    return output


def _phase_beep_and_sweep(antennas, arrayXYZ, t, freq, azimuth, elevation, beam_shape=1.0, sample_rate=dp_common.fS, direction=(0.0, 90.0)):
    """
    Worker function for phase_beam_shape that 'beep's (makes a simulated signals) and
    'sweep's (phases it appropriately).
    """
    
    # Convert from degrees to radian
    rAz = azimuth*numpy.pi/180.0
    rEl = elevation*numpy.pi/180.0
        
    # Unit vector for the current on-sky location
    currPos = numpy.array([numpy.cos(rEl)*numpy.sin(rAz), 
                    numpy.cos(rEl)*numpy.cos(rAz), 
                    numpy.sin(rEl)])
    # Stand response in this direction
    currResponse = beam_shape
    
    # Loop over stands to build the simulated signals
    signals = currResponse * numpy.exp(2j*numpy.pi*freq*t).astype(numpy.complex64)
    signals.shape = (1,signals.size)
    signals = numpy.repeat(signals, len(antennas), axis=0)
    delays = numpy.zeros((len(antennas),1), dtype=numpy.complex64)
    for i in xrange(len(antennas)):
        currDelay = antennas[i].cable.delay(freq) - numpy.dot(currPos, arrayXYZ[i,:]) / speedOfLight
        delays[i,0] = numpy.exp(-2j*numpy.pi*freq*currDelay)
    signals *= delays
        
    # Beamform with delay-and-sum and store the RMS result
    beam = phase_and_sum(antennas, signals, sample_rate=sample_rate, central_freq=freq, 
                         azimuth=direction[0], elevation=direction[1])
                        
    # Return
    sigHere = (numpy.abs(beam)**2).mean()
    return sigHere


def phase_beam_shape(antennas, sample_rate=dp_common.fS, central_freq=49.0e6, azimuth=0.0, elevation=90.0, progress=False):
    """
    Given a list of antennas, compute the on-sky response of the delay-and-sum
    scheme implemented in int_delay_and_sum.  A 360x90 numpy array spanning azimuth
    and elevation is returned.
    
    .. versionchanged:: 1.2.1
        Removed the 'disable_pool' keyword since recent optimztions to 
        the function have caused the multiprocessing.Pool feature to
        actually be slower than the signle-threaded version.
        
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.stations` module instead of a list of stand ID
        numbers.
    """
    
    # Get the stand delays and convert the delay times from seconds to samples
    delays = calc_delay(antennas, freq=central_freq, azimuth=azimuth, elevation=elevation)
    delays = numpy.round(delays*sample_rate).astype(numpy.int16)
    
    # Build up a base time array, load in the cable delays, and get the stand 
    # positions for geometric delay calculations.
    t = numpy.arange(0,1500)/sample_rate
    xyz = numpy.zeros((len(antennas),3))
    i = 0
    good = []
    for ant in antennas:
        if ant.stand.id <= 256:
            good.append(i)
            
        xyz[i,0] = ant.stand.x
        xyz[i,1] = ant.stand.y
        xyz[i,2] = ant.stand.z
        i += 1
        
    arrayX = xyz[good,0].mean()
    arrayY = xyz[good,1].mean()
    arrayZ = xyz[good,2].mean()
    arrayXYZ = xyz - numpy.array([arrayX, arrayY, arrayZ])
    
    # Load in the response of a single isolated stand
    standBeam = _load_stand_response(freq=central_freq)
    
    # Build up the beam shape over all azimuths and elevations
    beam_shape =  numpy.zeros((360,90))
    for az in xrange(360):
        rAz = az*numpy.pi/180.0
        for el in xrange(90):
            rEl = el*numpy.pi/180.0
            beam_shape[az,el] = standBeam.response(aipy.coord.azalt2top(numpy.concatenate([[rAz], [rEl]])))[0][0]
            
    # Build the output array and loop over all azimuths and elevations
    output = numpy.zeros((360,90))
    for az in xrange(360):
        rAz = az*numpy.pi/180.0
        for el in xrange(90):
            rEl = el*numpy.pi/180.0
            
            # Display the progress meter if the `progress' keyword is set to True.  The
            # progress meter displays a `.' every 2% complete and the percentages every
            # 10%.  At 100%, `Done' is displayed.
            if progress:
                fracDone = (az*90+el) / 32400.0 * 100
                if fracDone % 10 == 0 and round(fracDone, 2) != 100:
                    sys.stdout.write("%i%%" % fracDone)
                elif round(fracDone, 2) == 100:
                    sys.stdout.write("Done\n")
                elif round(fracDone, 3) % 2 == 0:
                    sys.stdout.write(".")
                else:
                    pass
                sys.stdout.flush()
                
            output[az,el] = _phase_beep_and_sweep(antennas, arrayXYZ, t, central_freq, az, el, beam_shape=beam_shape[az,el], sample_rate=sample_rate, direction=(azimuth, elevation))
            
    # Done
    return output


def circularize(x, y, iau=True):
    """
    Given a 1-D Numpy array of X polarization timeseries data and Y 
    polarization timeseries data, generate the two circular polarization.
    Returns a two-element tuple of L and R.
    
    .. versionchanged:: 1.0.1
        Added the 'iau' keyword to help define the convention 
        of positive Stokes V.  V = RCP - LCP.
        
    .. versionadded:: 1.0.0
    """
    
    lrSign = 1.0
    if not iau:
        lrSign = -1.0
        
    # Compute the circular terms
    l = (x + lrSign*1.0j*y) / numpy.sqrt(2)
    r = (x - lrSign*1.0j*y) / numpy.sqrt(2)
    
    # Done
    return l, r
