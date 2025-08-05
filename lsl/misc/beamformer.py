"""
Module to allow for post-acquisition delay-and-sum beamforming with integer 
sample delays for TBW time series data (int_delay_and_sum) and phase-and-sum 
beamforming for TBN time series data (delayAndSum).
"""

import os
import sys
import aipy
import ephem
import numpy as np
import concurrent.futures as cf

from astropy.constants import c as speedOfLight
from astropy.coordinates import Angle as AstroAngle

from lsl.common.data_access import DataAccess
from lsl.common import ndp as ndp_common

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.6'
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
    with DataAccess.open('antenna/beam-shape.npz', 'rb') as fh:
        dd = np.load(fh)
        coeffs = dd['coeffs']
        
        # Calculate how many harmonics are stored in the data set and reorder the data
        # to AIPY's liking
        deg = coeffs.shape[0]-1
        lmax = int((np.sqrt(1+8*coeffs.shape[1])-3)/2)
        beam_shapeDict = {}
        for i in range(deg+1):
            beam_shapeDict[i] = np.squeeze(coeffs[-1-i,:])
            
    # Build the beam object and done
    return aipy.amp.BeamAlm(np.array([freq/1e9]), lmax=lmax, mmax=lmax, deg=deg, nside=128, coeffs=beam_shapeDict)


def calc_delay(antennas, freq=49.0e6, azimuth=0.0, altitude=90.0):
    """
    Calculate the time delays for delay-and-sum beam forming a collection of 
    stands looking in at a particular azimuth and altitude (both in degrees).  
    A numpy array of the geometric + cable delays in seconds is returned.
    
    .. versionchanged:: 0.5.0
        Changed the computed array center to exclude stands #257 through #260
    """
    
    # Convert
    if isinstance(azimuth, AstroAngle):
        azimuth = azimuth.deg
    elif isinstance(azimuth, ephem.Angle):
        azimuth = azimuth * 180/np.pi
        
    if isinstance(altitude, AstroAngle):
        altitude = altitude.deg
    elif isinstance(altitude, ephem.Angle):
        altitude = altitude * 180/np.pi
        
    # Make sure the pointing coordinates make sense
    if altitude < 0 or altitude > 90:
        raise ValueError(f"Pointing altitude ({altitude:.2f} deg) is out of range [0, 90]")
    if azimuth < 0 or azimuth > 360:
        raise ValueError(f"Pointing azimuth ({azimuth:.2f} deg) is out of range [0, 360]")
        
    # Get the positions of the stands and compute the mean center of the array
    xyz = np.zeros((len(antennas),3))
    i = 0
    good = []
    for ant in antennas:
        if ant.stand.id <= 256:
            good.append(i)
            
        xyz[i,:] = ant.stand.xyz
        i += 1
        
    arrayX = xyz[good,0].mean()
    arrayY = xyz[good,1].mean()
    arrayZ = xyz[good,2].mean()
    
    # Build up a unit vector that points in the direction azimuth,altitude
    rAz = azimuth*np.pi/180.0
    rAlt = altitude*np.pi/180.0
    source = np.array([np.cos(rAlt)*np.sin(rAz), 
                       np.cos(rAlt)*np.cos(rAz), 
                       np.sin(rAlt)])
                        
    # Compute the stand positions relative to the average and loop over stands
    # to compute the time delays in seconds
    arrayXYZ = xyz - np.array([arrayX, arrayY, arrayZ])
    delays = np.zeros((len(antennas),))
    for i in list(range(len(antennas))):
        delays[i] = np.dot(source, arrayXYZ[i,:]) / speedOfLight
        
    # Get the cable delays for each stand and add that in as well
    for i in list(range(len(antennas))):
        delays[i] = antennas[i].cable.delay(freq) - delays[i]
        
    # Done
    return delays


def int_delay_and_sum(antennas, data, sample_rate=ndp_common.fS, freq=49e6, azimuth=0.0, altitude=90.0):
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
    delays = calc_delay(antennas, freq=freq, azimuth=azimuth, altitude=altitude)
    delays = np.round(delays*sample_rate).astype(np.int16)
    
    # Figure out the polarizations
    pols = np.array([ant.pol for ant in antennas])
    pol0 = np.where( pols == 0 )[0]
    pol1 = np.where( pols == 1 )[0]
    
    # Make the delays into something meaningful for the shifting of the data 
    # streams
    delays -= delays.min()
    
    # Delay and sum by looping over stands inside of looping over times
    output = np.zeros((2, (data.shape[1]-delays.max())), dtype=data.dtype)
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


def _int_beep_and_sweep(antennas, arrayXYZ, t, freq, azimuth, altitude, beam_shape=1.0, sample_rate=ndp_common.fS, direction=(0.0, 90.0)):
    """
    Worker function for int_beam_shape that 'beep's (makes a simulated signals) and
    'sweep's (delays it appropriately).
    """
    
    # Convert from degrees to radian
    rAz = azimuth*np.pi/180.0
    rAlt = altitude*np.pi/180.0
    
    # Unit vector for the currect on-sky location
    currPos = np.array([np.cos(rAlt)*np.sin(rAz), 
                        np.cos(rAlt)*np.cos(rAz), 
                        np.sin(rAlt)])
    # Stand response in this direction
    currResponse = beam_shape
    
    # Loop over stands to build the simulated singnals
    signals = currResponse * np.exp(2j*np.pi*freq*t).astype(np.complex64)
    signals.shape = (1,signals.size)
    signals = np.repeat(signals, len(antennas), axis=0)
    delays = np.zeros((len(antennas),1), dtype=np.complex64)
    for i in range(len(antennas)):
        currDelay = antennas[i].cable.delay(freq) - np.dot(currPos, arrayXYZ[i,:]) / speedOfLight
        delays[i,0] = np.exp(-2j*np.pi*freq*currDelay)
    signals *= delays
        
    # Beamform with delay-and-sum and store the RMS result
    beamHere = int_delay_and_sum(antennas, signals, sample_rate=sample_rate, freq=freq, azimuth=direction[0], altitude=direction[1])
    
    # Reduce the array dimensions
    beamHere = beamHere[0,:]
    
    # Return
    sigHere = (np.abs(beamHere)**2).mean()
    return sigHere


def int_beam_shape(antennas, sample_rate=ndp_common.fS, freq=49e6, azimuth=0.0, altitude=90.0, progress=False):
    """
    Given a list of antennas, compute the on-sky response of the delay-and-sum
    scheme implemented in int_delay_and_sum.  A 360x90 numpy array spanning azimuth
    and altitude is returned.
    
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.stations` module instead of a list of stand ID
        numbers.
        
    .. versionchanged:: 0.4.2
        Allowed for multiple polarization data to be delayed-and-summed correctly 
        and insured that a 2-D array is always returned (pol(s) by samples)
    
    .. versionchanged:: 3.0.0
        Dropped the 'disable_pool' keyword.
    """
    
    # Convert
    if isinstance(azimuth, AstroAngle):
        azimuth = azimuth.deg
    elif isinstance(azimuth, ephem.Angle):
        azimuth = azimuth * 180/np.pi
        
    if isinstance(altitude, AstroAngle):
        altitude = altitude.deg
    elif isinstance(altitude, ephem.Angle):
        altitude = altitude * 180/np.pi
        
    # Build up a base time array, load in the cable delays, and get the stand 
    # positions for geometric delay calculations.
    t = np.arange(0,1500)/sample_rate
    xyz = np.zeros((len(antennas),3))
    i = 0
    good = []
    for ant in antennas:
        if ant.stand.id <= 256:
            good.append(i)
            
        xyz[i,:] = ant.stand.xyz
        i += 1
        
    arrayX = xyz[good,0].mean()
    arrayY = xyz[good,1].mean()
    arrayZ = xyz[good,2].mean()
    arrayXYZ = xyz - np.array([arrayX, arrayY, arrayZ])
    
    # Load in the response of a single isolated stand
    standBeam = _load_stand_response(freq)
    
    # Build up the beam shape over all azimuths and altitudes
    beam_shape =  np.zeros((360,90))
    for az in list(range(360)):
        rAz = az*np.pi/180.0
        for el in list(range(90)):
            rAlt = el*np.pi/180.0
            beam_shape[az,el] = standBeam.response(aipy.coord.azalt2top(np.concatenate([[rAz], [rAlt]])))[0][0]
            
    # Build the output array and loop over all azimuths and altitudes
    output = np.zeros((360,90))
    with cf.ProcessPoolExecutor() as tpe:
        futures = {}
        for az in list(range(360)):
            for el in list(range(90)):
                task = tpe.submit(_int_beep_and_sweep,
                                  antennas, arrayXYZ, t, freq, az, el,
                                  beam_shape=beam_shape[az,el], sample_rate=sample_rate,
                                  direction=(azimuth, altitude))
                futures[task] = (az,el)
                
        i = 0
        for task in cf.as_completed(futures):
            # Display the progress meter if the `progress' keyword is set to True.  The
            # progress meter displays a `.' every 2% complete and the percentages every
            # 10%.  At 100%, `Done' is displayed.
            if progress:
                fracDone = i / 32400.0 * 100
                if fracDone % 10 == 0 and round(fracDone, 2) != 100:
                    sys.stdout.write("%i%%" % fracDone)
                elif round(fracDone, 2) == 100:
                    sys.stdout.write("Done\n")
                elif round(fracDone, 3) % 2 == 0:
                    sys.stdout.write(".")
                else:
                    pass
                sys.stdout.flush()
                
            az,el = futures[task]
            output[az,el] = task.result()
            i += 1
                
    # Done
    return output


def phase_and_sum(antennas, data, sample_rate=ndp_common.fS, central_freq=49.0e6, azimuth=0.0, altitude=90.0):
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
    delays = calc_delay(antennas, freq=central_freq, azimuth=azimuth, altitude=altitude)
    
    # Make the delays into something meaningful for the shifting of the data 
    # streams.  Then, get the beamforming coefficients (b^l_n (a la Steve's 
    # "Fun with TBN" memo).
    delays -= delays.min()
    bln = np.exp(2j*np.pi*central_freq*delays)
    
    # Figure out the polarizations
    pols = np.array([ant.pol for ant in antennas])
    stat = np.array([ant.combined_status for ant in antennas])
    pol0 = np.where( (pols == 0) & (stat == 33) )[0]
    pol1 = np.where( (pols == 1) & (stat == 33) )[0]
    
    # Loop over stands to compute the formed beam
    output = np.zeros((2, data.shape[1]), dtype=np.complex64)
    if len(pol0):
        output[0,:] = np.dot(bln[pol0], data[pol0,:]) / len(pol0)
    if len(pol1):
        output[1,:] = np.dot(bln[pol1], data[pol1,:]) / len(pol1)
        
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


def _phase_beep_and_sweep(antennas, arrayXYZ, t, freq, azimuth, altitude, beam_shape=1.0, sample_rate=ndp_common.fS, direction=(0.0, 90.0)):
    """
    Worker function for phase_beam_shape that 'beep's (makes a simulated signals) and
    'sweep's (phases it appropriately).
    """
    
    # Convert from degrees to radian
    rAz = azimuth*np.pi/180.0
    rAlt = altitude*np.pi/180.0
        
    # Unit vector for the current on-sky location
    currPos = np.array([np.cos(rAlt)*np.sin(rAz), 
                        np.cos(rAlt)*np.cos(rAz), 
                        np.sin(rAlt)])
    # Stand response in this direction
    currResponse = beam_shape
    
    # Loop over stands to build the simulated signals
    signals = currResponse * np.exp(2j*np.pi*freq*t).astype(np.complex64)
    signals.shape = (1,signals.size)
    signals = np.repeat(signals, len(antennas), axis=0)
    delays = np.zeros((len(antennas),1), dtype=np.complex64)
    for i in range(len(antennas)):
        currDelay = antennas[i].cable.delay(freq) - np.dot(currPos, arrayXYZ[i,:]) / speedOfLight
        delays[i,0] = np.exp(-2j*np.pi*freq*currDelay)
    signals *= delays
        
    # Beamform with delay-and-sum and store the RMS result
    beam = phase_and_sum(antennas, signals, sample_rate=sample_rate, central_freq=freq, 
                         azimuth=direction[0], altitude=direction[1])
                        
    # Return
    sigHere = (np.abs(beam)**2).mean()
    return sigHere


def phase_beam_shape(antennas, sample_rate=ndp_common.fS, central_freq=49.0e6, azimuth=0.0, altitude=90.0, progress=False):
    """
    Given a list of antennas, compute the on-sky response of the delay-and-sum
    scheme implemented in int_delay_and_sum.  A 360x90 numpy array spanning azimuth
    and altitude is returned.
    
    .. versionchanged:: 1.2.1
        Removed the 'disable_pool' keyword since recent optimztions to 
        the function have caused the multiprocessing.Pool feature to
        actually be slower than the single-threaded version.
        
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.stations` module instead of a list of stand ID
        numbers.
    """
    
    if isinstance(azimuth, AstroAngle):
        azimuth = azimuth.deg
    elif isinstance(azimuth, ephem.Angle):
        azimuth = azimuth * 180/np.pi
        
    if isinstance(altitude, AstroAngle):
        altitude = altitude.deg
    elif isinstance(altitude, ephem.Angle):
        altitude = altitude * 180/np.pi
        
    # Build up a base time array, load in the cable delays, and get the stand 
    # positions for geometric delay calculations.
    t = np.arange(0,1500)/sample_rate
    xyz = np.zeros((len(antennas),3))
    i = 0
    good = []
    for ant in antennas:
        if ant.stand.id <= 256:
            good.append(i)
            
        xyz[i,:] = ant.stand.xyz
        i += 1
        
    arrayX = xyz[good,0].mean()
    arrayY = xyz[good,1].mean()
    arrayZ = xyz[good,2].mean()
    arrayXYZ = xyz - np.array([arrayX, arrayY, arrayZ])
    
    # Load in the response of a single isolated stand
    standBeam = _load_stand_response(freq=central_freq)
    
    # Build up the beam shape over all azimuths and altitudes
    beam_shape =  np.zeros((360,90))
    for az in range(360):
        rAz = az*np.pi/180.0
        for el in range(90):
            rAlt = el*np.pi/180.0
            beam_shape[az,el] = standBeam.response(aipy.coord.azalt2top(np.concatenate([[rAz], [rAlt]])))[0][0]
            
    # Build the output array and loop over all azimuths and altitudes
    output = np.zeros((360,90))
    for az in range(360):
        for el in range(90):
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
                
            output[az,el] = _phase_beep_and_sweep(antennas, arrayXYZ, t, central_freq, az, el, beam_shape=beam_shape[az,el], sample_rate=sample_rate, direction=(azimuth, altitude))
            
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
    l = (x + lrSign*1.0j*y) / np.sqrt(2)
    r = (x - lrSign*1.0j*y) / np.sqrt(2)
    
    # Done
    return l, r
