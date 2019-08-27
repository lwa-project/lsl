# -*- coding: utf-8 -*-

"""
Classes and methods to model sky brightness and visibility.

.. versionchanged:: 1.2.0
    Removed the orignal SkyMap class that uses LFmap at 73.9 MHz
    Updated SkyMapGSM to use the actual GSM rather than a fit to the GSM
    Added a new SkyMapLFSM that uses the 5.1 degree resolution LFSM
"""

### This module handles the skymaps from the 74 MHz skymap.
### David Munton, ARL:UT Jan 2007

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
from numpy import pi, float32, log, exp, log10, sin, cos, arcsin, arccos, empty, arange, compress, clip, where, load

from scipy.interpolate import interp1d

from lsl import astro
from lsl.common.paths import DATA as dataPath

__version__   = '0.3'
__revision__ = '$Rev$'
__all__ = ['SkyMapGSM', 'SkyMapLFSM', 'ProjectedSkyMap']
__author__    = 'J. York'
__maintainer__ = 'Jayce Dowell'


### This code is the base class for the sky map. It takes as input a skymap file name and frequency to
### which the skymap corresponds.  It has the following methods:
###     _init_ - takes the array coordinate filename as an input argument.
###     normalize_power - Converts the skymap powers (in Kelvin radiated into 4 pi ster) into a power 
###                      seen at the antenna.
###     compute_total_power - Sums the power for all sources in the sky
###     ScaleSourcePowerstoFrequency - Scales the skymap from the base 73.8 MHz to the desired frequency.


class SkyMapGSM(object):
    """
    Extension of the SkyMap class to use the Global Sky Model.
    
    For more information on the Global Sky Model, see: http://space.mit.edu/~angelica/gsm/index.html
    
    .. note:: This class uses a slightly different interpolation method than
              the original GSM and introduces a few percent difference at 
              74 MHz.
    
    .. versionchanged:: 1.2.0
        Reworked the GSM model to use the actual GSM that has been 
        downsampled to 64 sides rather than the fit.
    """
    
    # Class data 
    degToRad = (pi/180.)             # Usual conversion factor
    
    _input = os.path.join(dataPath, 'skymap', 'gsm-408locked.npz')
    
    def __init__(self, filename=None, freq_MHz=73.9):
        """
        Initialize the SkyMapGSM object with an optional full file path to 
        the skymap file.
        """
        
        if filename is None:
            filename = self._input
            
        # Since we are using a pre-computed GSM which is a NPZ file, read it
        # in with numpy.load.
        dataDict = load(filename)
        
        # RA and dec. are stored in the dictionary as radians
        self.ra = dataDict['ra'].ravel() / self.degToRad
        self.dec = dataDict['dec'].ravel() / self.degToRad
        
        # Compute the temperature for the current frequency
        ## Load in the data
        freqs = dataDict['freqs']
        sigmas = dataDict['sigmas']
        comps = dataDict['comps']
        maps = dataDict['maps']
        ## Build the scale and spectral component interpolation functions so that we can
        ## move from the surveys to an arbitrary frequency.  This is done using cubic 
        ## interpolation across log(freq)
        scaleFunc = interp1d(log(freqs), log(sigmas), kind='cubic')
        compFuncs = []
        for i in xrange(comps.shape[1]):
            compFuncs.append( interp1d(log(freqs), comps[:,i], kind='cubic') )
        ## Actually create the realization by running the interplation and adding up the
        ## compnent maps
        output = maps[:,0]*0.0
        for i,compFunc in enumerate(compFuncs):
            output += compFunc(log(freq_MHz))*maps[:,i]
        output *= exp(scaleFunc(log(freq_MHz)))
        ## Save
        self._power = output
        
        # Close out the dictionary
        try:
            dataDict.close()
        except AttributeError:
            pass
            
    def normalize_power(self):
        """
        Compute the skymap power (total power radiated into 4 pi steradians) into 
        a power at antenna, based on pixel count.
        """
        
        return self._power
        
    def compute_total_power(self):
        """
        Compute and return the the total power from the sky.
        """
        
        if len(self._power) == 0:
            raise RuntimeError("%s contains no data" % type(self).__name__)
        return self.normalize_power().sum()


class SkyMapLFSM(SkyMapGSM):
    """
    Extension of the SkyMap class to use the Low Frequency Sky Model with 5.1
    degree resolution.
    
    For more information on the Low Frequency Sky Model, see: https://lda10g.alliance.unm.edu/LWA1LowFrequencySkySurvey/
    
    .. versionadded:: 1.2.0
    """
    
    _input = os.path.join(dataPath, 'skymap', 'lfsm-5.1deg.npz')
    
    def __init__(self, filename=None, freq_MHz=73.9):
        """
        Initialize the SkyMapLFSM object with an optional full file path to 
        the skymap file.
        """
        
        if filename is None:
            filename = self._input
            
        # Since we are using a pre-computed GSM which is a NPZ file, read it
        # in with numpy.load.
        dataDict = load(filename)
        
        # RA and dec. are stored in the dictionary as radians
        self.ra = dataDict['ra'].ravel() / self.degToRad
        self.dec = dataDict['dec'].ravel() / self.degToRad
        
        # Compute the temperature for the current frequency
        ## Load in the data
        freqs = dataDict['freqs']
        sigmas = dataDict['sigmas']
        comps = dataDict['comps']
        maps = dataDict['maps']
        ## Build the scale and spectral component interpolation functions so that we can
        ## move from the surveys to an arbitrary frequency.  This is done using spline
        ## interpolation for the scale factor and cubic for the 2-D structure.
        ## interpolation across log(freq)
        scaleFunc = interp1d(log(freqs), log(sigmas), kind='slinear')
        compFuncs = []
        for i in xrange(comps.shape[1]):
            compFuncs.append( interp1d(log(freqs), comps[:,i], kind='cubic') )
        ## Actually create the realization by running the interplation and adding up the
        ## compnent maps
        output = maps[:,0]*0.0
        for i,compFunc in enumerate(compFuncs):
            output += compFunc(log(freq_MHz))*maps[:,i]
        output *= exp(scaleFunc(log(freq_MHz)))
        ## Save
        self._power = output
        
        # Close out the dictionary
        try:
            dataDict.close()
        except AttributeError:
            pass


class ProjectedSkyMap(object):
    """
    The class for handling the model sky brightness maps over a particular site.
    This code is the base class for the sky map visible at a specific location. It 
    takes as input a skymap file name and frequency to which the skymap corresponds.
    It inherits from class SkyMap. It has the following methods:
    1. _init_ - takes the array coordinate filename as an input argument.
    2. get_direction_cosines - Computes the direction cosines 
    3. compute_visibile_power - Sums the power for all visible sources in 
    the sky.
    """
    
    def __init__(self, skyMapObject, lat, lon, utc_jd):
        """
        Initialize the skymap at input lat,lon (decimal degrees) and time (in 
        UTC julian day).
        """
        
        self.skyMapObject = skyMapObject
        
        assert  -90 <= lat <=  90, ValueError('lat = %g not in [-90,90]' % lat)
        assert -360 <= lon <= 360, ValueError('lon = %g not in [-360,360]' % lon)
        self.lat  = lat
        self.lon  = lon
        self.time = utc_jd
        
        # Compute the ra and dec locations to a altitude and azimuth. This requires a lat/lon and a time (UTC).
        # Alt and Az are compressed to show only the visible source.
        ## Extract the RA/dec values
        ra = self.skyMapObject.ra
        dec = self.skyMapObject.dec
        
        ## Compute the LST at the lat/lon (hours)
        lst = astro.get_local_sidereal_time(self.lon, self.time)
        
        ## RA -> HA (deg)
        ha = lst*15.0 - ra
        
        ## HA, dec, and lat in degrees to alt and az in radians
        sinAlt = sin(dec*pi/180)*sin(self.lat*pi/180) + cos(dec*pi/180)*cos(lat*pi/180)*cos(ha*pi/180)
        sinAlt = clip(sinAlt, -1.0, 1.0)
        alt = arcsin(sinAlt)
        cosAz = (sin(dec*pi/180)-sinAlt*sin(self.lat*pi/180))/(cos(alt)*cos(lat*pi/180))
        cosAz = clip(cosAz, -1.0, 1.0)
        az = arccos(cosAz)
        swap = where(sin(ha*pi/180) > 0)
        az[swap] = 2*pi-az[swap]
        
        ## Convert to alt and az to degrees
        self.alt = alt*180/pi
        self.az  = az*180/pi
        
        # Compress the az/alt so that only the visible sources are available. "
        visibleMask = self.alt > 0.0
        
        #Compress arrays to hide non-visible sources
        self.visibleAlt = compress(visibleMask, self.alt)
        self.visibleAz = compress(visibleMask, self.az)
        self.visiblePower = compress(visibleMask, self.skyMapObject._power)
        try:
            pixelSolidAngle = (2.0 * pi)**2/(self.skyMapObject.numPixelsX*self.skyMapObject.numPixelsY)
        except AttributeError:
            pixelSolidAngle = 2.5566346e-4
        fractionSolidAngle = (1.0/4.0*pi) * pixelSolidAngle
        
        # The cosine term is the projection of the receiving area onto the direction 
        # of the source
        normalizedPower = self.skyMapObject.normalize_power() * fractionSolidAngle  
        self.visibleNormalizedPower = compress(visibleMask, normalizedPower)
        #self.visibleNormalizedPower = self.visiblePower*cos(self.visibleAlt * self.skyMapObject.degToRad)
        self.visibleRa = compress(visibleMask, self.skyMapObject.ra)
        self.visibleDec = compress(visibleMask, self.skyMapObject.dec)
        
    def get_direction_cosines(self):
        """
        Compute the direction cosines and return the tuple of arrays (l,m,n).
        """
        
        altRad = self.visibleAlt*self.skyMapObject.degToRad
        azRad = self.visibleAz*self.skyMapObject.degToRad
        l = cos(altRad)*sin(azRad)
        m = cos(altRad)*cos(azRad)
        n = sin(altRad)
        return (l, m, n)
        
    def compute_visible_power(self):
        """
        Compute and return the the total power from visible portion of the sky.
        """
        
        if len(self.visibleNormalizedPower) == 0:
            raise RuntimeError("%s contains no data" % type(self).__name__)
        totalVisiblePower = sum(self.visibleNormalizedPower)
        return totalVisiblePower
