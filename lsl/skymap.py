"""
Classes and methods to model sky brightness and visibility.

.. versionchanged:: 1.2.0
    Removed the orignal SkyMap class that uses LFmap at 73.9 MHz
    Updated SkyMapGSM to use the actual GSM rather than a fit to the GSM
    Added a new SkyMapLFSM that uses the 5.1 degree resolution LFSM
"""

### This module handles the skymaps from the 74 MHz skymap.
### David Munton, ARL:UT Jan 2007

import os
import abc
import numpy as np

from scipy.interpolate import interp1d

from astropy.coordinates import ICRS, EarthLocation, AltAz
from astropy import units as astrounits
from astropy.time import Time as AstroTime

from lsl import astro
from lsl.common.data_access import DataAccess

from lsl.misc import telemetry
telemetry.track_module()


__version__   = '0.5'
__all__ = ['SkyMapBase', 'SkyMapBase', 'SkyMapGSM', 'SkyMapLFSM', 'ProjectedSkyMap']
__author__    = 'J. York'
__maintainer__ = 'Jayce Dowell'


class SkyMapBase(object):
    """
    This code is the base class for the sky map. It takes as input a skymap file
    name and frequency to which the skymap corresponds.  It has the following
    methods:
     * _init_ - takes the array coordinate filename as an input argument.
     * normalize_power - Converts the skymap powers (in Kelvin radiated into 4 pi
                         ster) into a power seen at the antenna.
     * compute_total_power - Sums the power for all sources in the sky
     * ScaleSourcePowerstoFrequency - Scales the skymap from the base 73.8 MHz to
                                      the desired frequency.
    """
    
    __metaclass__ = abc.ABCMeta
    
    # Class data 
    degToRad = (np.pi/180.)             # Usual conversion factor
    
    def __init__(self, filename=None, freq_MHz=73.9):
        self.filename = filename
        self.freq_MHz = freq_MHz
        self._power = None
        
        # Load in the coordinates and data
        self._load()
        
    @abc.abstractmethod
    def _load(self):
        """
        Load in the specified filename and populate the ra, dec, and _power
        attributes.  This will be over-ridden in the skymap-specific subclasses.
        """
        
        raise NotImplementedError
        
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
            raise RuntimeError(f"{type(self).__name__} contains no data")
        return self.normalize_power().sum()


class SkyMapGSM(SkyMapBase):
    """
    Extension of the SkyMapBase class to use the Global Sky Model.
    
    For more information on the Global Sky Model, see: http://space.mit.edu/~angelica/gsm/index.html
    
    .. note:: This class uses a slightly different interpolation method than
              the original GSM and introduces a few percent difference at 
              74 MHz.
    
    .. versionchanged:: 1.2.0
        Reworked the GSM model to use the actual GSM that has been 
        downsampled to 64 sides rather than the fit.
    """
    
    _input = os.path.join('skymap', 'gsm-408locked.npz')
    
    def __init__(self, filename=None, freq_MHz=73.9):
        """
        Initialize the SkyMapGSM object with an optional full file path to 
        the skymap file.
        """
        
        if filename is None:
            filename = self._input
        SkyMapBase.__init__(self, filename=filename, freq_MHz=freq_MHz)
        
    def _load(self):
        # Since we are using a pre-computed GSM which is a NPZ file, read it
        # in with numpy.load.
        with DataAccess.open(self.filename, 'rb') as fh:
            dataDict = np.load(fh)
            
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
            scaleFunc = interp1d(np.log(freqs), np.log(sigmas), kind='cubic')
            compFuncs = []
            for i in range(comps.shape[1]):
                compFuncs.append( interp1d(np.log(freqs), comps[:,i], kind='cubic') )
            ## Actually create the realization by running the interplation and adding up the
            ## compnent maps
            output = maps[:,0]*0.0
            for i,compFunc in enumerate(compFuncs):
                output += compFunc(np.log(self.freq_MHz))*maps[:,i]
            output *= np.exp(scaleFunc(np.log(self.freq_MHz)))
            ## Save
            self._power = output


class SkyMapLFSM(SkyMapGSM):
    """
    Extension of the SkyMapBase class to use the Low Frequency Sky Model with 5.1
    degree resolution.
    
    For more information on the Low Frequency Sky Model, see: https://lda10g.alliance.unm.edu/LWA1LowFrequencySkySurvey/
    
    .. versionadded:: 1.2.0
    """
    
    _input = os.path.join('skymap', 'lfsm-5.1deg.npz')
    
    def __init__(self, filename=None, freq_MHz=73.9):
        """
        Initialize the SkyMapLFSM object with an optional full file path to 
        the skymap file.
        """
        
        if filename is None:
            filename = self._input
        SkyMapBase.__init__(self, filename=filename, freq_MHz=freq_MHz)
        
    def _load(self):
        # Since we are using a pre-computed GSM which is a NPZ file, read it
        # in with numpy.load.
        with DataAccess.open(self.filename, 'rb') as fh:
            dataDict = np.load(fh)
            
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
            scaleFunc = interp1d(np.log(freqs), np.log(sigmas), kind='slinear')
            compFuncs = []
            for i in range(comps.shape[1]):
                compFuncs.append( interp1d(np.log(freqs), comps[:,i], kind='cubic') )
            ## Actually create the realization by running the interplation and adding up the
            ## compnent maps
            output = maps[:,0]*0.0
            for i,compFunc in enumerate(compFuncs):
                output += compFunc(np.log(self.freq_MHz))*maps[:,i]
            output *= np.exp(scaleFunc(np.log(self.freq_MHz)))
            ## Save
            self._power = output


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
    
    def __init__(self, skymap_object, lat, lon, utc_jd):
        """
        Initialize the skymap at input lat,lon (decimal degrees) and time (in 
        UTC julian day).
        """
        
        self.skymap_object = skymap_object
        
        assert  -90 <= lat <=  90, ValueError(f"lat = {lat:g} not in [-90,90]")
        assert -360 <= lon <= 360, ValueError(f"lon = {lon:g} not in [-360,360]")
        self.lat  = lat
        self.lon  = lon
        self.time = utc_jd
        
        # Compute the ra and dec locations to a altitude and azimuth. This requires a lat/lon and a time (UTC).
        # Alt and Az are compressed to show only the visible source.
        sc = ICRS(ra=self.skymap_object.ra*astrounits.deg, dec=self.skymap_object.dec*astrounits.deg)
        el = EarthLocation.from_geodetic(lon*astrounits.deg, lat*astrounits.deg, height=0*astrounits.m)
        aa = AltAz(location=el, obstime=AstroTime(utc_jd, format='jd', scale='utc'))
        tc = sc.transform_to(aa)
        self.alt = tc.alt.to('deg').value
        self.az = tc.az.to('deg').value
        
        # Compress the az/alt so that only the visible sources are available. "
        visibleMask = self.alt > 0.0
        
        #Compress arrays to hide non-visible sources
        self.visibleAlt = np.compress(visibleMask, self.alt)
        self.visibleAz = np.compress(visibleMask, self.az)
        self.visiblePower = np.compress(visibleMask, self.skymap_object._power)
        try:
            pixelSolidAngle = (2.0 * np.pi)**2/(self.skymap_object.numPixelsX*self.skymap_object.numPixelsY)
        except AttributeError:
            pixelSolidAngle = 2.5566346e-4
        fractionSolidAngle = (1.0/4.0*np.pi) * pixelSolidAngle
        
        # The cosine term is the projection of the receiving area onto the direction 
        # of the source
        normalizedPower = self.skymap_object.normalize_power() * fractionSolidAngle  
        self.visibleNormalizedPower = np.compress(visibleMask, normalizedPower)
        #self.visibleNormalizedPower = self.visiblePower*np.cos(self.visibleAlt * self.skymap_object.degToRad)
        self.visibleRa = np.compress(visibleMask, self.skymap_object.ra)
        self.visibleDec = np.compress(visibleMask, self.skymap_object.dec)
        
    def get_direction_cosines(self):
        """
        Compute the direction cosines and return the tuple of arrays (l,m,n).
        """
        
        altRad = self.visibleAlt*self.skymap_object.degToRad
        azRad = self.visibleAz*self.skymap_object.degToRad
        l = np.cos(altRad)*np.sin(azRad)
        m = np.cos(altRad)*np.cos(azRad)
        n = np.sin(altRad)
        return (l, m, n)
        
    def compute_visible_power(self):
        """
        Compute and return the the total power from visible portion of the sky.
        """
        
        if len(self.visibleNormalizedPower) == 0:
            raise RuntimeError("%s contains no data" % type(self).__name__)
        totalVisiblePower = sum(self.visibleNormalizedPower)
        return totalVisiblePower
