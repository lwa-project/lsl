"""
Python module for providing ephem-like FixedBody and planetary classes.
"""

from __future__ import print_function, division

import numpy
from functools import wraps

from skyfield import api
from skyfield.starlib import Star
from skyfield.toposlib import Topos
from skyfield.timelib import Time as SkyTime
from skyfield.units import Angle as SkyAngle, Distance as SkyDistance
ts = api.load.timescale()

from lsl._skyephem.angles import hours, degrees
from lsl._skyephem.dates import Date, J2000
from lsl._skyephem.observers import Observer
from lsl._skyephem.cache import load_planetary_ephemeris


__all__ = ['prepare_date_or_observer', 'Body', 'FixedBody', 'readdb',
           'Sun', 'Mercury', 'Venus', 'Moon', 'Mars', 'Jupiter', 'Saturn',
           'Uranus', 'Neptune']
           
           
_solar_system = load_planetary_ephemeris()



def prepare_date_or_observer(func):
    """
    Wrapper to make sure that methods that accept a Date, Observer, or None
    do the right thing.
    """
    
    @wraps(func)
    def wrapper(cls, date_or_observer=None):
        if date_or_observer is None:
            date_or_observer = ts.now()
        elif isinstance(date_or_observer, str):
            date_or_observer = Date(date_or_observer)
        elif isinstance(date_or_observer, float):
            date_or_observer = Date(date_or_observer)
        output = func(cls, date_or_observer=date_or_observer)
        return output
    return wrapper


class Body(object):
    """
    Base class for a celestial body that provides the compute() method.  This is
    sub-classed to generate both FixedBody and Planet.
    """
    
    name = ''
    
    def __init__(self, body=None):
        if body is None:
            ra = SkyAngle(degrees=0.0, preference='hours')
            dec = SkyAngle(degrees=0.0, preference='degrees')
            body = Star(ra=ra, dec=dec, epoch=J2000)
        self._body = body
        
    def __repr__(self):
        return "<%s.%s %s at 0x%x>" % (type(self).__module__, type(self).__name__, '"%s"' % self.name if self.name else 'None', id(self))
        
    @prepare_date_or_observer
    def compute(self, date_or_observer=None):
        """
        Given a Date or Observer instance or None, compute the location of
        the body.
        
        For a Date, this computes:
         * the astrometric position of the body - a_ra and a_dec
         
        For an Observer, this computes everything as in the case of a Date
        and also includes:
         * the apparent position of the body - ra and dec
         * the topocentric position of the body - az and alt
         
        If None is provided, the current time is used.
        """
        
        if isinstance(date_or_observer, SkyTime):
            t = date_or_observer
            obs = _solar_system['earth']
        else:
            t = date_or_observer.date
            obs = _solar_system['earth'] + date_or_observer._wgs84
        pos = obs.at(t).observe(self._body)
        
        ra, dec, _ = pos.radec()
        self.a_ra = hours(ra.radians, wrap=True)
        self.a_dec = degrees(dec.radians, wrap180=True)
        
        try:
            geo = pos.apparent()
            ra, dec, _ = geo.radec('date')
        except AttributeError as e:
            pass
        self.g_ra = hours(ra.radians, wrap=True)
        self.g_dec = degrees(dec.radians, wrap180=True)
        
        try:
            ra, dec, _ = pos.apparent().radec('date')
            self.ra = hours(ra.radians, wrap=True)
            self.dec = degrees(dec.radians, wrap180=True)
        except AttributeError as e:
            pass
            
        try:
            alt, az, _ = pos.apparent().altaz()
            #print(alt, az)
            self.az = degrees(az, wrap360=True)
            self.alt = degrees(alt, wrap180=True)
        except (AttributeError, ValueError) as e:
            pass
            
        self._rise_transit_set(date_or_observer)
        
    def _rise_transit_set(self, observer):
        """
        From:
            https://github.com/brandon-rhodes/pyephem/blob/master/libastro-3.7.7/riset.c
        """
        
        if not isinstance(observer, Observer):
            return False
            
        try:
            self.ra
        except AttributeError:
            return False
            
        # Load the values we need in radians
        ra = self.ra.radians
        dec = self.dec.radians
        lat = observer.lat.radians
        
        # Flip things around if we are south of the equator
        southern = (lat < 0)
        if southern:
            lat = -lat
            dec = -dec
            
        # Figure out how high it gets, if it is circumpolar, and if it never rises
        z = numpy.pi/2
        zmin = abs(dec - lat)
        zmax = numpy.pi - abs(dec + lat)
        if zmax <= z + 1e-9:
            circumpolar = True
        else:
            circumpolar = False
        if zmin >= z - 1e-9:
            neverup = True
        else:
            neverup = False
        transit_alt = zmin
            
        # Find the hour angle that it rises at
        cos_h = (numpy.cos(z)-numpy.sin(lat)*numpy.sin(dec))/(numpy.cos(lat)*numpy.cos(dec));
        if cos_h >= 1:
            rising_ha =  0.0
        elif cos_h <= -1:
            rising_ha = numpy.pi
        else:
            rising_ha = numpy.arccos(cos_h)
            
        # Find the setting azimuth
        xaz = numpy.sin(dec)*numpy.cos(lat)-numpy.cos(dec)*numpy.cos(rising_ha)*numpy.sin(lat);
        yaz = -numpy.cos(dec)*numpy.sin(rising_ha);
        if xaz == 0.:
	        if yaz > 0:
	            setting_az = numpy.pi/2
	        else:
	            setting_az = -numpy.pi/2
        else:
	        setting_az = numpy.arctan2(yaz, xaz)
        if southern:
	        setting_az = numpy.pi - setting_az
        rising_az = 2*numpy.pi - setting_az
	    
        # Find the rising and setting local sidereal times
        rising_lst = ra - rising_ha
        setting_lst = ra + rising_ha

        # Save it all
        ## Direct access
        self.circumpolar = circumpolar
        self.neverup = neverup
        self.transit_alt = degrees(transit_alt)
        self.rising_az = degrees(rising_az, wrap360=True)
        self.setting_az = degrees(setting_az, wrap360=True)
        ## Used later
        self._rising_lst = hours(rising_lst, wrap=True)
        self._setting_lst = hours(setting_lst, wrap=True)


class FixedBody(Body):
    """
    A celestial body in the ICRS frame, that can compute() its sky position.
    """
    
    name = ''
    __ra = SkyAngle(degrees=0.0, preference='hours')
    __dec = SkyAngle(degrees=0.0, preference='degrees')
    __epoch = J2000
    __pmra = 0.0    # mas/yr
    __pmdec = 0.0   # mas/yr
    
    def _parameter_callback(self):
        self._body = Star(ra=self.__ra, dec=self.__dec,
                          ra_mas_per_year=self.__pmra, dec_mas_per_year=self.__pmdec,
                          names=(self.name,), epoch=self.__epoch)
        
    @property
    def _ra(self):
        return self.__ra
    @_ra.setter
    def _ra(self, value):
        self.__ra = hours(value)
        self._parameter_callback()
        
    @property
    def _dec(self):
        return self.__dec
    @_dec.setter
    def _dec(self, value):
        self.__dec = degrees(value)
        self._parameter_callback()
        
    @property
    def _epoch(self):
        return self.__epoch
    @_epoch.setter
    def _epoch(self, value):
        self.__epoch = Date(value)
        self._parameter_callback()
        
    @property
    def _pmra(self):
        return self.__pmra
    @_pmra.setter
    def _pmra(self, value):
        self.__pmra = value
        self._parameter_callback()
    @property
    def _pmdec(self):
        return self.__pmdec
    @_pmdec.setter
    def _pmdec(self, value):
        self.__pmdec = value
        self._parameter_callback() 


def readdb(line):
    """
    Read in a line of text from an ephem database and return a FixedBody
    corresponding to that entry.
    """
    
    fields = line.split(',')
    name, type, ra, dec, _ = fields
    if type != 'f|J':
        raise ValueError("Can only deal with database entries that are 'f|J'")
    bdy = FixedBody()
    bdy.name = name
    bdy._ra = ra
    bdy._dec = dec
    return bdy


class Planet(Body):
    """
    A solar system body, that can compute() its sky position.
    """
    
    def __init__(self, name):
        Body.__init__(self)
        self.name = name
        try:
            self._body = _solar_system[name.lower()]
        except KeyError:
            self._body = _solar_system[name.lower()+" barycenter"]


class Sun(Planet):
    """
    A Planet instance representing the Sun.
    """
    
    def __init__(self):
        Planet.__init__(self, 'sun')


class Mercury(Planet):
    """
    A Planet instance representing Mercury.
    """
    
    def __init__(self):
        Planet.__init__(self, 'mercury')


class Venus(Planet):
    """
    A Planet instance representing Venus.
    """
    
    def __init__(self):
        Planet.__init__(self, 'venus')


class Moon(Planet):
    """
    A Planet instance representing Earth's moon
    """
    
    def __init__(self):
        Planet.__init__(self, 'moon')


class Mars(Planet):
    """
    A Planet instance representing Mars.
    """
    
    def __init__(self):
        Planet.__init__(self, 'mars')


class Jupiter(Planet):
    """
    A Planet instance representing Jupiter.
    """
    
    def __init__(self):
        Planet.__init__(self, 'jupiter')


class Saturn(Planet):
    """
    A Planet instance representing Saturn.
    """
    
    def __init__(self):
        Planet.__init__(self, 'saturn')


class Uranus(Planet):
    """
    A Planet instance representing Uranus.
    """
    
    def __init__(self):
        Planet.__init__(self, 'uranus')


class Neptune(Planet):
    """
    A Planet instance representing Neptune.
    """
    
    def __init__(self):
        Planet.__init__(self, 'neptune')
