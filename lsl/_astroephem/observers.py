from __future__ import print_function, division

import numpy
from functools import wraps
from scipy.optimize import minimize_scalar

from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, ICRS
import astropy.units as u

from angles import hours, degrees
from dates import Date


__all__ = ['CircumpolarError', 'NeverUpError', 'AlwaysUpError', 'Observer']


class CircumpolarError(ValueError):
    pass


class NeverUpError(CircumpolarError):
    pass


class AlwaysUpError(CircumpolarError):
    pass


def _location(djd, obs, bdy, value, rising):
    obs.date = Date(djd)
    bdy.compute(obs)
    diff = bdy.alt.rad
    if bdy.az.rad % (2*numpy.pi) <= numpy.pi:
        if not rising:
            diff += numpy.pi/2
    else:
        if rising:
            diff += numpy.pi/2
    diff = abs(diff - value) % (2*numpy.pi)
    return diff


def protect_date(func):
    @wraps(func)
    def wrapper(*args, **kwds):
        initial_date = args[0].date
        output = func(*args, **kwds)
        args[0].date = initial_date
        return output
    return wrapper


class Observer(object):
    """
    A location on earth for which positions are to be computed.
    
    An `Observer` instance allows you to compute the positions of
    celestial bodies as seen from a particular latitude and longitude on
    the Earth's surface.  The constructor takes no parameters; instead,
    set its attributes once you have created it.  Defaults:
    
    `date` - the moment the `Observer` is created
    `lat` - zero degrees latitude
    `lon` - zero degrees longitude
    `elevation` - 0 meters above sea level
    """
    
    def __init__(self):
        self.__lat = 0.0*u.deg
        self.__lon = 0.0*u.deg
        self.__elev = 0.0*u.m
        self.__date = Time.now()
        self._update()
        
    def _update(self):
        self._el = EarthLocation.from_geodetic(self.__lon, self.__lat, self.__elev)
        
    @property
    def lat(self):
        return degrees(self.__lat)
    @lat.setter
    def lat(self, value):
        self.__lat = degrees(value)
        self._update()
        
    @property
    def lon(self):
        return degrees(self.__lon)
    @lon.setter
    def lon(self, value):
        self.__lon = degrees(value)
        self._update()
    @property
    def long(self):
        return self.lon
    @long.setter
    def long(self, value):
        self.lon = value
        
    @property
    def elev(self):
        return self.__elev.to('m').value
    @elev.setter
    def elev(self, value):
        self.__elev = float(value)*u.m
        self._update()
    @property
    def elevation(self):
        return self.elev
    @elevation.setter
    def elevation(self, value):
        self.elev = value
        
    @property
    def date(self):
        return self.__date
    @date.setter
    def date(self, value):
        self.__date = Date(value)
        
    def as_astropy(self):
        return self._el, self.__date
        
    def sidereal_time(self):
        return self.__date.sidereal_time('apparent', longitude=self.__lon)
        
    def radec_of(self, az, alt):
        alt = degrees(alt)
        az = degrees(az)
        topo = AltAz(alt, az, obstime=self.__date, location=self._el)
        equ = topo.transform_to(ICRS())
        return hours(equ.ra.to('radian').value), degrees(equ.dec.to('radian').value)
        
    @protect_date
    def previous_transit(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, numpy.pi/2, True),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def next_transit(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, numpy.pi/2, True),
                              method='bounded',
                              bounds=(start, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def previous_antitransit(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, -numpy.pi/2, True),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def next_antitransit(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, -numpy.pi/2, True),
                              method='bounded',
                              bounds=(start, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def previous_rising(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, 0.0, True),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
        
    @protect_date
    def next_rising(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, 0.0, True),
                              method='bounded',
                              bounds=(start, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
        
    @protect_date
    def previous_setting(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, 0.0, False),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start),
                              options={'xatol': 1/86400.0})
        print(sol)
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
        
    @protect_date
    def next_setting(self, body, start=None):
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, 0.0, False),
                              method='bounded',
                              bounds=(start, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
