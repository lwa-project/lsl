"""
Python module to provide a ephem.Observer like class.
"""

from __future__ import print_function, division

import numpy
from functools import wraps
from scipy.optimize import minimize_scalar

from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, ICRS
import astropy.units as u

from lsl._astroephem.angles import hours, degrees
from lsl._astroephem.dates import Date


__all__ = ['CircumpolarError', 'NeverUpError', 'AlwaysUpError', 'Observer']


class CircumpolarError(ValueError):
    """
    Error class for when a body is circumpolar.
    """
    
    pass


class NeverUpError(CircumpolarError):
    """
    Error class for when a body never rises above the horizon for the given
    observer.
    """
    
    pass


class AlwaysUpError(CircumpolarError):
    """
    Error class for when a body never sets below the horizon for the given
    observer.
    """
    
    pass


def protect_date(func):
    """
    Wrapper to make sure that methods that search in time do not change the
    underlaying date of an Observer instance.
    """
    
    @wraps(func)
    def wrapper(*args, **kwds):
        initial_date = args[0].date
        output = func(*args, **kwds)
        args[0].date = initial_date
        return output
    return wrapper


def _location(djd, obs, bdy, value, selection):
    """
    Private function for helping to determine when a body is at a certain
    elevation.  selection controls whether the minimimzation of azimuth agnostic
    (0), in the east (1), or in the west (2).
    """
    
    obs.date = Date(djd)
    bdy.compute(obs)
    diff = bdy.alt.rad
    if bdy.az.rad % (2*numpy.pi) <= numpy.pi:
        if selection == 2:
            diff += numpy.pi/2
    else:
        if selection == 1:
            diff += numpy.pi/2
    diff = abs(diff - value) % (2*numpy.pi)
    return diff


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
    `horizon` - zero degrees
    """
    
    def __init__(self):
        self.__lat = 0.0*u.deg
        self.__lon = 0.0*u.deg
        self.__elev = 0.0*u.m
        self.__horz = 0.0*u.deg
        self.__date = Time.now()
        self._update()
        
    def _update(self):
        self._el = EarthLocation.from_geodetic(self.__lon, self.__lat, self.__elev)
        
    @property
    def lat(self):
        """
        The geodetic latitude of the Observer.  Positive is North.
        """
        
        return degrees(self.__lat)
    @lat.setter
    def lat(self, value):
        self.__lat = degrees(value)
        self._update()
        
    @property
    def lon(self):
        """
        The geodetic longtiude of the Observer.  Positive is East.
        """
        
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
        """
        The elevation of the Observer above the reference ellipsoid in m.
        """
        
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
    def horizon(self):
        """
        The effective horizon for the Observer.  Positive is above the true
        horizon.
        """
        
        return degrees(self.__horz)
    @horizon.setter
    def horizon(self, value):
        self.__horz = degrees(value)
        self._update()
        
    @property
    def date(self):
        """
        The date set for the Observer.
        """
        
        return self.__date
    @date.setter
    def date(self, value):
        self.__date = Date(value)
        
    def as_astropy(self):
        """
        Returns an astropy.coordinates.EarthLocation instance for the Observer.
        """
        
        return self._el, self.__date
        
    def sidereal_time(self):
        """
        Returns an Angle instance of the local apparent sidereal time for
        the Observer.
        """
        
        return hours(self.__date.sidereal_time('apparent', longitude=self.__lon))
        
    def radec_of(self, az, alt):
        """
        Given an azimuth and elevation as viewed by the Observer, return the
        right ascension and declination that they correspond to.  The RA and
        dec. values are returned as Angle instances in the ICRS frame.
        """
        
        alt = degrees(alt)
        az = degrees(az)
        topo = AltAz(az, alt, obstime=self.__date, location=self._el)
        equ = topo.transform_to(ICRS())
        return hours(equ.ra.to('radian').value), degrees(equ.dec.to('radian').value)
        
    @protect_date
    def previous_transit(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the previous transit
        time for the object.  If start is None, the time is relative to the
        Observer's current time.  Otherwise, it is relative to the time in
        start.  Returns an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, numpy.pi/2, 0),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start*1.0),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def next_transit(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the next transit time
        for the object.  If start is None, the time is relative to the
        Observer's current time.  Otherwise, it is relative to the time in
        start.  Returns an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, numpy.pi/2, 0),
                              method='bounded',
                              bounds=(start*1.0, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def previous_antitransit(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the previous anti-
        transit time for the object.  If start is None, the time is relative to
        the Observer's current time.  Otherwise, it is relative to the time in
        start.  Returns an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, -numpy.pi/2, 0),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start*1.0),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def next_antitransit(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the next anti-transit
        time for the object.  If start is None, the time is relative to the
        Observers current time.  Otherwise, it is relative to the time in start.
        Returns an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, -numpy.pi/2, 0),
                              method='bounded',
                              bounds=(start*1.0, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        return Date(sol.x)
        
    @protect_date
    def previous_rising(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the previous rise time
        for the object.  If start is None, the time is relative to the
        Observer's current time.  Otherwise, it is relative to the time in
        start.  Returns an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, self.__horz.to(u.rad).value, 1),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start*1.0),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
        
    @protect_date
    def next_rising(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the next rise time for
        the object.  If start is None, the time is relative to the Observer's
        current time.  Otherwise, it is relative to the time in start.  Returns
        an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, self.__horz.to(u.rad).value, 1),
                              method='bounded',
                              bounds=(start*1.0, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
        
    @protect_date
    def previous_setting(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the previous set time
        for the object.  If start is None, the time is relative to the
        Observer's current time.  Otherwise, it is relative to the time in
        start.  Returns an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, self.__horz.to(u.rad).value, 2),
                              method='bounded',
                              bounds=(start-u.sday.to(u.day), start*1.0),
                              options={'xatol': 1/86400.0})
        print(sol)
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
        
    @protect_date
    def next_setting(self, body, start=None):
        """
        Given a FixedBody or Planet instance, determine the next set time for
        the object.  If start is None, the time is relative to the Observer's
        current time.  Otherwise, it is relative to the time in start.  Returns
        an EphemTime instance.
        """
        
        if start is None:
            start = self.date
            
        sol = minimize_scalar(_location, args=(self, body, self.__horz.to(u.rad).value, 2),
                              method='bounded',
                              bounds=(start*1.0, start+u.sday.to(u.day)),
                              options={'xatol': 1/86400.0})
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise CircumpolarError()
        return Date(sol.x)
