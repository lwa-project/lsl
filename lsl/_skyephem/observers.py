"""
Python module to provide a ephem.Observer like class.
"""

from __future__ import print_function, division

import numpy
from functools import wraps
from scipy.optimize import minimize_scalar

from skyfield import api
from skyfield import almanac
from skyfield.toposlib import Topos
from skyfield.searchlib import find_discrete
from skyfield.units import Angle as SkyAngle, Distance as SkyDistance
ts = api.load.timescale()

from lsl._skyephem.angles import hours, degrees
from lsl._skyephem.dates import Date
from lsl._skyephem.cache import load_planetary_ephemeris


__all__ = ['CircumpolarError', 'NeverUpError', 'AlwaysUpError', 'Observer']



_solar_system = load_planetary_ephemeris()
_ter = _solar_system['earth']

# Sidereal day as a fraction of a day
_SIDEREAL_DAY = 0.9972695663194444


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
    diff = bdy.alt.radians
    if bdy.az.radians % (2*numpy.pi) <= numpy.pi:
        if selection == 2:
            diff += numpy.pi/2
    else:
        if selection == 1:
            diff += numpy.pi/2
    diff = abs(diff - value) % (2*numpy.pi)
    return diff


def _search(obs, bdy, value, selection, t0, t1, tol=0.001/86400.0, recursive=False):
    # Coarse grid to zero in on where we should look
    ts = numpy.linspace(t0-0.5/24, t1+0.5/24, 51)
    vs = [_location(t, obs, bdy, value, selection) for t in ts]
    best = numpy.argmin(vs)
    if best <= 1:
        best = -3
    if best >= ts.size - 2:
        best = 2
    t0_prime = ts[best] - (ts[1] - ts[0])
    t1_prime = ts[best] + (ts[1] - ts[0])
    
    # Fine search
    sol = minimize_scalar(_location, args=(obs, bdy, value, selection),
                          method='bounded',
                          bounds=(t0_prime, t1_prime),
                          options={'xatol': tol})
                          
    # Bounds check
    if recursive:
        if sol.x < t0:
            ## Outside the lower bound, shift up in time
            sol = _search(obs, bdy, value, selection, t0+_SIDEREAL_DAY, t1+_SIDEREAL_DAY, recursive=False)
        elif sol.x > t1:
            ## Outside the upper bound, shift down in time
            sol = _search(obs, bdy, value, selection, t0-_SIDEREAL_DAY, t1-_SIDEREAL_DAY, recursive=False)
            
    # Return the full optimization solution
    return sol


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
        self.__lat = SkyAngle(degrees=0.0)
        self.__lon = SkyAngle(degrees=0.0)
        self.__elev = SkyDistance(m=0.0)
        self.__horz = SkyAngle(degrees=0.0)
        self.__date = Date(ts.now())
        self._update()
        
    def _update(self):
        self._wgs84 = Topos(latitude=self.__lat, longitude=self.__lon, elevation_m=self.__elev.m)
        
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
        
        return self.__elev.m
    @elev.setter
    def elev(self, value):
        self.__elev = SkyDistance(m=float(value))
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
        
    def sidereal_time(self):
        """
        Returns an Angle instance of the local apparent sidereal time for
        the Observer.
        """
        
        return hours(self.__date.gast*numpy.pi/12 + self._wgs84.longitude.radians)
        
    def radec_of(self, az, alt):
        """
        Given an azimuth and elevation as viewed by the Observer, return the
        right ascension and declination that they correspond to.  The RA and
        dec. values are returned as Angle instances in the ICRS frame.
        """
        
        alt = degrees(alt)
        az = degrees(az)
        
        obs = _ter + self._wgs84
        pos = obs.at(self.__date).from_altaz(alt_degrees=alt.degrees, az_degrees=az.degrees)
        ra, dec, _ = pos.radec()
        return hours(ra), degrees(dec)
        
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
            
        sol = _search(self, body, numpy.pi/2, 0, start-1.0, start*1.0)
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
        
        sol = _search(self, body, numpy.pi/2, 0, start*1.0, start+1.0)
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
            
        sol = _search(self, body, -numpy.pi/2, 0, start-1.0, start*1.0)
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
            
        sol = _search(self, body, -numpy.pi/2, 0, start*1.0, start+1.0)
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
            
        sol = _search(self, body, self.__horz.radians, 1, start-1.0, start*1.0)
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise AlwaysUpError()
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
            
        sol = _search(self, body, self.__horz.radians, 1, start*1.0, start+1.0)
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise AlwaysUpError()
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
            
        sol = _search(self, body, self.__horz.radians, 2, start-1.0, start*1.0)
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise AlwaysUpError()
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
            
        sol = _search(self, body, self.__horz.radians, 2, start*1.0, start+1.0)
        if body.neverup:
            raise NeverUpError()
        elif body.circumpolar:
            raise AlwaysUpError()
        return Date(sol.x)
