from __future__ import print_function, division

import numpy

from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, ICRS
import astropy.units as u

from angles import hours, degrees
from dates import Date


__all__ = ['CircumpolarError', 'NeverUpError', 'AlwaysUpError', 'Observer']


class CircumpolarError(Exception):
    pass


class NeverUpError(CircumpolarError):
    pass


class AlwaysUpError(CircumpolarError):
    pass


def _piecewise_search(func, x):
    x = [x+v*u.minute for v in (-5, -3, -2, -1, 0, 1, 2, 3, 5)]
    f = numpy.array([func(v) for v in x])
    x = numpy.array([float(v) for v in x])
    left = numpy.polyfit(x[:4], f[:4], 1)
    right = numpy.polyfit(x[5:], f[5:], 1)
    xL = -left[1]/left[0]
    xR = -right[1]/right[0]
    return (xL+xR)/2.0


def _quadratic_search(func, x):
    x = [x+v*u.minute for v in (-20, -15, -10, -5, 0, 5, 10, 15, 20)]
    f = [func(v) for v in x]
    fit = numpy.polyfit([float(v) for v in x], f, 2)
    x1 = -fit[1] / (2*fit[0])
    return x1


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
        
    def _transit(self, delta, body, start=None):
        initial_date = self.__date
        if start is None:
            start = self.__date
            
        lst = self.sidereal_time()
        body.compute(self)
        if body.neverup:
            raise NeverUpError()
            
        diff = lst - body.ra
        t_transit = self.__date - TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        if diff.value > 0 and delta > 0:
            t_transit = t_transit + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_transit = t_transit - 1*u.sday
        def ninety(t):
            self.__date = Date(t)
            body.compute(self)
            return abs(90 - body.alt.to('deg').value)
        results = _quadratic_search(ninety, t_transit)
        t_transit = Date(results)
        
        self.__date = initial_date
        
        return t_transit
        
    def next_transit(self, body, start=None):
        return self._transit(1, body, start=start)
        
    def previous_transit(self, body, start=None):
        return self._transit(-1, body, start=start)
        
    def _antitransit(self, delta, body, start=None):
        initial_date = self.__date
        if start is None:
            start = self.__date
            
        lst = self.sidereal_time()
        body.compute(self)
        if body.neverup:
            raise NeverUpError()
            
        diff = lst - body.ra
        t_antitransit = self.__date - TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        t_antitransit = t_antitransit - 0.5*u.sday
        if diff.value > 0 and delta > 0:
            t_antitransit = t_antitransit + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_antitransit = t_antitransit - 1*u.sday
        def neg_ninety(t):
            self.__date = Date(t)
            body.compute(self)
            return body.alt.to('deg').value
        results = _quadratic_search(neg_ninety, t_antitransit)
        t_antitransit = Date(results)
        
        self.__date = initial_date
        
        return t_antitransit
        
    def next_antitransit(self, body, start=None):
        return self._antitransit(1, body, start=start)
        
    def previous_antitransit(self, body, start=None):
        return self._antitransit(-1, body, start=start)
        
    def _rise(self, delta, body, start=None):
        initial_date = self.__date
        if start is None:
            start = self.__date
            
        lst = self.sidereal_time()    
        body.compute(self)
        if body.neverup:
            raise NeverUpError()
        if body.circumpolar:
            raise AlwaysUpError()
            
        diff = lst - body._rising_lst
        t_rise = self.__date - TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        if diff.value > 0 and delta > 0:
            t_rise = t_rise + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_rise = t_rise - 1*u.sday
        def zero(t):
            self.__date = Date(t)
            body.compute(self)
            return abs(body.alt.to('arcsec').value)
        results = _piecewise_search(zero, t_rise)
        t_rise = Date(results)
           
        self.__date = initial_date
        
        return t_rise
        
    def next_rising(self, body, start=None):
        return self._rise(1, body, start=start)
        
    def previous_rising(self, body, start=None):
        return self._rise(-1, body, start=start)
        
    def _set(self, delta, body, start=None):
        initial_date = self.__date
        if start is None:
            start = self.__date
            
        lst = self.sidereal_time()   
        body.compute(self)
        if body.neverup:
            raise NeverUpError()
        if body.circumpolar:
            raise AlwaysUpError()
            
        diff = lst - body._setting_lst
        t_set = self.__date - TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        if diff.value > 0 and delta > 0:
            t_set = t_set + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_set = t_set - 1*u.sday
        def zero(t):
            self.__date = Date(t)
            body.compute(self)
            return abs(body.alt.to('arcsec').value)
        results = _piecewise_search(zero, t_set)
        t_set = Date(results)
        
        self.__date = initial_date
        
        return t_set
        
    def next_setting(self, body, start=None):
        return self._set(1, body, start=start)
        
    def previous_setting(self, body, start=None):
        return self._set(-1, body, start=start)