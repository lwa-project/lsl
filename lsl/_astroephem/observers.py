from __future__ import print_function, division

import numpy

from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, ICRS
import astropy.units as u

from scipy.optimize import fmin_tnc

from angles import degrees
from dates import Date, _DJD_OFFSET


__all__ = ['CircumpolarError', 'NeverUpError', 'AlwaysUpError', 'Observer']


class CircumpolarError(Exception):
    pass


class NeverUpError(CircumpolarError):
    pass


class AlwaysUpError(CircumpolarError):
    pass


class Observer(object):
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
        return self.__lat.to('radian').value
    @lat.setter
    def lat(self, value):
        self.__lat = degrees(value)
        self._update()
        
    @property
    def lon(self):
        return self.__lon.to('radian').value
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
        self.__elev = value*u.m
        self._update()
        
    @property
    def date(self):
        return self.__date
    @date.setter
    def date(self, value):
        self.__date = Date(value)
        
    def as_astropy(self):
        return self._el, self.__date
        
    def sidereal_time(self):
        return self.__date.sidereal_time('apparent', self.__lon)
        
    def radec_of(self, az, alt):
        alt = degrees(alt)
        az = degrees(az)
        topo = AltAz(alt, az, obstime=self.__date, location=self._el)
        equ = topo.transform_to(ICRS())
        return equ.ra.to('radian').value, equ.dec.to('radian').value 
        
    def _transit(self, delta, body, start=None):
        initial_date = self.__date
        if start is None:
            start = self.__date
            
        lst = self.sidereal_time()
        body.compute(self)
        if body.neverup:
            raise NeverUpError()
            
        diff = lst - body.ra*u.radian
        t_transit = self.__date + TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        if diff.value > 0 and delta > 0:
            t_transit = t_transit + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_transit = t_transit - 1*u.sday
            
        self.__date = initial_date
        
        return t_transit
        
    def next_transit(self, body, start=None):
        return self._transit(1, body, start=start)
        
    def previous_transit(self, body, start=None):
        return self._transit(-1, body, start=start)
        
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
        t_rise = self.__date + TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        if diff.value > 0 and delta > 0:
            t_rise = t_rise + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_rise = t_rise - 1*u.sday
            
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
        t_set = self.__date + TimeDelta(diff.to('radian').value*43200/numpy.pi*u.second, format='sec')
        if diff.value > 0 and delta > 0:
            t_set = t_set + 1*u.sday
        elif diff.value < 0 and delta < 0:
            t_set = t_set - 1*u.sday
            
        self.__date = initial_date
        
        return t_set
        
    def next_setting(self, body, start=None):
        return self._set(1, body, start=start)
        
    def previous_setting(self, body, start=None):
        return self._set(-1, body, start=start)