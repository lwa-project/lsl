from __future__ import print_function, division

import time
import calendar
import datetime

from astropy.time import Time

from functools import total_ordering 

__all__ = ['B1950', 'J2000', 'Date', 'now', 'localtime']


_DJD_OFFSET = 2415020.0


@total_ordering
class _FloatableTime(Time):
    def __str__(self):
        return "%i/%i/%i %02i:%02i:%06.3f" % self.tuple()
        
    def __float__(self):
        return self.jd - _DJD_OFFSET
        
    def __add__(self, other):
        if isinstance(other, (int, float)):
            return float(self) + other
        else:
            return Time.__add__(self, other)
            
    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return float(self) - other
        else:
            return Time.__sub__(self, other)
            
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return float(self)*other
        else:
            return Time.__mul__(self, other)
            
    def __neg__(self):
        return -float(self)
        
    def __eq__(self, other):
        if isinstance(other, (int, float)):
            return float(self) == other
        else:
            return Time.__eq__(self, other)
            
    def __lt__(self, other):
        if isinstance(other, (int, float)):
            return float(self) < float(other)
        else:
            return Time.__lt__(self, other)
           
    @property
    def djd(self):
        return float(self)
         
    def tuple(self):
        dt = self.datetime
        return (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second+dt.microsecond/1e6)
        
    def triple(self):
        y, m, d, _, _, _ = self.tuple()
        return (y, m, d)


B1950 = _FloatableTime(1950.0, format='byear', scale='utc')
J2000 = _FloatableTime(2000.0, format='jyear', scale='utc')


def Date(value):
    date = None
    if isinstance(value, Time):
        date = _FloatableTime(value)
    elif isinstance(value, (int, float)):
        value = value + _DJD_OFFSET
        date = _FloatableTime(value, format='jd', scale='utc')
    elif isinstance(value, str):
        value = value.replace('/', '-')
        date = _FloatableTime(value, format='iso', scale='utc')
    if date is None:
        raise ValueError("Cannot parse '%s'" % str(value))
    return date


def now():
    return Date(Time.now())


def localtime(date):
    timetuple = time.localtime(calendar.timegm(date.tuple()))
    return datetime.datetime(*timetuple[:7])
    
