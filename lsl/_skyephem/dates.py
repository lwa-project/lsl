"""
Python module providing an ephem-like interface to dates and times.
"""

from __future__ import print_function, division

import time
import calendar
from datetime import datetime
from functools import total_ordering

from skyfield import api
from skyfield.timelib import Time as SkyTime
ts = api.load.timescale()

from astropy.time import Time as AstroTime

from lsl._skyephem.config import PYEPHEM_REPR

__all__ = ['Date', 'B1950', 'J2000', 'now', 'localtime']


_DJD_OFFSET = 2415020.0


@total_ordering
class Date(SkyTime):
    """
    Base class for representing times in a way that behaves like ephem.Date.
    """
    
    def __init__(self, value):
        date = None
        if isinstance(value, SkyTime):
            date = value
        elif isinstance(value, AstroTime):
            date = ts.from_astropy(value)
        elif isinstance(value, (int, float)):
            value = value + _DJD_OFFSET
            value = AstroTime(value, format='jd', scale='utc')
            date = ts.from_astropy(value)
        elif isinstance(value, str):
            value = value.replace('/', '-')
            value = AstroTime(value, format='iso', scale='utc')
            date = ts.from_astropy(value)
        if date is None:
            raise ValueError("Cannot parse '%s'" % str(value))
            
        SkyTime.__init__(self, date.ts, date.whole, date.tt_fraction)
        
    @property
    def utc_jd(self):
        value = AstroTime(self.whole, self.tt_fraction, format='jd', scale='tt')
        return value.utc.jd
        
    def __repr__(self):
        if PYEPHEM_REPR:
            return str(float(self))
        else:
            return SkyTime.__repr__(self)
            
    def __str__(self):
        if PYEPHEM_REPR:
            return "%i/%i/%i %02i:%02i:%06.3f" % self.tuple()
        else:
            return SkyTime.__str__(self)
            
    def __float__(self):
        return float(self.utc_jd - _DJD_OFFSET)
        
    def __add__(self, other):
        if isinstance(other, (int, float)):
            return float(self) + other
        else:
            return SkyTime.__add__(self, other)
            
    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return float(self) - other
        else:
            return SkyTime.__sub__(self, other)
            
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return float(self) * other
        else:
            return SkyTime.__mul__(self, other)
            
    def __neg__(self):
        return -float(self)
        
    def __eq__(self, other):
        if isinstance(other, (int, float)):
            return float(self) == other
        else:
            return SkyTime.__eq__(self, other)
            
    def __lt__(self, other):
        if isinstance(other, (int, float)):
            return float(self) < float(other)
        else:
            return SkyTime.__lt__(self, other)
            
    def tuple(self):
        """
        Returns a six-element tuple representing the date and time.  The
        elements are:
         * year
         * month
         * day
         * hours
         * minutes
         * seconds
        """
        
        dt = self.utc_datetime()
        return (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second+dt.microsecond/1e6)
        
    def triple(self):
        """
        Returns a three-element tuple representing the date.  The elements are:
         * year
         * month
         * day
        """
        
        y, m, d, _, _, _ = self.tuple()
        return (y, m, d)


B1950 = Date(AstroTime(1950.0, format='byear', scale='utc'))
J2000 = Date(AstroTime(2000.0, format='jyear', scale='utc'))


def now():
    """
    Return an EphemTime instance representing the current time.
    """
    
    return Date(ts.now())


def localtime(date):
    """
    Given a EphemTime instance, return a Python datetime instance that
    represents in the input as a local time.
    """
    
    timetuple = time.localtime(calendar.timegm(date.tuple()))
    return datetime(*timetuple[:7])
