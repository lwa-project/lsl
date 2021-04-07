"""
Python module providing an ephem-like interface to dates and times.
"""

from __future__ import print_function, division

import time
import calendar
import datetime
from functools import total_ordering

from astropy.time import Time

from lsl._astroephem.config import PYEPHEM_REPR

__all__ = ['B1950', 'J2000', 'Date', 'now', 'localtime']


_DJD_OFFSET = 2415020.0


@total_ordering
class EphemTime(Time):
    """
    Base class for representing times in a way that behaves like ephem.date.
    """
    
    def __repr__(self):
        if PYEPHEM_REPR:
            return str(float(self))
        else:
            return Time.__repr__(self)
            
    def __str__(self):
        if PYEPHEM_REPR:
            return "%i/%i/%i %02i:%02i:%06.3f" % self.tuple()
        else:
            return Time.__str__(self)
            
    def __float__(self):
        return float(self.jd - _DJD_OFFSET)
        
    def __round__(self, ndigits=0):
        return round(float(self), ndigits)
        
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
        """
        Dublin Julian Date
        """
        
        return float(self)
         
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
        
        dt = self.datetime
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


B1950 = EphemTime(1950.0, format='byear', scale='utc')
J2000 = EphemTime(2000.0, format='jyear', scale='utc')


def Date(value):
    """
    Build an EphemTime instance representing a date and time.
    """
    
    date = None
    if isinstance(value, Time):
        date = EphemTime(value)
    elif isinstance(value, (int, float)):
        value = value + _DJD_OFFSET
        date = EphemTime(value, format='jd', scale='utc')
    elif isinstance(value, str):
        value = value.replace('/', '-')
        date = EphemTime(value, format='iso', scale='utc')
    if date is None:
        raise ValueError("Cannot parse '%s'" % str(value))
    return date


def now():
    """
    Return an EphemTime instance representing the current time.
    """
    
    return Date(Time.now())


def localtime(date):
    """
    Given a EphemTime instance, return a Python datetime instance that
    represents in the input as a local time.
    """
    
    timetuple = time.localtime(calendar.timegm(date.tuple()))
    return datetime.datetime(*timetuple[:7])
    
