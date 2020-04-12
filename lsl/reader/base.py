# -*- coding: utf-8 -*-

"""
Python module that contains the base FrameHeader, FramePayload, and Frame 
classes for all of the LSL readers.

.. versionadded:: 1.3.0
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range

import copy
import numpy
from textwrap import fill as tw_fill
from datetime import datetime, timedelta

from astropy.time import Time as AstroTime

from lsl.common import dp as dp_common
from lsl.astro import unix_to_utcjd, MJD_OFFSET


__version__ = '0.2'
__all__ = ['FrameHeaderBase', 'FramePayloadBase', 'FrameBase', 'FrameTimestamp']


def _build_repr(name, attrs=[]):
    name = '.'.join(name.split('.')[-2:])
    output = "<%s" % name
    first = True
    for key,value in attrs:
        output += "%s %s=%s" % (('' if first else ','), key, value)
        first = False
    output += ">"
    return output


class FrameHeaderBase(object):
    """
    Base class for all lsl.reader FrameHeader-type objects.
    """
    
    _header_attrs = []
    
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in self._header_attrs]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')


class FramePayloadBase(object):
    """
    Base class for all lsl.reader FramePayload-type objects.
    """
    
    _payload_attrs = []
    
    def __init__(self, data):
        self._data = data
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in self._payload_attrs]
        if self._data is not None:
            a.append(('dtype',str(self._data.dtype)))
            a.append(('shape',str(self._data.shape)))
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    @property
    def data(self):
        return self._data


class FrameBase(object):
    """
    Base class for all lsl.reader Frame-type objects.
    """
    
    _header_class = FrameHeaderBase
    _payload_class = FramePayloadBase
    
    def __init__(self, header=None, payload=None, valid=True):
        if header is None:
            self.header = self._header_class()
        else:
            if not isinstance(header, self._header_class):
                raise TypeError("Excepted header of type '%s' but found '%s'" % (self._header_class.__type__.__name__, header.__type__.__name__))
            self.header = header
            
        if payload is None:
            self.payload = self._payload_class()
        else:
            if not isinstance(payload, self._payload_class):
                raise TypeError("Excepted payload of type '%s' but found '%s'" % (self._payload_class.__type__.__name__, payload.__type__.__name__))
            self.payload = payload
            
        self.valid = valid
        
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [('header',repr(self.header).replace(',\n    ', ', ')),
             ('payload',repr(self.payload).replace(',\n    ', ', ')),
             ('valid', self.valid)]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')
        
    def __add__(self, y):
        """
        Add the data sections of two frames together or add a number 
        to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (type(self), int, float, complex, numpy.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        newFrame = copy.deepcopy(self)
        newFrame += y
        return newFrame
        
    def __iadd__(self, y):
        """
        In-place add the data sections of two frames together or add 
        a number to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (type(self), int, float, complex, numpy.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        try:
            self.payload._data += y.payload._data
        except AttributeError:
            self.payload._data += self.payload._data.dtype.type(y)
        return self
        
    def __mul__(self, y):
        """
        Multiple the data sections of two frames together or multiply 
        a number to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (type(self), int, float, complex, numpy.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        newFrame = copy.deepcopy(self)
        newFrame *= y
        return newFrame
            
    def __imul__(self, y):
        """
        In-place multiple the data sections of two frames together or 
        multiply a number to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (type(self), int, float, complex, numpy.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        try:
            self.payload._data *= y.payload._data
        except AttributeError:
            self.payload._data *= self.payload._data.dtype.type(y)
        return self
        
    def __eq__(self, y):
        """
        Check if the time tags of two frames are equal or if the time
        tag is equal to a particular value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, FrameTimestamp)):
            tY = y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
        if tX == tY:
            return True
        else:
            return False
            
    def __ne__(self, y):
        """
        Check if the time tags of two frames are not equal or if the time
        tag is not equal to a particular value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, FrameTimestamp)):
            tY = y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
        if tX != tY:
            return True
        else:
            return False
            
    def __gt__(self, y):
        """
        Check if the time tag of the first frame is greater than that of a
        second frame or if the time tag is greater than a particular value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, FrameTimestamp)):
            tY = y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
        if tX > tY:
            return True
        else:
            return False
            
    def __ge__(self, y):
        """
        Check if the time tag of the first frame is greater than or equal to 
        that of a second frame or if the time tag is greater than a particular 
        value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, FrameTimestamp)):
            tY = y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
        if tX >= tY:
            return True
        else:
            return False
            
    def __lt__(self, y):
        """
        Check if the time tag of the first frame is less than that of a
        second frame or if the time tag is greater than a particular value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, FrameTimestamp)):
            tY = y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
        if tX < tY:
            return True
        else:
            return False
            
    def __le__(self, y):
        """
        Check if the time tag of the first frame is less than or equal to 
        that of a second frame or if the time tag is greater than a particular 
        value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, FrameTimestamp)):
            tY = y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
        if tX <= tY:
            return True
        else:
            return False
            
    def __cmp__(self, y):
        """
        Compare two frames based on the time tags.  This is helpful for 
        sorting things.
        """
        
        tX = self.time
        if not isinstance(y, FrameBase):
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
        tY = y.time
           
        if tY > tX:
            return -1
        elif tX > tY:
            return 1
        else:
            return 0


class FrameTimestamp(object):
    """
    Class to represent the UNIX timestamp of a data frame as an integer 
    number of seconds and a fractional number of seconds.
    """
    
    def __init__(self, si=0, sf=0.0):
        if isinstance(si, float):
            sf = si - int(si)
            si = int(si)
        self._int = int(si)
        self._frac = float(sf)
        
    @classmethod
    def from_dp_timetag(cls, value, offset=0):
        """
        Create a new FrameTimestamp instance from a raw DP timetag with an optional
        offset.
        """
        
        tt = int(value) - offset
        s = tt // int(dp_common.fS)
        f = (tt - s*int(dp_common.fS)) / dp_common.fS
        return cls(s, f)
        
    @classmethod
    def from_mjd_mpm(cls, mjd, mpm):
        """
        Create a new FrameTimestamp from a MJD/MPM (milliseconds past midnight) pair.
        """
        
        imjd = int(mjd)
        fmjd = mjd - imjd
        mpm = mpm + int(fmjd*86400*1000)
        s =  mpm // 1000
        f = (mpm - s*1000) / 1000.0
        s = s + (imjd - 40587)*86400
        return cls(s, f)
        
    def __str__(self):
        dt = self.datetime
        return str(dt)
        
    def __repr__(self):
        return "<FrameTimestamp i=%i, f=%.9f>" % (self._int, self._frac)
        
    def __int__(self):
        return self._int
        
    def __float__(self):
        return self._int+self._frac
        
    def __getitem__(self, i):
        if i == 0:
            return self._int
        elif i == 1:
            return self._frac
        else:
            raise IndexError
            
    def __add__(self, other):
        try:
            oi, of = other[0], other[1]
        except (TypeError, IndexError):
            oi = int(other)
            of = other - oi
        _int = self._int + oi
        _frac = self._frac + of
        if _frac >= 1:
            _int += 1
            _frac -= 1
        return FrameTimestamp(_int, _frac)
        
    def __iadd_(self, other):
        try:
            oi, of = other[0], other[1]
        except (TypeError, IndexError):
            oi = int(other)
            of = other - oi
        self._int += oi
        self._frac += of
        if self._frac >= 1:
            self._int += 1
            self._frac -= 1
            
    def __sub__(self, other):
        try:
            oi, of = other[0], other[1]
        except (TypeError, IndexError):
            oi = int(other)
            of = other - oi
        _int = self._int - oi
        _frac = self._frac - of
        if _frac < 0:
            _int -= 1
            _frac += 1
        return _int+_frac
       
    def __isub__(self, other):
        try:
            oi, of = other[0], other[1]
        except (TypeError, IndexError):
            oi = int(other)
            of = other - oi
        self._int -= oi
        self._frac -= of
        if self._frac < 0:
            self._int -= 1
            self._frac += 1
            
    def __eq__(self, y):
        if isinstance(y, FrameTimestamp):
            return self._int == y._int and self._frac == y._frac
        elif isinstance(y, int):
            return self._int == y and self._frac == 0.0
        elif isinstance(y, float):
            return float(self) == float(y)
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
    def __ne__(self, y):
        return not (self == y)
            
    def __gt__(self, y):
        if isinstance(y, FrameTimestamp):
            return (self._int > y._int) or (self._int == y._int and self._frac > y._frac)
        elif isinstance(y, int):
            return self._int > y or (self._int == y and self._frac > 0.0)
        elif isinstance(y, float):
            return float(self) > y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
    def __ge__(self, y):
        return (self > y or self == y)
        
    def __lt__(self, y):
        if isinstance(y, FrameTimestamp):
            return (self._int < y._int) or (self._int == y._int and self._frac < y._frac)
        elif isinstance(y, int):
            return self._int < y
        elif isinstance(y, float):
            return float(self) < y
        else:
            raise TypeError("Unsupported type: '%s'" % type(y).__name__)
            
    def __le__(self, y):
        return (self < y or self == y)
        
    def __cmp__(self, y):
        if y > self:
            return -1
        elif self > y:
            return 1
        else:
            return 0
            
    @property
    def unix(self):
        """
        UNIX timestamp as a floating point value.
        """
        
        return self._int + self._frac
            
    @property
    def mjd(self):
        """
        MJD as a floating point value.
        """
        
        return unix_to_utcjd(self) - MJD_OFFSET
        
    @property
    def pulsar_mjd(self):
        """
        MJD as  three-element tuple of integer number of MJD days, fractional
        MJD day, and fractional seconds.
        """
        
        days = self._int // 86400
        frac = (self._int - days*86400) / 86400.0
        return (days + 40587, frac, self._frac)
        
    @property
    def dp_timetag(self):
        """
        Timestamp as a DP timetag (ticks of a 196 MHz clock since UTC midnight
        on January 1, 1970).
        """
        
        tt = self._int * int(dp_common.fS)
        tt = tt + int(self._frac*dp_common.fS)
        return tt
        
    @property
    def datetime(self):
        """
        Timestamp as a naive `datetime.datetime` instance.
        """
        
        s = self._int
        us = int(self._frac*1e6)
        if us >= 1000000:
            s += 1
            us -= 1000000
        dt = datetime.utcfromtimestamp(s)
        dt += timedelta(microseconds=us)
        return dt
        
    @property
    def astropy(self):
        """
        Timestamp as an `astropy.time.Time` instance.
        """
        
        return AstroTime(self._int, self._frac, format='unix', scale='utc')