"""
Python module that contains the base FrameHeader, FramePayload, and Frame 
classes for all of the LSL readers.

.. versionadded:: 2.0.0
"""

import copy
import pytz
import numpy as np
from functools import total_ordering
from textwrap import fill as tw_fill
from datetime import datetime, timedelta

from astropy.time import Time as AstroTime, TimeDelta as AstroDelta

from lsl.common import dp as dp_common
from lsl.astro import MJD_OFFSET


__version__ = '0.4'
__all__ = ['FrameHeaderBase', 'FramePayloadBase', 'FrameBase', 'FrameTimestamp', 'CI8']


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


@total_ordering
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
                raise TypeError(f"Excepted header of type '{self._header_class.__type__.__name__}' but found '{header.__type__.__name__}'")
            self.header = header
            
        if payload is None:
            self.payload = self._payload_class()
        else:
            if not isinstance(payload, self._payload_class):
                raise TypeError(f"Excepted payload of type '{self._payload_class.__type__.__name__}' but found '{payload.__type__.__name__}'")
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
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
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
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
        try:
            self.payload._data += y.payload._data
        except AttributeError:
            self.payload._data += self.payload._data.dtype.type(y)
        return self
        
    def __sub__(self, y):
        """
        Subtract the data sections of two frames or subtract a number 
        from every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
        newFrame = copy.deepcopy(self)
        newFrame -= y
        return newFrame
        
    def __isub__(self, y):
        """
        In-place subtract the data sections of two frames together or subtract 
        a number from every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
        try:
            self.payload._data -= y.payload._data
        except AttributeError:
            self.payload._data -= self.payload._data.dtype.type(y)
        return self
        
    def __mul__(self, y):
        """
        Multiple the data sections of two frames together or multiply 
        a number to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        try:
            self.payload._data *= y.payload._data
        except AttributeError:
            self.payload._data *= self.payload._data.dtype.type(y)
        return self
        
    def __floordiv__(self, y):
        """
        Divide the data sections of two frames together or divide 
        a number into every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        newFrame = copy.deepcopy(self)
        newFrame //= y
        return newFrame
        
    def __ifloordiv__(self, y):
        """
        In-place divide the data sections of two frames together or 
        divide a number into every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        try:
            self.payload._data //= y.payload._data
        except AttributeError:
            self.payload._data //= self.payload._data.dtype.type(y)
        return self
        
    def __truediv__(self, y):
        """
        Divide the data sections of two frames together or divide 
        a number into every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        newFrame = copy.deepcopy(self)
        newFrame /= y
        return newFrame
        
    def __itruediv__(self, y):
        """
        In-place divide the data sections of two frames together or 
        divide a number into every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
        if not isinstance(y, (FrameBase, int, float, complex, np.ndarray)):
            raise TypeError("Unsupported type '%s'" % type(y).__name__)
            
        try:
            self.payload._data /= y.payload._data
        except AttributeError:
            self.payload._data /= self.payload._data.dtype.type(y)
        return self
        
    def __div__(self, y):
        return self.__floordiv__(y)
        
    def __idiv__(self, y):
        return self.__ifloordiv__(y)
        
    def __eq__(self, y):
        """
        Check if the time tags of two frames are equal or if the time
        tag is equal to a particular value.
        """
        
        tX = self.time
        if isinstance(y, FrameBase):
            tY = y.time
        elif isinstance(y, (int, float, np.integer, np.floating, FrameTimestamp)):
            tY = y
        else:
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
        if tX == tY:
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
        elif isinstance(y, (int, float, np.integer, np.floating, FrameTimestamp)):
            tY = y
        else:
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
        if tX < tY:
            return True
        else:
            return False


@total_ordering
class FrameTimestamp(object):
    """
    Class to represent the UNIX timestamp of a data frame as an integer 
    number of seconds and a fractional number of seconds.
    """
    
    def __init__(self, si=0, sf=0.0):
        if isinstance(si, AstroTime):
            self._time = si
            if sf != 0.0:
                self._time += AstroDelta(sf, format='sec')
        else:
            if isinstance(si, (float, np.floating)):
                sf = sf + (si - int(si))
                si = int(si)
            # Make sure sf is [0.0, 1.0)
            if sf >= 1:
                sfi = int(sf)
                sff = sf - sfi
                si += sfi
                sf = sff
            elif sf < 0:
                sfi = int(sf) - 1
                sff = sf - sfi
                si += sfi
                sf = sff
            self._time = AstroTime(si, sf, format='unix', scale='utc')
            
    @classmethod
    def now(cls):
        """
        Create a new FrameTimestamp instance for the current time as determined
        from `time.time()`.
        """
        
        return cls(AstroTime.now())
        
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
        
        t = AstroTime(mjd, mpm//1000 / 86400, format='mjd', scale='utc')
        f = (mpm % 1000) / 1000.0
        return cls(t, f)
        
    @classmethod
    def from_pulsar_mjd(cls, mjd, mjd_frac, sec_frac):
        """
        Create a new FrameTimstamp from a three-element tuple of integer number 
        of MJD days, integer seconds since 0h on that day, and fractional seconds.
        """
        
        t = AstroTime(mjd, mjd_frac/86400, format='mjd', scale='utc')
        return cls(t, sec_frac)
        
    @classmethod
    def from_astropy(cls, t):
        """
        Create a new FrameTimestamp from a astropy.time.Time object.
        """
        
        return cls(t, 0.0)
        
    def __str__(self):
        dt = self._time.datetime
        return str(dt)
        
    def __repr__(self):
        t = self._time.unix
        return "<FrameTimestamp i=%i, f=%.9f (%.9f %.9f)>" % (int(t), t-int(t), self._time.jd1, self._time.jd2)
        
    def __format__(self, format_spec):
        if format_spec == '' or format_spec[-1] == 's':
            return format(str(self), format_spec)
        elif format_spec[-1] in ('e', 'E', 'f', 'F', 'g', 'G', 'n'):
            return format(float(self), format_spec)
        elif format_spec[-1] in ('d', 'n'):
            t = self._time.unix
            return format(int(t), format_spec)
        else:
            raise TypeError("unsupported format string passed to %s.__format__" % type(self).__name__)
            
    def __float__(self):
        return float(self._time.unix)
        
    def __getitem__(self, i):
        t = self._time.unix
        if i == 0:
            return int(t)
        elif i == 1:
            return t - int(t)
        else:
            raise IndexError
            
    def __add__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            t = self._time + AstroDelta(other, format='sec')
            return FrameTimestamp(t)
        elif isinstance(other, timedelta):
            t = self._time + AstroDelta(other, format='datetime')
            return FrameTimestamp(t)
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __iadd__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            self._time += AstroDelta(other, format='sec')
            return self
        elif isinstance(other, timedelta):
            self._time += AstroDelta(other, format='datetime')
            return self
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __sub__(self, other):
        if isinstance(other, FrameTimestamp):
            diff = self._time - other._time
            return diff.sec
        elif isinstance(other, (int, float, np.integer, np.floating)):
            t = self._time - AstroDelta(other, format='sec')
            return FrameTimestamp(t)
        elif isinstance(other, timedelta):
            t = self._time - AstroDelta(other, format='datetime')
            return FrameTimestamp(t)
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __isub__(self, other):
        if isinstance(other, (int, float, np.integer, np.floating)):
            self._time -= AstroDelta(other, format='sec')
            return self
        elif isinstance(other, timedelta):
            self._time -= AstroDelta(other, format='datetime')
            return self
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __eq__(self, y):
        if isinstance(y, FrameTimestamp):
            return self._time == y._time
        elif isinstance(y, (int, np.integer, float, np.floating)):
            return float(self) == y
        else:
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
    def __lt__(self, y):
        if isinstance(y, FrameTimestamp):
            return self._time < y._time
        elif isinstance(y, (int, np.integer, float, np.floating)):
            return float(self) < y
        else:
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
    @property
    def unix(self):
        """
        UNIX timestamp as a floating point value.
        """
        
        return self._time.unix
        
    @property
    def jd(self):
        """
        UTC JD as a floating point value.
        """
        
        return self._time.jd
        
    @property
    def mjd(self):
        """
        UTC MJD as a floating point value.
        """
        
        return self._time.mjd
        
    @property
    def pulsar_mjd(self):
        """
        UTC MJD as  three-element tuple of integer number of MJD days,
        integer number of seconds since 0h on that day, and fractional seconds.
        """
        
        j1i = int(self._time.jd1)
        j1f = self._time.jd1 - j1i
        j2i = int(self._time.jd2)
        j2f = self._time.jd2 - j2i
        
        mjd = j1i + j2i - MJD_OFFSET
        mjd_frac1 = j1f + (mjd - int(mjd))
        mjd_frac2 = j2f
        mjd = int(mjd)
        while mjd_frac1 >= 1.0:
            mjd += 1
            mjd_frac1 -= 1
        while mjd_frac2 >= 1.0:
            mjd += 1
            mjd_frac2 -= 1
        while mjd_frac1 < 0.0:
            mjd -= 1
            mjd_frac1 += 1
        while mjd_frac2 < 0.0:
            mjd -= 1
            mjd_frac2 += 1
            
        day_sec1 = mjd_frac1 * 86400
        day_sec2 = mjd_frac2 * 86400
        day_sec = int(day_sec1) + int(day_sec2)
        if day_sec > 86400:
            mjd += 1
            day_sec -= 86400
        sec_frac = (day_sec1 - int(day_sec1)) + (day_sec2 - int(day_sec2))
        return mjd, day_sec, sec_frac
        
    @property
    def tai_mjd(self):
        """
        TAI MJD as a floating point value.
        """
        
        return self._time.tai.mjd
        
    @property
    def dp_timetag(self):
        """
        Timestamp as a DP timetag (ticks of a 196 MHz clock since UTC midnight
        on January 1, 1970).
        """
        
        sec = int(self._time.unix)
        sec_frac = self._time - AstroTime(sec, format='unix', scale='utc')
        
        tt = sec * int(dp_common.fS)
        tt = tt + int(round(sec_frac.sec*dp_common.fS, 2))
        return tt
        
    @property
    def datetime(self):
        """
        Timestamp as a naive `datetime.datetime` instance in UTC.
        """
        
        return self._time.datetime
        
    @property
    def utc_datetime(self):
        """
        Timestamp as a time zone-aware datetime instance in UTC.
        """
        
        return pytz.utc.localize(self.datetime)
        
    @property
    def astropy(self):
        """
        Timestamp as an `astropy.time.Time` instance.
        """
        
        return self._time


CI8 = np.dtype([('re', 'i1'), ('im', 'i1')])
