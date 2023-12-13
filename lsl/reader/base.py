"""
Python module that contains the base FrameHeader, FramePayload, and Frame 
classes for all of the LSL readers.

.. versionadded:: 2.0.0
"""

import copy
import pytz
import numpy
from functools import total_ordering
from textwrap import fill as tw_fill
from datetime import datetime, timedelta

from astropy.time import Time as AstroTime

from lsl.common import dp as dp_common
from lsl.astro import unix_to_utcjd, MJD_OFFSET, unix_to_taimjd


__version__ = '0.3'
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        
        if not isinstance(y, (FrameBase, int, float, complex, numpy.ndarray)):
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
        elif isinstance(y, (int, float, numpy.integer, numpy.floating, FrameTimestamp)):
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
        elif isinstance(y, (int, float, numpy.integer, numpy.floating, FrameTimestamp)):
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
        if isinstance(si, (float, numpy.floating)):
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
        self._int = int(si)
        self._frac = float(sf)
        
    @classmethod
    def now(cls):
        """
        Create a new FrameTimestamp instance for the current time as determined
        from `time.time()`.
        """
        
        return cls(time.time())
        
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
        
    @classmethod
    def from_pulsar_mjd(cls, mjd, mjd_frac, sec_frac):
        """
        Create a new FrameTimstamp from a three-element tuple of integer number 
        of MJD days, fractional MJD day, and fractional seconds.
        """
        
        s = (mjd - 40587)*86400
        f = sec_frac
        s1 = int(mjd_frac * 86400)
        return cls(s+s1, f)
        
    def __str__(self):
        dt = self.datetime
        return str(dt)
        
    def __repr__(self):
        return "<FrameTimestamp i=%i, f=%.9f>" % (self._int, self._frac)
        
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
        if isinstance(other, (int, float, numpy.integer, numpy.floating)):
            oi = int(other)
            of = other - oi
            _int = self._int + oi
            _frac = self._frac + of
            if _frac >= 1:
                _int += 1
                _frac -= 1
            return FrameTimestamp(_int, _frac)
        elif isinstance(other, timedelta):
            oi = other.days*86400 + other.seconds
            of = other.microseconds/1e6
            _int = self._int + oi
            _frac = self._frac + of
            if _frac >= 1:
                _int += 1
                _frac -= 1
            return FrameTimestamp(_int, _frac)
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __iadd__(self, other):
        if isinstance(other, (int, float, numpy.integer, numpy.floating)):
            oi = int(other)
            of = other - oi
            self._int += oi
            self._frac += of
            if self._frac >= 1:
                self._int += 1
                self._frac -= 1
            return self
        elif isinstance(other, timedelta):
            oi = other.days*86400 + other.seconds
            of = other.microseconds/1e6
            self._int += oi
            self._frac += of
            if self._frac >= 1:
                self._int += 1
                self._frac -= 1
            return self
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __sub__(self, other):
        if isinstance(other, FrameTimestamp):
            oi, of = other[0], other[1]
            _int = self._int - oi
            _frac = self._frac - of
            if _frac < 0:
                _int -= 1
                _frac += 1
            return _int+_frac
        elif isinstance(other, (int, float, numpy.integer, numpy.floating)):
            oi = int(other)
            of = other - oi
            _int = self._int - oi
            _frac = self._frac - of
            if _frac < 0:
                _int -= 1
                _frac += 1
            return FrameTimestamp(_int, _frac)
        elif isinstance(other, timedelta):
            oi = other.days*86400 + other.seconds
            of = other.microseconds/1e6
            _int = self._int - oi
            _frac = self._frac - of
            if _frac < 0:
                _int -= 1
                _frac += 1
            return FrameTimestamp(_int, _frac)
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __isub__(self, other):
        if isinstance(other, (int, float, numpy.integer, numpy.floating)):
            oi = int(other)
            of = other - oi
            self._int -= oi
            self._frac -= of
            if self._frac < 0:
                self._int -= 1
                self._frac += 1
            return self
        elif isinstance(other, timedelta):
            oi = other.days*86400 + other.seconds
            of = other.microseconds/1e6
            self._int -= oi
            self._frac -= of
            if self._frac < 0:
                self._int -= 1
                self._frac += 1
            return self
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __eq__(self, y):
        if isinstance(y, FrameTimestamp):
            return self._int == y._int and self._frac == y._frac
        elif isinstance(y, (int, numpy.integer)):
            return self._int == y and self._frac == 0.0
        elif isinstance(y, (float, numpy.floating)):
            return float(self) == float(y)
        else:
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
    def __lt__(self, y):
        if isinstance(y, FrameTimestamp):
            return (self._int < y._int) or (self._int == y._int and self._frac < y._frac)
        elif isinstance(y, (int, numpy.integer)):
            return self._int < y
        elif isinstance(y, (float, numpy.floating)):
            return float(self) < y
        else:
            raise TypeError(f"Unsupported type: '{type(y).__name__}'")
            
    @property
    def unix(self):
        """
        UNIX timestamp as a floating point value.
        """
        
        return float(self)
        
    @property
    def jd(self):
        """
        UTC JD as a floating point value.
        """
        
        return unix_to_utcjd(self)
        
    @property
    def mjd(self):
        """
        UTC MJD as a floating point value.
        """
        
        return self.jd - MJD_OFFSET
        
    @property
    def pulsar_mjd(self):
        """
        UTC MJD as  three-element tuple of integer number of MJD days,
        fractional MJD day, and fractional seconds.
        """
        
        days = self._int // 86400
        frac = (self._int - days*86400) / 86400.0
        return (days + 40587, frac, self._frac)
        
    @property
    def tai_mjd(self):
        """
        TAI MJD as a floating point value.
        """
        
        return unix_to_taimjd(self)
        
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
        Timestamp as a naive `datetime.datetime` instance in UTC.
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
        
        return AstroTime(self._int, self._frac, format='unix', scale='utc')


CI8 = numpy.dtype([('re', 'i1'), ('im', 'i1')])
