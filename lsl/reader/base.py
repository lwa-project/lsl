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

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['FrameHeaderBase', 'FramePayloadBase', 'FrameBase']


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
    _header_attrs = []
    
    def __repr__(self):
        n = self.__class__.__module__+'.'+self.__class__.__name__
        a = [(attr,getattr(self, attr, None)) for attr in self._header_attrs]
        return tw_fill(_build_repr(n,a), subsequent_indent='    ')


class FramePayloadBase(object):
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


class FrameBase(object):
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
        
    def __str__(self):
        return "%s for %s at time %s s" % (self.__class__.__module__, self.__class__.__name__, self.id, self.time)
        
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
        elif isinstance(y, (int, float)):
            tY = (int(y), y%1)
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
        elif isinstance(y, (int, float)):
            tY = (int(y), y%1)
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
        elif isinstance(y, (int, float)):
            tY = (int(y), y%1)
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
        elif isinstance(y, (int, float)):
            tY = (int(y), y%1)
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
        elif isinstance(y, (int, float)):
            tY = (int(y), y%1)
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
        elif isinstance(y, (int, float)):
            tY = (int(y), y%1)
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
