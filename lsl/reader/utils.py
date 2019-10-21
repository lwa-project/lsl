# -*- coding: utf-8 -*-

"""
Python module that contains various helpers for the lsl.reader module.

.. versionadded:: 1.3.0
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
from bisect import bisect

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['FilePositionSaver', 'SplitFileWrapper']


class FilePositionSaver(object):
    """
    Simple class to save the current location in a file and
    return to that position when we are done with it.
    """
    
    def __init__(self, fh):
        self.fh = fh
        self.marker = self.fh.tell()
        
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, exc_tb):
        self.fh.seek(self.marker, 0)


class SplitFileWrapper(object):
    """
    Class to allow seamless access to a file that has been split into parts.
    """
    
    def __init__(self, filenames, sort=True):
        self._filenames = filenames
        if sort:
            self._filenames.sort()
            
        self._nfiles = len(self._filenames)
        self._sizes = [os.path.getsize(filename) for filename in self._filenames]
        self._offsets = [0,]
        for size in self._sizes[:-1]:
            self._offsets.append(self._offsets[-1]+size)
        self._total_size = sum(self._sizes)
        
        self.name = "%s+%i_others" % (self._filenames[0], self._nfiles-1)
        self.mode = 'rb'
        self._open_part(0)
        self.closed = False
        
    def __enter__(self):
	    return self
        
    def __exit__(self, exc_type, exc_value, exc_tb):
        self.close()
        
    def _open_part(self, index=None):
        try:
            self._fh.close()
        except AttributeError:
            pass
            
        if index is None:
            index = min([self._idx + 1, self._nfiles-1])
        self._idx = index
        self._fh = open(self._filenames[self._idx], self.mode)
        self._pos = self._offsets[self._idx]
        
    def close(self):
        if not self.closed:
            self.closed = True
            self._fh.close()
            del self._fh
            del self._pos
            
    def read(self, size):
        if self.closed:
            raise ValueError("I/O operation on a closed file")
            
        data = self._fh.read(size)
        self._pos += len(data)
        if len(data) < size and self._idx < self._nfiles-1:
            self._open_part()
            new_data = self._fh.read(size-len(data))
            self._pos += len(new_data)
            data += new_data
        return data
        
    @property
    def size(self):
        return self._total_size
        
    def tell(self):
        if self.closed:
            raise ValueError("I/O operation on a closed file")
            
        return self._pos
        
    def seek(self, pos, whence=0):
        if self.closed:
            raise ValueError("I/O operation on a closed file")
            
        if whence == 0:
            idx = bisect(self._offsets, pos) - 1
            self._open_part(idx)
            offset = pos - self._pos
            self._fh.seek(offset, 0)
            self._pos += offset
            
        elif whence == 1:
            offset = self._pos + pos
            self.seek(offset, 0)
            
        else:
            offset = self._total_size - pos
            self.seek(offset, 0)
