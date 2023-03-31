"""
Python module to reading in data from OVRO-LWA trigger voltage buffer dump files.
This module defines the following classes for storing the voltage buffer data
found in a file:

Frame
  object that contains all data associated with a particular TBF frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the TBF frame header object and
    * FramePayload   - the TBF frame data object.  
Combined, these two objects contain all of the information found in the 
original TBF frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import json
import numpy
import struct

from lsl.common import adp as adp_common
from lsl.reader.base import *
from lsl.reader._gofast import read_ovro
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.1'

class FrameHeader(FrameHeaderBase):
    _header_attrs = ['timetag', 'seq0', 'sync_time', 'pipeline_id', 'first_chan', 'nchan']
    
    def __init__(self, timetag=None, seq0=None, sync_time=None, pipeline_id=None, first_chan=None, nchan=None):
        self.timetag = timetag
        self.seq0 = seq0
        self.sync_time = sync_time
        self.pipeline_id = pipeline_id
        self.first_chan = first_chan
        self.nchan = nchan
        FrameHeaderBase.__init__(self)
        
    @property
    def time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance.
        """
        
        return FrameTimestamp.from_dp_timetag(self.sync_time*196000000 + self.timetag*8192)
        
    @property
    def channel_freqs(self):
        """
        Return a numpy.float32 array for the center frequencies, in Hz, of
        each channel in the data.
        """
        
        return (numpy.arange(self.nchan, dtype=numpy.float32)+self.first_chan) * 196e6 / 8192


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a OVRO-LWA
    triggered voltage buffer dump file.
    """
    
    _payload_attrs = []
    
    def __init__(self, fDomain=None):
        FramePayloadBase.__init__(self, fDomain)


class Frame(FrameBase):
    """
    Class that stores the information contained within a OVRO-LWA triggered
    voltage buffer dump file.  It's properties are FrameHeader and FramePayload
    objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    
    
    @property
    def channel_freqs(self):
        """
        Convenience wrapper for the Frame.FrameHeader.channel_freqs property.
        """
        
        return self.header.channel_freqs
        
    @property
    def time(self):
        """
        Convenience wrapper for the Frame.FramePayload.time property.
        """
        
        return self.header.time


def read_frame(filehandle, verbose=False):
    """
    Function to read in a single OVRO-LWA triggered voltage buffer dump file
    (header+data) and store the contents as a Frame object.
    """
    
    try:
        hsize, hblock_size = struct.unpack('<II', filehandle.read(8))
        header = json.loads(filehandle.read(hsize))
        nchan = header['nchan']
        nstand = header['nstand']
        npol = header['npol']
        filehandle.seek(hblock_size)
        ntime = os.path.getsize(filehandle.name) - filehandle.tell()
        ntime //= (nchan*nstand*npol*1) 
        
        data = read_ovro(filehandle, ntime, nchan, nstand, npol)
        
        newFrame = Frame()
        newFrame.header.timetag = header['time_tag']
        newFrame.header.seq0 = header['seq0']
        newFrame.header.sync_time = header['sync_time']
        newFrame.header.pipeline_id = header['pipeline_id']
        newFrame.header.first_chan = header['chan0']
        newFrame.header.nchan = nchan
        newFrame.payload._data = data
    except (OSError, gEOFError):
        raise EOFError
        
    return newFrame
