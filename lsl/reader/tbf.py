"""
Python module to reading in data from TBF files.  This module defines the 
following classes for storing the TBF data found in a file:

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

..versionchanged:: 2.1.3
    Added a new data_ci8 field to the FramePayload to store data more
    efficiently.

.. versionadded:: 1.2.0
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import numpy

from lsl.common import adp as adp_common
from lsl.reader.base import *
from lsl.reader._gofast import read_tbf, read_tbf_ci8
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'read_frame',
           'FRAME_SIZE', 'FRAME_CHANNEL_COUNT', 'get_frames_per_obs',
           'get_first_frame_count', 'get_channel_count', 'get_first_channel']

#: TBF packet size (header + payload)
FRAME_SIZE = 6168

#: Number of frequency channels in a TBF packet
FRAME_CHANNEL_COUNT = 12


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a TBF 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    _header_attrs = ['adp_id', 'frame_count', 'second_count', 'first_chan']
    
    def __init__(self, adp_id=None, frame_count=None, second_count=None, first_chan=None):
        self.adp_id = adp_id
        self.frame_count = frame_count
        self.second_count = second_count
        self.first_chan = first_chan
        FrameHeaderBase.__init__(self)
        
    @property
    def is_tbf(self):
        """
        Function to check if the data is really TBF.  Returns True if the 
        data is TBF, false otherwise.
        """
        
        if self.adp_id == 0x01:
            return True
        else:
            return False
            
    @property
    def channel_freqs(self):
        """
        Return a numpy.float32 array for the center frequencies, in Hz, of
        each channel in the data.
        """
        
        return (numpy.arange(FRAME_CHANNEL_COUNT, dtype=numpy.float32)+self.first_chan) * adp_common.fC


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a TBF
    frame.  Both fields listed in the DP ICD version H are stored.
    """
    
    _payload_attrs = ['timetag']
    
    def __init__(self, timetag=None, fDomain=None, fDomain_ci8=None):
        self.timetag = timetag
        FramePayloadBase.__init__(self, fDomain)
        
        if fDomain_ci8 is not None:
            self._data_ci8 = fDomain_ci8
            del self._data
            
    @property
    def time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance.
        """
        
        return FrameTimestamp.from_dp_timetag(self.timetag)
        
    @property
    def data(self):
        try:
            assert(self._data is not None)
            return self._data
        except (AttributeError, AssertionError):
            self._data = self._data_ci8[:,:,:,0] + 1j*self._data_ci8[:,:,:,1]
            self._data = self._data.astype(numpy.complex64)
            return self._data
            
    @property
    def data_ci8(self):
        """
        Read-only data stored as array of numpy.int8 values with an additional
        axis for the real/imaginary values.
        
        .. note:: This field is unaffected by mathematical operations performed
                  on its parent Frame.
        """
        
        return self._data_ci8


class Frame(FrameBase):
    """
    Class that stores the information contained within a single TBF 
    frame.  It's properties are FrameHeader and FramePayload objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    
    @property
    def is_tbf(self):
        """
        Convenience wrapper for the Frame.FrameHeader.is_tbf property.
        """
        
        return self.header.is_tbf
        
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
        
        return self.payload.time


def read_frame(filehandle, verbose=False):
    """
    Function to read in a single TBF frame (header+data) and store the 
    contents as a Frame object.
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_tbf(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
        
    return newFrame


def get_frames_per_obs(filehandle):
    """
    Find out how many frames are present per time stamp by examining the 
    first 1000 TBF records.  Return the number of frames per observation.
    """
    
    with FilePositionSaver(filehandle):
        # Build up the list-of-lists that store the index of the first frequency
        # channel in each frame.
        channels = []
        for i in range(1000):
            try:
                cFrame = read_frame(filehandle)
                if not cFrame.is_tbf:
                    continue
            except EOFError:
                break
            except SyncError:
                continue
                
            chan = cFrame.header.first_chan
            if chan not in channels:
                channels.append( chan )
                
    # Return the number of channels
    return len(channels)


def get_first_frame_count(filehandle):
    """
    Find and return the lowest frame count encountered in a TBF file.
    """
    
    # Find out how many frames there are per observation
    nFrames = get_frames_per_obs(filehandle)
    
    with FilePositionSaver(filehandle):
        firstFrameCount = 2**64-1
        freqs = []
        while len(freqs) < nFrames:
            cFrame = read_frame(filehandle)
            freq = cFrame.header.first_chan
            
            if freq not in freqs:
                freqs.append(freq)
            if cFrame.header.frame_count < firstFrameCount:
                firstFrameCount = cFrame.header.frame_count
                
    # Return the lowest frame number found
    return firstFrameCount


def get_channel_count(filehandle):
    """
    Find out the total number of channels that are present by examining 
    the first 1000 TBF records.  Return the number of channels found.
    """
    
    # Find out how many frames there are per observation
    nFrames = get_frames_per_obs(filehandle)
    
    # Convert to channels
    nChannels = nFrames * FRAME_CHANNEL_COUNT
    
    # Return the number of channels
    return nChannels


def get_first_channel(filehandle, frequency=False):
    """
    Find and return the lowest frequency channel in a TBF file.  If the 
    `frequency` keyword is True the returned value is in Hz.
    """
    
    # Find out how many frames there are per observation
    nFrames = get_frames_per_obs(filehandle)
    
    with FilePositionSaver(filehandle):
        # Find the lowest frequency channel
        freqs = []
        while len(freqs) < nFrames:
            cFrame = read_frame(filehandle)
            if frequency:
                freq = cFrame.channel_freqs[0]
            else:
                freq = cFrame.header.first_chan
                
            if freq not in freqs:
                freqs.append(freq)
                
    # Return the lowest frequency channel
    return min(freqs)
