"""
Python module to reading in data from TBX files.  TBX data are a complex
frequency-domain product that contains blocks of up to 16 channels from all
antennas in the array.  Each channel has a bandwidth of f\ :sub:`C` (25 kHz)
and there may be up to 224 different blocks of channels within a single
recording.  The stand ordering is based on the input into the digital system
rather than the stand number in the array.

This module defines the following classes for storing the TBX data found in a
file:

Frame
  object that contains all data associated with a particular TBX frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the TBX frame header object and
    * FramePayload   - the TBX frame data object.  
Combined, these two objects contain all of the information found in the 
original TBX frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.

.. versionadded:: 3.1.0
"""

import numpy as np

from lsl.common import adp as adp_common
from lsl.common import ndp as ndp_common
from lsl.reader.base import *
from lsl.reader._gofast import read_tbf, read_tbf_ci8
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.1'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'read_frame', 'read_frame_ci8',
           'FRAME_CHANNEL_COUNT', 'get_frame_size', 'get_frames_per_obs',
           'get_first_frame_count', 'get_channel_count', 'get_first_channel']

#: Number of frequency channels in a TBF packet
FRAME_CHANNEL_COUNT = 12


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a TBF 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    _header_attrs = ['frame_count', 'second_count', 'nstand', 'nchan', 'first_chan']
    
    def __init__(self, frame_count=None, second_count=None, first_chan=None, nstand=None, nchan=None):
        self.frame_count = frame_count
        self.second_count = second_count
        self.first_chan = first_chan
        self.nstand = nstand
        self.nchan = nchan
        FrameHeaderBase.__init__(self)
        
    @property
    def channel_freqs(self):
        """
        Return a numpy.float32 array for the center frequencies, in Hz, of
        each channel in the data.
        """
        
        fC = ndp_common.fC
        return (np.arange(FRAME_CHANNEL_COUNT, dtype=np.float32)+self.first_chan) * fC


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a TBF
    frame.  Both fields listed in the DP ICD version H are stored.
    """
    
    _payload_attrs = ['timetag']
    
    def __init__(self, timetag=None, fDomain=None):
        self.timetag = timetag
        FramePayloadBase.__init__(self, fDomain)
        
    @property
    def time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance.
        """
        
        return FrameTimestamp.from_dp_timetag(self.timetag)


class Frame(FrameBase):
    """
    Class that stores the information contained within a single TBF 
    frame.  It's properties are FrameHeader and FramePayload objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    
    @property
    def nstand(self):
        """
        Convenience wrapper for the Frame.FrameHeader.nstand property.
        """
        
        return self.header.nstand
        
    @property
    def nchan(self):
        """
        Convenience wrapper for the Frame.FrameHeader.nchan property.
        """
        
        return self.header.nchan
        
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
        newFrame = read_tbx(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell()
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
        
    return newFrame


def read_frame_ci8(filehandle, verbose=False):
    """
    Function to read in a single TBX frame (header+data) and store the 
    contents as a Frame object.
    
    .. note::
        This function differs from `read_frame` in that it returns a
        `lsl.reader.tbf.FramePayload` that contains a 4-D numpy.int8 array
        (channels by stands by polarizations by by real/complex) rather than a
        3-D numpy.complex64 array.
    
    .. versionadded:: 2.1.3
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_tbx_ci8(filehandle, Frame())
        newFrame.payload._data = newFrame.payload.data.view(CI8)
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
        
    return newFrame


def get_frame_size(filehandle):
    """
    Find out what the frame size in a file is at the current file location.
    Returns the frame size in bytes.
    """
    
    nPos = -1
    with FilePositionSaver(filehandle):
        for i in range(2500):
            try:
                cPos = filehandle.tell()
                cFrame = read_frame(filehandle)
                if not cFrame.is_tbf:
                    continue
                nPos = filehandle.tell()
                break
            except EOFError:
                break
            except SyncError:
                filehandle.seek(1, 1)
                
    if nPos < 0:
        raise RuntimeError("Unable to determine a frame size")
        
    return nPos - cPos


def get_frames_per_obs(filehandle):
    """
    Find out how many frames are present per time stamp by examining the 
    first 2500 TBF records.  Return the number of frames per observation.
    """
    
    with FilePositionSaver(filehandle):
        # Build up the list-of-lists that store the index of the first frequency
        # channel in each frame.
        channels = []
        for i in range(2500):
            try:
                cFrame = read_frame(filehandle)
                if not cFrame.is_tbf:
                    continue
            except EOFError:
                break
            except SyncError:
                filehandle.seek(0, 1)
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


def get_first_channel(filehandle, frequency=False, all_frames=False):
    """
    Find and return the lowest frequency channel in a TBF file.  If the 
    `frequency` keyword is True the returned value is in Hz.  If `all` is
    True then the lowest frequency in each unique TBF frame is returned as
    a list.
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
                
    freqs.sort()
    if all_frames:
        # Return all unique first frequency channels
        return freqs
    else:
        # Return the lowest frequency channel
        return freqs[0]
