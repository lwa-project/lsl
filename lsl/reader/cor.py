# -*- coding: utf-8 -*-

"""
Python module to reading in data from COR files.  This module defines the 
following classes for storing the COR data found in a file:

Frame
  object that contains all data associated with a particular COR frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the COR frame header object and
    * FramePayload   - the COR frame data object.  
Combined, these two objects contain all of the information found in the 
original COR frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.

.. versionchanged:: 1.2.1
    Updated for the switch over to 72 channels, complex64 data, and no
    data weights

.. versionadded:: 1.2.0
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import copy
import numpy

from lsl.common import adp as adp_common
from lsl.reader.base import *
from lsl.reader._gofast import NCHAN_COR
from lsl.reader._gofast import read_cor
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'read_frame', 'FRAME_SIZE', 'FRAME_CHANNEL_COUNT', 
           'get_frames_per_obs', 'get_channel_count', 'get_baseline_count']

FRAME_SIZE = 32 + NCHAN_COR*4*8
FRAME_CHANNEL_COUNT = NCHAN_COR


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a COR 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    _header_attrs = ['adp_id', 'frame_count', 'second_count', 'first_chan', 'gain']
    
    def __init__(self, adp_id=None, frame_count=None, second_count=None, first_chan=None, gain=None):
        self.adp_id = adp_id
        self.frame_count = frame_count
        self.second_count = second_count
        self.first_chan = first_chan
        self.gain = gain
        FrameHeaderBase.__init__(self)
        
    @property
    def is_cor(self):
        """
        Function to check if the data is really COR.  Returns True if the 
        data is COR, false otherwise.
        """
        
        if self.adp_id == 0x02:
            return True
        else:
            return False
            
    @property
    def channel_freqs(self):
        """
        Return a numpy.float32 array for the center frequencies, in Hz, of
        each channel in the data.
        """
        
        return (numpy.arange(NCHAN_COR, dtype=numpy.float32)+self.first_chan) * adp_common.fC


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a COR
    frame.
    """
    
    _payload_attrs = ['timetag', 'navg', 'stand0', 'stand1']
    
    def __init__(self, timetag=None, navg=None, stand0=None, stand1=None, vis=None):
        self.timetag = timetag
        self.navg = navg
        self.stand0 = stand0
        self.stand1 = stand1
        FramePayloadBase.__init__(self, vis)
        
    @property
    def id(self):
        """
        Return a tuple of the two stands that contribute the this frame.
        """
        
        return (self.stand0, self.stand1)
        
    @property
    def time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch as a two-
        element tuple.
        """
        
        seconds_i = self.timetag // int(adp_common.fS)
        seconds_f = (self.timetag % int(adp_common.fS)) / adp_common.fS
        
        return seconds_i, seconds_f
        
    @property
    def integration_time(self):
        """
        Return the integration time of the visibility in seconds.
        """
        
        return self.navg * adp_common.T2


class Frame(FrameBase):
    """
    Class that stores the information contained within a single COR 
    frame.  It's properties are FrameHeader and FramePayload objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    
    @property
    def is_cor(self):
        """
        Convenience wrapper for the Frame.FrameHeader.is_cor property.
        """
        
        return self.header.is_cor
        
    @property
    def channel_freqs(self):
        """
        Convenience wrapper for the Frame.FrameHeader.channel_freqs property.
        """
        
        return self.header.channel_freqs
        
    @property
    def gain(self):
        """
        Convenience wrapper for the Frame.FrameHeader.gain property.
        """

        return self.header.gain
        
    @property
    def time(self):
        """
        Convenience wrapper for the Frame.FramePayload.time property.
        """
        
        return self.payload.time
        
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FramePayload.id property.
        """
        
        return self.payload.id
        
    @property
    def integration_time(self):
        """
        Convenience wrapper for the Frame.FramePayload.integration_time
        property.
        """
        
        return self.payload.integration_time


def read_frame(filehandle, verbose=False):
    """
    Function to read in a single COR frame (header+data) and store the 
    contents as a Frame object.
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_cor(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
        
    return newFrame


def get_frames_per_obs(filehandle):
    """
    Find out how many frames are present per time stamp by examining the 
    first several COR records.  Return the number of frames per observation.
    """
    
    # Get the number of channels in the file
    nChan = get_channel_count(filehandle)
    nFrames = nChan / NCHAN_COR
    
    # Multiply by the number of baselines
    nFrames *= get_baseline_count(filehandle)
    
    # Return the number of channel/baseline pairs
    return nFrames


def get_channel_count(filehandle):
    """
    Find out the total number of channels that are present by examining 
    the first several COR records.  Return the number of channels found.
    """
    
    with FilePositionSaver(filehandle):
        # Build up the list-of-lists that store the index of the first frequency
        # channel in each frame.
        channels = []
        for i in range(64):
            try:
                cFrame = read_frame(filehandle)
            except:
                break
                
            chan = cFrame.header.first_chan
            if chan not in channels:
                channels.append( chan )
                
    # Return the number of channels
    return len(channels)*NCHAN_COR


def get_baseline_count(filehandle):
    """
    Find out the total number of baselines that are present by examining the 
    first several COR records.  Return the number of baselines found.
    """
    
    # This is fixed based on how ADP works
    nBaseline = 256*(256+1) / 2
    
    return nBaseline
