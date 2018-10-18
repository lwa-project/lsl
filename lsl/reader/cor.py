# -*- coding: utf-8 -*-

"""
Python module to reading in data from COR files.  This module defines the 
following classes for storing the COR data found in a file:

Frame
  object that contains all data associated with a particular COR frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the COR frame header object and
    * FrameData   - the COR frame data object.  
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

import copy
import numpy

from lsl.common import adp as adp_common
from lsl.reader._gofast import NCHAN_COR
from lsl.reader._gofast import readCOR
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'read_frame', 'FRAME_SIZE', 'FRAME_CHANNEL_COUNT', 
           'get_frames_per_obs', 'get_channel_count', 'get_baseline_count']

FRAME_SIZE = 32 + NCHAN_COR*4*8
FRAME_CHANNEL_COUNT = NCHAN_COR


class FrameHeader(object):
    """
    Class that stores the information found in the header of a COR 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    def __init__(self, adp_id=None, frame_count=None, second_count=None, first_chan=None, gain=None):
        self.adp_id = adp_id
        self.frame_count = frame_count
        self.second_count = second_count
        self.first_chan = first_chan
        self.gain = gain
        
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


class FrameData(object):
    """
    Class that stores the information found in the data section of a COR
    frame.
    """
    
    def __init__(self, timetag=None, navg=None, stand0=None, stand1=None, vis=None):
        self.timetag = timetag
        self.navg = navg
        self.stand0 = stand0
        self.stand1 = stand1
        self.vis = vis
        
    @property
    def id(self):
        """
        Return a tuple of the two stands that contribute the this frame.
        """
        
        return (self.stand0, self.stand1)
        
    def get_time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch.
        """
        
        seconds = self.timetag / adp_common.fS
        
        return seconds
        
    @property
    def integration_time(self):
        """
        Return the integration time of the visibility in seconds.
        """
        
        return self.navg * adp_common.T2


class Frame(object):
    """
    Class that stores the information contained within a single COR 
    frame.  It's properties are FrameHeader and FrameData objects.
    """
    
    def __init__(self, header=None, data=None):
        if header is None:
            self.header = FrameHeader()
        else:
            self.header = header
            
        if data is None:
            self.data = FrameData()
        else:
            self.data = data
            
        self.valid = True
        
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
        
    def get_time(self):
        """
        Convenience wrapper for the Frame.FrameData.get_time function.
        """
        
        return self.data.get_time()
        
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FrameData.id property.
        """
        
        return self.data.id
        
    @property
    def integration_time(self):
        """
        Convenience wrapper for the Frame.FrameData.integration_time
        property.
        """
        
        return self.data.integration_time
        
    def __add__(self, y):
        """
        Add the data sections of two frames together or add a number 
        to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
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
        
        try:
            self.data.vis += y.data.vis
        except AttributeError:
            self.data.vis += numpy.complex64(y)
        return self
        
    def __mul__(self, y):
        """
        Multiple the data sections of two frames together or multiply 
        a number to every element in the data section.
        
        .. note::
            In the case where a frame is given the weights are
            ignored.
        """
        
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
        
        try:
            self.data.vis *= y.data.vis
        except AttributeError:
            self.data.vis *= numpy.complex64(y)
        return self
            
    def __eq__(self, y):
        """
        Check if the time tags of two frames are equal or if the time
        tag is equal to a particular value.
        """
        
        tX = self.data.timetag
        try:
            tY = y.data.timetag
        except AttributeError:
            tY = y
        
        if tX == tY:
            return True
        else:
            return False
            
    def __ne__(self, y):
        """
        Check if the time tags of two frames are not equal or if the time
        tag is not equal to a particular value.
        """
        
        tX = self.data.timetag
        try:
            tY = y.data.timetag
        except AttributeError:
            tY = y
        
        if tX != tY:
            return True
        else:
            return False
            
    def __gt__(self, y):
        """
        Check if the time tag of the first frame is greater than that of a
        second frame or if the time tag is greater than a particular value.
        """
        
        tX = self.data.timetag
        try:
            tY = y.data.timetag
        except AttributeError:
            tY = y
        
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
        
        tX = self.data.timetag
        try:
            tY = y.data.timetag
        except AttributeError:
            tY = y
        
        if tX >= tY:
            return True
        else:
            return False
            
    def __lt__(self, y):
        """
        Check if the time tag of the first frame is less than that of a
        second frame or if the time tag is greater than a particular value.
        """
        
        tX = self.data.timetag
        try:
            tY = y.data.timetag
        except AttributeError:
            tY = y
        
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
        
        tX = self.data.timetag
        try:
            tY = y.data.timetag
        except AttributeError:
            tY = y
        
        if tX <= tY:
            return True
        else:
            return False
            
    def __cmp__(self, y):
        """
        Compare two frames based on the time tags.  This is helpful for 
        sorting things.
        """
        
        tX = self.data.timetag
        tY = y.data.timetag
        if tY > tX:
            return -1
        elif tX > tY:
            return 1
        else:
            return 0


def read_frame(filehandle, Verbose=False):
    """
    Function to read in a single COR frame (header+data) and store the 
    contents as a Frame object.
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = readCOR(filehandle, Frame())
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
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
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
            
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
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
