# -*- coding: utf-8 -*-

"""
Python module to reading in data from TBF files.  This module defines the 
following classes for storing the TBF data found in a file:

Frame
  object that contains all data associated with a particular TBF frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the TBF frame header object and
    * FrameData   - the TBF frame data object.  
Combined, these two objects contain all of the information found in the 
original TBF frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.

.. versionadded:: 1.2.0
"""

import copy
import numpy

from lsl.common import adp as adp_common
from lsl.reader._gofast import readTBF
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'read_frame', 'FRAME_SIZE', 'FRAME_CHANNEL_COUNT',
           'get_frames_per_obs', 'getFirstFrameCount', 'get_channel_count', 'get_first_channel', 
           '__version__', '__revision__', '__all__']

FRAME_SIZE = 6168
FRAME_CHANNEL_COUNT = 12


class FrameHeader(object):
    """
    Class that stores the information found in the header of a TBF 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    def __init__(self, adp_id=None, frame_count=None, second_count=None, first_chan=None):
        self.adp_id = adp_id
        self.frame_count = frame_count
        self.second_count = second_count
        self.first_chan = first_chan
        
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
        
        return (numpy.arange(12, dtype=numpy.float32)+self.first_chan) * adp_common.fC


class FrameData(object):
    """
    Class that stores the information found in the data section of a TBF
    frame.  Both fields listed in the DP ICD version H are stored.
    """
    
    def __init__(self, timetag=None, fDomain=None):
        self.timetag = timetag
        self.fDomain = fDomain
        
    def get_time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch.
        """
        
        seconds = self.timetag / adp_common.fS
        
        return seconds


class Frame(object):
    """
    Class that stores the information contained within a single TBF 
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
        
    def get_time(self):
        """
        Convenience wrapper for the Frame.FrameData.get_time function.
        """
        
        return self.data.get_time()
        
    def __add__(self, y):
        """
        Add the data sections of two frames together or add a number 
        to every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame += y
        return newFrame
        
    def __iadd__(self, y):
        """
        In-place add the data sections of two frames together or add 
        a number to every element in the data section.
        """
        
        try:
            self.data.fDomain += y.data.fDomain
        except AttributeError:
            self.data.fDomain += numpy.complex64(y)
        return self
        
    def __mul__(self, y):
        """
        Multiple the data sections of two frames together or multiply 
        a number to every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame *= y
        return newFrame
            
    def __imul__(self, y):
        """
        In-place multiple the data sections of two frames together or 
        multiply a number to every element in the data section.
        """
        
        try:
            self.data.fDomain *= y.data.fDomain
        except AttributeError:
            self.data.fDomain *= numpy.complex64(y)
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
    Function to read in a single TBF frame (header+data) and store the 
    contents as a Frame object.
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = readTBF(filehandle, Frame())
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
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Build up the list-of-lists that store the index of the first frequency
    # channel in each frame.
    channels = []
    for i in range(1000):
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
    return len(channels)


def getFirstFrameCount(filehandle):
    """
    Find and return the lowest frame count encountered in a TBF file.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Find out how many frames there are per observation
    nFrames = get_frames_per_obs(filehandle)
    
    firstFrameCount = 2**64-1
    freqs = []
    while len(freqs) < nFrames:
        cFrame = read_frame(filehandle)
        freq = cFrame.header.first_chan
            
        if freq not in freqs:
            freqs.append(freq)
        if cFrame.header.frame_count < firstFrameCount:
            firstFrameCount = cFrame.header.frame_count
            
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
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
    nChannels = nFrames * 12
    
    # Return the number of channels
    return nChannels


def get_first_channel(filehandle, frequency=False):
    """
    Find and return the lowest frequency channel in a TBF file.  If the 
    `frequency` keyword is True the returned value is in Hz.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Find out how many frames there are per observation
    nFrames = get_frames_per_obs(filehandle)
    
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
            
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Return the lowest frequency channel
    return min(freqs)
