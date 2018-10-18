# -*- coding: utf-8 -*-

# Python3 compatibility
from __future__ import division

"""
Python module for reading data in from TBN files.This module defines the 
following classes for storing the TBN data found in a file:

Frame
  object that contains all data associated with a particular TBN frame.  
  The primary constituents of each frame are:
    * FrameHeader - the TBN frame header object and
    * FrameData   - the TBN frame data object.
Combined, these two objects contain all of the information found in the 
original TBN frame.

ObservingBlock
object that stores a collection of Frames for all stands/polarizations for 
a particular time.

In addition to storing the data available in the frame, the Frame object also
has attributes for holding information about the gain, central frequency, and
filter code used for the observations.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  The readBlock
function reads in a (user-defined) number of TBN frames and returns a 
ObservingBlock object.

For describing the format of data in the file, two function are provided:

get_sample_rate
  read in the few frame of an open file handle and return the sampling rate 
  of the data

get_frames_per_obs
  read in the first several frames to see how many stands are found in the data.

..versionchanged:: 1.2.0
    Dropped support for ObservingBlock since the lsl.reader.buffer modules does
    a better job.

.. versionchanged:: 0.4.0
    Switched over from pure Python readers to the new C-base Go Fast! readers.

.. versionchanged:: 0.5.0
    Support for ECR 11 TBN header format change.
"""

import copy
import numpy

from lsl.common import dp as dp_common
from lsl.reader._gofast import readTBN
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError

__version__ = '0.8'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'read_frame', 
           'get_sample_rate', 'get_frames_per_obs', 'FRAME_SIZE', 'FILTER_CODES']

FRAME_SIZE = 1048

# List of filter codes and their corresponding sample rates in Hz
FILTER_CODES = {1:   1000, 2:   3125, 3:    6250, 4:    12500, 5: 25000, 6: 50000, 7: 100000}


class FrameHeader(object):
    """
    Class that stores the information found in the header of a TBW 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    
    .. versionchanged:: 0.5.0
        Added various attributes to retrieve the central frequnecy 
        and gain that are part of the ECR 11 changes.
    """

    def __init__(self, frame_count=None, tuning_word=None, tbn_id=None, gain=None):
        self.frame_count = frame_count
        self.tuning_word = tuning_word
        self.tbn_id = tbn_id
        self.gain = gain
        self.sample_rate = None
        
    @property
    def is_tbn(self):
        """
        Function to check if the data is really TBN and not TBW by examining
        the TBN ID field.  Returns True if the data is TBN, false otherwise.
        """

        mode = (self.tbn_id>>15)&1
        if mode == 0:
            return True
        else:
            return False
            
    @property
    def id(self):
        """
        Function to parse the TBN ID field and return a tuple of the stand 
        number and polarization.
        """

        if self.tbn_id&1023 % 2 == 0:
            stand = (self.tbn_id&1023) // 2
            pol = 1
        else:
            stand = (self.tbn_id&1023) // 2 + 1
            pol = 0
            
        return (stand, pol)
    
    @property
    def central_freq(self):
        """
        Convert the tuning word to a frequency in Hz.
        """

        return dp_common.fS * self.tuning_word / 2**32
        
    @property
    def filter_code(self):
        """
        Function to convert the sample rate in Hz to a filter code.
        """
        
        if self.sample_rate is None:
            return None
        else:
            sampleCodes = {}
            for key in FILTER_CODES:
                value = FILTER_CODES[key]
                sampleCodes[value] = key

            return sampleCodes[self.sample_rate]


class FrameData(object):
    """
    Class that stores the information found in the data section of a TBN
    frame.  Both fields listed in the DP ICD version H are stored.
    
    .. versionchanged:: 0.5.0
        Removed various attributes related to storing a central frequnecy 
        and gain that aren't needed with ECR 11.
    """

    def __init__(self, timetag=None, iq=None):
        self.timetag = timetag
        self.iq = iq

    def get_time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch.
        """
        
        return self.timetag / dp_common.fS


class Frame(object):
    """
    Class that stores the information contained within a single TBN 
    frame.  It's properties are FrameHeader and FrameData objects.
    
    .. versionchanged:: 0.5.0
        Removed various attributes related to storing a central frequnecy 
        and gain that aren't needed with ECR 11.
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
    def is_tbn(self):
        """
        Convenience wrapper for the Frame.FrameHeader.is_tbn property.
        """
        
        return self.header.is_tbn
        
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FrameHeader.id property.
        """
        
        return self.header.id
        
    def get_time(self):
        """
        Convenience wrapper for the Frame.FrameData.get_time function.
        """
        
        return self.data.get_time()
        
    @property
    def sample_rate(self):
        """
        Convenience wrapper for the Frame.FrameData.sample_rate property.
        """

        return self.header.sample_rate
        
    @sample_rate.setter
    def sample_rate(self, value):
        """
        Convenience wrapper for setting the Frame.FrameData.sample_rate property.
        """

        self.header.sample_rate = value
        
    @property
    def filter_code(self):
        """
        Convenience wrapper for the Frame.FrameData.filter_code property.
        """

        return self.header.filter_code
    
    @property
    def central_freq(self):
        """
        Convenience wrapper for the Frame.FrameHeader.central_freq property.
        """

        return self.header.central_freq
        
    @property
    def gain(self):
        """
        Convenience wrapper for the Frame.FrameHeader.gain property.
        """

        return self.header.gain
        
    @gain.setter
    def gain(self, value):
        """
        Convenience wrapper to set the Frame.FrameHader.gain property.
        """
        
        self.header.gain = value
        
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
            self.data.iq += y.data.iq
        except AttributeError:
            self.data.iq += numpy.complex64(y)
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
            self.data.iq *= y.data.iq
        except AttributeError:
            self.data.iq *= numpy.complex64(y)
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


def read_frame(filehandle, sample_rate=None, Verbose=False):
    """
    Function to read in a single TBN frame (header+data) and store the 
    contents as a Frame object.
    """

    # New Go Fast! (TM) method
    try:
        newFrame = readTBN(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
    
    if sample_rate is not None:
        newFrame.setsample_rate(sample_rate)
        
    return newFrame


def get_sample_rate(filehandle, nFrames=None, FilterCode=False):
    """
    Find out what the sampling rate/filter code is from consecutive sets of 
    observations.  By default, the rate in Hz is returned.  However, the 
    corresponding filter code can be returned instead by setting the FilterCode
    keyword to True.
    """

    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()

    if nFrames is None:
        nFrames = 520
    nFrames = 4*nFrames

    # Build up the list-of-lists that store ID codes and loop through 2,080
    # frames.  In each case, parse pull the TBN ID, extract the stand 
    # number, and append the stand number to the relevant polarization array 
    # if it is not already there.
    frames = {}
    for i in range(nFrames):
        try:
            cFrame = read_frame(filehandle)
        except EOFError:
            break
        except SyncError:
            continue
        
        stand, pol = cFrame.id
        key = 2*stand + pol
        try:
            frames[key].append(cFrame)
        except:
            frames[key] = [cFrame,]
            
    # Return to the place in the file where we started
    filehandle.seek(fhStart)

    # Any key with complete data will work for this, so pick the first key with two
    # valid frames
    keyCount = 0
    frame1 = None
    frame2 = None
    frameKeys = list(frames.keys())
    while frame1 is None and frame2 is None:
        validKey = frameKeys[keyCount]
        
        try:
            frame1 = frames[validKey][0]
        except IndexError:
            frame1 = None
            
        try:
            frame2 = frames[validKey][1]
        except IndexError:
            frame2 = None

        keyCount = keyCount + 1

    # Now that we have two valid frames that follow one another in time, load in their
    # time tags and calculate the sampling rate.  Since the time tags are based off f_S
    # @ 196 MSPS, and each frame contains 512 samples, the sampling rate is:
    #  f_S / <difference in time tags per 512 samples>
    time1 = frame1.data.timetag
    time2 = frame2.data.timetag
    rate = dp_common.fS / (abs( time2 - time1 ) / 512)

    if not FilterCode:
        return rate
    else:
        sampleCodes = {}
        for key in FILTER_CODES:
            value = FILTER_CODES[key]
            sampleCodes[value] = key

        return sampleCodes[rate]


def get_frames_per_obs(filehandle):
    """
    Find out how many frames are present per observation by examining 
    the first 2,080 TBN frames.  Return the number of frames per observations 
    as a two-	element tuple, one for each polarization.
    
    So many TBN frames are read in order to try to compensate for the inter-
    leaving of the packets from the various DP1 boards during the recording.
    
    .. note::
        Post-IOC it is probably simpler to adopt a value of the number of 
        frames per observation of 520 rather than try to find it from the
        file.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()

    # Build up the list-of-lists that store ID codes and loop through 600
    # frames.  In each case, parse pull the TBN ID, extract the stand 
    # number, and append the stand number to the relevant polarization array 
    # if it is not already there.
    idCodes = [[], []]
    maxX = 0
    maxY = 0
    for i in range(4*520):
        try:
            cFrame = read_frame(filehandle)
        except EOFError:
            break
        except SyncError:
            continue
        
        cID, cPol = cFrame.header.id
        if cID not in idCodes[cPol]:
            idCodes[cPol].append(cID)
        
        # Also, look at the actual IDs to try to figure out how many of
        # each there are in case the frames are out of order.  Will this
        # help in all cases?  Probably not but we should at least try it
        if cPol == 0:
            if cID > maxX:
                maxX = cID
        else:
            if cID > maxY:
                maxY = cID
            
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Compare idCodes sizes with maxX and maxY.  Load maxX and maxY with 
    # the larger of the two values.
    if maxX < len(idCodes[0]):
        maxX = len(idCodes[0])
    if maxY < len(idCodes[1]):
        maxY = len(idCodes[1])
    
    # Get the length of each beam list and return them as a tuple
    return (maxX, maxY)
