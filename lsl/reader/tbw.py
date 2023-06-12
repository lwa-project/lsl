"""
Python module to reading in data from both 12-bit and 4-bit TBW files.  TBW data
consist of a time series of real valued data sampled at f\ :sub:`S` (196 MHz)
from all antennas in the array.  The stand numbering is based on the input into
the digital system rather than the stand number in the array.  The data are
divided into packets that contain either 400 samples (12-bit data) or 1200
samples (4-bit data) per stand per polarization.

This module defines the following classes for storing the TBW data found in
a file:

Frame
  object that contains all data associated with a particular TBW frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the TBW frame header object and
    * FramePayload   - the TBW frame data object.  
Combined, these two objects contain all of the information found in the 
original TBW frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  read_frame 
is designed to work with both 4-bit and 12-bit observations.

For describing the format of data in the file, two function are provided:

get_data_bits
  read in the first frame of an open file handle and return whether or not 
  the data is 12 or 4-bit

get_frames_per_obs
  read in the first several frames to see how many stands are found in the 
  data.
.. note::
    This function is a little flaky on TBW data sets that have less 
    than a full complement or 12M (36M) samples.

.. versionchanged:: 0.4.0
    Switched over from pure Python readers to the new C-base Go Fast! readers.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
from lsl.reader.base import *
from lsl.reader._gofast import read_tbw
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.6'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'read_frame', 
           'FRAME_SIZE', 'get_data_bits', 'get_frames_per_obs']

#: TBW packet size (header + payload)
FRAME_SIZE = 1224


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a TBW 
    frame.  All three fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    _header_attrs = ['frame_count', 'second_count', 'tbw_id']
    
    def __init__(self, frame_count=None, second_count=None, tbw_id=None):
        self.frame_count = frame_count
        self.second_count = second_count
        self.tbw_id = tbw_id
        FrameHeaderBase.__init__(self)
        
    @property
    def is_tbw(self):
        """
        Function to check if the data is really TBW and not TBN by examining
        the TBW ID field.  Returns True if the data is TBW, false otherwise.
        """

        mode = (self.tbw_id>>15)&1
        if mode == 1:
            return True
        else:
            return False
    
    @property
    def id(self):
        """
        Function to parse the TBW ID field and return the stand number.
        """

        # Why &1023?  Well, from DP ICD revision H, it seems that the stand count 
        # only goes up 260.  So, channel numbers should range up to 520, which can
        # be represented as 10 bits or 1023.
        stand = self.tbw_id&1023

        return stand
    
    @property
    def data_bits(self):
        """
        Function to parse the TBW ID field and return the size of number of 
        bits that comprise the data.  12 is returned for 12-bit data, and 4 
        for 4-bit data.
        """

        bits = (self.tbw_id>>14)&1
        if bits == 0:
            dataBits = 12
        else:
            dataBits = 4

        return dataBits


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a TBW
    frame.  Both fields listed in the DP ICD version H are stored.
    """
    
    _payload_attrs = ['timetag']
    
    def __init__(self, timetag=None, xy=None):
        self.timetag = timetag
        FramePayloadBase.__init__(self, xy)
        
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
    Class that stores the information contained within a single TBW 
    frame.  It's properties are FrameHeader and FramePayload objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    
    @property
    def is_tbw(self):
        """
        Convenience wrapper for the Frame.FrameHeader.is_tbw property.
        """
        
        return self.header.is_tbw
        
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FrameHeader.id 
        property.
        """
        
        return self.header.id

    @property
    def data_bits(self):
        """
        Convenience wrapper for the Frame.FrameHeader.data_bits 
        property.
        """
        
        return self.header.data_bits
        
    @property
    def time(self):
        """
        Convenience wrapper for the Frame.FramePayload.time property.
        """
        
        return self.payload.time


def read_frame(filehandle, verbose=False):
    """
    Function to read in a single TBW frame (header+data) and store the 
    contents as a Frame object.  This function wraps readerHeader and 
    readData[(12)|4].
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_tbw(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError

    return newFrame


def get_data_bits(filehandle):
    """
    Find out the number of data bits used in the file be reading in the 
    first frame.
    """

    with FilePositionSaver(filehandle):
        # Read a frame
        cFrame = read_frame(filehandle)
        
    # Get the number of bits used to represent the data
    dataBits = cFrame.data_bits
    
    return dataBits


def get_frames_per_obs(filehandle):
    """
    Find out how many frames are present per observation by examining 
    the first frames for what would be 260 stands.  This is done by reading
    two frames and then skipping the next 30,000.
    
    .. note::
        Post-IOC it is probably simpler to adopt a value of the number of 
        frames per observation of 260 rather than try to find it from the
        file.
    """
    
    with FilePositionSaver(filehandle):
        idCodes = []
        for i in range(260):
            currentPosition = filehandle.tell()
            try:
                cFrame1 = read_frame(filehandle)
                cFrame2 = read_frame(filehandle)
            except EOFError:
                break
            except SyncError:
                continue
                
            cID = cFrame1.id
            if cID not in idCodes:
                idCodes.append(cID)
            cID = cFrame2.id
            if cID not in idCodes:
                idCodes.append(cID)
                
            # Junk 30,000 frames since that is how many frames there are per stand
            filehandle.seek(currentPosition+30000*FRAME_SIZE)
            
    # Get the length of the stand list and return
    return len(idCodes)
