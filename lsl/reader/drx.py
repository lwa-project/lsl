"""
Python module to read in DRX data.  DRX data consist of a time series of complex
data a variable sample rate of up to 19.6 MHz from the beamformer.  The data are divided into
packets of 4096 samples per beam per tuning per polarization.

This module defines the following classes for storing the DRX data found in a
file:

Frame
  object that contains all data associated with a particular DRX frame.  
  The primary constituents of each frame are:
    * FrameHeader - the DRX frame header object and
    * FramePayload   - the DRX frame data object.
Combined, these two objects contain all of the information found in the 
original DRX frame.

ObservingBlock
object that stores a collection of Frames for all beams/tunings/
polarizations for a particular time.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  The readBlock
function reads in a (user-defined) number of DRX frames and returns a 
ObservingBlock object.

For describing the format of data in the file, two function are provided:

get_beam_count
  read in the first few frames of an open file handle and return how many 
  beams are present in the file.

get_frames_per_obs
  read in the first several frames to see how many frames (tunings/polarizations)
  are associated with each beam.

..versionchanged:: 2.1.3
    Added a new read_frame_ci8 function that returns 8-bit+8-bit complex integers
    instead of numpy.complex64

..versionchanged:: 1.2.0
    Dropped support for ObservingBlock since the lsl.reader.buffer modules does
    a better job.

.. versionchanged:: 0.4.0
    Switched over from pure Python readers to the new C-base Go Fast! readers.
"""

from lsl.common import dp as dp_common
from lsl.reader.base import *
from lsl.reader._gofast import read_drx, read_drx_ci8
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.9'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'read_frame', 'read_frame_ci8',
           'get_sample_rate', 'get_beam_count', 'get_frames_per_obs', 'FRAME_SIZE',
           'FILTER_CODES']

#: DRX packet size (header + payload)
FRAME_SIZE = 4128

#: List of filter codes and their corresponding sample rates in Hz
FILTER_CODES = {1: 250000, 2: 500000, 3: 1000000, 4: 2000000, 5: 4900000, 6: 9800000, 7: 19600000}


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a DRX 
    frame.  All six fields listed in the DP ICD version H are stored as 
    well as the original binary header data.
    """
    
    _header_attrs = ['frame_count', 'drx_id', 'second_count', 'decimation', 'time_offset']
    
    def __init__(self, frame_count=None, drx_id=None, second_count=None, decimation=None, time_offset=None):
        self.frame_count = frame_count
        self.drx_id = drx_id
        self.second_count = second_count
        self.decimation = decimation
        self.time_offset = time_offset
        FrameHeaderBase.__init__(self)
        
    @property
    def is_drx(self):
        """
        Function to check if the data is really DRX and not DRX8 by examining
        the drx_id field.  Returns True if the data is DRX, false otherwise.
        """
        
        return (self.drx_id>>6)&1 == 0
        
    @property
    def id(self):
        """
        Parse the DRX ID into a tuple containing the beam (1 through
        4), tunning (1 and 2), and polarization (0 and 1).
        """
        
        beam = self.drx_id&7
        tune = (self.drx_id>>3)&7
        pol  = (self.drx_id>>7)&1

        return (beam, tune, pol)
    
    @property
    def sample_rate(self):
        """
        Return the sample rate of the data in samples/second.
        """
        
        sample_rate = dp_common.fS / self.decimation
        return sample_rate
        
    @property
    def filter_code(self):
        """
        Function to convert the sample rate in Hz to a filter code.
        """
        
        sampleCodes = {}
        for key,value in FILTER_CODES.items():
            sampleCodes[value] = key
            
        return sampleCodes[self.sample_rate]


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a DRX
    frame.  All three fields listed in the DP ICD version H are stored.
    """
    
    _payload_attrs = ['timetag', 'tuning_word', 'flags']
    
    def __init__(self, timetag=None, tuning_word=None, flags=None, iq=None):
        self.gain = None
        self.timetag = timetag
        self.tuning_word = tuning_word
        self.flags = flags
        FramePayloadBase.__init__(self, iq)
        
    @property
    def central_freq(self):
        """
        Function to set the central frequency of the DRX data in Hz.
        """

        return dp_common.word_to_freq(self.tuning_word)


class Frame(FrameBase):
    """
    Class that stores the information contained within a single DRX 
    frame.  It's properties are FrameHeader and FramePayload objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    gain = None

    @property
    def is_drx(self):
        """
        Convenience wrapper for the Frame.FrameHeader.is_drx property.
        """
        
        return self.header.is_drx
        
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FrameHeader.id 
        property.
        """
        
        return self.header.id

    @property
    def sample_rate(self):
        """
        Convenience wrapper for the Frame.FrameHeader.sample_rate 
        property.
        """
        
        return self.header.sample_rate
        
    @property
    def filter_code(self):
        """
        Convenience wrapper for the Frame.FrameHeader.filter_code property.
        """

        return self.header.filter_code
        
    @property
    def time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance.
        """
        
        return FrameTimestamp.from_dp_timetag(self.payload.timetag, offset=self.header.time_offset)
        
    @property
    def central_freq(self):
        """
        Convenience wrapper for the Frame.FramePayload.central_freq property.
        """

        return self.payload.central_freq


def read_frame(filehandle, gain=None, verbose=False):
    """
    Function to read in a single DRX frame (header+data) and store the 
    contents as a Frame object.  This function wraps readerHeader and 
    readData.
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_drx(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
    
    if gain is not None:
        newFrame.gain = gain
        
    return newFrame


def read_frame_ci8(filehandle, gain=None, verbose=False):
    """
    Function to read in a single DRX frame (header+data) and store the 
    contents as a Frame object.  This function wraps readerHeader and 
    readData.
    
    .. note::
        This function differs from `read_frame` in that it returns a
        `lsl.reader.drx.FramePayload` that contains a 2-D numpy.int8 array
        (samples by real/complex) rather than a 1-D numpy.complex64 array.
    
    .. versionadded:: 2.1.3
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_drx_ci8(filehandle, Frame())
        newFrame.payload._data = newFrame.payload.data.view(CI8)
    except gSyncError:
        mark = filehandle.tell() - FRAME_SIZE
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
    
    if gain is not None:
        newFrame.gain = gain
        
    return newFrame


def get_sample_rate(filehandle, nframes=None, filter_code=False):
    """
    Find out what the sampling rate/filter code is from a single observations.  
    By default, the rate in Hz is returned.  However, the corresponding filter 
    code can be returned instead by setting the FilterCode keyword to true.
    
    This function is included to make easier to write code for TBN analysis and 
    modify it for DRX data.
    """
    
    with FilePositionSaver(filehandle):
        # Read in one frame
        newFrame = read_frame(filehandle)
        
    if not filter_code:
        return newFrame.sample_rate
    else:
        return newFrame.filter_code


def get_beam_count(filehandle):
    """
    Find out how many beams are present by examining the first 32 DRX
    records.  Return the number of beams found.
    """
    
    with FilePositionSaver(filehandle):
        # Build up the list-of-lists that store ID codes and loop through 32
        # frames.  In each case, parse pull the DRX ID, extract the beam number, 
        # and append the DRX ID to the relevant beam array if it is not already 
        # there.
        beams = []
        for i in range(32):
            cFrame = read_frame(filehandle)
            
            cID = cFrame.header.drx_id
            beam = cID&7
            if beam not in beams:
                beams.append(beam)
                
    # Return the number of beams found
    return len(beams)


def get_frames_per_obs(filehandle):
    """
    Find out how many frames are present per beam by examining the first 
    32 DRX records.  Return the number of frames per observations as a four-
    element tuple, one for each beam.
    """
    
    with FilePositionSaver(filehandle):
        # Build up the list-of-lists that store ID codes and loop through 32
        # frames.  In each case, parse pull the DRX ID, extract the beam number, 
        # and append the DRX ID to the relevant beam array if it is not already 
        # there.
        idCodes = [[], [], [], []]
        for i in range(32):
            cFrame = read_frame(filehandle)
            
            cID = cFrame.header.drx_id
            beam = cID&7
            if cID not in idCodes[beam-1]:
                idCodes[beam-1].append(cID)
                
    # Get the length of each beam list and return them as a tuple
    return (len(idCodes[0]), len(idCodes[1]), len(idCodes[2]), len(idCodes[3]))
