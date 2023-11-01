"""
Python module for creating creating, validating, and writing simulated 
TBN frames to a file.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import numpy

from lsl.common.dp import fS
from lsl.reader import tbn
from lsl.reader.base import CI8

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
__all__ = ['SimFrame', 'frame_to_frame']


def frame_to_frame(tbn_frame):
    """
    Convert a :class:`lsl.reader.tbn.Frame` object to a raw DP TBN frame.
    """

    # The raw frame
    rawFrame = numpy.zeros(tbn.FRAME_SIZE, dtype=numpy.uint8)

    # Part 1: The header
    ## Sync. words (0xDEC0DE5C)
    rawFrame[0] = 0xDE  # 222
    rawFrame[1] = 0xC0  # 192
    rawFrame[2] = 0xDE  # 222
    rawFrame[3] = 0x5C  #  92
    ## Frame count
    rawFrame[5] = (tbn_frame.header.frame_count>>16) & 255
    rawFrame[6] = (tbn_frame.header.frame_count>>8) & 255
    rawFrame[7] = tbn_frame.header.frame_count & 255
    ## Tuning word
    rawFrame[8] = (tbn_frame.header.tuning_word>>24) & 255
    rawFrame[9] = (tbn_frame.header.tuning_word>>16) & 255
    rawFrame[10] = (tbn_frame.header.tuning_word>>8) & 255
    rawFrame[11] = tbn_frame.header.tuning_word & 255
    ## TBN ID
    rawFrame[12] = (tbn_frame.header.tbn_id>>8) & 255
    rawFrame[13] = tbn_frame.header.tbn_id & 255
    ## Gain
    rawFrame[14] = (tbn_frame.header.gain>>8) & 255
    rawFrame[15] = tbn_frame.header.gain & 255
    
    # Part 2: The data
    ## Time tag
    rawFrame[16] = (tbn_frame.payload.timetag>>56) & 255
    rawFrame[17] = (tbn_frame.payload.timetag>>48) & 255
    rawFrame[18] = (tbn_frame.payload.timetag>>40) & 255
    rawFrame[19] = (tbn_frame.payload.timetag>>32) & 255
    rawFrame[20] = (tbn_frame.payload.timetag>>24) & 255
    rawFrame[21] = (tbn_frame.payload.timetag>>16) & 255
    rawFrame[22] = (tbn_frame.payload.timetag>>8) & 255
    rawFrame[23] = tbn_frame.payload.timetag & 255
    ## Data
    if tbn_frame.payload.data.dtype == CI8:
        iq = tbn_frame.payload.data.view(numpy.int8).ravel().copy()
    else:
        iq = tbn_frame.payload.data
        iq.real
        ### Round and convert to unsigned integers
        iq = numpy.round(iq)
        if iq.dtype == numpy.complex128:
            iq = iq.astype(numpy.complex64)
        iq = iq.view(numpy.float32)
        iq = iq.astype(numpy.int8)
        
    rawFrame[24:] = iq
    
    return rawFrame


class SimFrame(tbn.Frame):
    """
    tbn.SimFrame extends the :class:`lsl.reader.tbn.Frame` object to yield a method 
    for easily creating DP ICD-compliant raw TBN frames.  Frames created with
    this method can be written to a file via the methods write_raw_frame() function.
    """

    def __init__(self, stand=None, pol=None, central_freq=None, gain=None, frame_count=None, obs_time=None, data=None):
        """
        Given a list of parameters, build a tbn.SimFrame object.  The parameters
        needed are:
         * stand id (>0 & <259)
         * polarization (0 for x, or 1 for y)
         * central frequency of tuning in (Hz)
         * TBN gain
         * which frame number to create
         * observation time in samples at fS since the epoch
         * 1-D numpy array representing the frame I/Q (complex) data
        
        Not all of these parameters are needed at initialization of the object and
        the values can be added later.

        .. versionchanged:: 0.3.4
            obs_time now in samples at fS, not seconds

        .. versionchanged:: 0.5.0
            Added support for ECR 11 TBN headers
        """
        
        tbn.Frame.__init__(self)
        self.stand = stand
        self.pol = pol
        self.freq = central_freq
        self.frame_count = frame_count
        self.gain = gain
        self.obs_time = obs_time
        self.data = data
        
    def _update(self):
        """
        Private function to use the object's parameter values to build up 
        a tbn.Frame-like object.
        """
        
        self.header.frame_count = self.frame_count
        self.header.tuning_word = int( round(self.freq/fS*2**32) )
        self.header.tbn_id = 2*(self.stand-1) + self.pol + 1
        self.header.gain = self.gain
        
        self.payload.timetag = self.obs_time
        self.payload._data = self.data
    
    def load_frame(self, tbn_frame):
        """
        Populate the a tbn.SimFrame object with a pre-made frame.
        """
        
        self.header = tbn_frame.header
        self.payload = tbn_frame.payload
        
        # Back-fill the class' fields to make sure the object is consistent
        ## Header
        self.stand = self.header.id[0]
        self.pol = self.header.id[1]
        self.freq = self.header.central_freq
        self.gain = self.header.gain
        self.frame_count = self.header.frame_count
        ## Data
        self.obs_time = self.payload.timetag
        self.data = self.payload.data
    
    def is_valid(self, raise_errors=False):
        """
        Check if simulated TBN frame is valid or not.  Valid frames return 
        True and invalid frames False.  If the 'raise_errors' keyword is set, 
        is_valid raises an error when a problem is encountered.
        """

        # Make sure we have the latest values
        self._update()

        stand, pol = self.id
        # Is the stand number reasonable?
        if stand == 0 or stand > 260:
            if raise_errors:
                raise ValueError("Invalid stand: '%s'" % stand)
            return False

        # Is the polarization reasonable?
        if pol not in [0, 1]:
            if raise_errors:
                raise ValueError("Invalid polarization: '%s'" % tune)
            return False

        # Is there data loaded into frame.payload.data?
        if self.payload.data is None:
            if raise_errors:
                raise ValueError("Invalid data payload: '%s'" % self.payload.data)
            return False

        # Does the data length make sense?
        if self.payload.data.shape[0] != 512:
            if raise_errors:
                raise ValueError("Invalid data length: %i", self.payload.data.shape[0])
            return False

        # Does the data type make sense?
        if self.payload.data.dtype != CI8 and self.payload.data.dtype.kind != 'c':
            if raise_errors:
                raise ValueError("Invalid data type: '%s'" % self.payload.data.dtype.kind)
            return False
            
        # If we made it this far, it's valid
        return True

    def create_raw_frame(self):
        """
        Re-express a simulated TBN frame as a numpy array of unsigned 8-bit 
        integers.  Returns a numpy array if the frame  is valid.  If the frame 
        is not ICD-compliant, a errors.baseSimError-type error is raised.
        """

        # Make sure we have the latest values
        self._update()

        self.is_valid(raise_errors=True)
        return frame_to_frame(self)

    def write_raw_frame(self, fh):
        """
        Write a simulated TBN frame to a filehandle if the frame is valid.
        If the frame is not ICD-compliant, a errors.baseSimError-type error 
        is raised.
        """

        rawFrame = self.create_raw_frame()
        rawFrame.tofile(fh)

    def __str__(self):
        if self.stand is None:
            return "Empty TBN SimFrame object"
        else:
            return "TBN SimFrame for stand %i, pol. %i @ time %i" % (self.stand, self.pol, self.obs_time)
