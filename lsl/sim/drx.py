# -*- coding: utf-8 -*-

"""
Python module for creating creating, validating, and writing simulated 
DRX frames to a file.
"""

import numpy

from lsl.common.dp import fS
from lsl.reader import drx
from lsl.sim.errors import invalidBeam, invalidTune, invalidPol, invalidDataSize, invalidDataType

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['SimFrame', 'frame2frame', '__version__', '__revision__', '__all__']


def frame2frame(drx_frame):
    """
    Convert a :class:`lsl.reader.drx.Frame` object to a raw DP DRX frame.
    """

    # The raw frame
    rawFrame = numpy.zeros(drx.FrameSize, dtype=numpy.uint8)

    # Part 1: The header
    ## Sync. words (0xDEC0DE5C)
    rawFrame[0] = 0xDE  # 222
    rawFrame[1] = 0xC0  # 192
    rawFrame[2] = 0xDE  # 222
    rawFrame[3] = 0x5C  #  92
    ## DRX ID
    rawFrame[4] = drx_frame.header.drxID
    ## Frame count
    rawFrame[5] = (drx_frame.header.frameCount>>16) & 255
    rawFrame[6] = (drx_frame.header.frameCount>>8) & 255
    rawFrame[7] = drx_frame.header.frameCount & 255
    ## Seconds count
    rawFrame[8] = (drx_frame.header.secondsCount>>24) & 255
    rawFrame[9] = (drx_frame.header.secondsCount>>16) & 255
    rawFrame[10] = (drx_frame.header.secondsCount>>8) & 255
    rawFrame[11] = drx_frame.header.secondsCount & 255
    ## Decimation
    rawFrame[12] = (drx_frame.header.decimation>>8) & 255
    rawFrame[13] = drx_frame.header.decimation & 255
    ## Time offset
    rawFrame[14] = (drx_frame.header.timeOffset>>8) & 255
    rawFrame[15] = drx_frame.header.timeOffset & 255

    # Part 2: The data
    ## Time tag
    rawFrame[16] = (drx_frame.data.timeTag>>56) & 255
    rawFrame[17] = (drx_frame.data.timeTag>>48) & 255
    rawFrame[18] = (drx_frame.data.timeTag>>40) & 255
    rawFrame[19] = (drx_frame.data.timeTag>>32) & 255
    rawFrame[20] = (drx_frame.data.timeTag>>24) & 255
    rawFrame[21] = (drx_frame.data.timeTag>>16) & 255
    rawFrame[22] = (drx_frame.data.timeTag>>8) & 255
    rawFrame[23] = drx_frame.data.timeTag & 255
    ## Flags
    rawFrame[24] = (drx_frame.data.flags>>56) & 255
    rawFrame[25] = (drx_frame.data.flags>>48) & 255
    rawFrame[26] = (drx_frame.data.flags>>40) & 255
    rawFrame[27] = (drx_frame.data.flags>>32) & 255
    rawFrame[28] = (drx_frame.data.flags>>24) & 255
    rawFrame[29] = (drx_frame.data.flags>>16) & 255
    rawFrame[30] = (drx_frame.data.flags>>8) & 255
    rawFrame[31] = drx_frame.data.flags & 255
    ## Data
    i = drx_frame.data.iq.real
    q = drx_frame.data.iq.imag
    ### Round, clip, and convert to unsigned integers
    i = i.round()
    i = i.clip(-8, 7)
    i = i.astype(numpy.int8)
    i += ((i & 8) << 1)
    q = q.round()
    q = q.clip(-8, 7)
    q = q.astype(numpy.int8)
    q += ((q & 8) << 1)
    
    rawFrame[32:] = (((i &  0xF) << 4) | (q & 0xF))
    
    return rawFrame


class SimFrame(drx.Frame):
    """
    drx.SimFrame extends the :class:`lsl.reader.drx.Frame` object to yield a method 
    for easily creating DP ICD-compliant raw DRX frames.  Frames created with
    this method can be written to a file via the methods write_raw_frame() function.
    """

    def __init__(self, beam=None, tune=None, pol=None, filter_code=None, time_offset=None, frame_count=None, obs_time=None, flags=None, iq=None):
        """
        Given a list of parameters, build a drx.SimFrame object.  The parameters
        needed are:
         * beam id (>0 & <5)
         * tunning (1 or 2)
         * polarization (0 for x, or 1 for y)
         * which filter code the data corresponds to (>0 & <8)
         * what time offset in units of f_S to use
         * which frame number to create
         * observation time in samples at fS since the epoch
         * what flags are set on the data
         * 1-D numpy array representing the frame I/Q (complex) data
        
        Not all of these parameters are needed at initialization of the object and
        the values can be added later.

        .. versionchanged: 0.3.4
            obs_time now in samples at fS, not seconds
        """
        
        self.beam = beam
        self.tune = tune
        self.pol = pol
        self.filter_code = filter_code
        self.time_offset = time_offset
        self.frame_count = frame_count
        self.secondsCount = 0
        self.obs_time = obs_time
        self.flags = flags
        self.iq = iq
        super(SimFrame, self).__init__()
        
    def __update(self):
        """
        Private function to use the object's parameter values to build up 
        a drx.Frame-like object.
        """
        
        self.header.frameCount = 0*self.frame_count
        self.header.secondsCount = 0*long(self.obs_time / fS)
        self.header.decimation = int(fS / drx.filterCodes[self.filter_code])
        self.header.timeOffset = self.time_offset
        self.header.drxID = (self.beam & 7) | ((self.tune & 7) << 3) | ((self.pol & 1) << 7)
        
        self.data.timeTag = self.obs_time
        self.data.flags = self.flags
        self.data.iq = self.iq
        
    def load_frame(self, drx_frame):
        """
        Populate the a drx.SimFrame object with a pre-made frame.
        """
        
        self.header = drx_frame.header
        self.data = drx_frame.data

        inverseCodes = {}
        for code in drx.filterCodes:
            rate = drx.filterCodes[code]
            inverseCodes[int(rate)] = code
        
        # Back-fill the class' fields to make sure the object is consistent
        ## Header
        self.beam = self.header.parseID()[0]
        self.tune = self.header.parseID()[1]
        self.pol = self.header.parseID()[2]
        self.frame_count = self.header.frameCount
        self.secondsCount = self.header.secondsCount
        self.filter_code = inverseCodes[int(fS / self.header.decimation)]
        self.time_offset = self.header.timeOffset
        ## Data
        self.obs_time = self.data.timeTag
        self.flags = self.data.flags
        self.iq = self.data.iq
    
    def is_valid(self, raiseErrors=False):
        """
        Check if simulated DRX frame is valid or not.  Valid frames return 
        True and invalid frames False.  If the 'raiseErrors' keyword is set, 
        is_valid raises an error when a problem is encountered.
        """

        # Make sure we have the latest values
        self.__update()

        # Is the time offset reasonable?
        if self.header.timeOffset >= fS:
            return False

        beam, tune, pol = self.parseID()
        # Is the beam number reasonable?
        if beam not in [1, 2, 3, 4]:
            if raiseErrors:
                raise invalidBeam()
            return False

        # Is the tunning number reasonable?
        if tune not in [1, 2]:
            if raiseErrors:
                raise invalidTune()
            return False

        # Is the polarization reasonable?
        if pol not in [0, 1]:
            if raiseErrors:
                raise invalidPol()
            return False

        # Is there data loaded into frame.data.iq?
        if self.data.iq is None:
            if raiseErrors:
                raise invalidDataSize()
            return False

        # Does the data length make sense?
        if self.data.iq.shape[0] != 4096:
            if raiseErrors:
                raise invalidDataSize()
            return False

        # Does the data type make sense?
        if self.data.iq.dtype.kind != 'c':
            if raiseErrors:
                raise invalidDataType()
            return False

        # If we made it this far, it's valid
        return True

    def create_raw_frame(self):
        """
        Re-express a simulated DRX frame as a numpy array of unsigned 8-bit 
        integers.  Returns a numpy array if the frame is valid.  If the frame 
        is not ICD-compliant, a errors.baseSimError-type error is raised.
        """

        # Make sure we have the latest values
        self.__update()

        self.is_valid(raiseErrors=True)
        return frame2frame(self)

    def write_raw_frame(self, fh):
        """
        Write a simulated DRX frame to a filehandle if the frame is valid.
        If the frame is not ICD-compliant, a errors.baseSimError-type error 
        is raised.
        """

        rawFrame = self.create_raw_frame()
        rawFrame.tofile(fh)

    def __str__(self):
        if self.beam is None:
            return "Empty DRX SimFrame object"
        else:
            return "DRX SimFrame for beam %i, tunning %i, pol. %i @ time %i" % (self.beam, self.tune, self.pol, self.obs_time)
