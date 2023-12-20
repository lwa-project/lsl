"""
Module to write VDIF frames.  The implementation of this module is similar
to that of lsl.sim.tbw in that the primary element defined in this module is
a Frame object which as attribute functions that can create a numpy 
representation of the raw frame and write that raw frame to an open file-
handle.

.. seealso::
    :mod:`lsl.sim.tbw`
"""

import numpy as np
import ephem

from lsl.common import dp as dp_common
import lsl.astro as astro

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.1'
__all__ = ['Frame',]

VDIF_EPOCH = ephem.Date('2000/01/01 00:00:00.00')
UNIX_EPOCH = ephem.Date('1970/01/01 00:00:00.00')


class Frame(object):
    """
    Object to create and write a VDIF (VLBI Data Interchange Format) 
    frame (version 1, June 26, 2009).
    """
    
    def __init__(self, stand=0, time=0, bits=8, data=None, sample_rate=dp_common.fS):
        self.stand = stand
        self.time = time
        self.bits = bits
        self.data = np.squeeze(data)
        self.sample_rate = sample_rate
        if data.dtype.kind == 'c':
            self.dataReal = False
        else:
            self.dataReal = True
            
        self.seconds = 0
        self.frame = 0
        self.epoch = 0
        
        # Validate
        if self.bits not in (1, 2, 4, 8):
            raise ValueError(f"Unsupported output bitdepth: {self.bits}")
            
        # Convert the time from UNIX to epoch and make the data ready to 
        # be written to the disk.
        self._data_adjust()
        self._set_epoch()
        
    def _data_adjust(self):
        """
        Adjust the data to the number of bits required for the frame 
        and convert to unsigned 32-bit integers.  Also, interlace complex
        arrays.
        """
        
        if self.dataReal:
            interlaced = self.data.astype(np.int32)
        else:
            interlaced = np.zeros(2*len(self.data), dtype=np.int32)
            interlaced[0::2] = self.data.real.astype(np.int32)
            interlaced[1::2] = self.data.imag.astype(np.int32)
        
        biased = interlaced + 2**(self.bits-1)
        biased = biased.astype(np.uint32) & (2**self.bits-1)
        
        self.data = biased
        
    def _set_epoch(self):
        """
        Convert a UNIX timestap in seconds since 1970 to seconds since
        January 1st, 2000 using ephem.Data objects.  Also, calculate the 
        frame number since the beginning of the second.
        """
        
        try:
            seconds_i, seconds_f = self.time
        except TypeError:
            seconds_i = int(self.time)
            seconds_f = self.time - seconds_i
                
        # UNIX time to seconds since DJD  = 0
        curEpoch = float(UNIX_EPOCH)*astro.SECS_IN_DAY + seconds_i
        self.epoch = 0
        # Seconds since the VDIF epoch
        epochSeconds = curEpoch - float(VDIF_EPOCH)*astro.SECS_IN_DAY
        # Integer seconds
        self.seconds = int(epochSeconds)
        
        # Compute the frames since the beginning of the second
        frame = (epochSeconds - self.seconds + seconds_f) * (self.sample_rate/(len(self.data)//(2-self.dataReal)))
        self.frame = int(round(frame))

    def create_raw_frame(self):
        """
        Using the data and information stored in the object, create a
        numpy array on unsigned 1-byte integers that represents the frame.
        """
        
        # Find out how many samples can be packed into a work and build
        # the corresponding numpy array.
        samplesPerWord = int(32 // self.bits)
        raw = np.zeros(32 + 4*len(self.data)//samplesPerWord, dtype=np.uint8)

        # Valid data, standard (not legacy) 32-bit header, and seconds since 
        # the 01/01/2000 epoch.
        raw[3] = (0 << 7) | (0 << 6) | ((self.seconds >> 24) & 63)
        raw[2] = (self.seconds >> 16) & 255
        raw[1] = (self.seconds >> 8) & 255
        raw[0] = self.seconds & 255

        # Reference epoch (0 == 01/01/2000) and frame count
        raw[7] = self.epoch & 63
        raw[6] = (self.frame >> 16) & 255
        raw[5] = (self.frame >> 8) & 255
        raw[4] = self.frame & 255

        # VDIF version number, number of channels (just 1), and data frame 
        # length in units to 8-bytes (8 raw array elements)
        raw[11] = (1 << 6) | (0 & 31)
        raw[10] = ((len(raw) // 8) >> 16) & 255
        raw[9] = ((len(raw) // 8) >> 8) & 255
        raw[8] = (len(raw) // 8) & 255

        # Data type, bits per sample, thread ID, and station ID
        # NB:  The thread ID is fixed at `0'
        if self.dataReal:
            raw[15] = (0 << 7) | ((self.bits-1 & 31) << 2) | ((0 >> 8) & 3)
        else:
            raw[15] = (1 << 7) | ((self.bits-1 & 31) << 2) | ((0 >> 8) & 3)
        raw[14] = 0 & 255
        raw[13] = (self.stand >> 8) & 255
        raw[12] = self.stand & 255
        
        # User-defined words 4 through 7 (array entries 16 to 31)
        ## NICT format so that we can put in the sample rate
        sample_rate_kHz = int(self.sample_rate / 1000.0)
        raw[19] = 0x01
        raw[18] = (0 << 7) | ((sample_rate_kHz >> 16) & 0xFF)
        raw[17] = (sample_rate_kHz >> 8) & 0xFF
        raw[16] = sample_rate_kHz & 0xFF
        raw[20:24] = 0
        raw[24:28] = 0
        raw[28:32] = 0
        
        # Data values
        for f in range(32, len(raw), 4):
            i = (f - 32) // 4
            word = 0
            for p in range(samplesPerWord):
                word |= (self.data[i*samplesPerWord + p] << (self.bits*p))
            raw[f+3] = (word >> 24) & 255
            raw[f+2] = (word >> 16) & 255
            raw[f+1] = (word >> 8) & 255
            raw[f+0] = word & 255
            
        return raw

    def write_raw_frame(self, fh):
        """
        Create a numpy representation of the VDIF frame and then write 
        it to the specified file handle.
        """
        
        rawFrame = self.create_raw_frame()
        rawFrame.tofile(fh)
