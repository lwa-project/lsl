"""
Python module to read in VDIF data.  This module defines the following 
classes for storing the VIDF data found in a file:

Frame
  object that contains all data associated with a particular DRX frame.  
  The primary constituents of each frame are:
    * FrameHeader - the VDIF frame header object and
    * FramePayload   - the VDIF frame data object.
Combined, these two objects contain all of the information found in the 
original VDIF frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  The readBlock
function reads in a (user-defined) number of VDIF frames and returns a 
ObservingBlock object.

For describing the format of data in the file, two function are provided:

get_thread_count
  read in the first few frames of an open file handle and return how many 
  threads are present in the file.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import copy
import warnings
from datetime import datetime

from lsl import astro
from lsl.common.mcs import datetime_to_mjdmpm
from lsl.reader.base import *
from lsl.reader._gofast import read_vdif
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.4'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'has_guppi_header', 'read_guppi_header', 
           'read_frame', 'get_frame_size', 'get_thread_count', 'get_frames_per_second', 
           'get_sample_rate']



def _crcc(data, length=48, mask=0o40003, cycle=16):
    """
    Compute a CRC checksum for the VLBA BCD time code stored in a Mark 5B 
    header.
    
    From:  mk5blib.c
    """
    
    state = 0
    for i in range(length):
        q = state & 1
        if ((data >> i) & 1) ^ q == 0:
            state &= -2
        else:
            state ^= mask
            state |= i
        state = (state >> 1) | (state & 1) << (cycle-1)
    return state


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a VDIF 
    frame.  Most fields in the VDIF version 1.1.1 header are stored.
    """
    
    _header_attrs = ['is_invalid', 'is_legacy', 'second_from_epoch',
                     'ref_epoch', 'frame_in_second',
                     'version', 'nchan', 'frame_length',
                     'is_complex', 'bits_per_sample', 'thread_id', 'station_id',
                     'extended_data_1', 'extended_data_2', 'extended_data_3', 'extended_data_4']
    
    def __init__(self, is_invalid=0, is_legacy=0, seconds_from_epoch=0, ref_epoch=0, frame_in_second=0, version=1, nchan=0, frame_length=0, is_complex='C', bits_per_sample=0, thread_id=0, station_id=0, extended_data_1=None, extended_data_2=None, extended_data_3=None, extended_data_4=None, sample_rate=0.0, central_freq=0.0):
        self.is_invalid = is_invalid
        self.is_legacy = is_legacy
        self.seconds_from_epoch = seconds_from_epoch
        
        self.ref_epoch = ref_epoch
        self.frame_in_second = frame_in_second
        
        self.version = version
        self.nchan = nchan
        self.frame_length = frame_length
        
        self.is_complex = is_complex
        self.bits_per_sample = bits_per_sample
        self.thread_id = thread_id
        self.station_id = station_id
        
        self.extended_data_1 = extended_data_1
        self.extended_data_2 = extended_data_2
        self.extended_data_3 = extended_data_3
        self.extended_data_4 = extended_data_4
        
        self.sample_rate = sample_rate
        self.central_freq = central_freq
        FrameHeaderBase.__init__(self)
        
    @property
    def time(self):
        """
        Function to convert the time tag to seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance.
        """
        
        # Get the reference epoch in the strange way that it is stored in VDIF 
        # and convert it to a MJD
        epochDT = datetime(2000+self.ref_epoch//2, (self.ref_epoch % 2)*6+1, 1, 0, 0, 0, 0)
        epochMJD, epochMPM = datetime_to_mjdmpm(epochDT)
        epochMJD = epochMJD + epochMPM/1000.0/86400.0
        
        # Get the frame MJD by adding the seconds_from_epoch value to the epoch
        frameMJD_i = epochMJD + self.seconds_from_epoch // 86400
        frameMJD_f = (self.seconds_from_epoch % 86400) / 86400.0
        frameMJD_s = 0.0
        
        if self.sample_rate == 0.0:
            # Try to get the sub-second time by parsing the extended user data
            try:
                ## Is there a sample rate to grab?
                eud = self.extended_user_data
                sample_rate = eud['sample_rate']
                sample_rate *= 1e6 if eud['sample_rate_units'] == 'MHz' else 1e3
                
                ## How many samples are in each frame?
                dataSize = self.frame_length*8 - 32 + 16*self.is_legacy		     # 8-byte chunks -> bytes - full header + legacy offset
                samplesPerWord = 32 // self.bits_per_sample					     # dimensionless
                nSamples = dataSize // 4 * samplesPerWord					     # bytes -> words -> data samples
                nSamples = nSamples / self.nchan / (2 if self.is_complex else 1) # data samples -> time samples
                
                ## What is the frame rate?
                frameRate = sample_rate // nSamples
                
                frameMJD_s += 1.0*self.frame_in_second/frameRate
            
            except KeyError:
                warnings.warn("Insufficient information to determine exact frame timestamp, time will be approximate", RuntimeWarning)
                
        else:
            # Use what we already have been told
            ## How many samples are in each frame?
            dataSize = self.frame_length*8 - 32 + 16*self.is_legacy		     # 8-byte chunks -> bytes - full header + legacy offset
            samplesPerWord = 32 // self.bits_per_sample				         # dimensionless
            nSamples = dataSize // 4 * samplesPerWord				         # bytes -> words -> samples
            nSamples = nSamples / self.nchan / (2 if self.is_complex else 1) # data samples -> time samples
        
            ## What is the frame rate?
            frameRate = self.sample_rate // nSamples
            
            frameMJD_s += 1.0*self.frame_in_second/frameRate
            
        # Convert from MJD to UNIX time
        if frameMJD_f > 1:
            frameMJD_i += 1
            frameMJD_f -= 1
            
        return FrameTimestamp.from_pulsar_mjd(frameMJD_i, frameMJD_f, frameMJD_s)
        
    @property
    def id(self):
        """
        Return a two-element tuple of the station ID and thread ID.
        
        .. note::
            The station ID is always returned as numeric.
        """
        return (self.station_id, self.thread_id)
        
    @property
    def extended_user_data(self):
        """
        Parse the extended user data if it was included with the reader.  
        The data is returned as a dictionary.
        """
        
        # Setup the output dictionary
        fields = {}
        
        # Is there anything to look at?
        if self.extended_data_1 is None or self.extended_data_2 is None or self.extended_data_3 is None or self.extended_data_4 is None:
            return fields
            
        # Extract the version
        edv = int((self.extended_data_1 >> 24) & 0xFF)
        
        # Parse accordingly
        if edv == 0x00:
            ## Empty
            pass
            
        elif edv == 0x01:
            ## NICT
            fields['sample_rate'] = int(self.extended_data_1 & (2**23-1))
            fields['sample_rate_units'] = 'MHz' if int((self.extended_data_1>>23) & 1) else 'kHz'
            fields['sync_pattern'] = self.extended_data_2
            fields['station_name'] = (self.extended_data_4 << 32) | self.extended_data_3
            
        elif edv == 0x02:
            ## ALMA
            fields['sync_word'] = int(self.extended_data_1 & 0xFFFF)
            fields['pic_status_word'] = self.extended_data_2
            fields['packet_serial_number'] = (self.extended_data_4 << 32) | self.extended_data_3
            
        elif edv == 0x03:
            ## NRAO
            fields['sample_rate'] = int(self.extended_data_1 & (2**23-1))
            fields['sample_rate_units'] = 'MHz' if int((self.extended_data_1 >> 23) & 1) else 'kHz'
            fields['sync_pattern'] = self.extended_data_2
            fields['tuning_word'] = self.extended_data_3
            fields['dbe_unit'] = int((self.extended_data_4 >> 24) & 0xF)
            fields['if_input_number'] = int((self.extended_data_4 >> 20) & 0xF)
            fields['subband'] = int((self.extended_data_4 >> 17) & 0x7)
            fields['electronic_sideband'] = 'USB' if (self.extended_data_4 >> 16) & 1 else 'LSB'
            mj = int((self.extended_data_4 >> 12) & 0xF)
            mn = int((self.extended_data_4 >>  8) & 0xF)
            pr = int(self.extended_data_4 & 0xFF)
            fields['version'] = '%i.%i-%02f' % (mj, mn, pr)
            
        elif edv == 0xAB:
            ## Haystack (which is really an embedded Mark 5B header)
            fields['sync_word'] = self.extended_data_1
            fields['years_from_2000'] = int((self.extended_data_2 >> 28) & 0xF)
            fields['user_specified_data'] = int((self.extended_data_2 >> 16) & 0xFFF)
            fields['data_from_internal_tvg'] = (self.extended_data_2 >> 15) & 1
            fields['frame_in_second'] = int(self.extended_data_2 & (2**14-1))
            j0 = int((self.extended_data_3 >> 28) & 0xF)
            j1 = int((self.extended_data_3 >> 24) & 0xF)
            j2 = int((self.extended_data_3 >> 20) & 0xF)
            s0 = int((self.extended_data_3 >> 16) & 0xF)
            s1 = int((self.extended_data_3 >> 12) & 0xF)
            s2 = int((self.extended_data_3 >>  8) & 0xF)
            s3 = int((self.extended_data_3 >>  4) & 0xF)
            s4 = int((self.extended_data_3 >>  0) & 0xF)
            f0 = int((self.extended_data_4 >> 28) & 0xF)
            f1 = int((self.extended_data_4 >> 24) & 0xF)
            f2 = int((self.extended_data_4 >> 20) & 0xF)
            f3 = int((self.extended_data_4 >> 16) & 0xF)
            f4 = int((self.extended_data_4 >>  8) & 0xF)
            crcf = int(self.extended_data_4 & 0xFFFF)
            crcf = (crcf & 0xF) << 8 | ((crcf >> 8) & 0xF)
            crcc = _crcc((self.extended_data_3 << 16) | ((self.extended_data_3 >> 16) & 0xFF), 48)
            fields['vlba_timecode'] = j0*1e7+j1*1e6+j2*1e5 + s0*1e4+s1*1e3+s2*1e2+s3*1e1+s4*1e0 + f0/1e1+f1/1e2+f2/1e3+f3/1e4+f4/1e5
            fields['vlba_timecode_value'] = True if crcc == crcf else False
            
        else:
            raise RuntimeError("Unknown extended user data version: %i" % edv)
            
        return fields


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a VDIF
    frame.
    
    .. note::
        Unlike the other readers in the :mod:`lsl.reader` module the
        data are stored as numpy.float32 values.
    """
    
    def __init__(self, data=None):
        FramePayloadBase.__init__(self, data)
        
    @property
    def data(self):
        return self._data


class Frame(FrameBase):
    """
    Class that stores the information contained within a single VDIF 
    frame.  It's properties are FrameHeader and FramePayload objects.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FrameHeader.id 
        property.
        """
        
        return self.header.id
        
    @property
    def extended_user_data(self):
        """
        Convenience wrapper for the Frame.FrameHeader.extended_user_data
        property.
        """
        
        return self.header.extended_user_data
        
    @property
    def time(self):
        """
        Convenience wrapper for the Frame.FrameHeader.time property.
        """
        
        return self.header.time
        
    @property
    def sample_rate(self):
        """
        Convenience wrapper for the Frame.FrameHeader.sample_rate property.
        """
        
        return self.header.sample_rate
        
    @property
    def central_freq(self):
        """
        Convenience wrapper for the Frame.FrameHeader.central_freq property.
        """
        
        return self.header.central_freq


def has_guppi_header(filehandle):
    """
    Determine if a VDIF file has a GUPPI header or not.
    
    .. versionadded:: 2.0.0
    """
    
    has_header = False
    with FilePositionSaver(filehandle):
        # Read in the first 16kB
        block = filehandle.read(16384)
        try:
            block = block.decode(encoding='ascii', errors='ignore')
        except AttributeError:
            pass
            
        if block.find('TELESCOP') != -1 \
           or block.find('END') != -1 \
           or block.find('CONTINUE') != -1:
            has_header = True
            
    return has_header


def read_guppi_header(filehandle):
    """
    Read in a GUPPI header at the start of a VDIF file from the VLA.  The 
    contents of the header are returned as a dictionary.
    """
    
    # Is there a GUPPI header?
    header = {}
    if not has_guppi_header(filehandle):
        warnings.warn("GUPPI header not found, returning an empty dictionary", RuntimeWarning)
        return header
        
    # Read in the GUPPI header
    while True:
        line = filehandle.read(80)
        try:
            line = line.decode(encoding='ascii', errors='ignore')
        except AttributeError:
            pass
            
        if line[:3] == 'END':
            break
        elif line[:8] == 'CONTINUE':
            junk, value2 = line.split(None, 1)
            value = "%s%s" % (value[:-1], value2[:-1])
        else:
            name, value = line.split('=', 1)
            name = name.strip()
        try:
            value = int(value, 10)
        except:
            try:
                value = float(value)
            except:
                value = value.strip().replace("'", '')
        header[name.strip()] = value
    header['OBSBW'] *= 1e6
    header['OBSFREQ'] *= 1e6
    
    return header


def read_frame(filehandle, sample_rate=0.0, central_freq=0.0, verbose=False):
    """
    Function to read in a single VDIF frame (header+data) and store the 
    contents as a Frame object.  This function wraps the _readerHeader and 
    _readData functions.
    """
    
    # New _vdif method
    try:
        newFrame = read_vdif(filehandle, Frame(), central_freq=central_freq, sample_rate=sample_rate)
    except gSyncError:
        mark = filehandle.tell()
        raise SyncError(type='VDIF', location=mark)
    except gEOFError:
        raise EOFError
        
    return newFrame


def get_frame_size(filehandle, nframes=None):
    """
    Find out what the frame size is in bytes from a single observation.
    """
    
    with FilePositionSaver(filehandle):
        # Read in one frame
        newFrame = read_frame(filehandle)
        
    return newFrame.header.frame_length*8


def get_thread_count(filehandle):
    """
    Find out how many thrads are present by examining the first 1024
    records.  Return the number of threads found.
    """
    
    # Get the frame size
    frame_size = get_frame_size(filehandle)
    
    with FilePositionSaver(filehandle):
        # Build up the list-of-lists that store ID codes and loop through 1024
        # frames.  In each case, parse pull the thread ID and append the thread 
        # ID to the relevant thread array if it is not already there.
        threads = []
        i = 0
        while i < 1024:
            try:
                cFrame = read_frame(filehandle)
            except SyncError:
                filehandle.seek(frame_size, 1)
                continue
            except EOFError:
                break
                
            cID = cFrame.header.thread_id
            if cID not in threads:
                threads.append(cID)
            i += 1
            
    # Return the number of threads found
    return len(threads)


def get_frames_per_second(filehandle):
    """
    Find out the number of frames per second in a file by watching how the 
    headers change.  Returns the number of frames in a second.
    """
    
    # Get the frame size
    frame_size = get_frame_size(filehandle)
    
    # Get the number of threads
    nThreads = get_thread_count(filehandle)
    
    with FilePositionSaver(filehandle):
        # Get the current second counts for all threads
        ref = {}
        i = 0
        while i < nThreads:
            try:
                cFrame = read_frame(filehandle)
            except SyncError:
                filehandle.seek(frame_size, 1)
                continue
            except EOFError:
                break
                
            cID = cFrame.header.thread_id
            cSC = cFrame.header.seconds_from_epoch
            ref[cID] = cSC
            i += 1
            
        # Read frames until we see a change in the second counter
        cur = {}
        fnd = []
        while True:
            ## Get a frame
            try:
                cFrame = read_frame(filehandle)
            except SyncError:
                filehandle.seek(frame_size, 1)
                continue
            except EOFError:
                break
                
            ## Pull out the relevant metadata
            cID = cFrame.header.thread_id
            cSC = cFrame.header.seconds_from_epoch
            cFC = cFrame.header.frame_in_second
            
            ## Figure out what to do with it
            if cSC == ref[cID]:
                ### Same second as the reference, save the frame number
                cur[cID] = cFC
            else:
                ### Different second than the reference, we've found something
                ref[cID] = cSC
                if cID not in fnd:
                    fnd.append( cID )
                    
            if len(fnd) == nThreads:
                break
                
    # Pull out the mode
    mode = {}
    for key in cur.keys():
        value = cur[key]
        try:
            mode[value] += 1
        except KeyError:
            mode[value] = 1
    best, bestValue = 0, 0
    for key in mode.keys():
        value = mode[key]
        if value > bestValue:
            best = key
            bestValue = value
            
    # Correct for a zero-based counter and return
    best += 1
    return best


def get_sample_rate(filehandle):
    """
    Find and return the sample rate in Hz by looking at how many frames 
    there are per second and how many samples there are in a frame.
    """
    
    # Get the number of frames per second
    nFramesSecond = get_frames_per_second(filehandle)
    
    with FilePositionSaver(filehandle):
        # Read in a frame
        cFrame = read_frame(filehandle)
        
    # Get the sample rate
    sample_rate = cFrame.payload.data.shape[-1] * nFramesSecond
    return float(sample_rate)
