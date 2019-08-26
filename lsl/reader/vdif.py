# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import division
import sys
if sys.version_info > (3,):
    xrange = range
    
"""
Python module to read in VDIF data.  This module defines the following 
classes for storing the VIDF data found in a file:

Frame
  object that contains all data associated with a particular DRX frame.  
  The primary constituents of each frame are:
    * FrameHeader - the VDIF frame header object and
    * FrameData   - the VDIF frame data object.
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

import copy
from datetime import datetime

from lsl import astro
from lsl.common.mcs import datetime_to_mjdmpm
from lsl.reader._gofast import readVDIF
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError


__version__ = '0.4'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'read_guppi_header', 'read_frame', 
           'get_frame_size', 'get_thread_count', 'get_frames_per_second', 'get_sample_rate']



def _crcc(data, length=48, mask=0o40003, cycle=16):
    """
    Compute a CRC checksum for the VLBA BCD time code stored in a Mark 5B 
    header.
    
    From:  mk5blib.c
    """
    
    state = 0
    for i in xrange(length):
        q = state & 1
        if ((data >> i) & 1) ^ q == 0:
            state &= -2
        else:
            state ^= mask
            state |= i
        state = (state >> 1) | (state & 1) << (cycle-1)
    return state


class FrameHeader(object):
    """
    Class that stores the information found in the header of a VDIF 
    frame.  Most fields in the VDIF version 1.1.1 header are stored.
    """
    
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
        
    @property
    def time(self):
        """
        Function to convert the time tag to seconds since the UNIX epoch as a two-
        element tuple.
        """
        
        # Get the reference epoch in the strange way that it is stored in VDIF 
        # and convert it to a MJD
        epochDT = datetime(2000+self.ref_epoch//2, (self.ref_epoch % 2)*6+1, 1, 0, 0, 0, 0)
        epochMJD, epochMPM = datetime_to_mjdmpm(epochDT)
        epochMJD = epochMJD + epochMPM/1000.0/86400.0
        
        # Get the frame MJD by adding the seconds_from_epoch value to the epoch
        frameMJD_i = epochMJD + self.seconds_from_epoch // 86400
        frameMJD_f = (self.seconds_from_epoch % 86400) / 86400.0
        
        if self.sample_rate == 0.0:
            # Try to get the sub-second time by parsing the extended user data
            try:
                ## Is there a sample rate to grab?
                eud = self.extended_user_data
                sample_rate = eud['sample_rate']
                sample_rate *= 1e6 if eud['sample_rate_units'] == 'MHz' else 1e3
            
                ## How many samples are in each frame?
                dataSize = self.frame_length*8 - 32 + 16*self.is_legacy		# 8-byte chunks -> bytes - full header + legacy offset
                samplesPerWord = 32 // self.bits_per_sample					# dimensionless
                nSamples = dataSize // 4 * samplesPerWord					# bytes -> words -> samples
            
                ## What is the frame rate?
                frameRate = sample_rate // nSamples
                
                frameMJD_f += 1.0*self.frame_in_second/frameRate/86400.0
            
            except KeyError:
                pass
                
        else:
            # Use what we already have been told
            ## How many samples are in each frame?
            dataSize = self.frame_length*8 - 32 + 16*self.is_legacy		# 8-byte chunks -> bytes - full header + legacy offset
            samplesPerWord = 32 // self.bits_per_sample				# dimensionless
            nSamples = dataSize // 4 * samplesPerWord				# bytes -> words -> samples
        
            ## What is the frame rate?
            frameRate = self.sample_rate // nSamples
            
            frameMJD_f += 1.0*self.frame_in_second/frameRate/86400.0

        # Convert from MJD to UNIX time
        if frameMJD_f > 1:
            frameMJD_i += 1
            frameMJD_f -= 1
        seconds_i = int(astro.utcjd_to_unix(frameMJD_i + astro.MJD_OFFSET))
        seconds_f = frameMJD_f * 86400.0
        
        return seconds_i, seconds_f
        
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


class FrameData(object):
    """
    Class that stores the information found in the data section of a VDIF
    frame.
    
    .. note::
        Unlike the other readers in the :mod:`lsl.reader` module the
        data are stored as numpy.float32 values.
    """
    
    def __init__(self, data=None):
        self.data = data


class Frame(object):
    """
    Class that stores the information contained within a single VDIF 
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
            self.data.data += y.data.data
        except AttributeError:
            self.data.data += y
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
            self.data.data *= y.data.data
        except AttributeError:
            self.data.data *= y
        return self
        
    #def __eq__(self, y):
    #	"""
    #	Check if the time tags of two frames are equal or if the time
    #	tag is equal to a particular value.
    #	"""
    #	
    #	tX = self.data.timetag
    #	try:
    #		tY = y.data.timetag
    #	except AttributeError:
    #		tY = y
    #	
    #	if tX == tY:
    #		return True
    #	else:
    #		return False
    #		
    #def __ne__(self, y):
    #	"""
    #	Check if the time tags of two frames are not equal or if the time
    #	tag is not equal to a particular value.
    #	"""
    #	
    #	tX = self.data.timetag
    #	try:
    #		tY = y.data.timetag
    #	except AttributeError:
    #		tY = y
    #	
    #	if tX != tY:
    #		return True
    #	else:
    #		return False
    #		
    #def __gt__(self, y):
    #	"""
    #	Check if the time tag of the first frame is greater than that of a
    #	second frame or if the time tag is greater than a particular value.
    #	"""
    #	
    #	tX = self.data.timetag
    #	try:
    #		tY = y.data.timetag
    #	except AttributeError:
    #		tY = y
    #	
    #	if tX > tY:
    #		return True
    #	else:
    #		return False
    #		
    #def __ge__(self, y):
    #	"""
    #	Check if the time tag of the first frame is greater than or equal to 
    #	that of a second frame or if the time tag is greater than a particular 
    #	value.
    #	"""
    #	
    #	tX = self.data.timetag
    #	try:
    #		tY = y.data.timetag
    #	except AttributeError:
    #		tY = y
    #	
    #	if tX >= tY:
    #		return True
    #	else:
    #		return False
    #		
    #def __lt__(self, y):
    #	"""
    #	Check if the time tag of the first frame is less than that of a
    #	second frame or if the time tag is greater than a particular value.
    #	"""
    #	
    #	tX = self.data.timetag
    #	try:
    #		tY = y.data.timetag
    #	except AttributeError:
    #		tY = y
    #	
    #	if tX < tY:
    #		return True
    #	else:
    #		return False
    #		
    #def __le__(self, y):
    #	"""
    #	Check if the time tag of the first frame is less than or equal to 
    #	that of a second frame or if the time tag is greater than a particular 
    #	value.
    #	"""
    #	
    #	tX = self.data.timetag
    #	try:
    #		tY = y.data.timetag
    #	except AttributeError:
    #		tY = y
    #	
    #	if tX <= tY:
    #		return True
    #	else:
    #		return False
    #		
    #def __cmp__(self, y):
    #	"""
    #	Compare two frames based on the time tags.  This is helpful for 
    #	sorting things.
    #	"""
    #	
    #	tX = self.data.timetag
    #	tY = y.data.timetag
    #	if tY > tX:
    #		return -1
    #	elif tX > tY:
    #		return 1
    #	else:
    #		return 0


def read_guppi_header(filehandle):
    """
    Read in a GUPPI header at the start of a VDIF file from the VLA.  The 
    contents of the header are returned as a dictionary.
    """
    
    # Read in the GUPPI header
    header = {}
    while True:
        line = filehandle.read(80)
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


def read_frame(filehandle, sample_rate=0.0, central_freq=0.0, Verbose=False):
    """
    Function to read in a single VDIF frame (header+data) and store the 
    contents as a Frame object.  This function wraps the _readerHeader and 
    _readData functions.
    """
    
    # New _vdif method
    try:
        newFrame = readVDIF(filehandle, Frame(), central_freq=central_freq, sample_rate=sample_rate)
    except gSyncError:
        mark = filehandle.tell()
        raise SyncError(location=mark)
    except gEOFError:
        raise EOFError
        
    return newFrame


def get_frame_size(filehandle, nFrames=None):
    """
    Find out what the frame size is in bytes from a single observation.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()

    # Read in one frame
    newFrame = read_frame(filehandle)
    
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    return newFrame.header.frame_length*8


def get_thread_count(filehandle):
    """
    Find out how many thrads are present by examining the first 1024
    records.  Return the number of threads found.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Get the frame size
    frame_size = get_frame_size(filehandle)
    
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
        
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Return the number of threads found
    return len(threads)


def get_frames_per_second(filehandle):
    """
    Find out the number of frames per second in a file by watching how the 
    headers change.  Returns the number of frames in a second.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Get the frame size
    frame_size = get_frame_size(filehandle)
    
    # Get the number of threads
    nThreads = get_thread_count(filehandle)
    
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
            
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Pull out the mode
    mode = {}
    for key,value in cur.iteritems():
        try:
            mode[value] += 1
        except KeyError:
            mode[value] = 1
    best, bestValue = 0, 0
    for key,value in mode.iteritems():
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
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Get the number of frames per second
    nFramesSecond = get_frames_per_second(filehandle)
    
    # Read in a frame
    cFrame = read_frame(filehandle)
    
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Get the sample rate
    sample_rate = cFrame.data.data.shape[-1] * nFramesSecond
    return float(sample_rate)
