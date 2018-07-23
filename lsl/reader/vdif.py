# -*- coding: utf-8 -*-

# Python3 compatiability
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

For reading in data, use the readFrame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  The readBlock
function reads in a (user-defined) number of VDIF frames and returns a 
ObservingBlock object.

For describing the format of data in the file, two function are provided:

getThreadCount
  read in the first few frames of an open file handle and return how many 
  threads are present in the file.
"""

import copy
from datetime import datetime

from lsl import astro
from lsl.common.mcs import datetime2mjdmpm
from lsl.reader._gofast import readVDIF
from lsl.reader._gofast import syncError as gsyncError
from lsl.reader._gofast import eofError as geofError
from lsl.reader.errors import syncError, eofError


__version__ = '0.4'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readGUPPIHeader', 'readFrame', 
           'getFrameSize', 'getThreadCount', 'getFramesPerSecond', 'getSampleRate', 
           '__version__', '__revision__', '__all__']



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
    
    def __init__(self, isInvalid=0, isLegacy=0, secondsFromEpoch=0, refEpoch=0, frameInSecond=0, version=1, nChan=0, frameLength=0, isComplex='C', bitsPerSample=0, threadID=0, stationID=0, extendedData1=None, extendedData2=None, extendedData3=None, extendedData4=None, sampleRate=0.0, centralFreq=0.0):
        self.isInvalid = isInvalid
        self.isLegacy = isLegacy
        self.secondsFromEpoch = secondsFromEpoch
        
        self.refEpoch = refEpoch
        self.frameInSecond = frameInSecond
        
        self.version = version
        self.nChan = nChan
        self.frameLength = frameLength
        
        self.isComplex = isComplex
        self.bitsPerSample = bitsPerSample
        self.threadID = threadID
        self.stationID = stationID
        
        self.extendedData1 = extendedData1
        self.extendedData2 = extendedData2
        self.extendedData3 = extendedData3
        self.extendedData4 = extendedData4
        
        self.sampleRate = sampleRate
        self.centralFreq = centralFreq
        
    def getTime(self):
        """
        Function to convert the time tag to seconds since the UNIX epoch.
        """
        
        # Get the reference epoch in the strange way that it is stored in VDIF 
        # and convert it to a MJD
        epochDT = datetime(2000+self.refEpoch/2, (self.refEpoch % 2)*6+1, 1, 0, 0, 0, 0)
        epochMJD, epochMPM = datetime2mjdmpm(epochDT)
        epochMJD = epochMJD + epochMPM/1000.0/86400.0
        
        # Get the frame MJD by adding the secondsFromEpoch value to the epoch
        frameMJD = epochMJD + self.secondsFromEpoch / 86400.0
        
        if self.sampleRate == 0.0:
            # Try to get the sub-second time by parsing the extended user data
            try:
                ## Is there a sample rate to grab?
                eud = self.parseExtendedUserData()
                sampleRate = eud['sampleRate']
                sampleRate *= 1e6 if eud['sampleRateUnits'] == 'MHz' else 1.0
            
                ## How many samples are in each frame?
                dataSize = self.frameLength*8 - 32 + 16*self.isLegacy		# 8-byte chunks -> bytes - full header + legacy offset
                samplesPerWord = 32 / self.bitsPerSample					# dimensionless
                nSamples = dataSize / 4 * samplesPerWord					# bytes -> words -> samples
            
                ## What is the frame rate?
                frameRate = sampleRate / nSamples
            
                frameMJD += 1.0*self.frameInSecond/frameRate/86400.0
            
            except KeyError:
                pass
                
        else:
            # Use what we already have been told
            ## How many samples are in each frame?
            dataSize = self.frameLength*8 - 32 + 16*self.isLegacy		# 8-byte chunks -> bytes - full header + legacy offset
            samplesPerWord = 32 / self.bitsPerSample				# dimensionless
            nSamples = dataSize / 4 * samplesPerWord				# bytes -> words -> samples
        
            ## What is the frame rate?
            frameRate = self.sampleRate / nSamples
            
            frameMJD += 1.0*self.frameInSecond/frameRate/86400.0

        # Convert from MJD to UNIX time
        seconds = astro.utcjd_to_unix(frameMJD + astro.MJD_OFFSET)
        
        return seconds
        
    def parseID(self):
        """
        Return a two-element tuple of the station ID and thread ID.
        
        .. note::
            The station ID is always returned as numeric.
        """
        return (self.stationID, self.threadID)
        
    def parseExtendedUserData(self):
        """
        Parse the extended user data if it was included with the reader.  
        The data is returned as a dictionary.
        """
        
        # Setup the output dictionary
        fields = {}
        
        # Is there anything to look at?
        if self.extendedData1 is None or self.extendedData2 is None or self.extendedData3 is None or self.extendedData4 is None:
            return fields
            
        # Extract the version
        edv = int((self.extendedData1 >> 24) & 0xFF)
        
        # Parse accordingly
        if edv == 0x00:
            ## Empty
            pass
            
        elif edv == 0x01:
            ## NICT
            fields['sampleRate'] = int(self.extendedData1 & (2**23-1))
            fields['sampleRateUnits'] = 'MHz' if int((self.extendedData1>>23) & 1) else 'kHz'
            fields['syncPattern'] = self.extendedData2
            fields['stationName'] = (self.extendedData4 << 32) | self.extendedData3
            
        elif edv == 0x02:
            ## ALMA
            fields['syncWord'] = int(self.extendedData1 & 0xFFFF)
            fields['picStatusWord'] = self.extendedData2
            fields['packetSerialNumber'] = (self.extendedData4 << 32) | self.extendedData3
            
        elif edv == 0x03:
            ## NRAO
            fields['sampleRate'] = int(self.extendedData1 & (2**23-1))
            fields['sampleRateUnits'] = 'MHz' if int((self.extendedData1 >> 23) & 1) else 'kHz'
            fields['syncPattern'] = self.extendedData2
            fields['tuningWord'] = self.extendedData3
            fields['dbeUnit'] = int((self.extendedData4 >> 24) & 0xF)
            fields['ifInputNumber'] = int((self.extendedData4 >> 20) & 0xF)
            fields['subBand'] = int((self.extendedData4 >> 17) & 0x7)
            fields['electronicSideBand'] = 'USB' if (self.extendedData4 >> 16) & 1 else 'LSB'
            mj = int((self.extendedData4 >> 12) & 0xF)
            mn = int((self.extendedData4 >>  8) & 0xF)
            pr = int(self.extendedData4 & 0xFF)
            fields['version'] = '%i.%i-%02f' % (mj, mn, pr)
            
        elif edv == 0xAB:
            ## Haystack (which is really an embedded Mark 5B header)
            fields['syncWord'] = self.extendedData1
            fields['yearsFrom2000'] = int((self.extendedData2 >> 28) & 0xF)
            fields['userSpecifiedData'] = int((self.extendedData2 >> 16) & 0xFFF)
            fields['dataFromInternalTVG'] = (self.extendedData2 >> 15) & 1
            fields['frameInSecond'] = int(self.extendedData2 & (2**14-1))
            j0 = int((self.extendedData3 >> 28) & 0xF)
            j1 = int((self.extendedData3 >> 24) & 0xF)
            j2 = int((self.extendedData3 >> 20) & 0xF)
            s0 = int((self.extendedData3 >> 16) & 0xF)
            s1 = int((self.extendedData3 >> 12) & 0xF)
            s2 = int((self.extendedData3 >>  8) & 0xF)
            s3 = int((self.extendedData3 >>  4) & 0xF)
            s4 = int((self.extendedData3 >>  0) & 0xF)
            f0 = int((self.extendedData4 >> 28) & 0xF)
            f1 = int((self.extendedData4 >> 24) & 0xF)
            f2 = int((self.extendedData4 >> 20) & 0xF)
            f3 = int((self.extendedData4 >> 16) & 0xF)
            f4 = int((self.extendedData4 >>  8) & 0xF)
            crcf = int(self.extendedData4 & 0xFFFF)
            crcf = (crcf & 0xF) << 8 | ((crcf >> 8) & 0xF)
            crcc = _crcc((self.extendedData3 << 16) | ((self.extendedData3 >> 16) & 0xFF), 48)
            fields['vlbaTimeCode'] = j0*1e7+j1*1e6+j2*1e5 + s0*1e4+s1*1e3+s2*1e2+s3*1e1+s4*1e0 + f0/1e1+f1/1e2+f2/1e3+f3/1e4+f4/1e5
            fields['vlbaTimeCodeValue'] = True if crcc == crcf else False
            
        else:
            raise RuntimeError("Unknown extended user data version: %i" % edv)
            
        return fields
        
    def getSampleRate(self):
        """
        Return the sample rate of the data in samples/second.
        """
        
        return self.sampleRate*1.0
        
    def getCentralFreq(self):
        """
        Function to get the central frequency of the VDIF data in Hz.
        """
        
        return self.centralFreq*1.0


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

    def parseID(self):
        """
        Convenience wrapper for the Frame.FrameHeader.parseID 
        function.
        """
        
        return self.header.parseID()
        
    def parseExtendedUserData(self):
        """
        Convenience wrapper for the Frame.FrameHeader.parseExtendedUserData
        function.
        """
        
        return self.header.parseExtendedUserData()
        
    def getTime(self):
        """
        Convenience wrapper for the Frame.FrameHeader.getTime function.
        """
        
        return self.header.getTime()
        
    def getSampleRate(self):
        """
        Convenience wrapper for the Frame.FrameHeader.getSampleRate function.
        """
        
        return self.header.getSampleRate()
        
    def getCentralFreq(self):
        """
        Convenience wrapper for the Frame.FrameHeader.getCentralFreq function.
        """
        
        return self.header.getCentralFreq()
        
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
    #	tX = self.data.timeTag
    #	try:
    #		tY = y.data.timeTag
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
    #	tX = self.data.timeTag
    #	try:
    #		tY = y.data.timeTag
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
    #	tX = self.data.timeTag
    #	try:
    #		tY = y.data.timeTag
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
    #	tX = self.data.timeTag
    #	try:
    #		tY = y.data.timeTag
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
    #	tX = self.data.timeTag
    #	try:
    #		tY = y.data.timeTag
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
    #	tX = self.data.timeTag
    #	try:
    #		tY = y.data.timeTag
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
    #	tX = self.data.timeTag
    #	tY = y.data.timeTag
    #	if tY > tX:
    #		return -1
    #	elif tX > tY:
    #		return 1
    #	else:
    #		return 0


def readGUPPIHeader(filehandle):
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


def readFrame(filehandle, sampleRate=0.0, centralFreq=0.0, Verbose=False):
    """
    Function to read in a single VDIF frame (header+data) and store the 
    contents as a Frame object.  This function wraps the _readerHeader and 
    _readData functions.
    """
    
    # New _vdif method
    try:
        newFrame = readVDIF(filehandle, Frame(), centralFreq=centralFreq, sampleRate=sampleRate)
    except gsyncError:
        mark = filehandle.tell()
        raise syncError(location=mark)
    except geofError:
        raise eofError
        
    return newFrame


def getFrameSize(filehandle, nFrames=None):
    """
    Find out what the frame size is in bytes from a single observation.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()

    # Read in one frame
    newFrame = readFrame(filehandle)
    
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    return newFrame.header.frameLength*8


def getThreadCount(filehandle):
    """
    Find out how many thrads are present by examining the first 1024
    records.  Return the number of threads found.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Get the frame size
    FrameSize = getFrameSize(filehandle)
    
    # Build up the list-of-lists that store ID codes and loop through 1024
    # frames.  In each case, parse pull the thread ID and append the thread 
    # ID to the relevant thread array if it is not already there.
    threads = []
    i = 0
    while i < 1024:
        try:
            cFrame = readFrame(filehandle)
        except syncError:
            filehandle.seek(FrameSize, 1)
            continue
        except eofError:
            break
            
        cID = cFrame.header.threadID
        if cID not in threads:
            threads.append(cID)
        i += 1
        
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Return the number of threads found
    return len(threads)


def getFramesPerSecond(filehandle):
    """
    Find out the number of frames per second in a file by watching how the 
    headers change.  Returns the number of frames in a second.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Get the frame size
    FrameSize = getFrameSize(filehandle)
    
    # Get the number of threads
    nThreads = getThreadCount(filehandle)
    
    # Get the current second counts for all threads
    ref = {}
    i = 0
    while i < nThreads:
        try:
            cFrame = readFrame(filehandle)
        except syncError:
            filehandle.seek(FrameSize, 1)
            continue
        except eofError:
            break
            
        cID = cFrame.header.threadID
        cSC = cFrame.header.secondsFromEpoch
        ref[cID] = cSC
        i += 1
        
    # Read frames until we see a change in the second counter
    cur = {}
    fnd = []
    while True:
        ## Get a frame
        try:
            cFrame = readFrame(filehandle)
        except syncError:
            filehandle.seek(FrameSize, 1)
            continue
        except eofError:
            break
            
        ## Pull out the relevant metadata
        cID = cFrame.header.threadID
        cSC = cFrame.header.secondsFromEpoch
        cFC = cFrame.header.frameInSecond
        
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


def getSampleRate(filehandle):
    """
    Find and return the sample rate in Hz by looking at how many frames 
    there are per second and how many samples there are in a frame.
    """
    
    # Save the current position in the file so we can return to that point
    fhStart = filehandle.tell()
    
    # Get the number of frames per second
    nFramesSecond = getFramesPerSecond(filehandle)
    
    # Read in a frame
    cFrame = readFrame(filehandle)
    
    # Return to the place in the file where we started
    filehandle.seek(fhStart)
    
    # Get the sample rate
    sampleRate = cFrame.data.data.shape[-1] * nFramesSecond
    return float(sampleRate)