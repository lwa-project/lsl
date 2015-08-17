# -*- coding: utf-8 -*-

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
function reads in a (user-defined) number of DRX frames and returns a 
ObservingBlock object.

For describing the format of data in the file, two function are provided:

getThreadCount
  read in the first few frames of an open file handle and return how many 
  threads are present in the file.
"""

import copy
import numpy
import struct
from datetime import datetime

from errors import *

from lsl import astro
from lsl.common.mcs import datetime2mjdmpm


__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readFrame', 'getThreadCount', '__version__', '__revision__', '__all__']



def _crcc(data, length=48, mask=040003, cycle=16):
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
	
	def __init__(self, isInvalid=0, isLegacy=0, secondsFromEpoch=0, refEpoch=0, frameInSecond=0, version=1, nChan=0, frameLength=0, isComplex='C', bitsPerSample=0, threadID=0, stationID=0, extendedData1=None, extendedData2=None, extendedData3=None, extendedData4=None):
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
		
	def getTime(self):
		"""
		Function to convert the time tag to seconds since the UNIX epoch.
		"""
		
		# Get the reference epoch in the strange way that it is stored in VDIF 
		# and convert it to a MJD
		epochDT = datetime(2000+self.refEpoch/2, (self.refEpoch % 2)*6+1, 1, 0, 0, 0, 0)
		epochMJD, epochMPM = datetime2mjdmpm(epochDT)
		
		# Get the frame MJD by adding the secondsFromEpoch value to the epoch
		frameMJD = epochMJD + self.secondsFromEpoch / 86400.0
		
		# Try to get the sub-second time by parsing the extended user data
		try:
			## Is there a sample rate to grab?
			eud = self.parseExtendedUserData()
			sampleRate = eud['sampleRate']
			sampleRate *= 1e6 if eud['sampleRateUnits'] == 'MHz' else 1.0
			
			## How many samples are in each frame?
			dataSize = self.frameLength*8 - 32 + 16*self.isLegacy		# 8-byte chunks -> bytes - full header + legacy offset
			samplesPerWord = 32 / self.bitsPerSample				# dimensionless
			nSamples = dataSize / 4 * samplesPerWord				# bytes -> words -> samples
			
			## What is the frame rate?
			frameRate = sampleRate / nSamples
			
			frameMJD += 1.0*self.frame/frameRate
			
		except KeyError:
			pass
			
		# Convert from MJD to UNIX time
		seconds = astro.utcjd_to_unix(frameMJD + astro.MJD_OFFSET)
		
		return seconds
		
	def parseID(self):
		"""
		Parse the thread ID into a...  thread ID.
		"""
		
		return self.threadID
		
	def parseExtendedUserData(self):
		"""
		Parse the extended user data if it was included with the reader.  
		The data is returned as a dictionary.
		"""
		
		# Setup the output dictionary
		fields = {}
		
		# Is there anything to look at?
		if self.extendedData1 is None or self.extendedData2 or self.extendedData3 is None or self.extendedData4 is None:
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
			fields['version'] = '%i.%i-%02f' % (mj,mn,pr)
			
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
			crcf = int(self.extendedData4 & 0xFFFF)
			crcf = (crcf & 0xF) << 8 | ((crcf >> 8) & 0xF)
			crcf = _crcc((self.extendedData3 << 16) | ((self.extendedData >> 16) & 0xFF), 48)
			fields['vlbaTimeCode'] = j0*1e7+j1*1e6+j2*1e5 + s0*1e4+s1*1e3+s2*1e2+s3*1e1+s4*1e0 + f0/1e1+f1/1e2+f2/1e3+f4/1e4
			fields['vlbaTimeCodeValue'] = True if crcc == crcf else False
			
		else:
			raise RuntimeError("Unknown extended user data version: %i" % edv)
			
		return fields


class FrameData(object):
	"""
	Class that stores the information found in the data section of a VDIF
	frame.
	"""
	
	def __init__(self, data=None):
		self.data = data


class Frame(object):
	"""
	Class that stores the information contained within a single DRX 
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
			self.data.iq += y
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
			self.data.iq *= y
		return self
			
	def __eq__(self, y):
		"""
		Check if the time tags of two frames are equal or if the time
		tag is equal to a particular value.
		"""
		
		tX = self.data.timeTag
		try:
			tY = y.data.timeTag
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
		
		tX = self.data.timeTag
		try:
			tY = y.data.timeTag
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
		
		tX = self.data.timeTag
		try:
			tY = y.data.timeTag
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
		
		tX = self.data.timeTag
		try:
			tY = y.data.timeTag
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
		
		tX = self.data.timeTag
		try:
			tY = y.data.timeTag
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
		
		tX = self.data.timeTag
		try:
			tY = y.data.timeTag
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
		
		tX = self.data.timeTag
		tY = y.data.timeTag
		if tY > tX:
			return -1
		elif tX > tY:
			return 1
		else:
			return 0


def _readHeader(filehandle):
	"""
	Function to read in a header block and parse it.  This function returns 
	a fully-populated FrameHeader instance.
	"""
	
	# Base header that is common to both standard and legacy formats
	baseHeaderData = filehandle.read(16)
	word0, word1, word2, word3 = struct.unpack('<4I', baseHeaderData)
	
	## First word:  invalid, legacy, and seconds from the reference epoch
	isInvalid = (word0 >> 31) & 1
	isLegacy = (word0 >> 30) & 1
	secondsFromEpoch = int(word0 & (2**30-1))
	
	## Second word: reference epoch and the frame number within the second
	refEpoch = int((word1>>24) & (2**6-1))
	frameInSecond = int(word1 & (2**24-1))
	
	## Third word: VIDF version, the number of channels in the frame, and 
	## the frame lenght in units of 8 bytes
	version = int((word2 >> 29) & 0x7)
	log2nChan = int((word2 >> 24) & (2**5-1))
	nChan = 2**log2nChan
	frameLength = int(word2 & (2**23-1))
	
	## Fourth word: Data type/size, thread ID, and station ID
	isComplex = int((word3 >> 31) & 1)
	bitsPerSample = int((word3 >> 26) & (2**5-1))
	threadID = int((word3 >> 16) & (2**10-1))
	stationID = word3 & 0xFFFF
	if int(stationID & 0xFF) < 48:
		### If ord(of the first character) is less than 48 then this is a numeric ID
		stationID = int(stationID)
	else:
		### Otherwise this is a string
		stationID = "%s%s" % (chr((stationID >> 8) & 0xFF), chr(stationID & 0xFF))
		
	# Legacy vs. Standard headers
	if isLegacy:
		## For legacy headers we don't have any extended user data
		extendedData1, extendedData2, extendedData3, extendedData4 = None, None, None, None
	else:
		## For standard headers there is an additional 16 bytes (4  words) 
		## that contain various other information.  The parseing of this 
		## data is EDV-specific and is offloaded to the FrameHeader class.
		extHeaderData = filehandle.read(16)
		extendedData1, extendedData2, extendedData3, extendedData4 = struct.unpack('<4I', extHeaderData)
		
	# Build the FrameHeader instance
	header = FrameHeader(isInvalid=isInvalid, isLegacy=isLegacy, secondsFromEpoch=secondsFromEpoch, 
					 refEpoch=refEpoch, frameInSecond=frameInSecond, 
					 version=version, nChan=nChan, frameLength=frameLength, 
					 isComplex=isComplex, bitsPerSample=bitsPerSample, threadID=threadID, stationID=stationID, 
					 extendedData1=extendedData1, extendedData2=extendedData2, extendedData3=extendedData3, extendedData4=extendedData4)
					 
	# Done
	return header


def _readData(filehandle, header):
	"""
	Function to read in a data block and parse it.  This function returns 
	a fully-populated FrameData instance.
	
	.. note::
		In order to correctly unpack the data the frame header needs to also
		by supplied.
	"""
	
	# Take the frame length, legacy header flag, and bits per sample and use them
	# to figure out how many samples are in the data section of the frame
	dataSize = header.frameLength*8 - 32 + 16*header.isLegacy		# 8-byte chunks -> bytes - full header + legacy offset
	samplesPerWord = 32 / header.bitsPerSample					# dimensionless
	nSamples = dataSize / 4 * samplesPerWord					# bytes -> words -> samples
	
	# Loop over the 4-bytes words in the data second to read in and parse the data
	samples = []
	for i in xrange(dataSize/4):
		## Read in the word and interpret it as an unsigned integer
		dataWord = filehandle.read(4)
		dataWord, = struct.unpack('>I', dataWord)
		
		## Loop over the samples in the word and unpack
		for j in xrange(samplesPerWord):
			### Unpack
			sample = int((dataWord>>(header.bitsPerSample*j)) & (2**header.bitsPerSample-1))
			### Converted to signed (if we have more than a bit) using the 
			### fixed binary offset scheme
			if header.bitsPerSample > 1:
				sample = sample - (2**(header.bitsPerSample-1) - 1)
			### Save
			samples.append( sample)
	samples = numpy.array(samples)
	
	# Re-arrange the data as needed.
	## Real vs. Complex
	if header.isComplex:
		### For complex data this involves separating I and Q and then converting 
		### the data to numpy.complex64
		data = samples[0::2] + 1j*samples[1::2]
		data = data.astype(numpy.complex64)
		
	else:
		### For real data we need to play some games to make sure that the data
		### are cast to the correct times
		if header.bitsPerSample > 16:
			dtype = numpy.int32
		elif header.bitsPerSample > 8:
			dtype = numpy.int16
		elif header.bitsPerSample > 1:
			dtype = numpy.int8
		else:
			dtype = numpy.bool
		data = samples.astype(dtype)
	## Multi-channel vs. single channel
	if header.nChan > 1:
		### For multi-channel data convert the data into a 2-D array with
		### channel number running over the first axis
		data.shape = (data.size/header.nChan, header.nChan)
		data = data.T
		
	# Build the FrameData instance
	payload = FrameData(data=data)
	
	# Done
	return payload


def readFrame(filehandle, Verbose=False):
	"""
	Function to read in a single VDIF frame (header+data) and store the 
	contents as a Frame object.  This function wraps the _readerHeader and 
	_readData functions.
	"""
	
	# Read in the header and data payload
	try:
		header = _readHeader(filehandle)
		payload = _readData(filehandle, header)
	except (OSError, struct.error):
		raise eofError()
		
	# Build up the Frame instance
	frame = Frame(header, payload)
	
	# Done
	return frame


def getThreadCount(filehandle):
	"""
	Find out how many thrads are present by examining the first 1024
	records.  Return the number of threads found.
	"""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Build up the list-of-lists that store ID codes and loop through 1024
	# frames.  In each case, parse pull the thread ID and append the thread 
	# ID to the relevant thread array if it is not already there.
	threads = []
	for i in range(1024):
		try:
			cFrame = readFrame(filehandle)
		except eofError:
			continue
			
		cID = cFrame.header.threadID
		if cID not in threads:
			threads.append(cID)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the number of threads found
	return len(threads)
	