# -*- coding: utf-8 -*-

"""Python module to reading in data from both 12-bit and 4-bit TBW files."""

import os
import sys
import time
import numpy
import array
import struct
import pyfits

from  lsl.common import dp as dp_common
from errors import *

__version__ = '0.4'
__revision__ = '$ Revision: 17 $'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readFrame', 'FrameSize', 'getDataBits', 'getFramesPerObs', '__version__', '__revision__', '__all__']

FrameSize = 1224


class FrameHeader(object):
	"""Class that stores the information found in the header of a TBW 
	frame.  All three fields listed in the DP IDC version H are stored as 
	well as the original binary header data."""

	def __init__(self, frameCount=None, secondsCount=None, tbwID=None,  raw=None):
		self.frameCount = frameCount
		self.secondsCount = secondsCount
		self.tbwID = tbwID
		self.raw = raw

	def isTBW(self):
		"""Function to check if the data is really TBW and not TBN by examining
		the TBW ID field.  Returns True if the data is TBW, false otherwise."""

		mode = (self.tbwID>>15)&1
		if mode == 1:
			return True
		else:
			return False

	def parseID(self):
		"""Function to parse the TBW ID field and return the stand number."""

		# Why &1023?  Well, from DP ICH revision H, it seems that the stand count 
		# only goes up 260.  So, channel numbers should range up to 520, which can
		# be represented as 10 bits or 1023.
		stand = self.tbwID&1023

		return stand

	def getDataBits(self):
		"""Function to parse the TBW ID field and return the size of number of 
		bits that comprise the data.  12 is returned for 12-bit data, and 4 
		for 4-bit data."""

		bits = (self.tbwID>>14)&1
		if bits == 0:
			dataBits = 12
		else:
			dataBits = 4

		return dataBits


class FrameData(object):
	"""Class that stores the information found in the data section of a TBW
	frame.  Both fields listed in the DP IDC version H are stored."""

	def __init__self(self, timeTag=None, samples=400, xy=None):
		self.timeTag = timeTag
		self.xy = xy

	def getTime(self):
		"""Function to convert the time tag from samples since station 
		midnight to seconds since station midnight.  This function needs to 
		dp_common module in order to work."""

		seconds = self.timeTag / dp_common.fS
		
		return seconds


class Frame(object):
	"""Class that stores the information contained within a single TBW 
	frame.  It's properties are FrameHeader and FrameData objects."""

	def __init__(self, header=FrameHeader(), data=FrameData()):
		self.header = header
		self.data = data
		self.valid = True

	def parseID(self):
		"""Convenience wrapper for the Frame.FrameHeader.parseID 
		function."""
		
		return self.header.parseID()

	def getDataBits(self):
		"""Convenience wrapper for the Frame.FrameHeader.getDataBits 
		function."""
		
		return self.header.getDataBits()

	def getTime(self):
		"""Convenience wrapper for the Frame.FrameData.getTime function."""
		
		return self.data.getTime()


def __readHeader(filehandle, Verbose=False):
	"""Private function to read in a TBW header.  Returns a FrameHeader object."""

	rawHeader = ''
	try:
		s = filehandle.read(4)
		rawHeader = rawHeader + s
		sync4, sync3, sync2, sync1 = struct.unpack(">BBBB", s)
		s = filehandle.read(4)
		rawHeader = rawHeader + s
		m5cID, frameCount3, frameCount2, frameCount1 = struct.unpack(">BBBB", s)
		frameCount = (long(frameCount3)<<16) | (long(frameCount2)<<8) | long(frameCount1)
		s = filehandle.read(4)
		rawHeader = rawHeader + s
		secondsCount = struct.unpack(">L", s)
		s = filehandle.read(4)
		rawHeader = rawHeader + s
		tbwID, junk = struct.unpack(">HH", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()

	if sync1 != 92 or sync2 != 222 or sync3 != 192 or sync4 != 222:
		raise syncError(sync1=sync1, sync2=sync2, sync3=sync3, sync4=sync4)

	newHeader = FrameHeader()
	newHeader.frameCount = frameCount
	newHeader.secondsCount = secondsCount[0]
	newHeader.tbwID = tbwID

	if Verbose:
		stand = newHeader.parseID()
		print "Header: ", tbwID, secondsCount, junk
		print "  Stand: ", stand

	return newHeader


def __readData12(filehandle):
	"""Private function to read in a TBW frame data section and unpack that data
	when is the 12-bit.  Returns a FrameData object."""

	try:
		s = filehandle.read(8)
		timeTag = struct.unpack(">Q", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()

	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=1200)
	rawData = rawData.astype(numpy.uint16)
	data = numpy.zeros((2,400), dtype=numpy.int16)
	if rawData.shape[0] < 3*data.shape[1]:
		raise numpyError()

	data[0,:] = (rawData[0::3]<<4) | ((rawData[1::3]>>4)&15)
	data[1,:] = ((rawData[1::3]&15)<<8) | rawData[2::3]
	# The data are signed, so apply the two-complement rule to generate 
	# the negative values
	negativeValues = numpy.where( data >= 2048 )
	data[negativeValues] -= 4096

	newData = FrameData()
	newData.timeTag = long(timeTag[0])
	newData.xy = data

	return newData
	
	
def __readData4(filehandle):
	"""Private function to read in a TBW frame data section and unpack that data
	when is the 4-bit.  Returns a FrameData object."""

	try:
		s = filehandle.read(8)
		timeTag = struct.unpack(">Q", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()
	
	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=1200)
	data = numpy.zeros((2,1200), dtype=numpy.int8)
	if rawData.shape[0] < data.shape[1]:
		raise numpyError()

	data[0,:] = (rawData>>4)&15
	data[1,:] = rawData&15
	# The data are signed, so apply the two-complement rule to generate 
	# the negative values
	negativeValues = numpy.where( data >= 8 )
	data[negativeValues] -= 16
	
	newData = FrameData()
	newData.timeTag = long(timeTag[0])
	newData.xy = data
	
	return newData


def readFrame(filehandle, Verbose=False):
	"""Function to read in a single TBW frame (header+data) and store the 
	contents as a Frame object.  This function wraps readerHeader and 
	readData[(12)|4]."""

	try:
		hdr = __readHeader(filehandle, Verbose=Verbose)
	except syncError, err:
		# Why?  If we run into a sync error here, then the following frame is invalid.  
		# Thus, we need to skip over this frame be advancing the file pointer 8+1200 B 
		currPos = filehandle.tell()
		frameEnd = currPos + FrameSize - 16
		filehandle.seek(frameEnd)
		raise err

	if hdr.getDataBits() == 12:
		dat = __readData12(filehandle)
	else:
		dat = __readData4(filehandle)
	
	newFrame = Frame()
	newFrame.header = hdr
	newFrame.data = dat

	return newFrame


def getDataBits(filehandle):
	"""Find out the number of data bits used in the file be reading in the 
	first frame."""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Go back to the beginning...
	filehandle.seek(0)

	cFrame = readFrame(filehandle)

	dataBits = cFrame.getDataBits()

	# Return to the place in the file where we started
	filehandle.seek(fhStart)

	return dataBits


def getFramesPerObs(filehandle):
	"""Find out how many frames are present per observation by examining 
	the first frames for what would be 256 stands.  This is done by reading
	two frames and then skipping the next 30,000."""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Go back to the beginning...
	filehandle.seek(0)

	idCodes = []
	for i in range(3*256):
		currentPosition = filehandle.tell()
		try:
			cFrame1 = readFrame(filehandle)
			cFrame2 = readFrame(filehandle)
		except eofError:
			break
		except syncError:
			continue
		except numpyError:
			break

		cID = cFrame1.parseID()
		if cID not in idCodes:
			idCodes.append(cID)
		cID = cFrame2.parseID()
		if cID not in idCodes:
			idCodes.append(cID)

		# Junk 30,000 frames since that is how many frames there are per stand
		filehandle.seek(currentPosition+30000*FrameSize)

	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Get the length of the stand list and return
	return len(idCodes)
