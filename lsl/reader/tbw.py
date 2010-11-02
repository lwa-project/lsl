#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python module to reading in data from both 12-bit and 4-bit TBW files."""

import os
import sys
import time
import numpy
import array
import struct
import pyfits

from  ..common import dp as dp_common
from errors import *

__version__ = '0.3'
__revision__ = '$ Revision: 14 $'
__all__ = ['TBWFrameHeader', 'TBWFrameData', 'TBWFrame', 'readTBWFrame', 'TBWFrameSize', 'getDataBits', 'getFramesPerObs', '__version__', '__revision__', '__all__']

TBWFrameSize = 1224


class TBWFrameHeader(object):
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


class TBWFrameData(object):
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


class TBWFrame(object):
	"""Class that stores the information contained within a single TBW 
	frame.  It's properties are TBWFrameHeader and TBWFrameData objects."""

	def __init__(self, header=TBWFrameHeader(), data=TBWFrameData()):
		self.header = header
		self.data = data

	def parseID(self):
		"""Convenience wrapper for the TBWFrame.TBWFrameHeader.parseID 
		function."""
		
		return self.header.parseID()

	def getDataBits(self):
		"""Convenience wrapper for the TBWFrame.TBWFrameHeader.getDataBits 
		function."""
		
		return self.header.getDataBits()

	def getTime(self):
		"""Convenience wrapper for the TBWFrame.TBWFrameData.getTime function."""
		
		return self.data.getTime()


def __readTBWHeader(filehandle, Verbose=False):
	"""Private function to read in a TBW header.  Returns a TBWFrameHeader object."""

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

	newHeader = TBWFrameHeader()
	newHeader.frameCount = frameCount
	newHeader.secondsCount = secondsCount[0]
	newHeader.tbwID = tbwID

	if Verbose:
		stand = newHeader.parseID()
		print "Header: ", tbwID, secondsCount, junk
		print "  Stand: ", stand

	return newHeader


def __readTBWData12(filehandle):
	"""Private function to read in a TBW frame data section and unpack that data
	when is the 12-bit.  Returns a TBWFrameData object."""

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
	if rawData.shape[0] < 3*data.shape[0]:
		raise numpyError()

	data[0,:] = (rawData[0::3]<<4) | ((rawData[1::3]>>4)&15)
	data[1,:] = ((rawData[1::3]&15)<<8) | rawData[2::3]
	# The data are signed, so apply the two-complement rule to generate 
	# the negative values
	negativeValues = numpy.where( data >= 2048 )
	data[negativeValues] -= 4096

	newData = TBWFrameData()
	newData.timeTag = long(timeTag[0])
	newData.xy = data

	return newData
	
	
def __readTBWData4(filehandle):
	"""Private function to read in a TBW frame data section and unpack that data
	when is the 4-bit.  Returns a TBWFrameData object."""

	try:
		s = filehandle.read(8)
		timeTag = struct.unpack(">Q", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()
	
	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=1200)
	data = numpy.zeros((2,1200), dtype=numpy.int8)
	if rawData.shape[0] < data.shape[0]:
		raise numpyError()

	data[0,:] = (rawData>>4)&15
	data[1,:] = rawData&15
	# The data are signed, so apply the two-complement rule to generate 
	# the negative values
	negativeValues = numpy.where( data >= 8 )
	data[negativeValues] -= 16
	
	newData = TBWFrameData()
	newData.timeTag = long(timeTag[0])
	newData.xy = data
	
	return newData


def readTBWFrame(filehandle, Verbose=False):
	"""Function to read in a single TBW frame (header+data) and store the 
	contents as a TBWFrame object.  This function wraps readerHeader and 
	readData[(12)|4]."""

	try:
		hdr = __readTBWHeader(filehandle, Verbose=Verbose)
	except syncError, err:
		# Why?  If we run into a sync error here, then the following frame is invalid.  
		# Thus, we need to skip over this frame be advancing the file pointer 8+1200 B 
		currPos = filehandle.tell()
		frameEnd = currPos + TBWFrameSize - 16
		filehandle.seek(frameEnd)
		raise err

	if hdr.getDataBits() == 12:
		dat = __readTBWData12(filehandle)
	else:
		dat = __readTBWData4(filehandle)
	
	newFrame = TBWFrame()
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

	cFrame = readTBWFrame(filehandle)

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
			cFrame1 = readTBWFrame(filehandle)
			cFrame2 = readTBWFrame(filehandle)
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
		filehandle.seek(currentPosition+30000*TBWFrameSize)

	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Get the length of the stand list and return
	return len(idCodes)


def main(args):
	from ..writer import tsfits
	import matplotlib.pyplot as plt

	# Determine the number of samples in the specified file
	nSamples = os.path.getsize(args[0]) / TBWFrameSize
	print "Samples in file: ", nSamples
	fh = open(args[0], "rb", buffering=TBWFrameSize)

	# Make sure that the data is TBW and determine the data length
	test = readTBWFrame(fh)
	print "TBW Data:  %s" % test.header.isTBW()
	if not test.header.isTBW():
		raise notTBWError()
	print "Data Length: %i bits" % test.getDataBits()
	if test.header.getDataBits() == 12:
		nData = 400
	else:
		nData = 1200
	fh.seek(0)

	# Due to the size of the FITS files being generated, the number of frames that 
	# can be read in is limited to 300,000, or 30,000 frames for 10 stands.  Getting
	# around this limit will require finding out how to do on-the-fly FITS binary 
	# table resizing.  
	nSamples = 900000

	tStart = time.time()

	# Create a new FITS file with the name 'tbw.fits'
	fitsFile = tsfits.TBW('tbw-tsfits-test.fits')

	# Read in the data and add it to the FITS file created above
	count = {}
	syncCount = 0
	masterCount = 0
	for i in range(nSamples):
		# Read in the next frame and anticipate any problems that could occur
		try:
			cFrame = readTBWFrame(fh, Verbose=False)
		except eofError:
			break
		except syncError:
			#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFrameSize-1)
			syncCount = syncCount + 1
			continue
		except numpyError:
			break

		stand = cFrame.header.parseID()
		if cFrame.header.frameCount % 10000 == 0:
			print "%2i  %14i  %6.3f  %5i  %5i" % (stand, cFrame.data.timeTag, cFrame.getTime(), cFrame.header.frameCount, cFrame.header.secondsCount)
		if stand not in count.keys():
			count[stand] = 0

		fitsFile.addStandData(cFrame)

		count[stand] = count[stand] + 1
		masterCount = masterCount + 1

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (masterCount, (tEnd-tStart), masterCount/(tEnd-tStart))

	fh.close()
	fitsFile.info()

	# Summary information about the file that was just read in
	print "Summary:"
	for stand in sorted(count.keys()):
		print "Stand: %2i, Frames: %5i" % (stand, count[stand])
	print "Sync Errors: %5i" % syncCount


if __name__ == "__main__":
	main(sys.argv[1:])
