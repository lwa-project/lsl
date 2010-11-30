# -*- coding: utf-8 -*-

"""Python module to read in DRX data.  This module defines the following 
classes for storing the DRX data found in a file:
  Frame          - object that contains all data associated with a particular DRX
                   frame.  The primary constituents of each frame are:
                     FrameHeader - the DRX frame header object and
                     FrameData   - the DRX frame data object.
                   Combined, these two objects contain all of the information found in
                   the original DRX frame.
  ObservingBlock - object that stores a collection of Frames for all beams/
                   tunnings/polarizations for a particular time.

The functions defined in this module fall into two class:
  1)  convert a frame in a file to a Frame object and
  2)  describe the format of the data in the file.

For reading in data, use the readFrame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  The readBlock
function reads in a (user-defined) number of DRX frames and returns a 
ObservingBlock object.

For describing the format of data in the file, two function are provided:
  getBeamCount -    read in the first few frames of an open file handle and 
                    return how many beams are present in the file.
  getFramesPerObs - read in the first several frames to see how many frames
                    (tunnings/polarizations) are associated with each beam.

Finally, there are also two experimental functions defined in this file for 
averaging DRX observations.  averageObservations performs a straight channel-
by-channel average of a collection of DRX frames, which averageObservations2
allows for both temporal and spatial averaging."""

import os
import sys
import math
import time
import numpy
import struct

from  lsl.common import dp as dp_common
from errors import *

__version__ = '0.4'
__revision__ = '$ Revision: 16 $'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'ObservingBlock', 'readFrame', 'readBlock', 'getBeamCount', 'getFramesPerObs', 'averageObservations', 'averageObservations2', 'FrameSize', 'filterCodes', '__version__', '__revision__', '__all__']

FrameSize = 4128

# List of filter codes and their corresponding sample rates in Hz
filterCodes = {1: 250000, 2: 500000, 3: 1000000, 4: 2000000, 5: 4000000, 6: 9800000, 7: 19600000}


class FrameHeader(object):
	"""Class that stores the information found in the header of a DRX 
	frame.  All six fields listed in the DP IDC version H are stored as 
	well as the original binary header data."""
	
	def __init__(self, frameCount=None, drxID=None, secondsCount=None, decimation=None, timeOffset=None, raw=None):
		self.frameCount = frameCount
		self.drxID = drxID
		self.secondsCount = secondsCount
		self.decimation = decimation
		self.timeOffset = timeOffset
		self.raw = raw
	
	def parseID(self):
		"""Parse the DRX ID into a tuple containing the beam (1 through
		4), tunning (1 and 2), and polarization (0 and 1)."""
		
		beam = self.drxID&7
		tune = (self.drxID>>3)&7
		pol  = (self.drxID>>7)&1

		return (beam, tune, pol)
	
	def getSampleRate(self):
		"""Return the sample rate of the data in samples/second."""
		
		sampleRate = dp_common.fS / self.decimation
		return sampleRate


class FrameData(object):
	"""Class that stores the information found in the data section of a DRX
	frame.  All three fields listed in the DP IDC version H are stored."""

	def __init__(self, timeTag=None, flags=None, iq=None):
		self.timeTag = timeTag
		self.flags = flags
		self.iq = iq

	def getTime(self):
		"""Function to convert the time tag from samples since station 
		midnight to seconds since station midnight.  This function needs 
		the dp_common module in order to work."""

		seconds = self.timeTag / dp_common.fS
		
		return seconds


class Frame(object):
	"""Class that stores the information contained within a single DRX 
	frame.  It's properties are FrameHeader and FrameData objects."""

	def __init__(self, header=FrameHeader(), data=FrameData()):
		self.header = header
		self.data = data
		self.valid = True

	def parseID(self):
		"""Convenience wrapper for the Frame.FrameHeader.parseID 
		function."""
		
		return self.header.parseID()

	def getSampleRate(self):
		"""Convenience wrapper for the Frame.FrameHeader.getSampleRate 
		function."""
		
		return self.header.getSampleRate()

	def getTime(self):
		"""Convenience wrapper for the Frame.FrameData.getTime function."""
		
		return self.data.getTime()
			
	def __add__(self, y):
		"""Add the data sections of two frames together or add a number 
		to every element in the data section."""
		
		try:
			return self.data.xy + y.self.data.xy
		except AttributeError:
			return self.data.xy + y
			
	def __iadd__(self, y):
		"""In-place add the data sections of two frames together or add 
		a number to every element in the data section."""
		
		try:
			self.data.xy += y.self.data.xy
		except AttributeError:
			self.data.xy += y
		
	def __mul__(self, y):
		"""Multiple the data sections of two frames together or multiply 
		a number to every element in the data section."""
		
		try:
			return self.data.xy * y.self.data.xy
		except AttributeError:
			return self.data.xy * y
			
	def __imul__(self, y):
		"""In-place multiple the data sections of two frames together or 
		multiply a number to every element in the data section."""
		
		try:
			self.data.xy *= y.self.data.xy
		except AttributeError:
			self.data.xy *= y
			
	def __eq__(self, y):
		"""Check if the time tags of two frames are equal or if the time
		tag is equal to a particular value."""
		
		tX = self.data.timeTag
		try:
			tY = self.data.timeTag
		except AttributeError:
			tY = y
		
		if tX == tY:
			return True
		else:
			return False
			
	def __ne__(self, y):
		"""Check if the time tags of two frames are not equal or if the time
		tag is not equal to a particular value."""
		
		tX = self.data.timeTag
		try:
			tY = self.data.timeTag
		except AttributeError:
			tY = y
		
		if tX != tY:
			return True
		else:
			return False
			
	def __gt__(self, y):
		"""Check if the time tag of the first frame is greater than that of a
		second frame or if the time tag is greater than a particular value."""
		
		tX = self.data.timeTag
		try:
			tY = self.data.timeTag
		except AttributeError:
			tY = y
		
		if tX > tY:
			return True
		else:
			return False
			
	def __ge__(self, y):
		"""Check if the time tag of the first frame is greater than or equal to 
		that of a second frame or if the time tag is greater than a particular 
		value."""
		
		tX = self.data.timeTag
		try:
			tY = self.data.timeTag
		except AttributeError:
			tY = y
		
		if tX >= tY:
			return True
		else:
			return False
			
	def __lt__(self, y):
		"""Check if the time tag of the first frame is less than that of a
		second frame or if the time tag is greater than a particular value."""
		
		tX = self.data.timeTag
		try:
			tY = self.data.timeTag
		except AttributeError:
			tY = y
		
		if tX < tY:
			return True
		else:
			return False
			
	def __le__(self, y):
		"""Check if the time tag of the first frame is less than or equal to 
		that of a second frame or if the time tag is greater than a particular 
		value."""
		
		tX = self.data.timeTag
		try:
			tY = self.data.timeTag
		except AttributeError:
			tY = y
		
		if tX <= tY:
			return True
		else:
			return False
			
	def __cmp__(self, y):
		"""Compare two frames based on the time tags.  This is helpful for 
		sorting things."""
		
		tX = self.data.timeTag
		tY = y.data.timeTag
		if tY > tX:
			return -1
		elif tX > tY:
			return 1
		else:
			return 0


class ObservingBlock(object):
	"""Class that stores all frames associates with a particular beam at a
	particular time."""

	def __init__(self, x1=Frame(), y1=Frame(), x2=Frame(), y2=Frame()):
		self.x1 = x1
		self.y1 = y1
		self.x2 = x2
		self.y2 = y2
		

def __readHeader(filehandle, Verbose=False):
	"""Private function to read in a DRX header.  Returns a FrameHeader object."""

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
		decimation, timeOffset = struct.unpack(">HH", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()

	if sync1 != 92 or sync2 != 222 or sync3 != 192 or sync4 != 222:
		raise syncError(sync1=sync1, sync2=sync2, sync3=sync3, sync4=sync4)

	drxID = m5cID
	newHeader = FrameHeader()
	newHeader.frameCount = frameCount
	newHeader.drxID = drxID
	newHeader.secondsCount = secondsCount
	newHeader.decimation = decimation
	newHeader.timeOffset = timeOffset
	newHeader.raw = rawHeader

	if Verbose:
		beam, tune, pol = newHeader.parseID()
		print "Header: ", drxID, decimation, timeOffset
		print "  Beam: ", beam
		print "  Tuning: ", tune
		print "  Polarization: ", pol

	return newHeader


def __readData(filehandle):
	"""Private function to read in a DRX frame data section.  Returns a 
	FrameData object."""

	try:
		s = filehandle.read(8)
		timeTag = struct.unpack(">Q", s)
		s = filehandle.read(8)
		flags = struct.unpack(">Q", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()
	
	# A truly excellent idea from Dan Wood
	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=4096)
	data = numpy.zeros(4096, dtype=numpy.complex_)
	if rawData.shape[0] < data.shape[0]:
		raise numpyError()
	
	data.real = (rawData>>4)&15
	data.imag = rawData&15

	# The data are signed, so apply the two-complement rule to generate 
	# the negative values
	negativeValues = numpy.where( data.real >= 8 )
	data.real[negativeValues] -= 16
	negativeValues = numpy.where( data.imag >= 8 )
	data.imag[negativeValues] -= 16
	
	newData = FrameData()
	newData.timeTag = timeTag[0]
	newData.flags = flags
	newData.iq = data

	return newData


def readFrame(filehandle, Verbose=False):
	"""Function to read in a single DRX frame (header+data) and store the 
	contents as a Frame object.  This function wraps readerHeader and 
	readData."""
	
	try:
		hdr = __readHeader(filehandle, Verbose=Verbose)
	except syncError, err:
		# Why?  If we run into a sync error here, then the following frame is invalid.  
		# Thus, we need to skip over this frame be advancing the file pointer 8+8+4096 B 
		currPos = filehandle.tell()
		frameEnd = currPos + FrameSize - 16
		filehandle.seek(frameEnd)
		raise err


	dat = __readData(filehandle)
	
	# Create the new frame object and return
	newFrame = Frame()
	newFrame.header = hdr
	newFrame.data = dat

	return newFrame


def readBlock(filehandle):
	"""Function to read in a single DRX block (four frames) and store the 
	contents as a ObservingBlock object.  This function wraps 
	readFrame."""
	
	try:
		x1 = readFrame(filehandle)
		y1 = readFrame(filehandle)
		x2 = readFrame(filehandle)
		y2 = readFrame(filehandle)
	except baseReaderError, err:
		raise err

	block = ObservingBlock(x1=x1, y1=y1, x2=x2, y2=y2)
	return block


def getBeamCount(filehandle):
	"""Find out how many beams are present by examining the first 16 DRX
	records.  Return the number of beams found."""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Go back to the beginning...
	filehandle.seek(0)

	# Build up the list-of-lists that store ID codes and loop through 32
	# frames.  In each case, parse pull the DRX ID, extract the beam number, 
	# and append the DRX ID to the relevant beam array if it is not already 
	# there.
	beams = []
	for i in range(16):
		cFrame = readFrame(filehandle)
		cID = cFrame.header.drxID
		beam = cID&7
		if beam not in beams:
			beams.append(beam)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)

	# Return the number of beams found
	return len(beams)


def getFramesPerObs(filehandle):
	"""Find out how many frames are present per beam by examining the first 
	16 DRX records.  Return the number of frames per observations as a four-
	element tuple, one for each beam."""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()

	# Go back to the beginning...
	filehandle.seek(0)
	
	# Build up the list-of-lists that store ID codes and loop through 32
	# frames.  In each case, parse pull the DRX ID, extract the beam number, 
	# and append the DRX ID to the relevant beam array if it is not already 
	# there.
	idCodes = [[], [], [], []]
	for i in range(16):
		cFrame = readFrame(filehandle)
		cID = cFrame.header.drxID
		beam = cID&7
		if cID not in idCodes[beam-1]:
			idCodes[beam-1].append(cID)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Get the length of each beam list and return them as a tuple
	return (len(idCodes[0]), len(idCodes[1]), len(idCodes[2]), len(idCodes[3]))
	

def averageObservations(Observations):
	"""Given a list of ObservingBlock objects, average the observations 
	together on a per tuning, per polarization basis.  A new ObsevingBlock 
	object is returned that contains the averaged data."""
	
	newBlock = Observations[0]
	tempx1 = newBlock.x1.data.iq * complex(0.0, 0.0)
	tempy1 = newBlock.y1.data.iq * complex(0.0, 0.0)
	tempx2 = newBlock.x2.data.iq * complex(0.0, 0.0)
	tempy2 = newBlock.y2.data.iq * complex(0.0, 0.0)
	
	obsCount = 0
	for Observation in Observations:
		tempx1 = tempx1 + Observation.x1.data.iq
		tempy1 = tempy1 + Observation.y1.data.iq
		tempx2 = tempx2 + Observation.x2.data.iq
		tempy2 = tempy2 + Observation.y2.data.iq
		obsCount = obsCount + 1
	
	tempx1 = tempx1 / float(obsCount)
	tempy1 = tempy1 / float(obsCount)
	tempx2 = tempx2 / float(obsCount)
	tempy2 = tempy2 / float(obsCount)
			
	newBlock.x1.data.iq = tempx1
	newBlock.y1.data.iq = tempy1
	newBlock.x2.data.iq = tempx2
	newBlock.y2.data.iq = tempy2
			
	return newBlock


def averageObservations2(Observations, timeAvg=1, chanAvg=1):
	"""Given a list of ObservingBlock objects, average the observations 
	together on a per tuning, per polarization basis.  A new ObsevingBlock 
	object is returned that contains the averaged data."""
	
	caBlocks = []
	if chanAvg > 1:
		for Observation in Observations:
			currBlock = Observation
			for i in range(len(Observation.x1.data.iq)/chanAvg):
				currBlock.x1.data.iq[i] = numpy.average( Observation.x1.data.iq[i*chanAvg:(i+1)*chanAvg] )
				currBlock.y1.data.iq[i] = numpy.average( Observation.y1.data.iq[i*chanAvg:(i+1)*chanAvg] )
				currBlock.x2.data.iq[i] = numpy.average( Observation.x2.data.iq[i*chanAvg:(i+1)*chanAvg] )
				currBlock.y2.data.iq[i] = numpy.average( Observation.y2.data.iq[i*chanAvg:(i+1)*chanAvg] )
			currBlock.x1.data.iq = currBlock.x1.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			currBlock.y1.data.iq = currBlock.y1.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			currBlock.x2.data.iq = currBlock.x2.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			currBlock.y2.data.iq = currBlock.y2.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			caBlocks.append(currBlock)
	else:
		caBlocks = Observations
			
	taBlocks = []
	if timeAvg > 1:
		for i in range(len(Observations)/timeAvg):
			taBlocks.append( averageObservations(caBlocks[i*timeAvg:(i+1)*timeAvg]) )
	else:
		taBlocks = caBlocks
			
	return taBlocks
