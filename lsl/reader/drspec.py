# -*- coding: utf-8 -*-

"""
Python module to read in DR spectrometer data.  This module defines the following 
classes for storing the spectra found in a file:

Frame
  object that contains all data associated with a particular spectrometer frame.  
  The primary constituents of each frame are:
    * FrameHeader - the spectrometer frame header object and
    * FrameData   - the spectral data object.
  Combined, these two objects contain all of the information found in the 
  original spectrometer frame.

The functions defined in this module fall into two class:
 1. convert a frame in a file to a Frame object and
 2. describe the format of the data in the file.

For reading in data, use the readFrame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object
For describing the format of data in the file, three function are provided:
  * getSampleRate - get the sample rate in the file
  * getFrameSize - get the total (header + data) frame size
  * getTransformSize - get the FFT length
"""

import copy
import numpy

from lsl.common import dp as dp_common
from lsl.reader.drx import filterCodes as drxFilterCodes
from _gofast import readDRSpec
from _gofast import syncError as gsyncError
from _gofast import eofError as geofError
from errors import *

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readFrame', 'getSampleRate', 'getFrameSize', 'getTransformSize', 'filterCodes', '__version__', '__revision__', '__all__']

# List of filter codes and their corresponding sample rates in Hz.  
# .. note::
#		These are just the DRX filter codes
filterCodes = drxFilterCodes


class FrameHeader(object):
	"""
	Class that stores the information found in the header of a DRX 
	frame.  All six fields listed in the DP ICD version H are stored as 
	well as the original binary header data.
	"""
	
	def __init__(self, beam=0, decimation=None, timeOffset=None, nInts=None):
		self.beam = beam
		self.decimation = decimation
		self.timeOffset = timeOffset
		
		if nInts is None:
			self.nInts = 0
		else:
			self.nInts = nInts
		
	def parseID(self):
		"""
		Return the beam the frame corresponds to.
		"""
		
		return self.beam
		
	def getSampleRate(self):
		"""
		Return the sample rate of the data in samples/second.
		"""
		
		sampleRate = dp_common.fS / self.decimation
		return sampleRate
		
	def getFilterCode(self):
		"""
		Function to convert the sample rate in Hz to a filter code.
		"""
		
		sampleCodes = {}
		for key,value in filterCodes.iteritems():
			sampleCodes[value] = key

		return sampleCodes[self.getSampleRate()]


class FrameData(object):
	"""
	Class that stores the information found in the data section of a DRX
	frame.  All three fields listed in the DP ICD version H are stored.
	"""

	def __init__(self, timeTag=None, tuningWords=None, fills=None, flags=None, X0=None, Y0=None, X1=None, Y1=None):
		self.timeTag = timeTag
		if tuningWords is None:
			self.tuningWords = [0, 0]
		else:
			self.tuningWords = tuningWords
		if fills is None:
			self.fills = [0, 0, 0, 0]
		else:
			self.fills = fills
		if flags is None:
			self.flags = [0, 0, 0, 0]
		else:
			self.flags = flags
		self.X0 = X0
		self.Y0 = Y0
		self.X1 = X1
		self.Y1 = Y1
		
	def getCentralFreq(self, which=None):
		"""
		Function to set the central frequency of the DRX data in Hz.
		"""

		if which is None:
			return [dp_common.fS * i / 2**32 for i in self.tuningWords]
		elif which == 1:
			return dp_common.fS * self.tuningWords[0] / 2**32
		elif which == 2:
			return dp_common.fS * self.tuningWords[1] / 2**32
			raise ValueError("Unknown tuning/polarization combination: '%i'" % which)


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

	def getSampleRate(self):
		"""
		Convenience wrapper for the Frame.FrameHeader.getSampleRate 
		function.
		"""
		
		return self.header.getSampleRate()
		
	def getFilterCode(self):
		"""
		Convenience wrapper for the Frame.FrameHeader.getFilterCode function.
		"""

		return self.header.getFilterCode()

	def getTime(self):
		"""
		Function to convert the time tag from samples since the UNIX epoch
		(UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch.
		"""

		seconds = (self.data.timeTag - self.header.timeOffset) / dp_common.fS
		
		return seconds
	
	def getCentralFreq(self, which=None):
		"""
		Convenience wrapper for the Frame.FrameData.getCentralFreq function.
		"""

		return self.data.getCentralFreq(which=which)

	def setGain(self, gain):
		"""
		Convenience wrapper for the Frame.FrameData.setGain function.
		"""

		self.data.setGain(gain)

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
			self.data.X0 += y.data.X0
			self.data.Y0 += y.data.Y0
			self.data.X1 += y.data.X1
			self.data.Y1 += y.data.Y1
		except AttributeError:
			self.data.X0 += y
			self.data.Y0 += y
			self.data.X1 += y
			self.data.Y1 += y
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
			self.data.X0 *= y.data.X0
			self.data.Y0 *= y.data.Y0
			self.data.X1 *= y.data.X1
			self.data.Y1 *= y.data.Y1
		except AttributeError:
			self.data.X0 *= y
			self.data.Y0 *= y
			self.data.X1 *= y
			self.data.Y1 *= y
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


def readFrame(filehandle, Gain=None, Verbose=False):
	"""
	Function to read in a single DR spectrometer frame (header+data) and 
	store the contents as a Frame object.
	"""
	
	# New Go Fast! (TM) method
	try:
		newFrame = readDRSpec(filehandle, Frame())
	except gsyncError:
		raise syncError
	except geofError:
		raise eofError
	
	if Gain is not None:
		newFrame.setGain(Gain)

	return newFrame



def getSampleRate(filehandle, nFrames=None, FilterCode=False):
	"""
	Find out what the sampling rate/filter code is from a single observations.  
	By default, the rate in Hz is returned.  However, the corresponding filter 
	code can be returned instead by setting the FilterCode keyword to true.
	"""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()

	# Read in one frame
	newFrame = readFrame(filehandle)
	
	# Return to the place in the file where we started
	filehandle.seek(fhStart)

	if not FilterCode:
		return newFrame.getSampleRate()
	else:
		return newFrame.getFilterCode()


def getFrameSize(filehandle):
	"""
	Find out what the frame size in a file is at the current file location.  Returns 
	the frame size in bytes.
	"""
	
	cPos = filehandle.tell()
	frame = readFrame(filehandle)
	nPos = filehandle.tell()
	
	FrameSize = nPos - cPos
	filehandle.seek(cPos)
	
	return FrameSize


def getTransformSize(filehandle):
	"""
	Find out what the transform size in a file is at the current file location.  
	"""
	
	cPos = filehandle.tell()
	frame = readFrame(filehandle)
	filehandle.seek(cPos)
	
	return frame.data.X0.size
	