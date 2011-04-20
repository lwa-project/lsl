# -*- coding: utf-8 -*-

"""Python module to reading in data from both 12-bit and 4-bit TBW files.  
This module defines the following classes for storing the TBW data found in
a file:

Frame
  object that contains all data associated with a particular TBW frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the TBW frame header object and
    * FrameData   - the TBW frame data object.  
  Combined, these two objects contain all of the information found in the 
  original TBW frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the readFrame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.  readFrame 
is designed to work with both 4-bit and 12-bit observations.

For describing the format of data in the file, two function are provided:

getDataBits
  read in the first frame of an open file handle and return whether or not 
  the data is 12 or 4-bit

getFramesPerObs
  read in the first several frames to see how many stands are found in the 
  data.
  .. note::

	This function is a little flaky on TBW data sets that have less 
	than a full complement or 12M (36M) samples.

.. versionchanged:: 0.4.0
	Switched over from pure Python readers to the new C-base Go Fast! readers.
"""

import copy
import numpy

from  lsl.common import dp as dp_common
from _gofast import readTBW
from _gofast import syncError as gsyncError
from _gofast import eofError as geofError
from errors import *

__version__ = '0.6'
__revision__ = '$ Revision: 23 $'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readFrame', 'FrameSize', 'getDataBits', 'getFramesPerObs', '__version__', '__revision__', '__all__']

FrameSize = 1224


class FrameHeader(object):
	"""Class that stores the information found in the header of a TBW 
	frame.  All three fields listed in the DP ICD version H are stored as 
	well as the original binary header data."""

	def __init__(self, frameCount=None, secondsCount=None, tbwID=None):
		self.frameCount = frameCount
		self.secondsCount = secondsCount
		self.tbwID = tbwID

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

		# Why &1023?  Well, from DP ICD revision H, it seems that the stand count 
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
	frame.  Both fields listed in the DP ICD version H are stored."""

	def __init__(self, timeTag=None, samples=400, xy=None):
		self.timeTag = timeTag
		self.xy = xy

	def getTime(self):
		"""Function to convert the time tag from samples since the UNIX epoch
		(UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch."""

		seconds = self.timeTag / dp_common.fS
		
		return seconds


class Frame(object):
	"""Class that stores the information contained within a single TBW 
	frame.  It's properties are FrameHeader and FrameData objects."""

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
			
	def __add__(self, y):
		"""Add the data sections of two frames together or add a number 
		to every element in the data section."""
		
		newFrame = copy.deepcopy(self)
		newFrame += y
		return newFrame
	
	def __iadd__(self, y):
		"""In-place add the data sections of two frames together or add 
		a number to every element in the data section."""
		
		try:
			self.data.xy += y.data.xy
		except AttributeError:
			self.data.xy += y
		return self
		
	def __mul__(self, y):
		"""Multiple the data sections of two frames together or multiply 
		a number to every element in the data section."""

		newFrame = copy.deepcopy(self)
		newFrame *= y		
		return newFrame
			
	def __imul__(self, y):
		"""In-place multiple the data sections of two frames together or 
		multiply a number to every element in the data section."""
		
		try:
			self.data.xy *= y.data.xy
		except AttributeError:
			self.data.xy *= y
		return self
			
	def __eq__(self, y):
		"""Check if the time tags of two frames are equal or if the time
		tag is equal to a particular value."""
		
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
		"""Check if the time tags of two frames are not equal or if the time
		tag is not equal to a particular value."""
		
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
		"""Check if the time tag of the first frame is greater than that of a
		second frame or if the time tag is greater than a particular value."""
		
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
		"""Check if the time tag of the first frame is greater than or equal to 
		that of a second frame or if the time tag is greater than a particular 
		value."""
		
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
		"""Check if the time tag of the first frame is less than that of a
		second frame or if the time tag is greater than a particular value."""
		
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
		"""Check if the time tag of the first frame is less than or equal to 
		that of a second frame or if the time tag is greater than a particular 
		value."""
		
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


def readFrame(filehandle, Verbose=False):
	"""Function to read in a single TBW frame (header+data) and store the 
	contents as a Frame object.  This function wraps readerHeader and 
	readData[(12)|4]."""
	
	# New Go Fast! (TM) method
	try:
		newFrame = readTBW(filehandle, Frame())
	except gsyncError:
		raise syncError
	except geofError:
		raise eofError

	return newFrame


def getDataBits(filehandle):
	"""Find out the number of data bits used in the file be reading in the 
	first frame."""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()

	# Read a frame
	cFrame = readFrame(filehandle)

	# Get the number of bits used to represent the data
	dataBits = cFrame.getDataBits()

	# Return to the place in the file where we started
	filehandle.seek(fhStart)

	return dataBits


def getFramesPerObs(filehandle):
	"""Find out how many frames are present per observation by examining 
	the first frames for what would be 260 stands.  This is done by reading
	two frames and then skipping the next 30,000.
	
	.. note::
		Post-IOC it is probably simpler to adopt a value of the number of 
		frames per observation of 260 rather than try to find it from the
		file.
	"""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()

	idCodes = []
	for i in range(260):
		currentPosition = filehandle.tell()
		try:
			cFrame1 = readFrame(filehandle)
			cFrame2 = readFrame(filehandle)
		except eofError:
			break
		except syncError:
			continue

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
