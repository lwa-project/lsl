# -*- coding: utf-8 -*-

"""Python module for creating creating, validating, and writing simulated 
TBW frames to a file."""

import numpy

from lsl.reader import tbw
from errors import *

__version__ = '0.2'
__revision__ = '$ Revision: 9 $'
__all__ = ['SimFrame', 'frame2frame', '__version__', '__revision__', '__all__']


def frame2frame(tbwFrame):
	"""Convert a tbw.Frame/tbw.SimFrame object to a raw DP TBW frame."""

	# The raw frame
	rawFrame = numpy.zeros(tbw.FrameSize, dtype=numpy.uint8)

	# Part 1: The header
	## Sync. words (0xDEC0DE5C)
	rawFrame[0] = 0xDE  # 222
	rawFrame[1] = 0xC0  # 192
	rawFrame[2] = 0xDE  # 222
	rawFrame[3] = 0x5C  #  92
	## Frame count
	rawFrame[5] = (tbwFrame.header.frameCount>>16) & 255
	rawFrame[6] = (tbwFrame.header.frameCount>>8) & 255
	rawFrame[7] = tbwFrame.header.frameCount & 255
	## Seconds count
	rawFrame[8] = (tbwFrame.header.secondsCount>>24) & 255
	rawFrame[9] = (tbwFrame.header.secondsCount>>16) & 255
	rawFrame[10] = (tbwFrame.header.secondsCount>>8) & 255
	rawFrame[11] = tbwFrame.header.secondsCount & 255
	## TBW ID
	rawFrame[12] = (tbwFrame.header.tbwID>>8) & 255
	rawFrame[13] = tbwFrame.header.tbwID & 255
	## NB: Next two bytes are unsigned
	
	# Part 2: The data
	## Time tag
	rawFrame[16] = (tbwFrame.data.timeTag>>56) & 255
	rawFrame[17] = (tbwFrame.data.timeTag>>48) & 255
	rawFrame[18] = (tbwFrame.data.timeTag>>40) & 255
	rawFrame[19] = (tbwFrame.data.timeTag>>32) & 255
	rawFrame[20] = (tbwFrame.data.timeTag>>24) & 255
	rawFrame[21] = (tbwFrame.data.timeTag>>16) & 255
	rawFrame[22] = (tbwFrame.data.timeTag>>8) & 255
	rawFrame[23] = tbwFrame.data.timeTag & 255
	## Data - common
	x = numpy.squeeze(tbwFrame.data.xy[0,:])
	y = numpy.squeeze(tbwFrame.data.xy[1,:])
	## Data - 12 bit
	if tbwFrame.getDataBits() == 12:
		### Round, clip, and convert to unsigned integers
		x = x.round()
		x = x.clip(-2048, 2047)
		x = x.astype(numpy.uint16)
		y = y.round()
		y = y.clip(-2048, 2047)
		y = y.astype(numpy.uint16)
		
		rawFrame[24::3] = (x >> 4) & 255
		rawFrame[25::3] = ((x & 15) << 4) | ((y >> 8) & 15)
		rawFrame[26::3] = y & 255
		
	## Data - 4 bit
	if tbwFrame.getDataBits() == 4:
		### Round, clip, and convert to unsigned integers
		x = x.round()
		x = x.clip(-8, 7)
		x = x.astype(numpy.uint8)
		y = y.round()
		y = y.clip(-8, 7)
		y = y.astype(numpy.uint8)

		rawFrame[24:] = (x << 4) | y
	
	return rawFrame


class SimFrame(tbw.Frame):
	"""tbw.SimFrame extends the lsl.reader.tbw.Frame object to yield a method 
	for easily creating DP ICD-compliant raw TBW frames.  Frames created with
	this method can be written to a file via the methods writeRawFrame() function."""

	def __init__(self, stand=None, frameCount=None, dataBits=12, obsTime=None, xy=None):
		"""Given a list of parameters, build a tbw.SimFrame object.  The parameters
		needed are:
		  + stand id (>0 & <259)
		  + which frame number to create
		  + dataBits (12 or 4)
		  + observation time in samples at fS since the epoch
		  + 2-D numpy array representing the frame data for both polarizations.
		Not all of these parameters are needed at initialization of the object and
		the values can be added later.

		.. versionchanged:: 0.3.4
			obsTime now in samples at fS, not seconds
		"""
		
		self.stand = stand
		self.frameCount = frameCount
		self.dataBits = dataBits
		self.obsTime = obsTime
		self.xy = xy
		
		self.header = tbw.FrameHeader()
		self.data = tbw.FrameData()
		
	def __update(self):
		"""Private function to use the object's parameter values to build up 
		a tbw.Frame-like object."""
		
		self.header.frameCount = self.frameCount
		self.header.secondsCount = int(self.obsTime)
		if self.dataBits == 12:
			self.header.tbwID = 32768 | self.stand
		else:
			self.header.tbwID = 32768 | 16384 | self.stand
		
		self.data.timeTag = self.obsTime
		self.data.xy = self.xy
		
	def loadFrame(self, tbwFrame):
		"""Populate the a tbw.SimFrame object with a pre-made frame."""
		
		self.header = tbwFrame.header
		self.data = tbwFrame.data
		
		# Back-fill the class' fields to make sure the object is consistent
		## Header
		self.stand = self.header.parseID()
		self.frameCount = self.header.frameCount
		self.dataBits = self.header.getDataBits()
		## Data
		self.obsTime = self.data.timeTag
		self.xy = self.data.xy
	
	def isValid(self, raiseErrors=False):
		"""Check if simulated TBW frame is valid or not.  Valid frames return 
		True and invalid frames False.  If the 'raiseErrors' keyword is set, 
		isValid() raises an error when a problem with the frame structure is 
		encountered."""

		# Make sure we have the latest values
		self.__update()

		stand = self.parseID()
		# Is the stand number reasonable?
		if stand == 0 or stand > 258:
			if raiseErrors:
				raise invalidStand()
			return False

		# Is there data loaded into frame.data.xy?
		if self.data.xy is None:
			if raiseErrors:
				raise invalidDataSize()
			return False

		# Does the data shape make sense?
		if self.data.xy.shape[0] != 2:
			if raiseErrors:
				raise invalidDataSize()
			return False

		# Does the data length match the data bits in the header?
		if self.getDataBits() == 12 and self.data.xy.shape[1] != 400:
			if raiseErrors:
				raise invalidDataSize()
			return False
		if self.getDataBits() == 4 and self.data.xy.shape[1] != 1200:
			if raiseErrors:
				raise invalidDataSize()
			return False

		# Does the data type make sense?
		if self.data.xy.dtype.kind == 'c':
			if raiseErrors:
				raise invalidDataType()
			return False

		# If we made it this far, it's valid
		return True

	def createRawFrame(self):
		"""Re-express a simulated TBW frame as a numpy array of unsigned 8-bit 
		integers.  Returns a numpy array if the frame  is valid.  If the frame 
		is not ICD-compliant, a errors.baseSimError-type error is raised."""

		# Make sure we have the latest values
		self.__update()

		self.isValid(raiseErrors=True)
		return frame2frame(self)

	def writeRawFrame(self, fh):
		"""Write a simulated TBW frame to a filehandle if the frame is valid.
		If the frame is not ICD-compliant, a errors.baseSimError-type error is 
		raised."""

		# Make sure we have the latest values
		self.__update()

		rawFrame = self.createRawFrame()
		rawFrame.tofile(fh)

	def __str__(self):
		if self.stand is None:
			return "Empty TBW SimFrame object"
		else:
			return "%i-bit TBW SimFrame for stand %i @ time %i" % (self.dataBits, self.stand, self.obsTime)
