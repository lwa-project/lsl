# -*- coding: utf-8 -*-

"""Python module for creating creating, validating, and writing simulated 
DRX frames to a file."""

import numpy

from lsl.common import dp as dp_common
from lsl.reader import drx
from errors import *

__version__ = '0.1'
__revision__ = '$ Revision: 6 $'
__all__ = ['SimFrame', 'frame2frame', '__version__', '__revision__', '__all__']


def frame2frame(drxFrame):
	"""Convert a drx.Frame object to a raw DP DRX frame."""

	# The raw frame
	rawFrame = numpy.zeros(drx.FrameSize, dtype=numpy.uint8)

	# Part 1: The header
	## Sync. words (0xDEC0DE5C)
	rawFrame[0] = 0xDE  # 222
	rawFrame[1] = 0xC0  # 192
	rawFrame[2] = 0xDE  # 222
	rawFrame[3] = 0x5C  #  92
	## DRX ID
	rawFrame[4] = drxFrame.header.drxID
	## Frame count
	rawFrame[5] = (drxFrame.header.frameCount>>16) & 255
	rawFrame[6] = (drxFrame.header.frameCount>>8) & 255
	rawFrame[7] = drxFrame.header.frameCount & 255
	## Seconds count
	rawFrame[8] = (drxFrame.header.timeTag>>24) & 255
	rawFrame[9] = (drxFrame.header.timeTag>>16) & 255
	rawFrame[10] = (drxFrame.header.timeTag>>8) & 255
	rawFrame[11] = drxFrame.header.timeTag & 255
	## Decimation
	rawFrame[12] = (drxFrame.decimation>>8) & 255
	rawFrame[13] = drxFrame.decimation & 255
	## Time offset
	rawFrame[14] = (drxFrame.timeOffset>>8) & 255
	rawFrame[15] = drxFrame.timeOffset & 255

	# Part 2: The data
	## Time tag
	rawFrame[16] = (drxFrame.data.timeTag>>56) & 255
	rawFrame[17] = (drxFrame.data.timeTag>>48) & 255
	rawFrame[18] = (drxFrame.data.timeTag>>40) & 255
	rawFrame[19] = (drxFrame.data.timeTag>>32) & 255
	rawFrame[20] = (drxFrame.data.timeTag>>24) & 255
	rawFrame[21] = (drxFrame.data.timeTag>>16) & 255
	rawFrame[22] = (drxFrame.data.timeTag>>8) & 255
	rawFrame[23] = drxFrame.data.timeTag & 255
	## Flags
	rawFrame[24] = (drxFrame.data.flags>>56) & 255
	rawFrame[25] = (drxFrame.data.flags>>48) & 255
	rawFrame[26] = (drxFrame.data.flags>>40) & 255
	rawFrame[27] = (drxFrame.data.flags>>32) & 255
	rawFrame[28] = (drxFrame.data.flags>>24) & 255
	rawFrame[29] = (drxFrame.data.flags>>16) & 255
	rawFrame[30] = (drxFrame.data.flags>>8) & 255
	rawFrame[31] = drxFrame.data.flags & 255
	## Data
	i = drxFrame.data.iq.real
	q = drxFrame.data.iq.imag
	### Round, clip, and convert to unsigned integers
	i = i.round()
	i = i.clip(-8, 7)
	i = i.astype(numpy.uint8)
	q = q.round()
	q = q.clip(-8, 7)
	q = q.astype(numpy.uint8)

	rawFrame[32:] = (i << 4) | q

	return rawFrame


class SimFrame(drx.Frame):
	"""drx.SimFrame extends the lsl.reader.tbn.Frame object to yield a method 
	for easily creating DP ICD-compliant raw DRX frames.  Frames created with
	this method can be written to a file via the methods writeRawFrame() function."""

	def __init__(self, beam=None, tune=None, pol=None, filterCode=None, timeOffset=None, frameCount=None, obsTime=None, flags=None, iq=None):
		"""Given a list of parameters, build a drx.SimFrame object.  The parameters
		needed are:
		  + beam id (>0 & <5)
		  + tunning (1 or 2)
		  + polarization (0 for x, or 1 for y)
		  + which filter code the data corresponds to (>0 & <8)
		  + what time offset in units of f_S to use
		  + which frame number to create
		  + observation time in seconds since the epoch
		  + what flags are set on the data
		  + 1-D numpy array representing the frame I/Q (complex) data
		Not all of these parameters are needed at initialization of the object and
		the values can be added later."""
		
		self.beam = beam
		self.tune = tune
		self.pol = pol
		self.filterCode = filterCode
		self.timeOffset = timeOffset
		self.frameCount = frameCount
		self.obsTime = obsTime
		self.flags = flags
		self.iq = iq
		
		self.header = drx.FrameHeader()
		self.data = drx.FrameData()
		
	def __update(self):
		"""Private function to use the object's parameter values to build up 
		a drx.Frame-like object."""
		
		self.header.frameCount = self.frameCount
		self.header.decimation = int(dp_common.fS / drx.filterCodes[self.filterCode])
		self.header.timeOffset = self.timeOffset
		self.header.drxID = (self.beam & 7) | ((self.tune & 7) << 3) | ((self.pol & 1) << 7)
		
		self.data.timeTag = long(self.obsTime * dp_common.fS)
		self.data.flags = self.flags
		self.data.iq = self.iq
		
	def loadFrame(self, drxFrame):
		"""Populate the a drx.SimFrame object with a pre-made frame."""
		
		self.header = drxFrame.header
		self.data = drxFrame.data
		
		# Back-fill the class' fields to make sure the object is consistent
		## Header
		self.beam = self.header.parseID()[0]
		self.tune = self.header.parseID()[1]
		self.pol = self.header.parseID()[2]
		self.frameCount = self.header.frameCount
		self.filterCode = int(dp_common.fS / self.header.decimations)
		self.timeOffset = self.header.timeOffset
		## Data
		self.obsTime = self.data.timeTag / dp_common.fS
		self.flags = self.data.flags
		self.iq = self.data.iq
	
	def isValid(self):
		"""Check if simulated DRX frame is valid or not.  Valid frames return 
		True and invalid frames False.  If the `raiseErrors' keyword is set, 
		isValid raises an error when a problem is encountered."""

		# Make sure we have the latest values
		self.__update()

		# Is the time offset reasonable?
		if self.header.timeOffset >= dp_common.fS:
			return False

		beam, tune, pol = self.parseID()
		# Is the beam number reasonable?
		if beam not in [1, 2, 3, 4]:
			if raiseErrors:
				raise invalidBeam()
			return False

		# Is the tunning number reasonable?
		if tune not in [1, 2]:
			if raiseErrors:
				raise invalidTune()
			return False

		# Is the polarization reasonable?
		if pol not in [0, 1]:
			if raiseErrors:
				raise invalidPol()
			return False

		# Is there data loaded into frame.data.iq?
		if self.data.iq is None:
			if raiseErrors:
				raise invalidDataSize()
			return False

		# Does the data length make sense?
		if self.data.iq.shape[0] != 4096:
			if raiseErrors:
				raise invalidDataSize()
			return False

		# Does the data type make sense?
		if self.data.iq.dtype.kind != 'c':
			if raiseErrors:
				raise invalidDataType()
			return False

		# If we made it this far, it's valid
		return True

	def createRawFrame(self):
		"""Re-express a simulated DRX frame as a numpy array of unsigned 8-bit 
		integers.  Returns a numpy array if the frame is valid.  If the frame 
		is not ICD-compliant, a errors.baseSimError-type error is raised."""

		# Make sure we have the latest values
		self.__update()

		self.isValid(raiseErrors=True)
		return frame2frame(self)

	def writeRawFrame(self, fh):
		"""Write a simulated DRX frame to a filehandle if the frame is valid.
		If the frame is not ICD-compliant, a errors.baseSimError-type error 
		is raised."""

		# Make sure we have the latest values
		self.__update()

		rawFrame = self.createRawFrame()
		rawFrame.tofile(fh)

	def __str__(self):
		if self.stand is None:
			return "Empty DRX SimFrame object"
		else:
			return "TBN SimFrame for beam %i, tunning %i, pol. %i @ time %i" % (self.beam, self.tune, self.pol, self.obsTime)
