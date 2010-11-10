# -*- coding: utf-8 -*-

"""Python module for creating creating, validating, and writing simulated 
TBW frames to a file."""

import numpy

from lsl.common import dp_common
from lsl.reader import tbw
from errors import *

__version__ = '0.1'
__revision__ = '$ Revision: 5 $'
__all__ = ['SimFrame', 'frame2frame', '__version__', '__revision__', '__all__']


def frame2frame(tbwFrame):
	"""Convert a tbw.Frame object to a raw DP TBW frame."""

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
	rawFrame[12] = (tbwFrame.tbwID>>8) & 255
	rawFrame[13] = tbwFrame.tbwID & 255
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
		x = x.axtype(numpy.uint16)
		y = y.round()
		y = y.clip(-2048, 2047)
		y = y.clip(numpy.uint16)

		rawFrame[24::3] = (x>>4) & 255
		rawFrame[25::3] = (x & 15) | ((y>>8) & 15)
		rawFrame[26::3] = y & 255
	## Data - 4 bit
	if tbwFrame.getDataBits() == 4:
		### Round, clip, and convert to unsigned integers
		x = x.round()
		x = x.clip(-8, 7)
		x = x.axtype(numpy.uint8)
		y = y.round()
		y = y.clip(-8, 7)
		y = y.clip(numpy.uint8)

		rawFrame[24:] = (x<<4) | y
	
	return rawFrame


class SimFrame(tbw.Frame):
	def isValid(self, raiseErrors=False):
		"""Check if simulated TBW frame is valid or not.  Valid frames return 
		True and invalid frames False.  If the `raiseErrors' keyword is set, 
		isValid raises an error when a problem is encountered."""

		# Is the frame actually a TBW frame?
		if not self.header.isTBW:
			if raiseErrors:
				raise invalidFrameType()
			return False

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
		integers.  Returns a numpy array if the frame  is valid, None otherwise."""

		self.isValid(raiseErrors=True)
		return frame2frame(self)

	def writeRawFrame(self, fh):
		"""Write a simulated TBW frame to a filehandle if the frame is valid."""

		self.isValid(raiseErrors=True)
		rawFrame = self.createRawFrame(self)
		rawFrame.tofile(fh)
