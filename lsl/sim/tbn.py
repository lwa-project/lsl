# -*- coding: utf-8 -*-

"""Python module for creating creating, validating, and writing simulated 
TBN frames to a file."""

import numpy

from lsl.common import dp_common
from lsl.reader import tbn

__version__ = '0.1'
__revision__ = '$ Revision: 5 $'
__all__ = ['SimFrame', 'frame2frame', '__version__', '__revision__', '__all__']


def frame2frame(tbnFrame):
	"""Convert a tbn.Frame object to a raw DP TBN frame."""

	# The raw frame
	rawFrame = numpy.zeros(tbn.FrameSize, dtype=numpy.uint8)

	# Part 1: The header
	## Sync. words (0xDEC0DE5C)
	rawFrame[0] = 0xDE  # 222
	rawFrame[1] = 0xC0  # 192
	rawFrame[2] = 0xDE  # 222
	rawFrame[3] = 0x5C  #  92
	## Frame count
	rawFrame[5] = (tbnFrame.header.frameCount>>16) & 255
	rawFrame[6] = (tbnFrame.header.frameCount>>8) & 255
	rawFrame[7] = tbnFrame.header.frameCount & 255
	## Seconds count
	rawFrame[8] = (tbnFrame.header.secondsCount>>24) & 255
	rawFrame[9] = (tbnFrame.header.secondsCount>>16) & 255
	rawFrame[10] = (tbnFrame.header.secondsCount>>8) & 255
	rawFrame[11] = tbnFrame.header.secondsCount & 255
	## TBN ID
	rawFrame[12] = (tbnFrame.tbnID>>8) & 255
	rawFrame[13] = tbnFrame.tbnID & 255
	## NB: Next two bytes are unsigned
	
	# Part 2: The data
	## Time tag
	rawFrame[16] = (tbnFrame.data.timeTag>>56) & 255
	rawFrame[17] = (tbnFrame.data.timeTag>>48) & 255
	rawFrame[18] = (tbnFrame.data.timeTag>>40) & 255
	rawFrame[19] = (tbnFrame.data.timeTag>>32) & 255
	rawFrame[20] = (tbnFrame.data.timeTag>>24) & 255
	rawFrame[21] = (tbnFrame.data.timeTag>>16) & 255
	rawFrame[22] = (tbnFrame.data.timeTag>>8) & 255
	rawFrame[23] = tbnFrame.data.timeTag & 255
	## Data
	i = tbnFrame.data.iq.real
	q = tbnFrame.data.iq.imag
	### Round, clip, and convert to unsigned integers
	i = i.round()
	i = i.clip(-128, 127)
	i = i.astype(numpy.uint8)
	q = q.round()
	q = q.clip(-128, 127)
	q = q.astype(numpy.uint8)
	
	rawFrame[24::2] = i
	rawFrame[23::2] = q
	
	return rawFrame


class SimFrame(tbn.Frame):
	def isValid(self, raiseErrors=False):
		"""Check if simulated TBN frame is valid or not.  Valid frames return 
		True and invalid frames False.  If the `raiseErrors' keyword is set, 
		isValid raises an error when a problem is encountered."""

		# Is the frame actually a TBN frame?
		if not self.header.isTBN:
			if raiseErrors:
				raise invalidFrameType()
			return False

		stand, pol = self.parseID()
		# Is the stand number reasonable?
		if stand == 0 or stand > 258:
			if raiseErrors:
				raise invalidStabd()
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
		if self.data.iq.shape[0] != 512:
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
		"""Re-express a simulated TBN frame as a numpy array of unsigned 8-bit 
		integers.  Returns a numpy array if the frame  is valid."""

		self.isValid(raiseErrors=True)
		return frame2frame(self)

	def writeRawFrame(self, fh):
		"""Write a simulated TBN frame to a filehandle if the frame is valid."""

		self.isValid(raiseErrors=True)
		rawFrame = self.createRawFrame(self)
		rawFrame.tofile(fh)
