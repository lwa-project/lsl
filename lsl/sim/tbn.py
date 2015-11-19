# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import print_function
import sys
if sys.version_info > (3,):
	xrange = range
	long = int

"""
Python module for creating creating, validating, and writing simulated 
TBN frames to a file.
"""

import numpy

from lsl.common.dp import fS
from lsl.reader import tbn
from errors import *

__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['SimFrame', 'frame2frame', '__version__', '__revision__', '__all__']


def frame2frame(tbnFrame):
	"""
	Convert a :class:`lsl.reader.tbn.Frame` object to a raw DP TBN frame.
	"""

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
	## Tuning word
	rawFrame[8] = (tbnFrame.header.tuningWord>>24) & 255
	rawFrame[9] = (tbnFrame.header.tuningWord>>16) & 255
	rawFrame[10] = (tbnFrame.header.tuningWord>>8) & 255
	rawFrame[11] = tbnFrame.header.tuningWord & 255
	## TBN ID
	rawFrame[12] = (tbnFrame.header.tbnID>>8) & 255
	rawFrame[13] = tbnFrame.header.tbnID & 255
	## Gain
	rawFrame[14] = (tbnFrame.header.gain>>8) & 255
	rawFrame[15] = tbnFrame.header.gain & 255
	
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
	rawFrame[25::2] = q
	
	return rawFrame


class SimFrame(tbn.Frame):
	"""
	tbn.SimFrame extends the :class:`lsl.reader.tbn.Frame` object to yield a method 
	for easily creating DP ICD-compliant raw TBN frames.  Frames created with
	this method can be written to a file via the methods writeRawFrame() function.
	"""

	def __init__(self, stand=None, pol=None, centralFreq=None, gain=None, frameCount=None, obsTime=None, iq=None):
		"""
		Given a list of parameters, build a tbn.SimFrame object.  The parameters
		needed are:
		  * stand id (>0 & <259)
		  * polarization (0 for x, or 1 for y)
		  * central frequency of tuning in (Hz)
		  * TBN gain
		  * which frame number to create
		  * observation time in samples at fS since the epoch
		  * 1-D numpy array representing the frame I/Q (complex) data
		  
		Not all of these parameters are needed at initialization of the object and
		the values can be added later.

		.. versionchanged:: 0.3.4
			obsTime now in samples at fS, not seconds

		.. versionchanged:: 0.5.0
			Added support for ECR 11 TBN headers
		"""
		
		self.stand = stand
		self.pol = pol
		self.freq = centralFreq
		self.gain = gain
		self.frameCount = frameCount
		self.obsTime = obsTime
		self.iq = iq
		
		self.header = tbn.FrameHeader()
		self.data = tbn.FrameData()
		
	def __update(self):
		"""
		Private function to use the object's parameter values to build up 
		a tbn.Frame-like object.
		"""
		
		self.header.frameCount = self.frameCount
		self.header.tuningWord = long( round(self.freq/fS*2**32) )
		self.header.tbnID = 2*(self.stand-1) + self.pol + 1
		self.header.gain = self.gain
		
		self.data.timeTag = self.obsTime
		self.data.iq = self.iq
	
	def loadFrame(self, tbnFrame):
		"""
		Populate the a tbn.SimFrame object with a pre-made frame.
		"""
		
		self.header = tbnFrame.header
		self.data = tbnFrame.data
		
		# Back-fill the class' fields to make sure the object is consistent
		## Header
		self.stand = self.header.parseID()[0]
		self.pol = self.header.parseID()[1]
		self.freq = self.header.getCentralFreq()
		self.gain = self.header.getGain()
		self.frameCount = self.header.frameCount
		## Data
		self.obsTime = self.data.timeTag
		self.iq = self.data.iq
	
	def isValid(self, raiseErrors=False):
		"""
		Check if simulated TBN frame is valid or not.  Valid frames return 
		True and invalid frames False.  If the 'raiseErrors' keyword is set, 
		isValid raises an error when a problem is encountered.
		"""

		# Make sure we have the latest values
		self.__update()

		stand, pol = self.parseID()
		# Is the stand number reasonable?
		if stand == 0 or stand > 258:
			if raiseErrors:
				raise invalidStand()
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
		"""
		Re-express a simulated TBN frame as a numpy array of unsigned 8-bit 
		integers.  Returns a numpy array if the frame  is valid.  If the frame 
		is not ICD-compliant, a errors.baseSimError-type error is raised.
		"""

		# Make sure we have the latest values
		self.__update()

		self.isValid(raiseErrors=True)
		return frame2frame(self)

	def writeRawFrame(self, fh):
		"""
		Write a simulated TBN frame to a filehandle if the frame is valid.
		If the frame is not ICD-compliant, a errors.baseSimError-type error 
		is raised.
		"""

		rawFrame = self.createRawFrame()
		rawFrame.tofile(fh)

	def __str__(self):
		if self.stand is None:
			return "Empty TBN SimFrame object"
		else:
			return "TBN SimFrame for stand %i, pol. %i @ time %i" % (self.stand, self.pol, self.obsTime)
