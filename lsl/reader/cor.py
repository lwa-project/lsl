# -*- coding: utf-8 -*-

"""
Python module to reading in data from COR files.  This module defines the 
following classes for storing the COR data found in a file:

Frame
  object that contains all data associated with a particular COR frame.  The 
  primary consituents of each frame are:
    * FrameHeader - the COR frame header object and
    * FrameData   - the COR frame data object.  
  Combined, these two objects contain all of the information found in the 
  original COR frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the readFrame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object.

.. versionadded:: 1.2.0
"""

import copy
import numpy

from  lsl.common import adp as adp_common
from _gofast import readCOR
from _gofast import syncError as gsyncError
from _gofast import eofError as geofError
from errors import *

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readFrame', 'FrameSize', 'getFramesPerObs', 'getChannelCount',
		 'getBaselineCount', '__version__', '__revision__', '__all__']

FrameSize = 4640


class FrameHeader(object):
	"""
	Class that stores the information found in the header of a COR 
	frame.  All three fields listed in the DP ICD version H are stored as 
	well as the original binary header data.
	"""
	
	def __init__(self, id=None, frameCount=None, secondsCount=None, firstChan=None, gain=None):
		self.id = id
		self.frameCount = frameCount
		self.secondsCount = secondsCount
		self.firstChan = firstChan
		self.gain = gain
		
	def isCOR(self):
		"""
		Function to check if the data is really COR.  Returns True if the 
		data is COR, false otherwise.
		"""
		
		if self.id == 0x02:
			return True
		else:
			return False
			
	def getChannelFreqs(self):
		"""
		Return a numpy.float32 array for the center frequencies, in Hz, of
		each channel in the data.
		"""
		
		return (numpy.arange(144, dtype=numpy.float32)+self.firstChan) * adp_common.fC
		
	def getGain(self):
		"""
		Get the current TBN gain for this frame.
		"""
		
		return self.gain


class FrameData(object):
	"""
	Class that stores the information found in the data section of a COR
	frame.
	"""
	
	def __init__(self, timeTag=None, nAvg=None, stand0=None, stand1=None, vis=None, wgt=None):
		self.timeTag = timeTag
		self.nAvg = nAvg
		self.stand0 = stand0
		self.stand1 = stand1
		self.vis = vis
		self.wgt = wgt
		
	def parseID(self):
		"""
		Return a tuple of the two stands that contribute the this frame.
		"""
		
		return (self.stand0, self.stand1)
		
	def getTime(self):
		"""
		Function to convert the time tag from samples since the UNIX epoch
		(UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch.
		"""
		
		seconds = self.timeTag / adp_common.fS
		
		return seconds
		
	def getIntegrationTime(self):
		"""
		Return the integration time of the visibility in seconds.
		"""
		
		return self.nAvg * adp_common.T2


class Frame(object):
	"""
	Class that stores the information contained within a single COR 
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
		
	def isCOR(self):
		"""
		Convenience wrapper for the Frame.FrameHeader.isCOR function.
		"""
		
		return self.header.isCOR()
		
	def getChannelFreqs(self):
		"""
		Convenience wrapper for the Frame.FrameHeader.getChannelFreqs function.
		"""
		
		return self.header.getChannelFreqs()
		
	def getGain(self):
		"""
		Convenience wrapper for the Frame.FrameHeader.getGain function.
		"""

		return self.header.getGain()
		
	def getTime(self):
		"""
		Convenience wrapper for the Frame.FrameData.getTime function.
		"""
		
		return self.data.getTime()
		
	def parseID(self):
		"""
		Convenience wrapper for the Frame.FrameData.parseID function.
		"""
		
		return self.data.parseID()
		
	def getIntegrationTime(self):
		"""
		Convenience wrapper for the Frame.FrameData.getIntegrationTime
		function.
		"""
		
		return self.data.getIntegrationTime()
		
	def __add__(self, y):
		"""
		Add the data sections of two frames together or add a number 
		to every element in the data section.
		
		.. note::
			In the case where a frame is given the weights are
			ignored.
		"""
		
		newFrame = copy.deepcopy(self)
		newFrame += y
		return newFrame
		
	def __iadd__(self, y):
		"""
		In-place add the data sections of two frames together or add 
		a number to every element in the data section.
		
		.. note::
			In the case where a frame is given the weights are
			ignored.
		"""
		
		try:
			self.data.vis += y.data.vis
		except AttributeError:
			self.data.vis += numpy.complex64(y)
		return self
		
	def __mul__(self, y):
		"""
		Multiple the data sections of two frames together or multiply 
		a number to every element in the data section.
		
		.. note::
			In the case where a frame is given the weights are
			ignored.
		"""
		
		newFrame = copy.deepcopy(self)
		newFrame *= y
		return newFrame
			
	def __imul__(self, y):
		"""
		In-place multiple the data sections of two frames together or 
		multiply a number to every element in the data section.
		
		.. note::
			In the case where a frame is given the weights are
			ignored.
		"""
		
		try:
			self.data.vis *= y.data.vis
		except AttributeError:
			self.data.vis *= numpy.complex64(y)
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


def readFrame(filehandle, Verbose=False):
	"""
	Function to read in a single COR frame (header+data) and store the 
	contents as a Frame object.
	"""
	
	# New Go Fast! (TM) method
	try:
		newFrame = readCOR(filehandle, Frame())
	except gsyncError:
		raise syncError
	except geofError:
		raise eofError
		
	return newFrame


def getFramesPerObs(filehandle):
	"""
	Find out how many frames are present per time stamp by examining the 
	first several thousand COR records.  Return the number of frames per 
	observation.
	"""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Build up the list-of-lists that store the index of the first frequency
	# channel in each frame.
	channelBaselinePairs = []
	for i in range(32896*2):
		try:
			cFrame = readFrame(filehandle)
		except:
			break
			
		chan = cFrame.header.firstChan
		baseline = cFrame.parseID()
		pair = (chan, baseline[0], baseline[1])
		if pair not in channelBaselinePairs:
			channelBaselinePairs.append( pair )
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the number of channel/baseline pairs
	return len(channelBaselinePairs)


def getChannelCount(filehandle):
	"""
	Find out the total number of channels that are present by examining 
	the first several thousand COR records.  Return the number of channels found.
	"""
	
	# Build up the list-of-lists that store the index of the first frequency
	# channel in each frame.
	channels = []
	for i in range(32896*2):
		try:
			cFrame = readFrame(filehandle)
		except:
			break
			
		chan = cFrame.header.firstChan
		if chan not in channels:
			channels.append( chan )
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the number of channels
	return len(channels)


def getBaselineCount(filehandle):
	"""
	Find out the total number of baselines that are present by examining the 
	first several thousand COR records.  Return the number of baselines found.
	observation.
	"""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Build up the list-of-lists that store the index of the first frequency
	# channel in each frame.
	baselines = []
	for i in range(32896*2):
		try:
			cFrame = readFrame(filehandle)
		except:
			break
			
		baseline = cFrame.parseID()
		if baseline not in baselines:
			baselines.append( baseline )
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the number of baselines
	return len(baselines)
