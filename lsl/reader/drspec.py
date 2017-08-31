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
  * getFFTsPerIntegration - get the number of FFT windows per integration
  * getTransformSize - get the FFT length
  * getIntegrationTime - get the integration time

.. note::
	This reader works with the most current version of the DR spectrometer data
	format as specified in the version 1.7.  To read data created with previous
	versions of the DR spectrometer, use LSL version 0.5.2.
	
.. versionchanged:: 1.0.1
	Added in new functions to help better describe the contents of a DR 
	spectrometer file.
"""

import copy
import numpy

from lsl.common import dp as dp_common
from lsl.reader.drx import filterCodes as drxFilterCodes
from lsl.reader._gofast import readDRSpec
from lsl.reader._gofast import syncError as gsyncError
from lsl.reader._gofast import eofError as geofError
from lsl.reader.errors import syncError, eofError

__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['FrameHeader', 'FrameData', 'Frame', 'readFrame', 'getDataProducts', 'containsLinearData',
		 'containsStokesData', 'getSampleRate', 'getFrameSize', 'getFFTsPerIntegration', 
		 'getTransformSize', 'getIntegrationTime', 'filterCodes', 
		 '__version__', '__revision__', '__all__']

# List of filter codes and their corresponding sample rates in Hz.  
# .. note::
#		These are just the DRX filter codes
filterCodes = drxFilterCodes


class FrameHeader(object):
	"""
	Class that stores the information found in the header of a DR spectrometer/DRX 
	frame.
	"""
	
	def __init__(self, beam=0, format=0, decimation=None, timeOffset=None, nInts=None):
		self.beam = beam
		self.format = format
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
		
	def getDataProducts(self):
		"""
		Return a list of data products contained in the file.
		
		.. versionadded:: 0.6.0
		"""
		
		products = []
		
		# Linear
		if self.format & 0x01:
			products.append('XX')
		if self.format & 0x02:
			products.append('XY')
		if self.format & 0x04:
			products.append('YX')
		if self.format & 0x08:
			products.append('YY')
			
		# Stokes
		if self.format & 0x10:
			products.append('I')
		if self.format & 0x20:
			products.append('Q')
		if self.format & 0x40:
			products.append('U')
		if self.format & 0x80:
			products.append('V')
		
		return products
		
	def containsLinearData(self):
		"""
		Return whether or not the frame contains linear polarization 
		products or not.
		
		.. versionadded:: 0.6.0
		"""
		
		if self.format < 0x10:
			return True
		else:
			return False
			
	def containsStokesData(self):
		"""
		Return whether or not the frame contains Stokes polarization
		parameters or not.
		
		.. versionadded:: 0.6.0
		"""
		
		if self.format < 0x10:
			return False
		else:
			return True
		
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
		
	def getFFTsPerIntegration(self):
		"""
		Return the number of FFT windows per integration.
		
		.. versionadded:: 1.0.1
		"""
		
		return self.nInts


class FrameData(object):
	"""
	Class that stores the information found in the data section of a DR spectrometer/
	DRX frame.
	
	.. versionchanged:: 0.5.3
		Added the saturations field to keep up with saturation counts.
		
	.. versionchanged:: 0.6.0
		The attributes that store the data are not defined until a frame is read in order
		to account for the fact that spectrometer files can store either linear or Stokes
		data.
	"""
	
	def __init__(self, timeTag=None, tuningWords=None, fills=None, errors=None, saturations=None):
		self.timeTag = timeTag
		if tuningWords is None:
			self.tuningWords = [0, 0]
		else:
			self.tuningWords = tuningWords
		if fills is None:
			self.fills = [0, 0, 0, 0]
		else:
			self.fills = fills
		if errors is None:
			self.errors = [0, 0, 0, 0]
		else:
			self.errors = errors
		if saturations is None:
			self.saturations = [0, 0, 0, 0]
		else:
			self.saturations = saturations
			
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
		else:
			raise ValueError("Unknown tuning/polarization combination: '%i'" % which)
			
	def setGain(self, gain):
		"""
		Function to set the gain of the DRX data.
		"""

		self.gain = gain


class Frame(object):
	"""
	Class that stores the information contained within a single DR spectrometer/
	DRX frame.  It's properties are FrameHeader and FrameDataLinear/FrameDataStokes
	objects.
	
	.. versionchanged:: 0.6.0
		By default the data contained with in a frame is normalized by the number of
		fills (header.fills parameter).  For data products that are a function of more
		than one primary input, i.e., XY* or I, the minimum fill of X and Y are used 
		for normalization.
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
		
	def getDataProducts(self):
		"""
		Convenience wrapper for the Frame.FrameHeder.getDataProducts
		function.
		"""
		
		return self.header.getDataProducts()
		
	def containsLinearData(self):
		"""
		Convenience wrapper for the Frame.FrameHeder.containsLinearData
		function.
		"""
		
		return self.header.containsLinearData()
		
	def containsStokesData(self):
		"""
		Convenience wrapper for the Frame.FrameHeder.containsStokesData
		function.
		"""
		
		return self.header.containsStokesData()
		
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
		
	def getFFTsPerIntegration(self):
		"""
		Conveinence wrapper for the Frame.FrameHeader.getFFTsPerIntegration 
		function.
		
		.. versionadded:: 1.0.1
		"""
		
		return self.header.getFFTsPerIntegration()
		
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
		
	def getTransformSize(self):
		"""
		Find out what the transform size is.
		
		.. versionadded:: 1.0.1
		"""
		
		p = self.getDataProducts()[0]
		return getattr(self.data, "%s0" % p, None).size
		
	def getIntegrationTime(self):
		"""
		Return the integration time for data in seconds.
		
		.. versionadded:: 1.0.1
		"""
		
		LFFT = self.getTransformSize()
		srate = self.getSampleRate()
		nInts = self.getFFTsPerIntegration()
		
		return nInts*LFFT/srate
		
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
		
		attrs = self.header.getDataProducts()
		
		for attrBase in attrs:
			for tuning in (0, 1):
				attr = "%s%i" % (attrBase, tuning)
				try:
					temp = getattr(self.data, attr, None) + getattr(y.data, attr, None)
				except TypeError:
					raise RuntimeError("Cannot add %s with %s" % (str(attrs), str(y.header.getDataProducts())))
				except AttributeError:
					temp = getattr(self.data, attr, None) + numpy.float32(y)
				setattr(self.data, attr, temp)
			
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
		
		attrs = self.header.getDataProducts()
		
		for attrBase in attrs:
			for tuning in (0, 1):
				attr = "%s%i" % (attrBase, tuning)
				try:
					temp = getattr(self.data, attr, None) * getattr(y.data, attr, None)
				except TypeError:
					raise RuntimeError("Cannot multiply %s with %s" % (str(attrs), str(y.header.getDataProducts())))
				except AttributeError:
					temp = getattr(self.data, attr, None) * numpy.float32(y)
				setattr(self.data, attr, temp)
			
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
	Function to read in a single DR spectrometer/DRX frame (header+data) and 
	store the contents as a Frame object.
	"""
	
	# New Go Fast! (TM) method
	try:
		newFrame = readDRSpec(filehandle, Frame())
	except gsyncError:
		mark = filehandle.tell()
		raise syncError(location=mark)
	except geofError:
		raise eofError
		
	if Gain is not None:
		newFrame.setGain(Gain)
		
	return newFrame


def getDataProducts(filehandle):
	"""
	Find out the data products contained in the file by looking at a frame.
	"""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Read in one frame
	newFrame = readFrame(filehandle)
	
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the data products
	return newFrame.header.getDataProducts()


def containsLinearData(filehandle):
	"""
	Find out if the file contains linear polarization products or not.
	"""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Read in one frame
	newFrame = readFrame(filehandle)
	
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the verdict
	return newFrame.header.containsLinearData()


def containsStokesData(filehandle):
	"""
	Find out if the file contains Stokes parameters or not.
	"""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Read in one frame
	newFrame = readFrame(filehandle)
	
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Return the verdict
	return newFrame.header.containsStokesData()


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
	Find out what the frame size in a file is at the current file location.
	Returns the frame size in bytes.
	"""
	
	cPos = filehandle.tell()
	frame = readFrame(filehandle)
	nPos = filehandle.tell()
	
	FrameSize = nPos - cPos
	filehandle.seek(cPos)
	
	return FrameSize


def getFFTsPerIntegration(filehandle):
	"""
	Find out what the number of FFT windows per integration is at the 
	current file location.
	
	.. versionadded:: 1.0.1
	"""
	
	cPos = filehandle.tell()
	frame = readFrame(filehandle)
	filehandle.seek(cPos)
	
	return frame.getFFTsPerIntegration()


def getTransformSize(filehandle):
	"""
	Find out what the transform size in a file is at the current file 
	location.  
	"""
	
	cPos = filehandle.tell()
	frame = readFrame(filehandle)
	filehandle.seek(cPos)
	
	return frame.getTransformSize()


def getIntegrationTime(filehandle):
	"""
	Find out what the integration time is at the current file location.
	"""
	
	cPos = filehandle.tell()
	frame = readFrame(filehandle)
	filehandle.seek(cPos)
	
	return frame.getIntegrationTime()
