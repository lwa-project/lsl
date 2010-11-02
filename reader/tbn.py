#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python module for reading data in from TBN files."""

import os
import sys
import time
import numpy
import array
import struct

from  ..common import dp as dp_common
from errors import *

__version__ = '0.3'
__revision__ = '$ Revision: 5 $'
__all__ = ['TBNFrameHeader', 'TBNFrameData', 'TBNFrame', 'TBNObservingBlock', 'readTBNFrame', 'readTBNBlock', 'getSampleRate', 'getFramesPerObs', 'TBNFrameSize', 'filterCodes', '__version__', '__revision__', '__all__']

TBNFrameSize = 1048

# List of filter codes and their corresponding sample rates in Hz
filterCodes = {1: 1000, 2: 3125, 3: 6250, 4: 12500, 5: 25000, 6: 50000, 7: 100000}


class TBNFrameHeader(object):
	"""Class that stores the information found in the header of a TBW 
	frame.  All three fields listed in the DP IDC version H are stored as 
	well as the original binary header data."""

	def __init__(self, frameCount=None, secondsCount=None, tbnID=None,  raw=None):
		self.frameCount = frameCount
		self.secondsCount = secondsCount
		self.tbnID = tbnID
		self.raw = raw
		
	def isTBN(self):
		"""Function to check if the data is really TBN and not TBW by examining
		the TBN ID field.  Returns True if the data is TBN, false otherwise."""

		mode = (self.tbnID>>15)&1
		if mode == 0:
			return True
		else:
			return False

	def parseID(self):
		"""Function to parse the TBN ID field and return a tuple of the stand 
		number and polarization."""

		if self.tbnID&1023 % 2 == 0:
			stand = (self.tbnID&1023) / 2
			pol = 1
		else:
			stand = (self.tbnID&1023) / 2 + 1
			pol = 0
			
		return (stand, pol)


class TBNFrameData(object):
	"""Class that stores the information found in the data section of a TBN
	frame.  Both fields listed in the DP IDC version H are stored."""

	def __init__self(self, timeTag=None, iq=None):
		self.sampleRate = None
		self.centralFreq = None
		self.gain
		self.timeTag = timeTag
		self.iq = iq

	def getTime(self):
		"""Function to convert the time tag from samples since station midnight to
		seconds since station midnight.  This function needs to dp_common module 
		in order to work."""
		
		return self.timeTag / dp_common.fS

	def getFilterCode(self):
		"""Function to convert the sample rate in Hz to a filter code."""
		if self.sampleRate is None:
			return None
		else:
			sampleCodes = {}
			for key,value in filterCodes.iteritems():
				sampleCodes[value] = key

			return sampleCodes[self.sampleRate]

	def setSampleRate(self, sampleRate):
		"""Function to set the sample rate of the TBN data in Hz."""

		self.sampleRate = sampleRate

	def setCentralFreq(self, centralFreq):
		"""Function to set the central frequency of the TBN data in Hz."""

		self.centralFreq = centralFreq

	def setGain(self, gain):
		"""Function to set the gain of the TBN data."""

		self.gain = gain


class TBNFrame(object):
	"""Class that stores the information contained within a single TBN 
	frame.  It's properties are TBNFrameHeader and TBNFrameData objects."""

	def __init__(self, header=TBNFrameHeader(), data=TBNFrameData()):
		self.header = header
		self.data = data

	def parseID(self):
		"""Convenience wrapper for the TBNFrame.TBNFrameHeader.parseID function."""
		
		return self.header.parseID()

	def getTime(self):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.getTime function."""
		
		return self.data.getTime()

	def getFilterCode(self):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.getFilterCode function."""

		return self.data.getFilterCode()

	def setSampleRate(self, sampleRate):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.setSampleRate function."""

		self.data.setSampleRate(sampleRate)

	def setCentralFreq(self, centralFreq):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.setCentralFreq function."""

		self.data.setCentralFreq(centralFreq)

	def setGain(self, gain):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.setGain function."""

		self.data.setGain(gain)


class TBNObservingBlock(object):
	def __init__(self, x=[], y=[]):
		self.x = x
		self.y = y

	def getTime(self):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.getTime function."""
		
		return self.x[0].data.getTime()

	def getFilterCode(self):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.getFilterCode function."""

		return self.x[0].data.getFilterCode()

	def setSampleRate(self, sampleRate):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.setSampleRate function."""

		for i in len(self.x):
			self.x[i].data.setSampleRate(sampleRate)
		for i in len(self.y):
			self.y[i].data.setSampleRate(sampleRate)

	def setCentralFreq(self, centralFreq):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.setCentralFreq function."""

		for i in len(self.x):
			self.x[i].data.setCentralFreq(centralFreq)
		for i in len(self.y):
			self.y[i].data.setCentralFreq(centralFreq)

	def setGain(self, gain):
		"""Convenience wrapper for the TBNFrame.TBNFrameData.setGain function."""

		for i in len(self.x):
			self.x[i].data.setGain(gain)
		for i in len(self.y):
			self.y[i].data.setGain(gain)


def __readTBNHeader(filehandle, Verbose=False):
	"""Private function to read in a TBN header.  Returns a TBNFrameHeader object."""

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
		tbnID, junk = struct.unpack(">HH", s)
	except IOError:
		raise eofError()
	except struct.error:
		raise eofError()

	if sync1 != 92 or sync2 != 222 or sync3 != 192 or sync4 != 222:
		raise syncError(sync1=sync1, sync2=sync2, sync3=sync3, sync4=sync4)

	newHeader = TBNFrameHeader()
	newHeader.frameCount = frameCount
	newHeader.secondsCount = secondsCount[0]
	newHeader.tbnID = tbnID
	newHeader.raw = rawHeader

	if Verbose:
		stand, pol = newHeader.parseID()
		print "Header: ", tbnID
		print "  Stand: ", stand
		print "  Polarization: ", pol

	return newHeader


def __readTBNData(filehandle):
	"""Private function to read in a TBN frame data section.  Returns a 
	TBNFrameData object."""

	try:
		s = filehandle.read(8)
		timeTag = struct.unpack(">Q", s)
	except IOError:
		raise eofError()

	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=1024)
	rawData = rawData.astype(numpy.int8)
	data = numpy.zeros(512, dtype=numpy.csingle)
	if rawData.shape[0] < 2*data.shape[0]:
		raise numpyError()

	# The data are signed, so apply the two-complement rule to generate 
	# the negative values
	negativeValues = numpy.where( rawData >= 128 )
	rawData[negativeValues] -= 256

	data.real = rawData[0::2]
	data.imag = rawData[1::2]
	#print data.real.min(), data.real.max(), data.imag.min(), data.imag.max()
	
	newData = TBNFrameData()
	newData.timeTag = timeTag[0]
	newData.iq = data

	return newData


def readTBNFrame(filehandle, SampleRate=None, CentralFreq=None, Gain=None, Verbose=False):
	"""Function to read in a single TBN frame (header+data) and store the 
	contents as a TBNFrame object.  This function wraps readerHeader and 
	readData."""

	try:
		hdr = __readTBNHeader(filehandle, Verbose=Verbose)
	except syncError, err:
		# Why?  If we run into a sync error here, then the following frame is invalid.  
		# Thus, we need to skip over this frame be advancing the file pointer 8+1000 B 
		currPos = filehandle.tell()
		frameEnd = currPos + TBNFrameSize - 16
		filehandle.seek(frameEnd)
		raise err

	dat = __readTBNData(filehandle)
	
	newFrame = TBNFrame()
	newFrame.header = hdr
	newFrame.data = dat
	
	newFrame.setSampleRate(SampleRate)
	newFrame.setCentralFreq(CentralFreq)
	newFrame.setGain(Gain)

	return newFrame


def readTBNBlock(filehandle, nFrames=512, SampleRate=None, CentralFreq=None, Gain=None, Verbose=False):

	frames = []
	for i in range(0, nFrames):
		try:
			frame = readTBNFrame(filehandle, SampleRate=SampleRate, CentralFreq=CentralFreq, Gain=Gain, Verbose=Verbose)
		except baseReaderError:
			frame = None
		frames.append( frame )

	return TBNObservingBlock(x=frames[0::2], y=frames[1::2])
	

def getSampleRate(filehandle, nFrames=None, FilterCode=False):
	"""Find out what the sampling rate/filter code is from consecutive sets of 
	observations.  By default, the rate in Hz is returned.  However, the 
	corresponding filter code can be returned instead by setting the FilterCode
	keyword to true."""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Go back to the beginning...
	filehandle.seek(0)

	if nFrames is None:
		nFrames = 512
	nFrames = 2*nFrames

	# Build up the list-of-lists that store ID codes and loop through 512
	# frames.  In each case, parse pull the TBN ID, extract the stand 
	# number, and append the stand number to the relevant polarization array 
	# if it is not already there.
	frames = {}
	for i in range(nFrames):
		try:
			cFrame = readTBNFrame(filehandle)
		except eofError:
			break
		except syncError:
			continue
		except numpyError:
			break
		
		stand, pol = cFrame.parseID()
		key = 2*stand + pol
		if key not in frames.keys():
			frames[key] = []
		frames[key].append(cFrame)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)

	# Any key with complete data will work for this, so pick the first key with two
	# valid frames
	keyCount = 0
	frame1 = None
	frame2 = None
	while frame1 is None and frame2 is None:
		validKey = (frames.keys())[keyCount]
		frame1 = frames[validKey][0]
		frame2 = frames[validKey][1]

		keyCount = keyCount + 1

	# Now that we have two valid frames that follow one another in time, load in their
	# time tags and calculate the sampling rate.  Since the time tags are based off f_S
	# @ 196 MSPS, and each frame contains 512 samples, the sampleing rate is:
	#  f_S / <difference in time tags per 512 samples>
	time1 = frame1.data.timeTag
	time2 = frame2.data.timeTag
	rate = dp_common.fS / (abs( time2 - time1 ) / 512)

	if not FilterCode:
		return rate
	else:
		sampleCodes = {}
		for key,value in filterCodes.iteritems():
			sampleCodes[value] = key

		return sampleCodes[rate]


def getFramesPerObs(filehandle):
	"""Find out how many frames are present per observation by examining 
	the first 512 TBN frames.  Return the number of frames per observations 
	as a two-	element tuple, one for each polarization."""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Go back to the beginning...
	filehandle.seek(0)

	# Build up the list-of-lists that store ID codes and loop through 512
	# frames.  In each case, parse pull the TBN ID, extract the stand 
	# number, and append the stand number to the relevant polarization array 
	# if it is not already there.
	idCodes = [[], []]
	for i in range(512):
		try:
			cFrame = readTBNFrame(filehandle)
		except eofError:
			break
		except syncError:
			continue
		except numpyError:
			break
		
		cID, cPol = cFrame.header.parseID()
		if cID not in idCodes[cPol]:
			idCodes[cPol].append(cID)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Get the length of each beam list and return them as a tuple
	return (len(idCodes[0]), len(idCodes[1]))


def main(args):
	from ..writer import tsfits
	import matplotlib.pyplot as plt

	nSamples = os.path.getsize(args[0]) / TBNFrameSize
	print "Samples in file: ", nSamples
	fh = open(args[0], "rb", buffering=TBNFrameSize)
	
	test = readTBNFrame(fh)
	print "TBN Data:  %s" % test.header.isTBN()
	if not test.header.isTBN():
		raise notTBNError()
	fh.seek(0)

	nFpO = getFramesPerObs(fh)
	print "Samples per observations: %i in x pol., %i in y pol." % nFpO
	nFpO = nFpO[0] + nFpO[1]

	sampleRate = getSampleRate(fh, nFrames=nFpO)
	print "Filter code is: %i" % getSampleRate(fh, nFrames=nFpO, FilterCode=True)
	print "Sampling rate is: %i Hz" % sampleRate

	tStart = time.time()

	# Create a new FITS file with the name 'tbw.fits'
	fitsFile = tsfits.TBN('tbn-tsfits-test.fits')

	nSamples = 340000

	blocks = []
	count = {}
	syncCount = 0
	masterCount = 0
	for i in range(nSamples):
		try:
			cFrame = readTBNBlock(fh, nFrames=nFpO, SampleRate=sampleRate)
		except eofError:
			break
		except syncError:
			syncCount = syncCount + 1
			continue
		except numpyError:
			break
		blocks.append(cFrame)

		for frame in blocks[-1].x:
			if frame is None:
				print "sync error"
				continue
			stand, pol = frame.parseID()
			if stand not in count.keys():
				count[stand] = 0
			count[stand] = count[stand] + 1

		for xFrame, yFrame in zip(cFrame.x, cFrame.y):
			if xFrame is not None:
				fitsFile.addStandData(xFrame)
			if yFrame is not None:
				fitsFile.addStandData(yFrame)

		masterCount = masterCount + 1
		#print cFrame.header.parseID(),cFrame.data.timeTag

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (nFpO*nSamples, (tEnd-tStart), nFpO*nSamples/(tEnd-tStart))

	fh.close()
	fitsFile.info()

	# Summary information about the file that was just read in
	print "Summary:"
	for stand in sorted(count.keys()):
		print "Stand: %2i, Frames: %5i" % (stand, count[stand])
	print "Sync Errors: %5i" % syncCount

	#i = 0
	#x = numpy.empty((nSamples,256))
	#y = numpy.empty((nSamples,256))
	#for block in blocks:
		#for j in range(256):
			#tX = block.x[j].data.iq
			#tX = tX - tX.mean()
			#tY = block.y[j].data.iq
			#tY = tY - tY.mean()
			#x[i,j] = (tX*tX.conj()).real.mean()
			#y[i,j] = (tY*tY.conj()).real.mean()
		#i = i + 1
	
	#fig = plt.figure()
	#ax1 = fig.add_subplot(311)
	#ax1.imshow(numpy.log10(x)*10.0, origin='lower', label='X')
	
	#ax2 = fig.add_subplot(312)
	#ax2.imshow(numpy.log10(x)*10.0, origin='lower', label='Y')

	#ax3 = fig.add_subplot(313)
	#ax3.plot(x[:,0])
	##for i in range(nFpO/2):
		##ax1.plot(numpy.log10(x[i,:])*10.0, label='X')
		##ax1.plot(numpy.log10(y[i,:])*10.0, label='Y')
	##ax1.legend(loc=0)
	#plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])

