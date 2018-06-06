# -*- coding: utf-8 -*-
"""
LWA Development Primitives - A set of utilities that should make developing 
new analysis software easier.  These functions wrap the nitty gritty of the 
file reading and unpacking behind Python objects.

Data format objects included are:
  * TBWFile
  * TBNFile
  * DRXFile
  * DRSpecFile
  * TBFFile
  
Also included are the LWA1DataFile, LWASVDataFile, and LWADataFile functions 
that take a filename and try to determine the correct data format object to
use.

.. versionchanged:: 1.2.0
	Added support for LWA-SV ADP data
"""

import os
import numpy
import warnings
from scipy.stats import norm
from collections import deque

from lsl.common.dp import fS
from lsl.common.adp import fC
from lsl.reader import tbw, tbn, drx, drspec, tbf, cor, errors
from lsl.reader.buffer import TBNFrameBuffer, DRXFrameBuffer

__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['TBWFile', 'TBNFile', 'DRXFile', 'DRSpecFile', 'TBFFile', 'LWA1DataFile', 
		 'LWASVDataFile', 'LWADataFile', 
		 '__version__', '__revision__', '__all__']


class LDPFileBase(object):
	"""
	Class to make it easy to interface with raw LWA1 data files and DR spectrometer
	data files.
	"""
	
	def __init__(self, filename=None, fh=None, ignoreTimeTagErrors=False):
		# Make sure that we are given either a filename or an open file handle
		if filename is None and fh is None:
			raise RuntimeError("Must specify either a filename or open file instance")
			
		# Store a valid file handle and mark the object as ready
		if fh is None:
			self.filename = filename
			self.fh = open(filename, 'rb')
		else:
			self.filename = fh.name
			if fh.mode.find('b') == -1:
				fh.close()
				fh = open(self.filename, 'rb')
			self.fh = fh
			
		# Set whether or not reading errors are fatal
		self.ignoreTimeTagErrors = ignoreTimeTagErrors
		
		# Ready the file
		self._readyFile()
		
		# Describe the contents of the file
		self.description = {}
		self._describeFile()
		
	def __getattr__(self, name):
		## Try to access the attribute as a real attribute
		try:
			return super(LDPFileBase, self).__getattr__(name)
		except AttributeError:
			pass
			
		## Try to access the attribute via the 'getInfo' method
		try:
			return self.getInfo(name)
		except ValueError:
			raise AttributeError("'%s' object has no attribute '%s'" % (type(self).__name__, name))
			
	def _readyFile(self):
		"""
		Method for finding the start of valid data.  This will be over-
		ridden in the format-specific subclasses.
		"""
		
		raise NotImplementedError
		
	def _describeFile(self):
		"""
		Method for describing the contents of a file using.  This will 
		be over-ridden in the format-specific subclasses.
		"""
		
		raise NotImplementedError
		
	def getInfo(self, key=None):
		"""
		Retrieve metadata about the file.  This will return a dictionary 
		of values if no key is specified.
		"""
		
		if key is None:
			return self.description
		else:
			try:
				return self.description[key]
			except KeyError:
				raise ValueError("Unknown key '%s'" % key)
				
	def getRemainingFrameCount(self):
		"""
		Return the number of frames left in the file.
		"""
		
		return (self.description['size'] - self.fh.tell()) / self.description['FrameSize']
		
	def reset(self):
		"""
		Reset the file to the beginning.
		"""
		
		self.fh.seek(0)
		
		# Ready the file
		self._readyFile()
		
		# Reset any buffers
		if getattr(self, "buffer", None) is not None:
			self.buffer.flush()
			
		# Reset the timetag checker
		if getattr(self, "_timetag", None) is not None:
			self._timetag = None
			
		# Describe the contents of the file
		self.description = {}
		self._describeFile()
		
	def close(self):
		"""
		Close the file.
		"""
		
		self.fh.close()
		
	def offset(self, *args, **kwds):
		"""
		Offset into the data.
		"""
		
		raise NotImplementedError
		
	def readFrame(self):
		"""
		Read a single frame from the data.
		"""
		
		raise NotImplementedError
		
	def read(self, *args, **kwds):
		"""
		Read a certain amount of time from the data.
		"""
		
		raise NotImplementedError
		
	def estimateLevels(self, *args, **kwds):
		"""
		Estimate the standard deviation of the data.
		"""
		
		raise NotImplementedError


class TBWFile(LDPFileBase):
	"""
	Class to make it easy to interface with a TBW file.  Methods defined for this class are:
	  * getInfo - Get information about the file's contents
	  * getRemainingFrameCount - Get the number of frames remaining in the file
	  * readFrame - Read and return a single `lsl.reader.tbw.Frame` instance
	  * read - Read in the capture and return it as a numpy array
	"""
	
	def _readyFile(self):
		"""
		Find the start of valid TBW data.  This function:
		1) Aligns on the first valid Mark 5C frame and
		2) Skips over any TBN frames at the beginning of the file.
		"""
		
		# Align on the start of a Mark5C packet
		while True:
			try:
				junkFrame = tbw.readFrame(self.fh)
				break
			except errors.syncError:
				self.fh.seek(-tbw.FrameSize+1, 1)
				
		# Jump over any TBN data in the file
		while not junkFrame.header.isTBW():
			try:
				junkFrame = tbw.readFrame(self.fh)
			except errors.syncError:
				## If we reached this then we are probably in an old TBW file that has
				## a bunch of TBN frames at the beginning.  We need to seek backwards,
				## realign on the sync word, and read forwards again.
				
				## Jump back a TBW frame
				self.fh.seek(-tbw.FrameSize, 1)
				
				## Find the sync word again
				while True:
					try:
						junkFrame = tbn.readFrame(self.fh)
						break
					except errors.syncError:
						self.fh.seek(-tbn.FrameSize+1, 1)
						
				## Find the end of the TBN data
				while True:
					try:
						junkFrame = tbn.readFrame(self.fh)
					except errors.syncError:
						break
				self.fh.seek(-2*tbn.FrameSize, 1)
				junkFrame = tbw.readFrame(self.fh)
		self.fh.seek(-tbw.FrameSize, 1)
		
		return True
		
	def _describeFile(self):
		"""
		Describe the TBW file.
		"""
		
		junkFrame = self.readFrame()
		self.fh.seek(-tbw.FrameSize, 1)
		
		# Basic file information
		filesize = os.fstat(self.fh.fileno()).st_size
		nFramesFile = filesize / tbw.FrameSize
		srate = 196e6
		bits = junkFrame.getDataBits()
		start = junkFrame.getTime()
		startRaw = junkFrame.data.timeTag
		
		# Trick to figure out how many antennas are in a file and the "real" 
		# start time.  For details of why this needs to be done, see the read()
		# function below.
		idsFound = []
		timesFound = []
		filePosRef = self.fh.tell()
		while True:
			try:
				for i in xrange(26):
					frame = tbw.readFrame(self.fh)
					while not frame.header.isTBW():
						frame = tbw.readFrame(self.fh)
					stand = frame.parseID()
					if stand not in idsFound:
						idsFound.append(stand)
					if frame.header.frameCount < 1000:
						timesFound.append( (frame.header.frameCount-1, frame.data.timeTag) )
				self.fh.seek(tbw.FrameSize*(30000-26), 1)
			except:
				break
		self.fh.seek(filePosRef)
		
		# What is that start time again?
		startTimeTag = None
		for fc,tt in timesFound:
			tt = tt - fc*(1200 if bits == 4 else 400)
			if startTimeTag is None or tt < startTimeTag:
				startTimeTag = tt
		start = startTimeTag / fS
		startRaw = startTimeTag
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 'FrameSize': tbw.FrameSize,
						'sampleRate': srate, 'dataBits': bits, 'nAntenna': 2*len(idsFound), 
						'tStart': start, 'tStartSamples': startRaw}
						
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.tbw.Frame` instance.
		"""
		
		frame = tbw.readFrame(self.fh)
		while not frame.header.isTBW():
			frame = tbw.readFrame(self.fh)
			
		return frame
		
	def read(self, duration=None, timeInSamples=False):
		"""
		Read and return the entire TBW capture.  This function returns 
		a three-element tuple with elements of:
		  0) the actual duration of data read in, 
		  1) the time tag for the first sample, and
		  2) a 2-D Numpy array of data.
		  
		The time tag is returned as seconds since the UNIX epoch by default.
		However, the time tags can be returns as samples at fS if the 
		timeInSamples keyword is set.
		
		The sorting order of the output data array is by 
		digitizer number - 1.
		
		.. note::
			Setting the 'duration' keyword has no effect on the read 
			process because the entire capture is always read in.
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Get the data frame size
		dataSize = 400
		if self.description['dataBits'] == 4:
			dataSize = 1200
			
		# Find out how many frames to work with at a time
		nFrames = int(30000)
		
		# Initialize the time variables
		# Explination:
		#   This is needed to work out what the "real" start time is of the 
		#   capture due to buffering in the data recorder.  What can happen 
		#   is that the last ~4 MB of a previous capture could be stuck in 
		#   the data recoder's buffer and that the buffer won't get dumped 
		#   until the next capture is launch.  Thus, you can end up in a 
		#   situation where the first few valid TBW frames in a file are from 
		#   the previous capture.
		#   
		#   To get around this we use the frame count-correction time tag of 
		#   the lowest frame number found.  This skips over the trailing edge of 
		#   the previous capture (which should have a high frame count) while
		#   allowing the code to deal with files that may be missing the first
		#   frame from the first board to send a frame.
		setTime = None
		setTimeRef = 1000
		
		# Initialize the output data array
		data = numpy.zeros((self.description['nAntenna'], nFrames*dataSize), dtype=numpy.int16)
		
		# Read in the next frame and anticipate any problems that could occur
		i = 0
		while i < ((self.description['nAntenna']/2)*nFrames):
			try:
				cFrame = tbw.readFrame(self.fh)
			except errors.eofError:
				break
			except errors.syncError:
				continue
				
			if not cFrame.header.isTBW():
				continue
				
			stand = cFrame.header.parseID()
			aStandX = 2*(stand-1) + 0
			aStandY = 2*(stand-1) + 1
			
			if cFrame.header.frameCount < setTimeRef:
				newSetTime = cFrame.data.timeTag - (cFrame.header.frameCount-1)*dataSize
				if setTime is None or cFrame.data.timeTag < setTime:
					setTime = newSetTime
					setTimeRef = cFrame.header.frameCount
					
			try:
				cnt = cFrame.header.frameCount - 1
				data[aStandX, cnt*dataSize:(cnt+1)*dataSize] = cFrame.data.xy[0,:]
				data[aStandY, cnt*dataSize:(cnt+1)*dataSize] = cFrame.data.xy[1,:]
				
				i += 1
			except ValueError:
				pass
				
		# Deal with the time if we don't want it in samples
		if not timeInSamples:
			setTime = setTime / fS
			
		# Calculate the duration
		duration = data.shape[1]/self.getInfo('sampleRate')
		
		return duration, setTime, data


class TBNFile(LDPFileBase):
	"""
	Class to make it easy to interface with a TBN file.  Methods defined for this class are:
	  * getInfo - Get information about the file's contents
	  * getRemainingFrameCount - Get the number of frames remaining in the file
	  * offset - Offset a specified number of seconds into the file
	  * readFrame - Read and return a single `lsl.reader.tbn.Frame` instance
	  * read - Read a chunk of data in and return it as a numpy array
	  * estimateLevels - Estimate the n-sigma level for the absolute value of the voltages 
	"""
	
	def _readyFile(self):
		"""
		Given an open file handle, find the start of valid TBN data.  This
		function:
		1) Aligns on the first valid Mark 5C frame and
		2) Skips over any TBW frames at the beginning of the file.
		"""
		
		# Align on the start of a Mark5C packet
		while True:
			try:
				junkFrame = tbn.readFrame(self.fh)
				break
			except errors.syncError:
				self.fh.seek(-tbn.FrameSize+1, 1)
				
		# Jump over any TBN data in the file
		while not junkFrame.header.isTBN():
			junkFrame = tbn.readFrame(self.fh)
		self.fh.seek(-tbn.FrameSize, 1)
		
		return True
		
	def _describeFile(self):
		"""
		Describe the TBN file and initialize the frame circular buffer.
		"""
		
		filesize = os.fstat(self.fh.fileno()).st_size
		nFramesFile = filesize / tbn.FrameSize
		framesPerObsX, framesPerObsY = tbn.getFramesPerObs(self.fh)
		srate =  tbn.getSampleRate(self.fh, nFrames=((framesPerObsX+framesPerObsY)*3))
		bits = 8
		
		junkFrame = self.readFrame()
		self.fh.seek(-tbn.FrameSize, 1)
		tuning1 = junkFrame.getCentralFreq()
		start = junkFrame.getTime()
		startRaw = junkFrame.data.timeTag
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 'FrameSize': tbn.FrameSize,
						'nAntenna': framesPerObsX+framesPerObsY, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start, 'tStartSamples': startRaw, 'freq1': tuning1}
						
		# Initialize the buffer as part of the description process
		pols = []
		if framesPerObsX != 0:
			pols.append(0)
		if framesPerObsY != 0:
			pols.append(1)
		nAntenna = framesPerObsX + framesPerObsY
		
		self.buffer = TBNFrameBuffer(stands=range(1,nAntenna/len(pols)+1), pols=pols)
		
	def offset(self, offset):
		"""
		Offset a specified number of seconds in an open TBN file.  This function 
		returns the exact offset time.
		
		.. note::
			The offset provided by this function is relatively crude due to the
			structure of TBN files.
		"""
		
		frameOffset = int(offset * self.description['sampleRate'] / 512 * self.description['nAntenna'])
		frameOffset = int(1.0 * frameOffset / self.description['nAntenna']) * self.description['nAntenna']
		self.fh.seek(frameOffset*tbn.FrameSize)
		
		# Update the file metadata
		self._describeFile()
		
		# Reset the buffer
		self.buffer.flush()
		
		# Reset the timetag checker
		self._timetag = None
		
		return 1.0 * frameOffset / self.description['nAntenna'] * 512 / self.description['sampleRate']
		
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.tbn.Frame` instance.
		"""
		
		# Reset the buffer
		if getattr(self, "buffer", None) is not None:
			self.buffer.flush()
			
		# Reset the timetag checker
		self._timetag = None
		
		return tbn.readFrame(self.fh)
		
	def read(self, duration, timeInSamples=False):
		"""
		Read in a chunk (in seconds) of TBN data.  This function returns 
		a three-element tuple with elements of:
		  0) the actual duration of data read in, 
		  1) the time tag for the first sample, and
		  2) a 2-D Numpy array of data.
		  
		The time tag is returned as seconds since the UNIX epoch by default.
		However, the time tags can be returns as samples at fS if the 
		timeInSamples keyword is set.
		
		The sorting order of the output data array is by 
		digitizer number - 1.
		""" 
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Covert the sample rate to an expected timetag skip
		timetagSkip = int(512 / self.description['sampleRate'] * fS)
		
		# Setup the counter variables:  frame count and time tag count
		if getattr(self, "_timetag", None) is None:
			self._timetag = 0
			junkFrame = tbn.readFrame(self.fh)
			self._timetag = junkFrame.data.timeTag - timetagSkip
			self.fh.seek(-tbn.FrameSize, 1)
			
		# Find out how many frames to read in
		frameCount = int(round(1.0 * duration * self.description['sampleRate'] / 512))
		frameCount = frameCount if frameCount else 1
		duration = frameCount * 512 / self.description['sampleRate']
		
		nFrameSets = 0
		eofFound = False
		setTime = None
		count = [0 for i in xrange(self.description['nAntenna'])]
		data = numpy.zeros((self.description['nAntenna'], frameCount*512), dtype=numpy.complex64)
		while True:
			if eofFound or nFrameSets == frameCount:
				break
			
			cFrames = deque()
			for i in xrange(self.description['nAntenna']/2):
				try:
					cFrames.append( tbn.readFrame(self.fh, Verbose=False) )
				except errors.eofError:
					eofFound = True
					self.buffer.append(cFrames)
					break
				except errors.syncError:
					continue
				
			self.buffer.append(cFrames)
			cFrames = self.buffer.get()
			
			# Continue adding frames if nothing comes out.
			if cFrames is None:
				continue
				
			# If something comes out, add it to the data array
			cTimetag = cFrames[0].data.timeTag
			if cTimetag != self._timetag+timetagSkip:
				actStep = cTimetag - self._timetag
				if self.ignoreTimeTagErrors:
					warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
				else:
					raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
			self._timetag = cFrames[0].data.timeTag
			
			for cFrame in cFrames:
				stand,pol = cFrame.header.parseID()
				aStand = 2*(stand-1)+pol
				
				if setTime is None:
					if timeInSamples:
						setTime = cFrame.data.timeTag
					else:
						setTime = cFrame.getTime()
						
				data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
				count[aStand] += 1
			nFrameSets += 1
			
		# If we've hit the end of the file and haven't read in enough frames, 
		# flush the buffer
		if eofFound or nFrameSets != frameCount:
			for cFrames in self.buffer.flush():
				cTimetag = cFrames[0].data.timeTag
				if cTimetag != self._timetag+timetagSkip:
					actStep = cTimetag - self._timetag
					if self.ignoreTimeTagErrors:
						warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
					else:
						raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
				self._timetag = cFrames[0].data.timeTag
				
				for cFrame in cFrames:
					stand,pol = cFrame.header.parseID()
					aStand = 2*(stand-1)+pol
					
					if setTime is None:
						if timeInSamples:
							setTime = cFrame.data.timeTag
						else:
							setTime = cFrame.getTime()
						
					data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
					count[aStand] += 1
				nFrameSets += 1
				
				if nFrameSets == frameCount:
					break
					
		# Adjust the duration to account for all of the things that could 
		# have gone wrong while reading the data
		duration = nFrameSets * 512 / self.description['sampleRate']
		
		return duration, setTime, data
		
	def estimateLevels(self, nFrames=100, Sigma=5.0):
		"""
		Estimate the n-sigma level for the absolute value of the voltages.  
		Returns a list with indicies that are the digitizer numbers minus one.
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Go!
		count = {}
		for i in xrange(self.description['nAntenna']):
			count[i] = 0
		data = numpy.zeros((self.description['nAntenna'], nFrames*512))
		for i in xrange(nFrames):
			for j in xrange(self.description['nAntenna']):
				# Read in the next frame and anticipate any problems that could occur
				try:
					cFrame = tbn.readFrame(self.fh, Verbose=False)
				except errors.eofError:
					break
				except errors.syncError:
					continue
					
				s,p = cFrame.parseID()
				aStand = 2*(s-1) + p
				
				data[aStand, count[aStand]*512:(count[aStand]+1)*512] = numpy.abs( cFrame.data.iq )
				count[aStand] +=  1
		self.fh.seek(-tbn.FrameSize*self.description['nAntenna']*nFrames, 1)
		
		# Statistics
		rv = norm()
		frac = rv.cdf(Sigma) - rv.cdf(-Sigma)
		index = int(round(data.shape[1]*frac))
		if index == data.shape[1]:
			index = data.shape[1] - 1
		
		levels = [0 for i in xrange(self.description['nAntenna'])]
		for i in xrange(self.description['nAntenna']):
			data2 = sorted(data[i,:])
			levels[i] = data2[index]
		
		return levels


class DRXFile(LDPFileBase):
	"""
	Class to make it easy to interface with a DRX file.  Methods defined for this class are:
	  * getInfo - Get information about the file's contents
	  * getRemainingFrameCount - Get the number of frames remaining in the file
	  * offset - Offset a specified number of seconds into the file
	  * readFrame - Read and return a single `lsl.reader.drx.Frame` instance
	  * read - Read a chunk of data in and return it as a numpy array
	  * estimateLevels - Estimate the n-sigma level for the absolute value of the voltages 
	 """
	
	def _readyFile(self):
		"""
		Given an open file handle, find the start of valid DRX data.  This function:
		1) aligns on the first valid Mark 5C frame and
		2) skips over frames with a decimation of zero. 
		3) aligns the tuning/polarization timetags
		"""
		
		# Align on the start of a Mark5C packet...
		while True:
			try:
				junkFrame = drx.readFrame(self.fh)
				try:
					# ... that has a valid decimation
					srate = junkFrame.getSampleRate()
					break
				except ZeroDivisionError:
					pass
			except errors.syncError:
				self.fh.seek(-drx.FrameSize+1, 1)
				
		self.fh.seek(-drx.FrameSize, 1)
		
		# Line up the time tags for the various tunings/polarizations
		ids = []
		timeTags = []
		for i in xrange(32):
			junkFrame = drx.readFrame(self.fh)
			b,t,p = junkFrame.parseID()
			id = (t,p)
			if id not in ids:
				ids.append(id)
			timeTags.append(junkFrame.data.timeTag)
		self.fh.seek(-32*drx.FrameSize, 1)
		
		return True
		
	def _describeFile(self):
		"""
		Describe the DRX file.
		"""
		
		filesize = os.fstat(self.fh.fileno()).st_size
		nFramesFile = filesize / drx.FrameSize
		beams = drx.getBeamCount(self.fh)
		tunepols = drx.getFramesPerObs(self.fh)
		tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
		beampols = tunepol
		bits = 4
		
		beams = []
		tunes = []
		pols = []
		tuning1 = 0.0
		tuning2 = 0.0
		for i in xrange(4):
			junkFrame = self.readFrame()
			b,t,p = junkFrame.parseID()
			srate = junkFrame.getSampleRate()
			if b not in beams:
				beams.append(b)
			if t not in tunes:
				tunes.append(t)
			if p not in pols:
				pols.append(p)
				
			if t == 1:
				tuning1 = junkFrame.getCentralFreq()
			else:
				tuning2 = junkFrame.getCentralFreq()
				
			if i == 0:
				start = junkFrame.getTime()
				startRaw = junkFrame.data.timeTag - junkFrame.header.timeOffset
		self.fh.seek(-drx.FrameSize*4, 1)
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 'FrameSize': drx.FrameSize,
						'beampols': beampols, 'beam': b, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start, 'tStartSamples': startRaw, 'freq1': tuning1, 'freq2': tuning2}
						
		# Initialize the buffer as part of the description process
		self.buffer = DRXFrameBuffer(beams=beams, tunes=tunes, pols=pols)
		
	def offset(self, offset):
		"""
		Offset a specified number of seconds in an open DRX file.  This function 
		returns the exact offset time.
		"""
		
		junkFrame = drx.readFrame(self.fh)
		self.fh.seek(-drx.FrameSize, 1)
		
		# Get the initial time, sample rate, and beampols
		t0 = junkFrame.getTime()
		sampleRate = junkFrame.getSampleRate()
		beampols = drx.getFramesPerObs(self.fh)
		beampols = reduce(int.__add__, beampols)
		
		# Offset in frames for beampols beam/tuning/pol. sets
		ioffset = int(offset * sampleRate / 4096 * beampols)
		ioffset = int(1.0 * ioffset / beampols) * beampols
		self.fh.seek(ioffset*drx.FrameSize, 1)
		
		# Iterate on the offsets until we reach the right point in the file.  This
		# is needed to deal with files that start with only one tuning and/or a 
		# different sample rate.  
		while True:
			junkFrame = drx.readFrame(self.fh)
			self.fh.seek(-drx.FrameSize, 1)
			
			## Figure out where in the file we are and what the current tuning/sample 
			## rate is
			t1 = junkFrame.getTime()
			sampleRate = junkFrame.getSampleRate()
			beampols = drx.getFramesPerObs(self.fh)
			beampols = reduce(int.__add__, beampols)
			
			## See how far off the current frame is from the target
			tDiff = t1 - (t0 + offset)
			
			## Half that to come up with a new seek parameter
			tCorr   = -tDiff / 2.0
			cOffset = int(tCorr * sampleRate / 4096 * beampols)
			cOffset = int(1.0 * cOffset / beampols) * beampols
			ioffset += cOffset
			
			## If the offset is zero, we are done.  Otherwise, apply the offset
			## and check the location in the file again/
			if cOffset is 0:
				break
			try:
				self.fh.seek(cOffset*drx.FrameSize, 1)
			except IOError:
				warnings.warn("Could not find the correct offset, giving up", RuntimeWarning)
				break
				
		# Update the file metadata
		self._describeFile()
		
		# Reset the buffer
		self.buffer.flush()
		
		# Zero out the time tag checker
		self._timetag = None
		
		return t1 - t0
		
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.drx.Frame` instance.
		"""
		
		# Reset the buffer
		if getattr(self, "buffer", None) is not None:
			self.buffer.flush()
			
		# Zero out the time tag checker
		self._timetagSkip = None
		self._timetag = None
		
		return drx.readFrame(self.fh)
		
	def read(self, duration, timeInSamples=False):
		"""
		Given an open DRX file and an amount of data to read in in seconds, read 
		in the data and return a three-element tuple of the actual duration read 
		in, the time for the first sample, and the data as numpy 
		array.
		
		..note::
			This function always returns a 2-D array with the first dimension
			holding four elements.  These elements contain, in order:
			  * Tuning 1, polarization X
			  * Tuning 1, polarization Y
			  * Tuning 2, polarization X
			  * Tuning 2, polarization Y
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Covert the sample rate to an expected timetag skip
		if getattr(self, "_timetagSkip", None) is None:
			self._timetagSkip = int(4096 / self.description['sampleRate'] * fS)
		
		# Setup the counter variables:  frame count and time tag count
		if getattr(self, "_timetag", None) is None:
			self._timetag = {0:0, 1:0, 2:0, 3:0}
			
		# Find out how many frames to read in
		frameCount = int(round(1.0 * duration * self.description['sampleRate'] / 4096))
		frameCount = frameCount if frameCount else 1
		duration = frameCount * 4096 / self.description['sampleRate']
		
		# Setup the output arrays
		setTime = None
		data = numpy.zeros((4,frameCount*4096), dtype=numpy.complex64)
		
		# Go!
		nFrameSets = 0
		eofFound = False
		count = {0:0, 1:0, 2:0, 3:0}
		while True:
			if eofFound or nFrameSets == frameCount:
				break
				
			cFrames = deque()
			for i in xrange(self.description['beampols']):
				try:
					cFrames.append( drx.readFrame(self.fh, Verbose=False) )
				except errors.eofError:
					eofFound = True
					self.buffer.append(cFrames)
					break
				except errors.syncError:
					continue
					
			self.buffer.append(cFrames)
			cFrames = self.buffer.get()
			
			# Continue adding frames if nothing comes out.
			if cFrames is None:
				continue
				
			# If something comes out, add it to the data array
			for cFrame in cFrames:
				b,t,p = cFrame.parseID()
				aStand = 2*(t-1) + p
				cTimetag = cFrame.data.timeTag
				if self._timetag[aStand] == 0:
					self._timetag[aStand] = cTimetag - self._timetagSkip
				if cTimetag != self._timetag[aStand]+self._timetagSkip:
					actStep = cTimetag - self._timetag[aStand]
					if self.ignoreTimeTagErrors:
						warnings.warn("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep), RuntimeWarning)
					else:
						raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep))
						
				if setTime is None:
					if timeInSamples:
						setTime = cFrame.data.timeTag - cFrame.header.timeOffset
					else:
						setTime = cFrame.getTime()
						
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
				count[aStand] +=  1
				self._timetag[aStand] = cTimetag
			nFrameSets += 1
			
		# If we've hit the end of the file and haven't read in enough frames, 
		# flush the buffer
		if eofFound or nFrameSets != frameCount:
			for cFrames in self.buffer.flush():
				for cFrame in cFrames:
					b,t,p = cFrame.parseID()
					aStand = 2*(t-1) + p
					cTimetag = cFrame.data.timeTag
					if self._timetag[aStand] == 0:
						self._timetag[aStand] = cTimetag - self._timetagSkip
					if cTimetag != self._timetag[aStand]+self._timetagSkip:
						actStep = cTimetag - self._timetag[aStand]
						if self.ignoreTimeTagErrors:
							warnings.warn("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep), RuntimeWarning)
						else:
							raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep))
							
					if setTime is None:
						if timeInSamples:
							setTime = cFrame.data.timeTag - cFrame.header.timeOffset
						else:
							setTime = cFrame.getTime()
							
					data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
					count[aStand] +=  1
					self._timetag[aStand] = cTimetag
				nFrameSets += 1
				
				if nFrameSets == frameCount:
					break
					
		# Adjust the duration to account for all of the things that could 
		# have gone wrong while reading the data
		duration = nFrameSets * 4096 / self.description['sampleRate']
			
		return duration, setTime, data
		
	def estimateLevels(self, nFrames=100, Sigma=5.0):
		"""
		Estimate the n-sigma level for the absolute value of the voltages.  
		Returns a list with indicies corresponding to:
		0)  Tuning 1, X pol.
		1)  Tuning 1, Y pol.
		2)  Tuning 2, X pol.
		3)  Tuning 2, Y pol.
		
		..note::
			The returned list always has four items, regardless of whether 
			or not the input DRX file has one or two tunings.
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
		
		# Sample the data
		count = {0:0, 1:0, 2:0, 3:0}
		data = numpy.zeros((4, nFrames*4096))
		for i in xrange(nFrames):
			for j in xrange(self.description['beampols']):
				# Read in the next frame and anticipate any problems that could occur
				try:
					cFrame = drx.readFrame(self.fh, Verbose=False)
				except errors.eofError:
					break
				except errors.syncError:
					continue
					
				b,t,p = cFrame.parseID()
				aStand = 2*(t-1) + p
				
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = numpy.abs( cFrame.data.iq )
				count[aStand] +=  1
		self.fh.seek(-drx.FrameSize*self.description['beampols']*nFrames, 1)
		
		# Statistics
		rv = norm()
		frac = rv.cdf(Sigma) - rv.cdf(-Sigma)
		index = int(round(data.shape[1]*frac))
		if index == data.shape[1]:
			index = data.shape[1] - 1
		
		levels = [0, 0, 0, 0]
		for i in xrange(4):
			data2 = sorted(data[i,:])
			levels[i] = data2[index]
			
		return levels


class DRSpecFile(LDPFileBase):
	"""
	Class to make it easy to interface with a DR Spectrometer file.  
	Methods defined for this class are:
	  * getInfo - Get information about the file's contents
	  * getRemainingFrameCount - Get the number of frames remaining in the file
	  * offset - Offset a specified number of seconds into the file
	  * readFrame - Read and return a single `lsl.reader.drspec.Frame` instance
	  * read - Read a chunk of data in and return it as a numpy array
	 """
	
	def _readyFile(self):
		"""
		Ready the DRSpec file.
		"""
		
		return True
		
	def _describeFile(self):
		"""
		Describe the DRSpec file.
		"""
		
		filesize = os.fstat(self.fh.fileno()).st_size
		FrameSize = drspec.getFrameSize(self.fh)
		nFramesFile = filesize / FrameSize
		LFFT = drspec.getTransformSize(self.fh)
		junkFrame = drspec.readFrame(self.fh)
		self.fh.seek(-FrameSize, 1)
		
		bits = 32
		beam = junkFrame.parseID()
		beampols = 4
		srate = junkFrame.getSampleRate()
		nInt = junkFrame.header.nInts
		tInt = nInt*LFFT/srate
		start = junkFrame.getTime()
		startRaw = junkFrame.data.timeTag - junkFrame.header.timeOffset
		tuning1, tuning2 = junkFrame.getCentralFreq()
		prod = junkFrame.getDataProducts()
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 'FrameSize': FrameSize, 
						'beampols': beampols, 'beam': beam, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start, 'tStartSamples': startRaw, 'freq1': tuning1, 'freq2': tuning2, 
						'nInt': nInt, 'tInt': tInt, 'LFFT': LFFT, 
						'nProducts': len(prod), 'dataProducts': prod}
						
	def offset(self, offset):
		"""
		Offset a specified number of seconds in an open DR spectrometer file.  This 
		function returns the exact offset time.
		"""
		
		# Gather some basic information and read in the first frame
		junkFrame = drspec.readFrame(self.fh)
		self.fh.seek(-self.description['FrameSize'], 1)
		
		# Get the initial time, sample rate, and integration time
		t0 = junkFrame.getTime()
		
		# Offset in frames for beampols beam/tuning/pol. sets
		ioffset = int(round(offset / self.description['tInt']))
		self.fh.seek(ioffset*self.description['FrameSize'], 1)
		
		# Iterate on the offsets until we reach the right point in the file.  This
		# is needed to deal with files that start with only one tuning and/or a 
		# different sample rate.  
		while True:
			junkFrame = drspec.readFrame(self.fh)
			self.fh.seek(-self.description['FrameSize'], 1)
			
			## Figure out where in the file we are and what the current tuning/sample 
			## rate is
			t1 = junkFrame.getTime()
			sampleRate = junkFrame.getSampleRate()
			LFFT = junkFrame.getTransformSize()
			tInt = junkFrame.header.nInts*LFFT/sampleRate
			
			## See how far off the current frame is from the target
			tDiff = t1 - (t0 + offset)
			
			## Half that to come up with a new seek parameter
			tCorr   = -tDiff / 2.0
			cOffset = int(round(tCorr / tInt))
			ioffset += cOffset
			
			## If the offset is zero, we are done.  Otherwise, apply the offset
			## and check the location in the file again/
			if cOffset is 0:
				break
			self.fh.seek(cOffset*self.description['FrameSize'], 1)
			
		# Update the file metadata
		self._describeFile()
		
		# Zero out the timetag checker
		self._timetag = None
		
		return t1 - t0
		
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.drspec.Frame` instance.
		"""
		
		# Update the timetag checker
		self._timetag = None
		
		return drspec.readFrame(self.fh)
		
	def read(self, duration, timeInSamples=False):
		"""
		Given an open DR spectrometer file and an amount of data read in in 
		seconds, read in the data and return a three-element tuple of the actual 
		duration read in, the times at the beginning of each stream, and the 
		data as numpy array.
		
		..note::
			This function always returns a 3-D array with the first dimension
			indexing over data product, the second over time and the third over
			frequency channel.
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Covert the sample rate to an expected timetag skip
		timetagSkip = self.description['tInt']
		
		# Setup the counter variables:  frame count and time tag count
		count = 0
		if getattr(self, "_timetag", None) is None:
			self._timetag = 0
			junkFrame = drspec.readFrame(self.fh)
			self._timetag = junkFrame.getTime() - timetagSkip
			self.fh.seek(-self.description['FrameSize'], 1)
		
		# Find out how many frames to read in
		frameCount = int(round(1.0 * duration / self.description['tInt']))
		frameCount = frameCount if frameCount else 1
		duration = frameCount * self.description['tInt']
		
		# Setup the output arrays
		data = numpy.zeros((2*self.description['nProducts'],frameCount,self.description['LFFT']), dtype=numpy.float32)
		
		# Go!
		nFrameSets = 0
		setTime = None
		for i in xrange(frameCount):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drspec.readFrame(self.fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				continue
				
			cTimetag = cFrame.getTime()
			if cTimetag > self._timetag + 1.001*timetagSkip:
				actStep = cTimetag - self._timetag
				if self.ignoreTimeTagErrors:
					warnings.warn("Invalid timetag skip encountered, expected %i but found %i" % (timetagSkip, actStep), RuntimeWarning)
				else:
					raise RuntimeError("Invalid timetag skip encountered, expected %i but found %i" % (timetagSkip, actStep))
					
			if setTime is None:
				if timeInSamples:
					setTime = cFrame.data.timeTag - cFrame.header.timeOffset
				else:
					setTime = cFrame.getTime()
					
			for j,p in enumerate(self.description['dataProducts']):
				data[j+0,                             count, :] = getattr(cFrame.data, '%s0' % p, None)
				data[j+self.description['nProducts'], count, :] = getattr(cFrame.data, '%s1' % p, None)
			count +=  1
			self._timetag = cTimetag
			nFrameSets += 1
			
		# Adjust the duration to account for all of the things that could 
		# have gone wrong while reading the data
		duration = nFrameSets * self.description['tInt']
		
		return duration, setTime, data


def LWA1DataFile(filename=None, fh=None, ignoreTimeTagErrors=False):
	"""
	Wrapper around the various LWA1-related classes defined here that takes
	a file, determines the data type, and initializes and returns the 
	appropriate LDP class.
	"""
	
	# Open the file as appropriate
	if fh is None:
		fh = open(filename, 'rb')
	else:
		filename = fh.name
		if fh.mode.find('b') == -1:
			fh.close()
			fh = open(filename, 'rb')
			
	# Read a bit of data to try to find the right type
	for mode in (drx, tbn, tbw, drspec):
		## Set if we find a valid frame marker
		foundMatch = False
		## Set if we can read more than one valid successfully
		foundMode = False
		
		## Sort out the frame size.  This is tricky because DR spectrometer files
		## have frames of different sizes depending on the mode
		if mode == drspec:
			try:
				mfs = drspec.getFrameSize(fh)
			except:
				mfs = 0
		else:
			mfs = mode.FrameSize
			
		## Loop over the frame size to try and find what looks like valid data.  If
		## is is found, set 'foundMatch' to True.
		for i in xrange(mfs):
			try:
				junkFrame = mode.readFrame(fh)
				foundMatch = True
				break
			except errors.eofError:
				break
			except errors.syncError:
				fh.seek(-mfs+1, 1)
				
		## Did we strike upon a valid frame?
		if foundMatch:
			### Is so, we now need to try and read more frames to make sure we have 
			### the correct type of file
			fh.seek(-mfs, 1)
			
			try:
				for i in xrange(2):
					junkFrame = mode.readFrame(fh)
				foundMode = True
			except errors.eofError:
				break
			except errors.syncError:
				### Reset for the next mode...
				fh.seek(0)
		else:
			### Reset for the next mode...
			fh.seek(0)
			
		## Did we read more than one valid frame?
		if foundMode:
			break
			
	# There is an ambiguity that can arise for TBW data such that it *looks* 
	# like TBN.  If the identified mode is TBN, skip halfway into the file and 
	# verify that it is still TBN.  We also need to catch the LWA-SV DRX vs.
	# TBF ambiguity since we could have been given an LWA-SV file by accident
	if mode in (tbn, drx):
		## Sort out the frame size
		omfs = mode.FrameSize
		
		## Seek half-way in
		nFrames = os.path.getsize(filename)/omfs
		fh.seek(nFrames/2*omfs)
		
		## Read a bit of data to try to find the right type
		for mode in (tbn, tbw, drx):
			### Set if we find a valid frame marker
			foundMatch = False
			### Set if we can read more than one valid successfully
			foundMode = False
			
			### Sort out the frame size.
			mfs = mode.FrameSize
			
			### Loop over the frame size to try and find what looks like valid data.  If
			### is is found, set 'foundMatch' to True.
			for i in xrange(mfs):
				try:
					junkFrame = mode.readFrame(fh)
					foundMatch = True
					break
				except errors.eofError:
					break
				except errors.syncError:
					fh.seek(-mfs+1, 1)
					
			### Did we strike upon a valid frame?
			if foundMatch:
				#### Is so, we now need to try and read more frames to make sure we have 
				#### the correct type of file
				fh.seek(-mfs, 1)
				
				try:
					for i in xrange(4):
						junkFrame = mode.readFrame(fh)
					foundMode = True
				except errors.syncError:
					#### Reset for the next mode...
					fh.seek(nFrames/2*omfs)
			else:
				#### Reset for the next mode...
				fh.seek(nFrames/2*omfs)
				
			### Did we read more than one valid frame?
			if foundMode:
				break
				
	fh.close()
	
	# Raise an error if nothing is found
	if not foundMode:
		raise RuntimeError("File '%s' does not appear to be a valid LWA1 data file" % filename)
		
	# Otherwise, build and return the correct LDPFileBase sub-class
	if mode == drx:
		ldpInstance = DRXFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
	elif mode == tbn:
		ldpInstance = TBNFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
	elif mode == tbw:
		ldpInstance = TBWFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
	else:
		ldpInstance = DRSpecFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
		
	# Done
	return ldpInstance


class TBFFile(LDPFileBase):
	"""
	Class to make it easy to interface with a TBF file.  Methods defined for this class are:
	  * getInfo - Get information about the file's contents
	  * getRemainingFrameCount - Get the number of frames remaining in the file
	  * readFrame - Read and return a single `lsl.reader.tbw.Frame` instance
	  * read - Read in the capture and return it as a numpy array
	"""
	
	def _readyFile(self):
		"""
		Find the start of valid TBF data.  This function:
		1) Aligns on the first valid Mark 5C frame.
		"""
		
		# Align on the start of a Mark5C packet
		while True:
			try:
				junkFrame = tbf.readFrame(self.fh)
				break
			except errors.syncError:
				self.fh.seek(-tbf.FrameSize+1, 1)
				
		# Skip over any DRX frames the start of the file
		i = 0
		while True:
			try:
				junkFrame = tbf.readFrame(self.fh)
				break
			except errors.syncError:
				i += 1
				self.fh.seek(-tbf.FrameSize+drx.FrameSize, 1)
		if i == 0:
			self.fh.seek(-tbf.FrameSize, 1)
		self.fh.seek(-tbf.FrameSize, 1)
		
		return True
		
	def _describeFile(self):
		"""
		Describe the TBF file.
		"""
		
		# Read in frame
		junkFrame = tbf.readFrame(self.fh)
		self.fh.seek(-tbf.FrameSize, 1)
		
		# Basic file information
		filesize = os.fstat(self.fh.fileno()).st_size
		nFramesFile = filesize / tbf.FrameSize
		srate = fC
		bits = 4
		nFramesPerObs = tbf.getFramesPerObs(self.fh)
		nChan = tbf.getChannelCount(self.fh)
		
		# Pre-load the channel mapper and find the first frame
		self.mapper = []
		marker = self.fh.tell()
		firstFrameCount = 2**64-1
		while len(self.mapper) < nChan/12:
			cFrame = tbf.readFrame(self.fh)
			if cFrame.header.firstChan not in self.mapper:
				self.mapper.append( cFrame.header.firstChan )
			if cFrame.header.frameCount < firstFrameCount:
				firstFrameCount = cFrame.header.frameCount
				start = junkFrame.getTime()
				startRaw = junkFrame.data.timeTag
		self.fh.seek(marker)
		self.mapper.sort()
		
		# Calculate the frequencies
		freq = numpy.zeros(nChan)
		for i,c in enumerate(self.mapper):
			freq[i*12:i*12+12] = c + numpy.arange(12)
		freq *= fC
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 'FrameSize': tbf.FrameSize,
						'firstFrameCount': firstFrameCount, 'sampleRate': srate, 'dataBits': bits, 
						'nAntenna': 512, 'nChan': nChan, 'freq1': freq, 'tStart': start, 
						'tStartSamples': startRaw}
						
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.tbw.Frame` instance.
		"""
		
		frame = tbf.readFrame(self.fh)
		while not frame.header.isTBF():
			frame = tbf.readFrame(self.fh)
			
		return frame
		
	def read(self, duration=None, timeInSamples=False):
		"""
		Read and return the entire TBW capture.  This function returns 
		a three-element tuple with elements of:
		  0) the actual duration of data read in, 
		  1) the time tag for the first sample, and
		  2) a 3-D Numpy array of data.
		  
		The time tag is returned as seconds since the UNIX epoch by default.
		However, the time tags can be returns as samples at fS if the 
		timeInSamples keyword is set.
		
		The sorting order of the output data array is by 
		digitizer number - 1.
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Find out how many frames to work with at a time
		framesPerObs = self.description['nChan']/12
		frameCount = int(round(1.0 * duration * self.description['sampleRate'] * framesPerObs))
		frameCount = frameCount if frameCount else 1
		duration = frameCount / framesPerObs / self.description['sampleRate']
		firstFrameCount = self.description['firstFrameCount']
		
		# Initialize the output data array
		data = numpy.zeros((self.description['nAntenna'], self.description['nChan'], frameCount/framesPerObs), dtype=numpy.complex64)
		
		# Read in the next frame and anticipate any problems that could occur
		i = 0
		while i < frameCount:
			try:
				cFrame = tbf.readFrame(self.fh)
				i += 1
			except errors.eofError:
				break
			except errors.syncError:
				continue
				
			if not cFrame.header.isTBF():
				continue
				
			if cFrame.header.frameCount == firstFrameCount:
				if timeInSamples:
					setTime = cFrame.data.timeTag
				else:
					setTime = cFrame.getTime()
					
			firstChan = cFrame.header.firstChan
			cnt = cFrame.header.frameCount - firstFrameCount
			
			subData = cFrame.data.fDomain
			subData.shape = (12,512)
			subData = subData.T
			
			aStand = self.mapper.index(firstChan)
			try:
				data[:,aStand*12:(aStand+1)*12,cnt] = subData
			except IndexError:
				pass
				
		# Calculate the duration
		duration = data.shape[2] / self.getInfo('sampleRate')
		
		return duration, setTime, data


def LWASVDataFile(filename=None, fh=None, ignoreTimeTagErrors=False):
	"""
	Wrapper around the various LWA-SV-related classes defined here that takes
	a file, determines the data type, and initializes and returns the 
	appropriate LDP class.
	"""
	
	# Open the file as appropriate
	if fh is None:
		fh = open(filename, 'rb')
	else:
		filename = fh.name
		if fh.mode.find('b') == -1:
			fh.close()
			fh = open(filename, 'rb')
			
	# Read a bit of data to try to find the right type
	for mode in (drx, tbn, tbf, cor, drspec):
		## Set if we find a valid frame marker
		foundMatch = False
		## Set if we can read more than one valid successfully
		foundMode = False
		
		## Sort out the frame size.  This is tricky because DR spectrometer files
		## have frames of different sizes depending on the mode
		if mode == drspec:
			try:
				mfs = drspec.getFrameSize(fh)
			except:
				mfs = 0
		else:
			mfs = mode.FrameSize
			
		## Loop over the frame size to try and find what looks like valid data.  If
		## is is found, set 'foundMatch' to True.
		for i in xrange(mfs):
			try:
				junkFrame = mode.readFrame(fh)
				foundMatch = True
				break
			except errors.eofError:
				break
			except errors.syncError:
				fh.seek(-mfs+1, 1)
				
		## Did we strike upon a valid frame?
		if foundMatch:
			### Is so, we now need to try and read more frames to make sure we have 
			### the correct type of file
			fh.seek(-mfs, 1)
			
			try:
				for i in xrange(2):
					junkFrame = mode.readFrame(fh)
				foundMode = True
			except errors.eofError:
				break
			except errors.syncError:
				### Reset for the next mode...
				fh.seek(0)
		else:
			### Reset for the next mode...
			fh.seek(0)
			
		## Did we read more than one valid frame?
		if foundMode:
			break
			
	# There is an ambiguity that can arise for TBF data such that it *looks* 
	# like DRX.  If the identified mode is DRX, skip halfway into the file and 
	# verify that it is still DRX.   We also need to catch the LWA1 TBN vs.
	# TBW ambiguity since we could have been given an LWA1 file by accident.
	if mode in (drx, tbn):
		## Sort out the frame size
		omfs = mode.FrameSize
		
		## Seek half-way in
		nFrames = os.path.getsize(filename)/omfs
		fh.seek(nFrames/2*omfs)
		
		## Read a bit of data to try to find the right type
		for mode in (tbn, drx, tbf):
			### Set if we find a valid frame marker
			foundMatch = False
			### Set if we can read more than one valid successfully
			foundMode = False
			
			### Sort out the frame size.
			mfs = mode.FrameSize
			
			### Loop over the frame size to try and find what looks like valid data.  If
			### is is found, set 'foundMatch' to True.
			for i in xrange(mfs):
				try:
					junkFrame = mode.readFrame(fh)
					foundMatch = True
					break
				except errors.eofError:
					break
				except errors.syncError:
					fh.seek(-mfs+1, 1)
					
			### Did we strike upon a valid frame?
			if foundMatch:
				#### Is so, we now need to try and read more frames to make sure we have 
				#### the correct type of file
				fh.seek(-mfs, 1)
				
				try:
					for i in xrange(4):
						junkFrame = mode.readFrame(fh)
					foundMode = True
				except errors.syncError:
					#### Reset for the next mode...
					fh.seek(nFrames/2*omfs)
			else:
				#### Reset for the next mode...
				fh.seek(nFrames/2*omfs)
				
			### Did we read more than one valid frame?
			if foundMode:
				break
				
	fh.close()
	
	# Raise an error if nothing is found
	if not foundMode:
		raise RuntimeError("File '%s' does not appear to be a valid LWA-SV data file" % filename)
		
	# Otherwise, build and return the correct LDPFileBase sub-class
	if mode == drx:
		ldpInstance = DRXFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
	elif mode == tbn:
		ldpInstance = TBNFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
	elif mode == tbf:
		ldpInstance = TBFFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
	else:
		ldpInstance = DRSpecFile(filename, ignoreTimeTagErrors=ignoreTimeTagErrors)
		
	# Done
	return ldpInstance


def LWADataFile(filename=None, fh=None, ignoreTimeTagErrors=False):
	"""
	Wrapper around the various classes defined here that takes a file, 
	determines the data type, and initializes and returns the appropriate
	LDP class.
	"""
	
	found = False
	
	# LWA-1?
	if not found:
		try:
			ldpInstance = LWA1DataFile(filename=filename, fh=fh, ignoreTimeTagErrors=ignoreTimeTagErrors)
			found = True
		except RuntimeError:
			pass
			
	# LWA-SV?
	if not found:
		try:
			ldpInstance = LWASVDataFile(filename=filename, fh=fh, ignoreTimeTagErrors=ignoreTimeTagErrors)
			found = True
		except RuntimeError:
			pass
			
	# Failed?
	if not found:
		raise RuntimeError("File '%s' does not appear to be a valid LWA1 or LWA-SV data file" % filename)
		
	return ldpInstance
