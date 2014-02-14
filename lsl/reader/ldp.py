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
  
Also included is the LWA1DataFile function that take a filename and tries to determine the 
correct data format object to use.
"""

import os
import numpy
from scipy.stats import norm
from collections import deque

from lsl.common.dp import fS
from lsl.reader import tbw, tbn, drx, drspec, errors
from lsl.reader.buffer import TBNFrameBuffer, DRXFrameBuffer

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['TBWFile', 'TBNFile', 'DRXFile', 'DRSpecFile', 'LWA1DataFile', '__version__', '__revision__', '__all__']


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
		
	def _readyFile(self):
		"""
		Method for finding the start of valid data.  This will be over-
		ridden in the format-specific subclasses.
		"""
		
		pass
		
	def _describeFile(self):
		"""
		Method for describing the contents of a file using.  This will 
		be over-ridden in the format-specific subclasses.
		"""
		
		pass
		
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
				raise ValueError("Unknown key '%d'" % key)
				
	def reset(self):
		"""
		Reset the file to the beginning.
		"""
		
		self.fh.seek(0)
		
		# Ready the file
		self._readyFile()
		
		# Describe the contents of the file
		self.description = {}
		self._describeFile()
		
	def close(self):
		self.fh.close()


class TBWFile(LDPFileBase):
	"""
	Class to make it easy to interface with a TBW file.  Method defined for this class are:
	  * getInfo - Get information about the file's contents
	  * readFrame - Read and return a single `lsl.reader.tbw.Frame` instance
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
			junkFrame = tbw.readFrame(self.fh)
		self.fh.seek(-tbw.FrameSize, 1)
		
		return True
		
	def _describeFile(self):
		"""
		Describe the TBW file.
		"""
		
		junkFrame = self.readFrame()
		self.fh.seek(-tbw.FrameSize, 1)
		
		filesize = os.fstat(self.fh.fileno()).st_size
		nFramesFile = filesize / tbw.FrameSize
		srate = 196e6
		bits = junkFrame.getDataBits()
		start = junkFrame.getTime()
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start}
						
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.tbw.Frame` instance.
		"""
		
		return tbw.readFrame(self.fh)


class TBNFile(LDPFileBase):
	"""
	Class to make it easy to interface with a TBN file.  Method defined for this class are:
	  * getInfo - Get information about the file's contents
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
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 
						'nAntenna': framesPerObsX+framesPerObsY, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start, 'freq1': tuning1}
						
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
		
		return 1.0 * frameOffset / self.description['nAntenna'] * 512 / self.description['sampleRate']
		
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.tbn.Frame` instance.
		"""
		
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
		""" 
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
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
			for j in xrange(nAntenna):
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
		self.fh.seek(-tbn.FrameSize*nAntenna*nFrames, 1)
		
		# Statistics
		rv = norm()
		frac = rv.cdf(Sigma) - rv.cdf(-Sigma)
		index = int(round(data.shape[1]*frac))
		if index == data.shape[1]:
			index = data.shape[1] - 1
		
		levels = [0 for i in xrange(nAntenna)]
		for i in xrange(nAntenna):
			data2 = sorted(data[i,:])
			levels[i] = data2[index]
		
		return levels


class DRXFile(LDPFileBase):
	"""
	Class to make it easy to interface with a DRX file.  Method defined for this class are:
	  * getInfo - Get information about the file's contents
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
		for i in xrange(16):
			junkFrame = drx.readFrame(self.fh)
			b,t,p = junkFrame.parseID()
			id = (t,p)
			if id not in ids:
				ids.append(id)
			timeTags.append(junkFrame.data.timeTag)
		self.fh.seek(-16*drx.FrameSize, 1)
		
		if len(ids) == 4:
			i = 0
			while (timeTags[i+0] != timeTags[i+1]) or (timeTags[i+0] != timeTags[i+2]) or (timeTags[i+0] != timeTags[i+3]):
				i += 1
				self.fh.seek(drx.FrameSize, 1)
		else:
			i = 0
			while (timeTags[i+0] != timeTags[i+1]):
				i += 1
				self.fh.seek(drx.FrameSize, 1)
				
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
				tuning1 = junkFrame.getCentralFreq()
				
			if i == 0:
				start = junkFrame.getTime()
		self.fh.seek(-drx.FrameSize*4, 1)
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 
						'beampols': beampols, 'beam': b, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start, 'freq1': tuning1, 'freq2': tuning2}
						
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
		offset = int(offset * sampleRate / 4096 * beampols)
		offset = int(1.0 * offset / beampols) * beampols
		self.fh.seek(offset*drx.FrameSize, 1)
		
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
			offset += cOffset
			
			## If the offset is zero, we are done.  Otherwise, apply the offset
			## and check the location in the file again/
			if cOffset is 0:
				break
			self.fh.seek(cOffset*drx.FrameSize, 1)
			
		self.description['beampols'] = beampols
		self.description['sampleRate'] = sampleRate
		
		return t1 - t0
		
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.drx.Frame` instance.
		"""
		
		return drx.readFrame(self.fh)
		
	def read(self, duration, timeInSamples=False):
		"""
		Given an open DRX file and an amount of data to read in in seconds, read 
		in the data and return a three-element tuple of the actual duration read 
		in, the time for the first sample, and the data as numpy 
		array.
		
		..note::
			This function always returns a 2-D array with the first dimension
			holding four elements.
		"""
		
		# Make sure there is file left to read
		if self.fh.tell() == os.fstat(self.fh.fileno()).st_size:
			raise errors.eofError()
			
		# Covert the sample rate to an expected timetag skip
		timetagSkip = int(4096 / self.description['sampleRate'] * fS)
		
		# Setup the counter variables:  frame count and time tag count
		timetag = {0:0, 1:0, 2:0, 3:0}
		for i in xrange(self.description['beampols']):
			junkFrame = drx.readFrame(self.fh)
			b,t,p = junkFrame.parseID()
			aStand = 2*(t-1) + p
			timetag[aStand] = junkFrame.data.timeTag - timetagSkip
		self.fh.seek(-drx.FrameSize*self.description['beampols'], 1)
		
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
				if not self.ignoreTimeTagErrors:
					cTimetag = cFrame.data.timeTag
					if cTimetag != timetag[aStand]+timetagSkip:
						actStep = cTimetag - timetag[aStand]
						raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (timetagSkip, t, p, actStep))
						
				if setTime is None:
					if timeInSamples:
						setTime = cFrame.data.timeTag - cFrame.header.timeOffset
					else:
						setTime = cFrame.getTime()
						
				data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
				count[aStand] +=  1
				timetag[aStand] = cTimetag
			nFrameSets += 1
			
		# If we've hit the end of the file and haven't read in enough frames, 
		# flush the buffer
		if eofFound or nFrameSets != frameCount:
			for cFrames in self.buffer.flush():
				for cFrame in cFrames:
					b,t,p = cFrame.parseID()
					aStand = 2*(t-1) + p
					if not self.ignoreTimeTagErrors:
						cTimetag = cFrame.data.timeTag
						if cTimetag != timetag[aStand]+timetagSkip:
							actStep = cTimetag - timetag[aStand]
							raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (timetagSkip, t, p, actStep))
							
					if setTime is None:
						if timeInSamples:
							setTime = cFrame.data.timeTag - cFrame.header.timeOffset
						else:
							setTime = cFrame.getTime()
							
					data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
					count[aStand] +=  1
					timetag[aStand] = cTimetag
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
		self.fh.seek(-drx.FrameSize*beampols*nFrames, 1)
		
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
		tuning1, tuning2 = junkFrame.getCentralFreq()
		prod = junkFrame.getDataProducts()
		
		self.description = {'size': filesize, 'nFrames': nFramesFile, 'FrameSize': FrameSize, 
						'beampols': beampols, 'beam': beam, 
						'sampleRate': srate, 'dataBits': bits, 
						'tStart': start, 'freq1': tuning1, 'freq2': tuning2, 
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
		offset = int(round(offset / self.description['tInt']))
		self.fh.seek(offset*FrameSize, 1)
		
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
			offset += cOffset
			
			## If the offset is zero, we are done.  Otherwise, apply the offset
			## and check the location in the file again/
			if cOffset is 0:
				break
			self.fh.seek(cOffset*self.description['FrameSize'], 1)
			
		return t1 - t0
		
	def readFrame(self):
		"""
		Read and return a single `lsl.reader.drspec.Frame` instance.
		"""
		
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
		timetag = 0
		junkFrame = drspec.readFrame(self.fh)
		timetag = junkFrame.getTime() - timetagSkip
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
				
			if not self.ignoreTimeTagErrors:
				cTimetag = cFrame.getTime()
				if cTimetag > timetag + 1.001*timetagSkip:
					actStep = cTimetag - timetag
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
			timetag = cTimetag
			nFrameSets += 1
			
		# Adjust the duration to account for all of the things that could 
		# have gone wrong while reading the data
		duration = nFrameSets * self.description['tInt']
		
		return duration, setTime, data


def LWA1DataFile(filename=None, fh=None, ignoreTimeTagErrors=False):
	"""
	Wrapper around the various classes defined here that takes a file, 
	determines the data type, and initializes and returns the appropriate
	LDP class.
	"""
	
	# Open the file as appropriate
	if fh is None:
		fh = open(filename, 'rb')
	else:
		filename = fh.name
		if fh.mode.find('b') == -1:
			fh.close()
			fh = open(self.filename, 'rb')
			
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
			except errors.syncError:
				### Reset for the next mode...
				fh.seek(0)
		else:
			### Reset for the next mode...
			fh.seek(0)
			
		## Did we read more than one valid frame?
		if foundMode:
			break
	fh.close()
	
	# Raise an error if nothing is found
	if not foundMode:
		raise RuntimeError("File '%s' does not appear to be a valid LWA1 data file" % filename)
		
	# Otherwise, build and return the correct LDPFileBase sub-class
	if mode == drx:
		ldpInstance = DRXFile(filename)
	elif mode == tbn:
		ldpInstance = TBNFile(filename)
	elif mode == tbw:
		ldpInstance = TBWFile(filename)
	else:
		ldpInstance = DRSpecFile(filename)
		
	# Done
	return ldpInstance