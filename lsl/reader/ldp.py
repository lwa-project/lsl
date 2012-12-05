"""
LWA Development Primitives - A set of utilities that should make developing 
new analysis software easier.  These functions wrap the nitty gritty of the 
file reading and unpacking behind bite sized funtions.

Functions included are:
  * tbwReadyFile() - Seek in an open filehandle to find the first valid 
                     TBW frame
                     
  * tbnReadyFile() - Seek in an open filehandle to find the first valid 
                     TBN frame
  * tbnOffsetFile() -  Offset by a specified number of seconds into a TBN file
  * tbnInitializeBuffer() - Create an approprite lsl.reader.buffer.TBNFrameBuffer
                            instance for a given TBN file
  * tbnReadFrames() - Read in a certain amount of data from a TBN filehandle
  * tbnEstimateLevels() - Estimate the sqrt(I^2+Q^2) value corresponding to a 
                          certain sigma value (assuming Gaussian noise) for each
                          antenna
  
  * drxReadyFile() - Seek in an open filehandle to find the first valid DRX frame
  * drxOffsetFile () - Offset by a specified number of seconds into a DRX file
  * drxReadFrames() -  Read in a certain amount of data from a DRX filehandle
  * drxEstimateLevels() - Estimate the sqrt(I^2+Q^2) value corresponding to a 
                          certain sigma value (assuming Gaussian noise) for each
                          tuning/polarization
  
The ReadFrames functions return three-element tuples with elements:
  0) Actual amount of data read in
  1) Timestamps for the first data sample (in seconds or ticks)
  2) The data as a NumPy array
"""

import os
import numpy
from scipy.stats import norm
from collections import deque

from lsl.common.dp import fS
from lsl.reader import tbw, tbn, drx, drspec, errors
from lsl.reader.buffer import TBNFrameBuffer

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['tbwReadyFile', 
		 'tbnReadyFile', 'tbnOffsetFile', 'tbnInitializeBuffer', 'tbnReadFrames', 'tbnEstimateLevels', 
		 'drxReadyFile', 'drxOffsetFile', 'drxReadFrames', 'drxEstimateLevels', 
		 'drspecOffsetFile', 'drspecReadFrames', 
		 '__version__', '__revision__', '__all__']


def tbwReadyFile(fh):
	"""
	Given an open file handle, find the start of valid TBW data.  This 
	function:
	  1) Aligns on the first valid Mark 5C frame and
	  2) Skips over any TBN frames at the beginning of the file.
	"""
	
	# Align on the start of a Mark5C packet
	while True:
		try:
			junkFrame = tbw.readFrame(fh)
		except errors.syncError:
			fh.seek(-tbw.FrameSize+1, 1)
			
	# Jump over any TBN data in the file
	while not junkFrame.isTBW():
		junkFrame = tbw.readFrame(fh)
	fh.seek(-tbw.FrameSize, 1)
	
	return True


def tbnReadyFile(fh):
	"""
	Given an open file handle, find the start of valid TBN data.  This
	function:
	  1) Aligns on the first valid Mark 5C frame and
	  2) Skips over any TBW frames at the beginning of the file.
	"""
	
	# Align on the start of a Mark5C packet
	while True:
		try:
			junkFrame = tbn.readFrame(fh)
		except errors.syncError:
			fh.seek(-tbn.FrameSize+1, 1)
			
	# Jump over any TBN data in the file
	while not junkFrame.isTBN():
		junkFrame = tbn.readFrame(fh)
	fh.seek(-tbn.FrameSize, 1)
	
	return True


def tbnOffsetFile(fh, offset, antpols=520):
	"""
	Offset a specified number of seconds in an open TBN file.  This function 
	returns the exact offset time.
	
	.. note::
		The offset provided by this function is relatively crude due to the
		structure of TBN files.
	"""
	
	sampleRate = tbn.getSampleRate(fh)
	frameOffset = int(offset * sampleRate / 512 * antpols)
	frameOffset = int(1.0 * frameOffset / antpols) * antpols
	fh.seek(frameOffset*tbn.FrameSize)
	
	return 1.0 * frameOffset / antpols * 512 / sampleRate


def tbnInitializeBuffer(fh):
	"""
	Create a new lsl.reader.buffer.TBNFrameBuffer instance tailored to the 
	contents of the given file.
	"""
	
	# Get the number of antennas to read in
	nX, nY = tbn.getFramesPerObs(fh)
	nAntenna = nX + nY
	
	pols = []
	if nX != 0:
		pols.append(0)
	if nY != 0:
		pols.append(1)
		
	buffer = TBNFrameBuffer(stands=range(1,nAntenna/len(pols)+1), pols=pols)
	return buffer


def tbnReadFrames(fh, duration, buffer=None, timeInSamples=False):
	"""
	Read in a chunk (in seconds) of TBN data from an open filehandle.  This
	function returns a three-element tuple with elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 2-D Numpy array of data.

	The time tag is returned as seconds sicne the UNIX epoch by default.
	However, the time tags can be returns as samples at fS if the 
	timeInSamples keyword is set.
	""" 
	
	# Make sure there is file left to read
	if fh.tell() == os.fstat(fh.fileno()).st_size:
		raise errors.eofError()
	
	# Get the number of antennas to read in
	nX, nY = tbn.getFramesPerObs(fh)
	nAntenna = nX + nY
	
	if buffer is None:
		# Create a new reordering buffer
		buffer = tbnInitializeBuffer(fh)
		
	# Get the sample rate for the TBN data
	sampleRate = tbn.getSampleRate(fh, nFrames=nAntenna*3)
	
	# Find out how many frames to read in
	frameCount = int(round(1.0 * duration * sampleRate / 512))
	duration = frameCount * 512 / sampleRate
	
	nFrameSets = 0
	eofFound = False
	setTime = None
	count = [0 for i in xrange(nAntenna)]
	data = numpy.zeros((nAntenna, frameCount*512), dtype=numpy.complex64)
	while True:
		if eofFound or nFrameSets == frameCount:
			break
		
		cFrames = deque()
		for i in xrange(nAntenna/2):
			try:
				cFrames.append( tbn.readFrame(fh, Verbose=False) )
			except errors.eofError:
				eofFound = True
				break
			except errors.syncError:
				continue
			
		buffer.append(cFrames)
		cFrames = buffer.get()
		
		# Continue adding frames if nothing comes out.
		if cFrames is None:
			continue
			
		# If something comes out, add it to the data array
		for cFrame in cFrames:
			stand,pol = cFrame.header.parseID()
			aStand = 2*(stand-1)+pol
			
			if setTime is None:
				setTime = cFrame.getTime()
				
			data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
			count[aStand] += 1
		nFrameSets += 1
		
	# If we've hit the end of the file and haven't read in enough frames, 
	# flush the buffer
	if eofFound or nFrameSets != frameCount:
		for cFrames in buffer.flush():
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
	# have gone wrong while readin the data
	duration = nFrameSets * 512 / sampleRate
	
	return duration, setTime, data


def tbnEstimateLevels(fh, nFrames=100, Sigma=5.0):
	"""
	Estimate the n-sigma level for the absolute value of the voltages.  
	Returns a list with indicies that are the digitizer numbers minus one.
	"""
	
	# Make sure there is file left to read
	if fh.tell() == os.fstat(fh.fileno()).st_size:
		raise errors.eofError()
	
	# Get the number of antennas to read in
	nX, nY = tbn.getFramesPerObs(fh)
	nAntenna = nX + nY
	
	# Go!
	count = {}
	for i in xrange(nAntenna):
		count[i] = 0
	data = numpy.zeros((nAntenna, nFrames*512))
	for i in xrange(nFrames):
		for j in xrange(nAntenna):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbn.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				continue
				
			s,p = cFrame.parseID()
			aStand = 2*(s-1) + p
			
			data[aStand, count[aStand]*512:(count[aStand]+1)*512] = numpy.abs( cFrame.data.iq )
			count[aStand] +=  1
	fh.seek(-tbn.FrameSize*nAntenna*nFrames, 1)
	
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


def drxReadyFile(fh):
	"""
	Given an open file handle, find the start of valid DRX data.  This function:
	 1) aligns on the first valid Mark 5C frame and
         2) skips over frames with a decimation of zero. 
	"""
	
	# Align on the start of a Mark5C packet...
	while True:
		try:
			junkFrame = drx.readFrame(fh)
			try:
				# ... that has a valid decimation
				srate = junkFrame.getSampleRate()
				break
			except ZeroDivisionError:
				pass
		except errors.syncError:
			fh.seek(-drx.FrameSize+1, 1)
			
	fh.seek(-drx.FrameSize, 1)
	
	return True


def drxOffsetFile(fh, offset):
	"""
	Offset a specified number of seconds in an open DRX file.  This function 
	returns the exact offset time.
	"""
	
	junkFrame = drx.readFrame(fh)
	fh.seek(-drx.FrameSize, 1)
	
	# Get the initial time, sample rate, and beampols
	t0 = junkFrame.getTime()
	sampleRate = junkFrame.getSampleRate()
	beampols = drx.getFramesPerObs(fh)
	beampols = reduce(int.__add__, beampols)
	
	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(offset * sampleRate / 4096 * beampols)
	offset = int(1.0 * offset / beampols) * beampols
	fh.seek(offset*drx.FrameSize, 1)
	
	# Iterate on the offsets until we reach the right point in the file.  This
	# is needed to deal with files that start with only one tuning and/or a 
	# different sample rate.  
	while True:
		junkFrame = drx.readFrame(fh)
		fh.seek(-drx.FrameSize, 1)
		
		## Figure out where in the file we are and what the current tuning/sample 
		## rate is
		t1 = junkFrame.getTime()
		sampleRate = junkFrame.getSampleRate()
		beampols = drx.getFramesPerObs(fh)
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
		fh.seek(cOffset*drx.FrameSize, 1)
	
	return t1 - t0


def drxReadFrames(fh, duration, timeInSamples=False):
	"""
	Given an open DRX file and an amount of data to read in in seconds, read 
	in the data and return a three-element tuple of the actual duration read 
	in, the times at the beginning of each stream, and the data as numpy 
	array.
	
	..note::
		This function always returns a 2-D array with the first dimension
		holding four elements.
	"""
	
	# Make sure there is file left to read
	if fh.tell() == os.fstat(fh.fileno()).st_size:
		raise errors.eofError()
	
	# Sample the data
	junkFrame = drx.readFrame(fh)
	fh.seek(-drx.FrameSize, 1)
	
	# Figure out the sample rate and how many beam/tuning/pols are in the DRX file
	sampleRate = junkFrame.getSampleRate()
	beampols = drx.getFramesPerObs(fh)
	beampols = reduce(int.__add__, beampols)
	
	# Covert the sample rate to an expected timetag skip
	timetagSkip = int(4096 / sampleRate * fS)
	
	# Setup the counter variables:  frame count and time tag count
	count = {0:0, 1:0, 2:0, 3:0}
	timetag = {0:0, 1:0, 2:0, 3:0}
	for i in xrange(beampols):
		junkFrame = drx.readFrame(fh)
		b,t,p = junkFrame.parseID()
		aStand = 2*(t-1) + p
		timetag[aStand] = junkFrame.data.timeTag - timetagSkip
	fh.seek(-drx.FrameSize*beampols, 1)
	
	# Find out how many frames to read in
	frameCount = int(round(1.0 * duration * sampleRate / 4096))
	duration = frameCount * 4096 / sampleRate
	
	# Setup the output arrays
	if timeInSamples:
		times = numpy.zeros(4, dtype=numpy.int64) - 1
	else:
		times = numpy.zeros(4, dtype=numpy.float64) - 1
	data = numpy.zeros((4,frameCount*4096), dtype=numpy.complex64)
	
	# Go!
	for i in xrange(frameCount*beampols):
		# Read in the next frame and anticipate any problems that could occur
		try:
			cFrame = drx.readFrame(fh, Verbose=False)
		except errors.eofError:
			break
		except errors.syncError:
			continue

		b,t,p = cFrame.parseID()
		aStand = 2*(t-1) + p
		cTimetag = cFrame.data.timeTag
		if cTimetag != timetag[aStand]+timetagSkip:
			actStep = cTimetag - timetag[aStand]
			raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (timetagSkip, t, p, actStep))
			
		if times[aStand] < 0:
			if timeInSamples:
				times[aStand] = cFrame.data.timeTag - cFrame.header.timeOffset
			else:
				times[aStand] = cFrame.getTime()
		
		data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
		count[aStand] +=  1
		timetag[aStand] = cTimetag
		
	return duration, times, data


def drxEstimateLevels(fh, nFrames=100, Sigma=5.0):
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
	if fh.tell() == os.fstat(fh.fileno()).st_size:
		raise errors.eofError()
	
	# Sample the data
	beampols = drx.getFramesPerObs(fh)
	beampols = reduce(int.__add__, beampols)
	
	count = {0:0, 1:0, 2:0, 3:0}
	data = numpy.zeros((4, nFrames*4096))
	for i in xrange(nFrames):
		for j in xrange(beampols):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				continue
				
			b,t,p = cFrame.parseID()
			aStand = 2*(t-1) + p
			
			data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = numpy.abs( cFrame.data.iq )
			count[aStand] +=  1
	fh.seek(-drx.FrameSize*beampols*nFrames, 1)
	
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


def drspecOffsetFile(fh, offset):
	"""
	Offset a specified number of seconds in an open DR spectrometer file.  This 
	function returns the exact offset time.
	"""
	
	# Gather some basic information and read in the first frame
	FrameSize = drspec.getFrameSize(fh)
	LFFT = drspec.getTransformSize(fh)
	junkFrame = drspec.readFrame(fh)
	fh.seek(-FrameSize, 1)
	
	# Get the initial time, sample rate, and integration time
	t0 = junkFrame.getTime()
	sampleRate = junkFrame.getSampleRate()
	tInt = junkFrame.header.nInts*LFFT/sampleRate
	
	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(round(offset / tInt))
	fh.seek(offset*FrameSize, 1)
	
	# Iterate on the offsets until we reach the right point in the file.  This
	# is needed to deal with files that start with only one tuning and/or a 
	# different sample rate.  
	while True:
		junkFrame = drspec.readFrame(fh)
		fh.seek(-drx.FrameSize, 1)
		
		## Figure out where in the file we are and what the current tuning/sample 
		## rate is
		t1 = junkFrame.getTime()
		sampleRate = junkFrame.getSampleRate()
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
		fh.seek(cOffset*FrameSize, 1)
	
	return t1 - t0


def drspecReadFrames(fh, duration, timeInSamples=False):
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
	if fh.tell() == os.fstat(fh.fileno()).st_size:
		raise errors.eofError()
	
	# Sample the data
	FrameSize = drspec.getFrameSize(fh)
	LFFT = drspec.getTransformSize(fh)
	junkFrame = drspec.readFrame(fh)
	fh.seek(-FrameSize, 1)
	
	# Figure out the sample rate, integration time, and number of data products
	sampleRate = junkFrame.getSampleRate()
	tInt = junkFrame.header.nInts*LFFT/sampleRate
	dataProducts = junkFrame.getDataProducts()
	nProducts = len(dataProducts)
	
	# Covert the sample rate to an expected timetag skip
	timetagSkip = int(tInt * fS)
	
	# Setup the counter variables:  frame count and time tag count
	count = 0
	timetag = 0
	junkFrame = drspec.readFrame(fh)
	timetag = junkFrame.data.timeTag - timetagSkip
	fh.seek(-FrameSize, 1)
	
	# Find out how many frames to read in
	frameCount = int(round(1.0 * duration / tInt))
	duration = frameCount * tInt
	
	# Setup the output arrays
	if timeInSamples:
		times = numpy.zeros(2*nProducts, dtype=numpy.int64) - 1
	else:
		times = numpy.zeros(2*nProducts, dtype=numpy.float64) - 1
	data = numpy.zeros((2*nProducts,frameCount,LFFT), dtype=numpy.complex64)
	
	# Go!
	for i in xrange(frameCount):
		# Read in the next frame and anticipate any problems that could occur
		try:
			cFrame = drspec.readFrame(fh, Verbose=False)
		except errors.eofError:
			break
		except errors.syncError:
			continue
			
		cTimetag = cFrame.data.timeTag
		if cTimetag != timetag[aStand]+timetagSkip:
			actStep = cTimetag - timetag[aStand]
			raise RuntimeError("Invalid timetag skip encountered, expected %i but found %i" % (timetagSkip, actStep))
			
		if times[aStand] < 0:
			if timeInSamples:
				times[:] = cFrame.data.timeTag - cFrame.header.timeOffset
			else:
				times[:] = cFrame.getTime()
				
		for j,p in enumerate(dataProducts):
			data[j+0,         count, :] = getattr(cFrame.data, '%s0' % p, None)
			data[j+nProducts, count, :] = getattr(cFrame.data, '%s1' % p, None)
		count +=  1
		timetag = cTimetag
		
	return duration, times, data