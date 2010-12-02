# -*- coding: utf-8 -*-

"""Buffer for dealing with out-of-order/missing frames."""

import copy

__version__ = '0.3'
__revision__ = '$ Revision: 5 $'
__all__ = ['FrameBuffer', 'TBNFrameBuffer', 'DRXFrameBuffer', 'TBWFrameBuffer', '__version__', '__revision__', '__all__']


def cmpStands(x, y):
	"""Function to compare two frames and sort by stand/beam number."""
	
	# Parse if frame IDs to extract the stand/beam, tunning, and polarization
	# information (where appropriate)
	idsX = x.parseID()
	idsY = y.parseID()
	try:
		len1 = len(idsX)
		if len1 == 2:
			sX = 2*(idsX[0]-1) + idsX[1]
			sY = 2*(idsY[0]-1) + idsY[1]
		else:
			sX = 4*(idsX[0]-1) + 2*(idsX[1]-1) + idsX[2]
			sY = 4*(idsY[0]-1) + 2*(idsY[1]-1) + idsY[2]
	except TypeError:
		sX = idsX
		sY = idsY
		
	# Do the comparison
	if sY > sX:
		return -1
	elif sX > sY:
		return 1
	else:
		return 0
	

class FrameBuffer(object):
	"""Frame buffer for re-ordering TBN and DRX frames in time order.  
	This class is filled with frames and a returns a frame list when 
	either:
	  * the buffer contains all frames for a given time, or
	  * another has started to be filled and has more than 50 frames in 
	    it.

	The buffer also keeps track of what has already been read out so 
	that tardy frames are just dropped.  For buffers that are read out,
	missing frames are replaced with frames filled with zeros."""

	def __init__(self, stands=[], pols=[], beams=[], tunes=[]):
		"""Initialize the buffer with a list of:
		  * TBN
		      * list of stands
		      * list of pols
		  * DRX
		      * list of beams
		      * list of tunnings
		      * list of pols
		By doing this, we should be able to keep up with when the buffer 
		is full and to help figure out which stands are missing."""

		# The buffer itself
		self.buffer = {}
		self.done = []
		
		# Buffer statistics
		self.full = 0		# Number of times a full buffer was emptied
		self.partial = 0	# Number of times a partial buffer was emptied
		self.dropped = 0	# Number of late frames dropped
		self.missing = 0	# Number of missing frames
		
		# Information describing the obserations
		self.stands = stands
		self.pols = pols
		self.beams = beams
		self.tunes = tunes
		
		# Figure out how many frames fill the buffer and the list of all
		# possible frames in the data set
		calcFramesOut = self.__calcFrames()
		self.nFrames = calcFramesOut[0]
		self.possibleFrames = calcFramesOut[1]
		
		# When to dump the dump partial packets
		self.dump = 50
		
	def __calcFrames(self):
		"""Calculate the maximum number of frames that we expect from 
		the setup of the observations and a list of tuples that describe
		all of the possible stand/beam/pol/tune combination."""
		
		# TBN Mode
		if len(self.beams) == 0:
			nFrames = len(self.stands)*len(self.pols)
				
			frameList = []
			for stand in self.stands:
				for pol in self.pols:
					frameList.append((stand,pol))
				
		# DRX mode
		else:
			nFrames = len(self.beams)*len(self.tunes)*len(self.pols)
			
			frameList = []
			for beam in self.beams:
				for tune in self.tunes:
					for pol in self.pols:
						frameList.append((beam, tune, pol))
			
		return (nFrames, frameList)

	def append(self, frame):
		"""Append a new frame to the buffer with the appropriate time
		tag.  True is returned if the frame was added to the buffer and
		False if the frame was dropped."""

		# Make sure that it is not in the `done' list.  If it is,
		# disgaurd the frame and make a note of it.
		if frame.data.timeTag in self.done:
			self.dropped = self.dropped + 1
			return False

		# If that time tag hasn't been done yet, add it to the 
		# buffer in the correct place.
		if frame.data.timeTag not in self.buffer.keys():
			self.buffer[frame.data.timeTag] = []
		self.buffer[frame.data.timeTag].append(frame)
		return True

	def get(self):
		"""Return a list of frames that consitute a 'full' buffer.  
		Afterwards, delete that buffer."""

		# Get the current status of the buffer
		nKeys, minKey, minCount, maxKey, maxCount = self.__bufferStatus()
		if nKeys == 0:
			output = None
		elif nKeys == 1:
			if maxCount == self.nFrames:
				self.done.append(maxKey)
				self.full = self.full + 1
				
				output = copy.deepcopy(self.buffer[maxKey])
				del(self.buffer[maxKey])
			else:
				output = None
		else:
			if maxCount == self.nFrames:
				self.done.append(maxKey)
				self.full = self.full + 1
				
				output = copy.deepcopy(self.buffer[maxKey])
				del(self.buffer[maxKey])
			elif minCount == self.dump:
				self.partial = self.partial + 1
				self.missing = self.missing + (len(self.buffer[maxKey]) - self.nFrames)
				self.done.append(maxKey)
				
				output = copy.deepcopy(self.buffer[maxKey])
				for frame in self.__missingList(maxKey):
					output.append( self.__createFrame(key, frame) )
				
				del(self.buffer[maxKey])
				print "%i is missing Frames" % maxKey
			else:
				output = None
		
		# Sort and return
		if output is None:
			return output
		else:
			return sorted(output, cmp=cmpStands)

	def __bufferStatus(self):
		"""Return the state of the buffer in the form of a tuple.  
		The values are:
		  * number of buffer keys, 
		  * buffer key with minimum length, 
		  * length of the shortest buffer, 
		  * buffer key with maximum length, and
		  * length of the longest buffer.
		If the buffer is empty, min/maxKeys are set to None."""

		# Build the list of keys to look at
		keyList = self.buffer.keys()

		# Setup the basic variables
		nKeys = len(keyList)
		minKey = None
		minCount = self.nFrames+1
		maxKey = None
		maxCount = 0
		
		# Loop to fill in the min/max values
		keyList = self.buffer.keys()
		for key, count in zip(keyList, [len(self.buffer[key]) for key in keyList]):
			if count < minCount:
				minKey = key
				minCount = count
			if count > maxCount:
				maxKey = key
				maxCount = count
				
		return (nKeys, minKey, minCount, maxKey, maxCount)

	def __missingList(self, key):
		"""Create a list of tuples of missing frame information."""
		
		# Find out what frames we have
		frameList = []
		for frame in self.buffer[key]:
			frameList.append(frame.parseID())
			
		# Compare the existing list with the possible list stored in the 
		# FrameBuffer object to build a list of missing frames.
		missingList = []
		for frame in self.frameList:
			if frame not in frameList:
				missingList.append(frame)
				
		return missingList

	def __createFill(self, key, frameParameters):
		"""Create a 'fill' frame of zeros using an existing good
		packet as a template."""

		fillFrame = copy.deepcopy(self.buffer[key][0])
		
		# Get out the frame parameters
		if len(frameParameters) == 2:
			stand, pol = frameParameters
		else:
			bean, tune, pol = frameParameters
		
		# Fix-up the header
		try:
			fillFrame.header.tbnID = 2*(stand-1) + pol + 1
		except AttributeError:
			pass
		try:
			fillFrame.header.drxID = (beam & 7) | ((tune & 7) << 3) | ((pol & 1) << 7)
		except AttributeError:
			pass
		fillFrame.header.frameCount = 0
		
		# Zero the data for the fill packet
		fillFrame.data.iq *= 0
		
		# Invalidate the frame
		fillFrame.valid = False
		
		return fillFrame
		
	def status(self):
		"""Print out the status of the buffer."""
		
		nf = 0
		for key in self.buffer.keys():
			nf = nf + len(self.buffer[key])
			
		outString = ''
		outString = '\n'.join([outString, "Current buffer level:  %i frames" % nf])
		outString = '\n'.join([outString, "Buffer dumps:  %i full / %i partial" % (self.filled, self.partial)])
		outString = '\n'.join([outString, "--"])
		outString = '\n'.join([outString, "Missing frames:  %i" % self.missing])
		outString = '\n'.join([outString, "Dropped frames:  %i" % self.dropped])

		print outString


class TBNFrameBuffer(FrameBuffer):
	"""A sub-type of FrameBuffer specifically for dealing with TBN frames.
	See FrameBuffer for a description of how the buffering is implemented."""
	
	def __init__(self, stands=[], pols=[0, 1]):
		super(TBNFrameBuffer, self).__init__(stands=stands, pols=pols)


class DRXFrameBuffer(FrameBuffer):
	"""A sub-type of FrameBuffer specifically for dealing with DRX frames.
	See FrameBuffer for a description of how the buffering is implemented."""
	
	def __init__(self, beams=[], tunes=[1, 2], pols=[1, 2]):
		super(DRXFrameBuffer, self).__init__(beams=beams, tunes=tunes, pols=pols)


class TBWFrameBuffer(object):
	"""Frame buffer for re-ordering TBW frames in time order.  
	This class is filled with frames and a returns a frame list when 
	either:
	  * the buffer contains all frames for a given time, or
	  * another has started to be filled and has more than 5 frames in 
	    it.

	The buffer also keeps track of what has already been read out so 
	that tardy frames are just dropped.  For buffers that are read out,
	missing frames are replaced with frames filled with zeros.
	
	TBW cannot be treated in the same way as TBN and DRX because of how 
	the packets are sent out by DP, i.e., in stands pairs rather than by 
	time."""

	def __init__(self, stands=[], nFrames=None):
		"""Initialize the buffer with a list of:
		  * TBW:
		      * list of stands
		By doing this, we should be able to keep up with when the buffer 
		is full and to help figure out which stands are missing."""

		# The buffer itself
		self.buffer = {}
		self.done = []
		
		# Buffer statistics
		self.full = 0		# Number of times a full buffer was emptied
		self.partial = 0	# Number of times a partial buffer was emptied
		self.dropped = 0	# Number of late frames dropped
		self.missing = 0	# Number of missing frames
		
		# Information describing the obserations
		self.stands = stands
		
		# Figure out how many frames fill the buffer and the list of all
		# possible frames in the data set
		self.nFrames = nFrames
		self.possibleFrames = self.stands
		
		# When to dump the dump partial packets
		self.dump = 50

	def append(self, frame):
		"""Append a new frame to the buffer with the appropriate time
		tag.  True is returned if the frame was added to the buffer and
		False if the frame was dropped."""

		# Make sure that it is not in the `done' list.  If it is,
		# disgaurd the frame and make a note of it.
		if frame.parseID() in self.done:
			self.dropped = self.dropped + 1
			return False

		# If that time tag hasn't been done yet, add it to the 
		# buffer in the correct place.
		if frame.parseID() not in self.buffer.keys():
			self.buffer[frame.parseID()] = []
		self.buffer[frame.parseID()].append(frame)
		return True

	def get(self):
		"""Return a list of frames that consitute a 'full' buffer.  
		Afterwards, delete that buffer."""
		
		# Get the current status of the buffer
		nKeys, minKey, minCount, maxKey, maxCount = self.__bufferStatus()
		if nKeys <= 1:
			output = None
		elif nKeys == 2:
			if maxCount == self.nFrames:
				self.done.append(maxKey)
				self.full = self.full + 1
				
				output = copy.deepcopy(self.buffer[maxKey])
				del(self.buffer[maxKey])
			else:
				output = None
		else:
			if maxCount == self.nFrames:
				self.done.append(maxKey)
				self.full = self.full + 1
				
				output = copy.deepcopy(self.buffer[maxKey])
				del(self.buffer[maxKey])
			elif minCount == self.dump:
				self.partial = self.partial + 1
				self.missing = self.missing + (len(self.buffer[maxKey]) - self.nFrames)
				self.done.append(maxKey)
				
				output = copy.deepcopy(self.buffer[maxKey])
				del(self.buffer[maxKey])
				print "%i is missing Frames" % maxKey
			else:
				output = None
		# Sort and return
		if output is None:
			return output
		else:
			return sorted(output)

	def __bufferStatus(self):
		"""Return the state of the buffer in the form of a tuple.  
		The values are:
		  * number of buffer keys, 
		  * buffer key with minimum length, 
		  * length of the shortest buffer, 
		  * buffer key with maximum length, and
		  * length of the longest buffer.
		If the buffer is empty, min/maxKeys are set to None."""

		# Build the list of keys to look at

		keyList = self.buffer.keys()
		# Setup the basic variables
		nKeys = len(keyList)
		minKey = None
		minCount = self.nFrames+1
		maxKey = None
		maxCount = 0
		
		# Loop to fill in the min/max values
		keyList = self.buffer.keys()
		for key, count in zip(keyList, [len(self.buffer[key]) for key in keyList]):
			if count < minCount:
				minKey = key
				minCount = count
			if count > maxCount:
				maxKey = key
				maxCount = count
				
		return (nKeys, minKey, minCount, maxKey, maxCount)
		
	def status(self):
		"""Print out the status of the buffer."""
		
		nf = 0
		for key in self.buffer.keys():
			nf = nf + len(self.buffer[key])
			
		outString = ''
		outString = '\n'.join([outString, "Current buffer level:  %i frames" % nf])
		outString = '\n'.join([outString, "Buffer dumps:  %i full / %i partial" % (self.filled, self.partial)])
		outString = '\n'.join([outString, "--"])
		outString = '\n'.join([outString, "Missing frames:  %i" % self.missing])
		outString = '\n'.join([outString, "Dropped frames:  %i" % self.dropped])

		print outString


	# Stuff that may work someday.  The trick is figuring out which 
	# time stamps are missing from a sequence:
	"""
	From get:
		for frame in self.__missingList(maxKey):
			output.append( self.__createFrame(key, frame) )
	
	def __missingList(self, key):
		Create a list of tuples of missing frame information.
		
		frameList = []
		for frame in self.buffer[key]:
			frameList.append(frame.parseID())
			

	def __createFill(self, key, frameParameters):
		Create a 'fill' frame of zeros using an existing good
		packet as a template.

		fillFrame = copy.deepcopy(self.buffer[key][0])
		
		# Fix-up the header
		fillFrame.header.frameCount = 0
		
		# Zero the data for the fill packet
		fillFrame.data.timeTag = frameParameters
		fillFrame.data.xy *= 0
		
		# Invalidate the frame
		fillFrame.valid = False
		
		return fillFrame
	"""