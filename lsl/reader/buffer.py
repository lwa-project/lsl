# -*- coding: utf-8 -*-

"""
Buffer for dealing with out-of-order/missing frames.

.. versionchanged:: 0.5
	Removed support for DRX FrameBuffers since they shouldn't be needed.
	
.. versionchanged:: 0.6
	Removed support for TBW FrameBuffers since they didn't really work.
"""

import copy
from collections import deque
try:
	from collections import OrderedDict
except ImportError:
	from lsl.misc.OrderedDict import OrderedDict


__version__ = '0.7'
__revision__ = '$Rev$'
__all__ = ['FrameBuffer', 'TBNFrameBuffer', '__version__', '__revision__', '__all__']


def _cmpStands(x, y):
	"""
	Function to compare two frames and sort by stand/beam number.  This 
	should work TBW, TBN, and DRX.
	"""
	
	# Parse if frame IDs to extract the stand/beam, tunning, and polarization
	# information (where appropriate)
	idsX = x.parseID()
	idsY = y.parseID()

	# Do a try...except block to catch TBW vs. TBN/DRX
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
	"""
	Frame buffer for re-ordering TBW and TBN frames in time order.  
	This class is filled with frames and a returns a frame list when 
	the 'nSegments' starts filling.  In that case, the oldest segment 
	is returned.

	The buffer also keeps track of what has already been read out so 
	that tardy frames are just dropped.  For buffers that are read out,
	missing frames are replaced with frames filled with zeros.
	
	.. note::
		Due to the nature of the buffer, it is possible that there will
		still be 'nSegments'-1 segements in the buffer that are either
		full or partially full.  This can be retrieved using the buffer's 
		'flush()' function.
		
	.. note::
		DRX should _not_ need a ring buffer since the four beam outputs
		are, effectively, on different networks.
	"""

	def __init__(self, mode='TBN', stands=[], pols=[], samples=12000000, bits=12, nSegments=6, ReorderFrames=False):
		"""
		Initialize the buffer with a list of:
		  * TBW:
		      * list of stands
		      * number of data samples (default of 12,000,000)
		  * TBN
		      * list of stands
		      * list of pols
		By doing this, we should be able to keep up with when the buffer 
		is full and to help figure out which stands are missing.
		"""

		# Input validation
		if mode.upper() not in ('TBW', 'TBN'):
			raise RuntimeError("Invalid observing mode '%s'" % mode)
		
		if mode.upper() == 'TBW':
			if bits == 12:
				if samples > 12000000:
					raise RuntimeError('Invalid number of samples for 12-bit TBW')
			else:
				if samples > 36000000:
					raise RuntimeError('Invalid number of samples for 4-bit TBW')
				
		else:
			for pol in pols:
				if pol not in (0, 1):
					raise RuntimeError("Invalid polarization '%i'" % pol)

		# The buffer itself
		self.nSegments = nSegments
		self.buffer = OrderedDict()
		self.done = deque([0,], maxlen=self.nSegments)
		
		# Buffer statistics
		self.full = 0		# Number of times a full buffer was emptied
		self.partial = 0	# Number of times a partial buffer was emptied
		self.dropped = 0	# Number of late frames dropped
		self.missing = 0	# Number of missing frames
		
		# Information describing the observations
		self.mode = mode.upper()
		self.stands = stands
		self.pols = pols
		self.bits = bits
		self.samples = samples
		
		# If we should reorder the returned frames by stand/pol or not
		self.reorder = ReorderFrames
		
		# Figure out how many frames fill the buffer and the list of all
		# possible frames in the data set
		calcFramesOut = self.__calcFrames()
		self.nFrames = calcFramesOut[0]
		self.possibleFrames = calcFramesOut[1]
		
	def __calcFrames(self):
		"""
		Calculate the maximum number of frames that we expect from 
		the setup of the observations and a list of tuples that describes
		all of the possible stand/pol combination.
		"""
		
		nFrames = 0
		frameList = []
		if self.mode == 'TBW':
			# TBW mode
			if self.bits == 12:
				nFrames = self.samples / 400
			else:
				nFrames = self.samples / 1200
				
			for stand in self.stands:
				frameList.append(stand)
			
		else:
			# TBN mode
			nFrames = len(self.stands)*len(self.pols)
				
			for stand in self.stands:
				for pol in self.pols:
					frameList.append((stand,pol))
			
		return (nFrames, frameList)
		
	def figureOfMerit(self, frame):
		"""
		Figure of merit for storing/sorting frames in the ring buffer.
		
		This will be overridden by sub-classes of FrameBuffer.
		"""
		
		pass

	def append(self, frames):
		"""
		Append a new frame to the buffer with the appropriate time tag.  
		True is returned if the frame was added to the buffer and False if 
		the frame was dropped because it belongs to a buffer that has 
		already been returned.
		"""

		# Convert input to a deque (if needed)
		typeName = type(frames).__name__
		if typeName == 'deque':
			pass
		elif typeName == 'list':
			frames = deque(frames)
		else:
			frames = deque([frames,])

		# Loop over frames
		while True:
			try:
				frame = frames.popleft()
			except IndexError:
				break

			# Make sure that it is not in the `done' list.  If it is,
			# disgaurd the frame and make a note of it.
			fom = self.figureOfMerit(frame)
			if fom < self.done[-1]:
				self.dropped += 1
				continue

			# If that time tag hasn't been done yet, add it to the 
			# buffer in the correct place.
			try:
				self.buffer[fom].append(frame)
			except KeyError:
				self.buffer[fom] = deque()
				self.buffer[fom].append(frame)
				
		return True

	def get(self):
		"""
		Return a list of frames that consitute a 'full' buffer.  Afterwards, 
		delete that buffer and mark it as closed so that any missing frames that
		are recieved late are dropped.  If none of the buffers are ready to be 
		dumped, None is returned.
		"""

		# Get the current status of the buffer
		keys = self.buffer.keys()
		
		# If the ring is full, dump the oldest
		if len(keys) < self.nSegments:
			return None
		
		oldestKey = keys[0]
		oldestCount = len(self.buffer[oldestKey])
		
		if oldestCount == self.nFrames:
			self.full = self.full + 1
			
			output = self.buffer[oldestKey]
		else:
			self.partial = self.partial + 1
			self.missing = self.missing + (self.nFrames - len(self.buffer[oldestKey]))
			
			output = self.buffer[oldestKey]
			for frame in self.__missingList(oldestKey):
				output.append( self.__createFill(oldestKey, frame) )
		
		del(self.buffer[oldestKey])
		self.done.append(oldestKey)
		
		# Sort and return
		if self.reorder:
			output = list(output)
			output.sort(cmp=_cmpStands)
		return output
		
	def flush(self):
		"""
		Return a list of lists containing all remaining frames in the 
		buffer from buffers that are considered 'full'.  Afterwards, 
		delete all buffers.  This is useful for emptying the buffer after
		reading in all of the data.
		
		.. note::
			It is possible for this function to return list of packets that
			are mostly invalid.
		"""
		
		keys = self.buffer.keys()
		
		output = []
		for key in keys:
			keySize = len(self.buffer[key])
			
			if keySize == self.nFrames:
				self.full = self.full + 1
				
				output2 = self.buffer[key]
			
			else:
				self.partial = self.partial + 1
				self.missing = self.missing + (self.nFrames - len(self.buffer[key]))
				
				output2 = self.buffer[key]
				for frame in self.__missingList(key):
					output2.append( self.__createFill(key, frame) )
			
			if self.reorder:
				output2 = list(output2)
				output2.sort(cmp=_cmpStands)
			output.append( output2 )
			del(self.buffer[key])
			self.done.append(key)
				
		return output

	def __missingList(self, key):
		"""
		Create a list of tuples of missing frame information.
		"""
		
		# Find out what frames we have
		frameList = []
		for frame in self.buffer[key]:
			frameList.append(frame.parseID())
			
		# Compare the existing list with the possible list stored in the 
		# FrameBuffer object to build a list of missing frames.
		missingList = []
		for frame in self.possibleFrames:
			if frame not in frameList:
				missingList.append(frame)
				
		return missingList

	def __createFill(self, key, frameParameters):
		"""
		Create a 'fill' frame of zeros using an existing good
		packet as a template.
		"""

		# Get a template based on the first frame for the current buffer
		fillFrame = copy.deepcopy(self.buffer[key][0])
		
		# Get out the frame parameters
		if self.mode == 'TBW':
			stand = frameParameters
		else:
			stand, pol = frameParameters
		
		# Fix-up the header
		if self.mode == 'TBW':
			fillFrame.header.tbwID = (fillFrame.header.tbwID&64512)|stand
		else:
			fillFrame.header.tbnID = 2*(stand-1) + pol + 1
		
		# Zero the data for the fill packet
		if self.mode == 'TBW':
			fillFrame.data.xy *= 0
		else:
			fillFrame.data.iq *= 0
		
		# Invalidate the frame
		fillFrame.valid = False
		
		return fillFrame
		
	def status(self):
		"""
		Print out the status of the buffer.  This contains information about:
		  1.  The current buffer fill level
		  2. The numer of full and partial buffer dumps preformed
		  3. The number of missing frames that fill packets needed to be created
		     for
		  4. The number of frames that arrived too late to be incorporated into 
		     one of the ring buffers
		"""
		
		nf = 0
		for key in self.buffer.keys():
			nf = nf + len(self.buffer[key])
			
		outString = ''
		outString = '\n'.join([outString, "Current buffer level:  %i frames" % nf])
		outString = '\n'.join([outString, "Buffer dumps:  %i full / %i partial" % (self.full, self.partial)])
		outString = '\n'.join([outString, "--"])
		outString = '\n'.join([outString, "Missing frames:  %i" % self.missing])
		outString = '\n'.join([outString, "Dropped frames:  %i" % self.dropped])

		print outString


class TBNFrameBuffer(FrameBuffer):
	"""
	A sub-type of FrameBuffer specifically for dealing with TBN frames.
	See :class:`lsl.reader.buffer.FrameBuffer` for a description of how the 
	buffering is implemented.
	
	Keywords:
	stands
	  list of stand to expect packets for
	  
	pols
	  list of polarizations to expect packets for
	  
	nSegments
	  number of ring segments to use for the buffer (default is 20)
	  
	ReorderFrames
	  whether or not to reorder frames returned by get() or flush() by 
	  stand/polarization (default is False)
	  
	The number of segements in the ring can be converted to a buffer time in 
	seconds:
	
	+----------+------------------------------------------------+
	|          |                TBN Filter Code                 |
	| Segments +------+------+------+------+------+------+------+
	|          |  1   |  2   |  3   |  4   |  5   |  6   |  7   |
	+----------+------+------+------+------+------+------+------+
	|    10    |  5.1 |  1.6 |  0.8 |  0.4 |  0.2 |  0.1 | 0.05 |
	+----------+------+------+------+------+------+------+------+
	|    20    | 10.2 |  3.3 |  1.6 |  0.8 |  0.4 |  0.2 | 0.10 |
	+----------+------+------+------+------+------+------+------+
	|    30    | 15.4 |  4.9 |  2.5 |  1.2 |  0.6 |  0.3 | 0.15 |
	+----------+------+------+------+------+------+------+------+
	|    40    | 20.5 |  6.6 |  3.3 |  1.6 |  0.8 |  0.4 | 0.20 |
	+----------+------+------+------+------+------+------+------+
	|    50    | 25.6 |  8.2 |  4.1 |  2.0 |  1.0 |  0.5 | 0.26 |
	+----------+------+------+------+------+------+------+------+
	|   100    | 51.2 | 16.4 |  8.2 |  4.1 |  2.0 |  1.0 | 0.51 |
	+----------+------+------+------+------+------+------+------+
	
	"""
	
	def __init__(self, stands=[], pols=[0, 1], nSegments=20, ReorderFrames=False):
		super(TBNFrameBuffer, self).__init__(mode='TBN', stands=stands, pols=pols, nSegments=nSegments, ReorderFrames=ReorderFrames)
		
	def figureOfMerit(self, frame):
		"""
		Figure of merit for sorting frames.  For TBN it is:
		    <frame count of frame>
		"""
		
		return frame.data.timeTag
		#return frame.header.frameCount
