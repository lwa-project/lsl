import copy

class FrameBuffer(object):
	"""Frame buffer for re-ordering TBW/TBN/DRX frames in time order.  
	This class is filled with frames and a returns a frame list when 
	either:
	  + the buffer contains all frames for a given time, or
	  + another has started to be filled and has more than 5 frames in 
	    it.
	The buffer also keeps track of what has already been read out so 
	that tardy frames are just dropped.  For buffers that are read out,
	missing frames are replaced with frames filled with zeros."""

	def __init__(self, frameObject=None, nFrames=20):
		"""Initialize the buffer with what type of frame you have and
		how many frame constitute a full buffer."""

		self.frameObject = frameObject
		self.buffer = {}
		self.nFrames = nFrames

	def append(self, frame):
		"""Append a new frame to the buffer with the appropriate time
		tag."""

		if frame.data.timeTag not in self.buffer.keys():
			self.buffer[frame.data.timeTag] = []
		self.buffer[frame.data.timeTag].append(frame)

	def get(self, frame):
		"""Return a list of frames that consitute a 'full' buffer.  
		Afterwards, delete that buffer."""
		
		

	def __createFill(self, key=None, id=0, pol=0, tune=None):
		"""Create a 'fill' frame of zeros using an existing good
		packet as a template."""

		fillFrame = copy.deepcopy(self.buffer[key][0])
		
		# Fix-up the header
		try:
			fillFrame.tbwID = id
		except AttributeError:
			pass
		try:
			fillFrame.tbnID = 2*id + pol
		except AttributeError:
			pass
		try:
			fillFrame.drxID = 4*id + 2*tune + pol
		except AttributeError:
			pass
		
		# Zero the data for the fill packet
		try:
			exampleFrame.data.iq *= 0
		except AttributeError:
			exampleFrame.data.xy *= 0
		
		return fillFrame

	def __fillLevel(self):
		fillLevel = {}
		for key in self.buffer.keys():
			fillLevel[key] = len(self.buffer[key])

		return fillLevel

