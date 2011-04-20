# -*- coding: utf-8 -*-

"""Module to provide direct access to the files stored on a DRSU.

.. warning::
	Direct access to DRSU is currently an experimental feature that has only
	been tested to a limited extent.  The direct access module also only 
	compiles on Linux systems and may not work with 32-bit Linux installation.
"""


import os
try:
	import _drsu
except ImportError:
	pass
from lsl.reader import errors

__version__ = '0.2'
__revision__ = '$ Revision: 5 $'
__all__ = ['File', 'listFiles', 'getFileByName', 'shepherdedReadFrame', '__version__', '__revision__', '__all__']

class File(object):
	"""Object to provide direct access to a file stored on a DRSU.  The File 
	object supports:
	  * open
	  * seek
	  * tell
	  * read
	  * close
	"""
	
	def __init__(self, device, name, size, ctime=-1, mtime=-1, atime=-1, mode=None):
		self.device = device
		self.name = name
		self.size = size
		self.mode = mode
		
		self.ctime = ctime
		self.mtime = mtime
		self.atime = atime
		
		self.start = 0
		self.stop = 0
		self.chunkSize = 0
		
		self.mjd = -1
		self.mpm = -1
		
		self.fh = None
		self.bytesRead = 0
		
	def open(self):
		"""Open the 'file' on the device and ready the File object for 
		reading."""
		
		# Open the device and set the self.fh attribute
		self.fh = open(self.device, 'rb', buffering=self.chunkSize)
		
		# Move to the start position of the file
		self.fh.seek(self.start)
		
		# Set the bytesRead attribute to zero
		self.bytesRead = 0
		
	def seek(self, offset, whence=0):
		"""Seek in the file.  All three 'whence' modes are supported."""
		
		if whence == 0:
			# From the start of the file (start + offset)
			if offset < 0:
				raise IOError(22, 'Invalid argument')
			
			self.fh.seek(self.start + offset, 0)
			self.bytesRead = offset
			
		elif whence == 1:
			# From the current position
			self.fh.seek(offset, 1)
			self.bytesRead += offset
			
		else:
			# From the end of the file (start + size + offset)
			self.fh.seek(self.start + self.size + offset, 0)
			self.bytesRead = self.start + self.size + offset
	
	def tell(self):
		"""Return the current position in the file by using the internal
		'bytesRead' attribute value."""
		
		return self.bytesRead
	
	def read(self, size=0):
		"""Read in a specified amount of data."""
		
		if self.bytesRead > self.size:
			return ''
		elif self.bytesRead + size > self.size:
			size = self.size - self.bytesRead
		else:
			pass
		self.bytesRead += size
		return self.fh.read(size)
		
	def close(self):
		"""Close the file object."""
		
		self.fh.close()
		self.fh = None
		
	def getFilehandle(self):
		"""Return the underlying file object associated with an open file.
		None is returned if the file isn't open."""
		
		return self.fh


def listFiles(device):
	"""Function to return a list of File instances describing the files on a 
	the specified DRSU device.
	
	.. note::
		Currently, the user needs to have read/write privileges to the device in
		question.  This typically means running the script calling this function
		via `sudo`.
	"""
	
	try:
		return  _drsu.listFiles(device, File)
	except NameError:
		if os.uname()[0] != 'Linux':
			raise RuntimeError("Direct DRSU access is not supported on non-linux OSes")


def getFileByName(device, filename):
	"""Function to return a File instance corresponding to the specified filename 
	on the provided device.  None is returned if the filename does not exists on
	the device."""
	
	try:
		fileList = _drsu.listFiles(device, File)
	except NameError:
		if os.uname()[0] != 'Linux':
			raise RuntimeError("Direct DRSU access is not supported on non-linux OSes")
		
	for drsuFile in fileList:
		if drsuFile.name == filename:
			return drsuFile
	
	return None


def shepherdedReadFrame(File, reader):
	"""Given a (open) File object and a reader module (e.g., :mod:`lsl.reader.tbn`), 
	read in a single Frame instance.  This function wraps the reader's readFrame 
	method and 'shepherds' the reading such that the file size specified by the 
	File.size attribute is enforced on the reading process.  This ensures that data 
	beyond the 'official' end of the file are not read in.
	
	.. note::
		Currently, the user needs to have read/write privileges to the device in
		question.  This typically means running the script calling this function
		via `sudo`.
	"""
	
	# Make sure that the file has actually be opened.  We could do this 
	# here but people should really open their own files so that they 
	# close their own files.
	if File.fh is None:
		raise IOError('File has not be opened for reading')
	
	# Make sure we haven't read past the 'official' end of the file
	if File.bytesRead + reader.FrameSize > File.size:
		raise errors.eofError()
	
	# Read in the frame, update the File's bytesRead attribute.
	#
	# Note:  even if the reader encounters an exception, the File's bytesRead
	# value will be updated since the reader has read all of the data and 
	# advanced the pointer
	File.bytesRead += reader.FrameSize
	frame = reader.readFrame(File.fh)
	
	# Return
	return frame
