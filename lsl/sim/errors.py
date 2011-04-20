# -*- coding: utf-8 -*-

"""
Module that contains the error classes for the DRX, TBN, and TBW simulated
data writers.  These errors are currently meant to deal with things that can 
invalidate a frame or set of frames.
"""

import numpy

__version__ = '0.1'
__revision__ = '$ Revision: 1 $'
__all__ = ['baseSimError', 'invalidStand', 'invalidPol', 'invalidBeam', 'invalidTune', 'invalidDataSize', 'invalidDataType', 'listErrorCodes', 'MinErrorNo', 'MaxErrorNo', '__version__', '__revision__', '__all__']


MinErrorNo = 1
MaxErrorNo = 6


class baseSimError(IOError):
	"""
	Base class for problems that can invalidate a data frame.
	"""

	def __init__(self, strerror, errno='-1'):
		self.errno = errno
		self.strerror = strerror
		self.filename = None
		self.args = (errno, strerror)

	def __str__(self):
		return "%s" % self.strerror


class invalidStand(baseSimError):
	"""
	Extension to the base class for dealing with frames with out-of-range 
	stand numbers.  The error code is 1.
	"""

	def __init__(self):
		self.errno = 1
		self.strerror = 'Stand number is out of range (==0 or >258)'
		self.filename = None
		self.args = (self.errno, self.strerror)


class invalidPol(baseSimError):
	"""
	Extension to the base class for dealing with frames with out-of-range 
	polarization numbers.  The error code is 2.
	"""

	def __init__(self):
		self.errno = 2
		self.strerror = 'Polarization is out of range (!=0 and !=1)'
		self.filename = None
		self.args = (self.errno, self.strerror)


class invalidBeam(baseSimError):
	"""
	Extension to the base class for dealing with frames with out-of-range 
	DRX beam numbers.  The error code is 3.
	"""

	def __init__(self):
		self.errno = 3
		self.strerror = 'Beam is out of range (==0 or >4)'
		self.filename = None
		self.args = (self.errno, self.strerror)


class invalidTune(baseSimError):
	"""
	Extension to the base class for dealing with frames with out-of-range 
	DRX tunning numbers.  The error code is 4.
	"""

	def __init__(self):
		self.errno = 4
		self.strerror = 'Tunning is out of range (!=1 and !=2)'
		self.filename = None
		self.args = (self.errno, self.strerror)


class invalidDataSize(baseSimError):
	"""
	Extension to the base class for dealing with frames with data sections that
	have the wrong array size.  The error code is 5.
	"""

	def __init__(self):
		self.errno = 5
		self.strerror = 'Data array size is not consistent with frame type'
		self.filename = None
		self.args = (self.errno, self.strerror)


class invalidDataType(baseSimError):
	"""
	Extension to the base class for dealing with frames with data sections that
	have the wrong type (real vs complex).  The error code is 6.
	"""

	def __init__(self):
		self.errno = 6
		self.strerror = 'Data array has the wrong general kind'
		self.filename = None
		self.args = (self.errno, self.strerror)


def listErrorCodes(errno=None):
	"""
	Function to provide a list of errors defined in this file.  It 
	alternatively takes an error code using the 'errno' keyword and returns its
	description.
	"""

	if errno is None:
		for i in range(MinErrorNo, (MaxErrorNo+1)):
			listErrorCodes(errno=i)
	else:
		if errno == 1:
			print "1: Stand number is out of range (==0 or >258)"
		elif errno == 2:
			print "2: Polarization is out of range (!=0 and !=1)"
		elif errno == 3:
			print "3: Beam is out of range (==0 or >4)"
		elif errno == 4:
			print "4: Tunning is out of range (!=1 and !=2)"
		elif errno == 5:
			print "5: Data array size is not consistent with frame type"
		elif errno == 6:
			print "6: Data array has the wrong general kind"
		else:
			print "Unknown error code '%i'" % errno
