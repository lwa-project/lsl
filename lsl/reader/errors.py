# -*- coding: utf-8 -*-

"""Module that contains the error classes for the DRX, TBN, and TBW readers.  
These errors are currently meant to deal with file I/O problems."""

__version__ = '0.1'
__revision__ = '$ Revision: 0 $'
__all__ = ['baseReaderError', 'eofError', 'numpyError', 'syncError', 'notTBNError', 'notTBWError', 'listErrorCodes', 'MinErrorNo', 'MaxErrorNo', '__version__', '__revision__', '__all__']


MinErrorNo = 1
MaxErrorNo = 5


class baseReaderError(IOError):
	"""Base class for file I/O problem during numpy.fromfile calls and out-of-
sync Mark5C headers."""

	def __init__(self, strerror, errno='-1'):
		self.errno = errno
		self.strerror = strerror
		self.filename = None
		self.args = (errno, strerror)

	def __str__(self):
		return "%s" % self.strerror


class eofError(baseReaderError):
	"""Extension to the base class for dealing with EOF errors.  The error code
	is 1."""

	def __init__(self):
		self.errno = 1
		self.strerror = 'End of file encountered during filehandle read'
		self.filename = None
		self.args = (self.errno, self.strerror)


class numpyError(baseReaderError):
	"""Extension to the base class for dealing with errors that occur when 
	numpy.fromfile tried to read more data can exists.  This is a specialized form
	of EOF error.  The error code is 2."""

	def __init__(self):
		self.errno = 2
		self.strerror = 'End of file encountered during numpy.fromfile call'
		self.filename = None
		self.args = (self.errno, self.strerror)


class syncError(baseReaderError):
	"""Extension to the base class for dealing with Mark 5C header sync word 
	problems.  If the sync word doesn't match what is expected.  The error code 
	is 3."""

	def __init__(self, sync1=None, sync2=None, sync3=None, sync4=None):
		self.errno = 3
		self.strerror = 'Mark 5C sync word differs from expected'
		self.filename = None
		self.args = (self.errno, self.strerror)
		self.syncWord = (sync1, sync2, sync3, sync4)

	def __str__(self):
		if self.syncWord[0] is None:
			return self.strerror
		else:
			return "%s:  %2X %2X %2X %2X" % (self.strerror, self.syncWord[0], self.syncWord[1], self.syncWord[2], self.syncWord[3])


class notTBNError(baseReaderError):
	"""Extenstion to the base class for dealing with trying to read in TBW data 
	with a TBN reader.  The error code is 4."""

	def __init__(self):
		self.errno = 4
		self.strerror = 'Data appears to be TBW, not TBN as expected'
		self.filename = None
		self.args = (self.errno, self.strerror)


class notTBWError(baseReaderError):
	"""Extenstion to the base class for dealing with trying to read in TBN data 
	with a TBW reader.  The error code is 5."""

	def __init__(self):
		self.errno = 5
		self.strerror = 'Data appears to be TBN, not TBW as expected'
		self.filename = None
		self.args = (self.errno, self.strerror)


def listErrorCodes(errno=None):
	"""Function to provide a list of errors defined in this file.  It 
	alternatively takes an error code using the 'errno' keyword and returns its
	description."""

	if errno is None:
		for i in range(MinErrorNo, (MaxErrorNo+1)):
			listErrorCodes(errno=i)
	else:
		if errno == 1:
			print "1: End of file encountered during filehandle read"
		elif errno == 2:
			print "2: End of file encountered during numpy.fromfile call"
		elif errno == 3:
			print "3: Mark 5C sync word differs from expected"
		elif errno == 4:
			print "4: Data appears to be TBW, not TBN as expected"
		elif errno == 5:
			print "5: Data appears to be TBN, not TBW as expected"
		else:
			print "Unknown error code '%i'" % errno


		