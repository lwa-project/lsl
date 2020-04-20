"""
Module that contains the error classes for the DRX, TBN, and TBW readers.  
These errors are currently meant to deal with file I/O problems.

.. versionchanged::1.1.4
    Removed numpyError and re-enumerated
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['BaseReaderError', 'EOFError', 'SyncError', 'notTBNError', 'notTBWError', 'list_error_codes', 
           'MinErrorNo', 'MaxErrorNo']


MinErrorNo = 1
MaxErrorNo = 4


class BaseReaderError(IOError):
    """
    Base class for file I/O problem during numpy.fromfile calls and out-of-
    sync Mark5C headers.
    """

    def __init__(self, strerror, errno='-1'):
        self.errno = errno
        self.strerror = strerror
        self.filename = None
        self.args = (errno, strerror)
        IOError.__init__(self)

    def __str__(self):
        return "%s" % self.strerror


class EOFError(BaseReaderError):
    """
    Extension to the base class for dealing with EOF errors.  The error code
    is 1.
    """

    def __init__(self):
        self.errno = 1
        self.strerror = 'End of file encountered during filehandle read'
        self.filename = None
        self.args = (self.errno, self.strerror)
        BaseReaderError.__init__(self, self.strerror, errno=self.errno)


class SyncError(BaseReaderError):
    """
    Extension to the base class for dealing with Mark 5C header sync word 
    problems.  If the sync word doesn't match what is expected.  The error code 
    is 3.
    """

    def __init__(self, location=None, sync1=None, sync2=None, sync3=None, sync4=None):
        self.errno = 2
        self.strerror = 'Mark 5C sync word differs from expected'
        self.filename = None
        self.args = (self.errno, self.strerror)
        self.location = location
        self.syncWord = (sync1, sync2, sync3, sync4)
        BaseReaderError.__init__(self, self.strerror, errno=self.errno)

    def __str__(self):
        output = self.strerror
        if self.location is not None:
            output = '%s at byte %i' % (output, self.location)
        if self.syncWord[0] is not None:
            output = '%s: %02X %02X %02X %02X' % (output, self.syncWord[0], self.syncWord[1], self.syncWord[2], self.syncWord[3])
        return output


class notTBNError(BaseReaderError):
    """
    Extenstion to the base class for dealing with trying to read in TBW data 
    with a TBN reader.  The error code is 4.
    """

    def __init__(self):
        self.errno = 3
        self.strerror = 'Data appears to be TBW, not TBN as expected'
        self.filename = None
        self.args = (self.errno, self.strerror)
        BaseReaderError.__init__(self, self.strerror, errno=self.errno)


class notTBWError(BaseReaderError):
    """
    Extenstion to the base class for dealing with trying to read in TBN data 
    with a TBW reader.  The error code is 5.
    """

    def __init__(self):
        self.errno = 4
        self.strerror = 'Data appears to be TBN, not TBW as expected'
        self.filename = None
        self.args = (self.errno, self.strerror)
        BaseReaderError.__init__(self, self.strerror, errno=self.errno)


def list_error_codes(errno=None):
    """
    Function to provide a list of errors defined in this file.  It 
    alternatively takes an error code using the 'errno' keyword and returns its
    description.
    """

    if errno is None:
        for i in range(MinErrorNo, (MaxErrorNo+1)):
            list_error_codes(errno=i)
    else:
        if errno == 1:
            print("1: End of file encountered during filehandle read")
        elif errno == 2:
            print("2: Mark 5C sync word differs from expected")
        elif errno == 3:
            print("3: Data appears to be TBW, not TBN as expected")
        elif errno == 4:
            print("4: Data appears to be TBN, not TBW as expected")
        else:
            print("Unknown error code '%i'" % errno)
