"""
Module that contains the error classes for the DRX, TBN, and TBW readers.  
These errors are currently meant to deal with file I/O problems.

.. versionchanged::1.1.4
    Removed numpyError and re-enumerated
"""

import logging

from lsl.logger import LSL_LOGGER

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
__all__ = ['BaseReaderError', 'EOFError', 'SyncError', 'list_error_codes', 
           'MinErrorNo', 'MaxErrorNo']


MinErrorNo = 1
MaxErrorNo = 2


class BaseReaderError(IOError):
    """
    Base class for file I/O problem during numpy.fromfile calls and out-of-
    sync Mark5C headers.
    """

    def __init__(self, strerror, errno='-1', filename=None):
        IOError.__init__(self, errno, strerror, filename)
        
    def __str__(self):
        return str(self.strerror)


class EOFError(BaseReaderError):
    """
    Extension to the base class for dealing with EOF errors.  The error code
    is 1.
    """

    def __init__(self):
        errno = 1
        strerror = 'End of file encountered during filehandle read'
        BaseReaderError.__init__(self, strerror, errno=errno)


class SyncError(BaseReaderError):
    """
    Extension to the base class for dealing with Mark 5C header sync word 
    problems.  If the sync word doesn't match what is expected.  The error code 
    is 2.
    """

    def __init__(self, type='Mark5C', location=None, sync1=None, sync2=None, sync3=None, sync4=None):
        errno = 2
        strerror = '%s sync word differs from expected' % type
        BaseReaderError.__init__(self, strerror, errno=errno)
        self.location = location
        self.syncWord = (sync1, sync2, sync3, sync4)
        
    def __str__(self):
        output = self.strerror
        if self.location is not None:
            output = f"{output} at byte {self.location}"
        if self.syncWord[0] is not None:
            output = f"{output}: {self.syncWord[0]:02X} {self.syncWord[1]:02X} {self.syncWord[2]:02X} {self.syncWord[3]:02X}"
        return output


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
            LSL_LOGGER.info("1: End of file encountered during filehandle read")
        elif errno == 2:
            print("2: Sync word differs from expected")
            LSL_LOGGER.info("2: Sync word differs from expected")
        else:
            print(f"Unknown error code '{errno}'")
            LSL_LOGGER.info(f"Unknown error code '{errno}'")
