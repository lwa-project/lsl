"""
Python module to read in DR spectrometer data.  DR Spectrometer data contain DRX
data that has been transformed to the Fourier domain, converted to power, and
integrated on-the-fly by the data recorder.  These data can have various
integration times, channel counts, and polarization products stored.

This module defines the following classes for storing the spectra found in a
file:

Frame
  object that contains all data associated with a particular spectrometer frame.  
  The primary constituents of each frame are:
    * FrameHeader - the spectrometer frame header object and
    * FramePayload   - the spectral data object.
Combined, these two objects contain all of the information found in the 
original spectrometer frame.

The functions defined in this module fall into two class:
  1. convert a frame in a file to a Frame object and
  2. describe the format of the data in the file.

For reading in data, use the read_frame function.  It takes a python file-
handle as an input and returns a fully-filled Frame object
For describing the format of data in the file, three function are provided:
  * get_sample_rate - get the sample rate in the file
  * getFRAME_SIZE - get the total (header + data) frame size
  * get_ffts_per_integration - get the number of FFT windows per integration
  * get_transform_size - get the FFT length
  * get_integration_time - get the integration time

.. note::
    This reader works with the most current version of the DR spectrometer data
    format as specified in the version 1.7.  To read data created with previous
    versions of the DR spectrometer, use LSL version 0.5.2.
    
.. versionchanged:: 1.0.1
    Added in new functions to help better describe the contents of a DR 
    spectrometer file.
"""

import copy
import numpy

from lsl.common import dp as dp_common
from lsl.reader.base import *
from lsl.reader.drx import FILTER_CODES as drx_FILTER_CODES
from lsl.reader._gofast import read_drspec
from lsl.reader._gofast import SyncError as gSyncError
from lsl.reader._gofast import EOFError as gEOFError
from lsl.reader.errors import SyncError, EOFError
from lsl.reader.utils import FilePositionSaver

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
__all__ = ['FrameHeader', 'FramePayload', 'Frame', 'read_frame', 'get_data_products', 'is_linear',
           'is_stokes', 'get_sample_rate', 'get_frame_size', 'get_ffts_per_integration', 
           'get_transform_size', 'get_integration_time', 'FILTER_CODES']

#: List of filter codes and their corresponding sample rates in Hz.  
#: .. note::
#:		These are just the DRX filter codes
FILTER_CODES = drx_FILTER_CODES


class FrameHeader(FrameHeaderBase):
    """
    Class that stores the information found in the header of a DR spectrometer/DRX 
    frame.
    """
    
    _header_attrs = ['beam', 'format', 'decimation', 'time_offset', 'nints']
    
    def __init__(self, beam=0, format=0, decimation=None, time_offset=None, nints=None):
        self.beam = beam
        self.format = format
        self.decimation = decimation
        self.time_offset = time_offset
        
        if nints is None:
            self.nints = 0
        else:
            self.nints = nints
            
        FrameHeaderBase.__init__(self)
        
    @property
    def id(self):
        """
        Return the beam the frame corresponds to.
        """
        
        return self.beam
        
    @property
    def data_products(self):
        """
        Return a list of data products contained in the file.
        
        .. versionadded:: 0.6.0
        """
        
        products = []
        
        # Linear
        if self.format & 0x01:
            products.append('XX')
        if self.format & 0x02:
            products.append('XY_real')
        if self.format & 0x04:
            products.append('XY_imag')
        if self.format & 0x08:
            products.append('YY')
            
        # Stokes
        if self.format & 0x10:
            products.append('I')
        if self.format & 0x20:
            products.append('Q')
        if self.format & 0x40:
            products.append('U')
        if self.format & 0x80:
            products.append('V')
        
        return products
        
    @property
    def is_linear(self):
        """
        Return whether or not the frame contains linear polarization 
        products or not.
        
        .. versionadded:: 0.6.0
        """
        
        if self.format < 0x10:
            return True
        else:
            return False
            
    @property
    def is_stokes(self):
        """
        Return whether or not the frame contains Stokes polarization
        parameters or not.
        
        .. versionadded:: 0.6.0
        """
        
        if self.format < 0x10:
            return False
        else:
            return True
        
    @property
    def sample_rate(self):
        """
        Return the sample rate of the data in samples/second.
        """
        
        sample_rate = dp_common.fS / self.decimation
        return sample_rate
        
    @property
    def filter_code(self):
        """
        Function to convert the sample rate in Hz to a filter code.
        """
        
        sampleCodes = {}
        for key,value in FILTER_CODES.items():
            sampleCodes[value] = key
            
        return sampleCodes[self.sample_rate]
        
    @property
    def ffts_per_integration(self):
        """
        Return the number of FFT windows per integration.
        
        .. versionadded:: 1.0.1
        """
        
        return self.nints


class FramePayload(FramePayloadBase):
    """
    Class that stores the information found in the data section of a DR spectrometer/
    DRX frame.
    
    .. versionchanged:: 0.5.3
        Added the saturations field to keep up with saturation counts.
        
    .. versionchanged:: 0.6.0
        The attributes that store the data are not defined until a frame is read in order
        to account for the fact that spectrometer files can store either linear or Stokes
        data.
    """
    
    _payload_attrs = ['timetag', 'tuning_words', 'fills', 'errors', 'saturations']
    
    def __init__(self, timetag=None, tuning_words=None, fills=None, errors=None, saturations=None):
        self.timetag = timetag
        if tuning_words is None:
            self.tuning_words = [0, 0]
        else:
            self.tuning_words = tuning_words
        if fills is None:
            self.fills = [0, 0, 0, 0]
        else:
            self.fills = fills
        if errors is None:
            self.errors = [0, 0, 0, 0]
        else:
            self.errors = errors
        if saturations is None:
            self.saturations = [0, 0, 0, 0]
        else:
            self.saturations = saturations
        FramePayloadBase.__init__(self, None)
        
    @property
    def data(self):
        packed = []
        dtype = []
        for p in ('XX', 'XY_real', 'XY_imag', 'YY', 'I', 'Q', 'U', 'V'):
            for t in (0, 1):
                sub = getattr(self, "%s%i" % (p, t), None)
                if sub is not None:
                    packed.append(sub)
                    dtype.append(("%s%i" % (p, t), packed[-1].dtype))
        return numpy.rec.array(packed, dtype=dtype)
        
    @property
    def central_freq(self):
        """
        Function to set the central frequency of the DRX data in Hz.
        """
        
        return [dp_common.fS * i / 2**32 for i in self.tuning_words]


class Frame(FrameBase):
    """
    Class that stores the information contained within a single DR spectrometer/
    DRX frame.  It's properties are FrameHeader and FramePayloadLinear/FramePayloadStokes
    objects.
    
    .. versionchanged:: 0.6.0
        By default the data contained with in a frame is normalized by the number of
        fills (header.fills parameter).  For data products that are a function of more
        than one primary input, i.e., real(XY*) or I, the minimum fill of X and Y are 
        used for normalization.
    """
    
    _header_class = FrameHeader
    _payload_class = FramePayload
    gain = None
    
    @property
    def id(self):
        """
        Convenience wrapper for the Frame.FrameHeader.id 
        property.
        """
        
        return self.header.id
        
    @property
    def data_products(self):
        """
        Convenience wrapper for the Frame.FrameHeder.data_products
        property.
        """
        
        return self.header.data_products
        
    @property
    def is_linear(self):
        """
        Convenience wrapper for the Frame.FrameHeder.is_linear
        property.
        """
        
        return self.header.is_linear
        
    @property
    def is_stokes(self):
        """
        Convenience wrapper for the Frame.FrameHeder.is_stokes
        property.
        """
        
        return self.header.is_stokes
        
    @property
    def sample_rate(self):
        """
        Convenience wrapper for the Frame.FrameHeader.sample_rate 
        property.
        """
        
        return self.header.sample_rate
        
    @property
    def filter_code(self):
        """
        Convenience wrapper for the Frame.FrameHeader.filter_code property.
        """
        
        return self.header.filter_code
        
    @property
    def ffts_per_integration(self):
        """
        Conveinence wrapper for the Frame.FrameHeader.ffts_per_integration 
        property.
        
        .. versionadded:: 1.0.1
        """
        
        return self.header.ffts_per_integration
        
    @property
    def time(self):
        """
        Function to convert the time tag from samples since the UNIX epoch
        (UTC 1970-01-01 00:00:00) to seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance.
        """
        
        return FrameTimestamp.from_dp_timetag(self.payload.timetag, offset=self.header.time_offset)
        
    @property
    def central_freq(self):
        """
        Convenience wrapper for the Frame.FramePayload.central_freq property.
        """
        
        return self.payload.central_freq
        
    @property
    def transform_size(self):
        """
        Find out what the transform size is.
        
        .. versionadded:: 1.0.1
        """
        
        p = self.data_products[0]
        return getattr(self.payload, "%s0" % p, None).size
        
    @property
    def integration_time(self):
        """
        Return the integration time for data in seconds.
        
        .. versionadded:: 1.0.1
        """
        
        LFFT = self.transform_size
        srate = self.sample_rate
        nints = self.ffts_per_integration
        
        return nints*LFFT/srate
        
    def __add__(self, y):
        """
        Add the data sections of two frames together or add a number 
        to every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame += y
        return newFrame
        
    def __iadd__(self, y):
        """
        In-place add the data sections of two frames together or add 
        a number to every element in the data section.
        """
        
        for attrBase in self.header.data_products:
            for tuning in (0, 1):
                attr = "%s%i" % (attrBase, tuning)
                try:
                    temp = getattr(self.payload, attr, None) + getattr(y.payload, attr, None)
                except TypeError:
                    raise RuntimeError(f"Cannot add {str(attrs)} with {str(y.header.get_data_products())}")
                except AttributeError:
                    temp = getattr(self.payload, attr, None) + numpy.float32(y)
                setattr(self.payload, attr, temp)
            
        return self
        
    def __sub__(self, y):
        """
        Subtract the data sections of two frames or subtract a number 
        from every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame -= y
        return newFrame
        
    def __isub__(self, y):
        """
        In-place subtract the data sections of two frames or subtract 
        a number from every element in the data section.
        """
        
        for attrBase in self.header.data_products:
            for tuning in (0, 1):
                attr = "%s%i" % (attrBase, tuning)
                try:
                    temp = getattr(self.payload, attr, None) - getattr(y.payload, attr, None)
                except TypeError:
                    raise RuntimeError(f"Cannot add {str(attrs)} with {str(y.header.get_data_products())}")
                except AttributeError:
                    temp = getattr(self.payload, attr, None) - numpy.float32(y)
                setattr(self.payload, attr, temp)
            
        return self
        
    def __mul__(self, y):
        """
        Multiple the data sections of two frames together or multiply 
        a number to every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame *= y
        return newFrame
            
    def __imul__(self, y):
        """
        In-place multiple the data sections of two frames together or 
        multiply a number to every element in the data section.
        """
        
        for attrBase in self.header.data_products:
            for tuning in (0, 1):
                attr = "%s%i" % (attrBase, tuning)
                try:
                    temp = getattr(self.payload, attr, None) * getattr(y.payload, attr, None)
                except TypeError:
                    raise RuntimeError(f"Cannot multiply {str(attrs)} with {str(y.header.get_data_products())}")
                except AttributeError:
                    temp = getattr(self.payload, attr, None) * numpy.float32(y)
                setattr(self.payload, attr, temp)
                
        return self
        
    def __floordiv__(self, y):
        """
        Divide the data sections of two frames or divide
        a number from every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame //= y
        return newFrame
            
    def __ifloordiv__(self, y):
        """
        In-place divide the data sections of two frames or 
        divide a number from every element in the data section.
        """
        
        for attrBase in self.header.data_products:
            for tuning in (0, 1):
                attr = "%s%i" % (attrBase, tuning)
                try:
                    temp = getattr(self.payload, attr, None) // getattr(y.payload, attr, None)
                except TypeError:
                    raise RuntimeError(f"Cannot multiply {str(attrs)} with {str(y.header.get_data_products())}")
                except AttributeError:
                    temp = getattr(self.payload, attr, None) // numpy.float32(y)
                setattr(self.payload, attr, temp)
                
        return self
        
    def __truediv__(self, y):
        """
        Divide the data sections of two frames or divide
        a number from every element in the data section.
        """
        
        newFrame = copy.deepcopy(self)
        newFrame /= y
        return newFrame
            
    def __itruediv__(self, y):
        """
        In-place divide the data sections of two frames or 
        divide a number from every element in the data section.
        """
        
        for attrBase in self.header.data_products:
            for tuning in (0, 1):
                attr = "%s%i" % (attrBase, tuning)
                try:
                    temp = getattr(self.payload, attr, None) / getattr(y.payload, attr, None)
                except TypeError:
                    raise RuntimeError(f"Cannot multiply {str(attrs)} with {str(y.header.get_data_products())}")
                except AttributeError:
                    temp = getattr(self.payload, attr, None) / numpy.float32(y)
                setattr(self.payload, attr, temp)
                
        return self


def read_frame(filehandle, gain=None, verbose=False):
    """
    Function to read in a single DR spectrometer/DRX frame (header+data) and 
    store the contents as a Frame object.
    """
    
    # New Go Fast! (TM) method
    try:
        newFrame = read_drspec(filehandle, Frame())
    except gSyncError:
        mark = filehandle.tell()
        raise SyncError(type='DRSpectrometer', location=mark)
    except gEOFError:
        raise EOFError
        
    if gain is not None:
        newFrame.gain = gain
        
    return newFrame


def get_data_products(filehandle):
    """
    Find out the data products contained in the file by looking at a frame.
    """
    
    with FilePositionSaver(filehandle):
        # Read in one frame
        newFrame = read_frame(filehandle)
        
    # Return the data products
    return newFrame.header.data_products


def is_linear(filehandle):
    """
    Find out if the file contains linear polarization products or not.
    """
    
    with FilePositionSaver(filehandle):
        # Read in one frame
        newFrame = read_frame(filehandle)
        
    # Return the verdict
    return newFrame.header.is_linear


def is_stokes(filehandle):
    """
    Find out if the file contains Stokes parameters or not.
    """
    
    with FilePositionSaver(filehandle):
        # Read in one frame
        newFrame = read_frame(filehandle)
        
    # Return the verdict
    return newFrame.header.is_stokes


def get_sample_rate(filehandle, nframes=None, filter_code=False):
    """
    Find out what the sampling rate/filter code is from a single observations.  
    By default, the rate in Hz is returned.  However, the corresponding filter 
    code can be returned instead by setting the FilterCode keyword to true.
    """
    
    with FilePositionSaver(filehandle):
        # Read in one frame
        newFrame = read_frame(filehandle)
        
    if not filter_code:
        return newFrame.sample_rate
    else:
        return newFrame.filter_code


def get_frame_size(filehandle):
    """
    Find out what the frame size in a file is at the current file location.
    Returns the frame size in bytes.
    """
    
    with FilePositionSaver(filehandle):
        cPos = filehandle.tell()
        read_frame(filehandle)
        nPos = filehandle.tell()
        
    return nPos - cPos


def get_ffts_per_integration(filehandle):
    """
    Find out what the number of FFT windows per integration is at the 
    current file location.
    
    .. versionadded:: 1.0.1
    """
    
    with FilePositionSaver(filehandle):
        frame = read_frame(filehandle)
        
    return frame.ffts_per_integration


def get_transform_size(filehandle):
    """
    Find out what the transform size in a file is at the current file 
    location.  
    """
    
    with FilePositionSaver(filehandle):
        frame = read_frame(filehandle)
        
    return frame.transform_size


def get_integration_time(filehandle):
    """
    Find out what the integration time is at the current file location.
    """
    
    with FilePositionSaver(filehandle):
        frame = read_frame(filehandle)
        
    return frame.integration_time
