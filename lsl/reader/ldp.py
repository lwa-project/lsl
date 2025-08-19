"""
LWA Development Primitives - A set of utilities that should make developing 
new analysis software easier.  These functions wrap the nitty gritty of the 
file reading and unpacking behind Python objects.

Data format objects included are:
  * TBWFile
  * TBNFile
  * DRXFile
  * DRX8File
  * DRSpecFile
  * TBFFile
  * CORFILE

Also included are the LWA1DataFile, LWASVDataFile, LWANADataFile, and LWADataFile
functions that take a filename and try to determine the correct data format
object to use.

.. versionchanged:: 1.2.0
    Added support for LWA-SV ADP data
    
.. versionchanged:: 3.0.4
    Added support for DRX8 data
"""

import os
import abc
import copy
import numpy as np
import warnings
from textwrap import fill as tw_fill
from scipy.stats import norm
from collections import deque, defaultdict

from lsl.common.ndp import fS, fC
from lsl.reader import drx, drx8, drspec, cor, tbx, errors
from lsl.reader.buffer import DRXFrameBuffer, DRX8FrameBuffer, CORFrameBuffer, TBXFrameBuffer
from lsl.reader.utils import *
from lsl.reader.base import FrameTimestamp, CI8
from lsl.common.color import colorfy

from lsl.config import LSL_CONFIG
LDP_CONFIG = LSL_CONFIG.view('ldp')

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.7'
__all__ = ['DRXFile', 'DRX8File', 'DRSpecFile', 'TBXFile', 'LWADataFile']


class _LDPFileRegistry(object):
    """
    Class to keep track of which files are open so that we can close them out
    when we exit.
    
    This concept/framework/class is borrowed from PyTables:
        https://github.com/PyTables/PyTables/blob/master/tables/file.py
    """
    
    def __init__(self):
        self._name_mapping = defaultdict(set)
        self._handlers = set()
        
    @property
    def filenames(self):
        return list(self._name_mapping.keys())
        
    @property
    def handlers(self):
        return self._handlers
        
    def __len__(self):
        return len(self._handlers)
        
    def __contains__(self, filename):
        return filename in self.filenames
        
    def add(self, handler):
        self._name_mapping[handler.filename].add(handler)
        self._handlers.add(handler)
        
    def remove(self, handler):
        filename = handler.filename
        self._name_mapping[filename].remove(handler)
        if not self._name_mapping[filename]:
            del self._name_mapping[filename]
        self._handlers.remove(handler)
        
    def close_all(self):
        handlers = list(self.handlers)  # make a copy
        for handler in handlers:
            handler.close()


_open_ldp_files = _LDPFileRegistry()


class LDPFileBase(object):
    """
    Class to make it easy to interface with raw LWA1 data files and DR spectrometer
    data files.
    """
    
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, filename=None, fh=None, ignore_timetag_errors=False, buffering=-1):
        # Make sure that we are given either a filename or an open file handle
        if filename is None and fh is None:
            raise RuntimeError("Must specify either a filename or open file instance")
            
        # Store a valid file handle and mark the object as ready
        if fh is None:
            self.filename = filename
            self.fh = open(filename, 'rb', buffering)
        else:
            self.filename = fh.name
            if not isinstance(fh, SplitFileWrapper):
                if fh.mode.find('b') == -1:
                    fh.close()
                    fh = open(self.filename, 'rb', buffering)
            self.fh = fh
        _open_ldp_files.add(self)
        
        # Set whether or not reading errors are fatal
        self.ignore_timetag_errors = ignore_timetag_errors
        
        # Ready the file
        self._ready_file()
        
        # Describe the contents of the file
        self.description = {}
        self._describe_file()
        
    def __enter__(self):
        return self
        
    def __exit__(self, type, value, tb):
        self.close()
        
    def __getattr__(self, name):
        ## Try to access the attribute as a real attribute
        try:
            return super(LDPFileBase, self).__getattr__(name)
        except AttributeError:
            pass
            
        ## Try to access the attribute via the 'get_info' method
        try:
            return self.get_info(name)
        except ValueError:
            raise AttributeError("'%s' object has no attribute '%s'" % (type(self).__name__, name))
            
    def __str__(self):
        return f"{type(self).__name__} @ {self.filename}"
        
    def __repr__(self):
        n = self.__class__.__name__
        a = [('filename',repr(self.filename)),
             ('ignore_timetag_errors',self.ignore_timetag_errors)]
        a.extend([(attr,self.get_info(attr)) for attr in self.description])
        output = "<%s" % n
        first = True
        for key,value in a:
            output += "%s %s=%s" % (('' if first else ','), key, value)
            first = False
        output += ">"
        return tw_fill(output, subsequent_indent='    ')
        
    @abc.abstractmethod
    def _ready_file(self):
        """
        Method for finding the start of valid data.  This will be over-
        ridden in the format-specific subclasses.
        """
        
        raise NotImplementedError
        
    @abc.abstractmethod
    def _describe_file(self):
        """
        Method for describing the contents of a file using.  This will 
        be over-ridden in the format-specific subclasses.
        """
        
        raise NotImplementedError
        
    def get_info(self, key=None):
        """
        Retrieve metadata about the file.  This will return a dictionary 
        of values if no key is specified.
        """
        
        if key is None:
            return self.description
        else:
            try:
                return self.description[key]
            except KeyError:
                raise ValueError(f"Unknown key '{key}'")
                
    @property
    def info(self):
        """
        Return a dictionary containing metadata about the file.  Equivalent to
        calling get_info(None).
        """
        
        return self.description
        
    def get_remaining_frame_count(self):
        """
        Return the number of frames left in the file.
        """
        
        return (self.description['size'] - self.fh.tell()) // self.description['frame_size']
        
    @property
    def nframe_remaining(self):
        """
        Alternate method of accessing the number of frames remaining.
        """
        
        return self.get_remaining_frame_count()
        
    def reset(self):
        """
        Reset the file to the beginning.
        """
        
        self.fh.seek(0)
        
        # Ready the file
        self._ready_file()
        
        # Reset any buffers
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        if hasattr(self, "_timetag"):
            self._timetag = None
            
        # Describe the contents of the file
        self.description = {}
        self._describe_file()
        
    def close(self):
        """
        Close the file.
        """
        
        self.fh.close()
        _open_ldp_files.remove(self)
        
    def offset(self, offset):
        """
        Offset into the data.
        """
        
        raise NotImplementedError
        
    @abc.abstractmethod
    def read_frame(self):
        """
        Read a single frame from the data.
        """
        
        raise NotImplementedError
        
    @abc.abstractmethod
    def read(self, duration, time_in_samples=False):
        """
        Read a certain amount of time from the data.
        """
        
        raise NotImplementedError
        
    def read_sequence(self, duration, time_in_samples=False):
        """
        Return a generator that yields the results of the read() method until 
        the end of the file is reached.
        """
        
        while True:
            try:
                output = self.read(duration, time_in_samples=time_in_samples)
                yield output
            except errors.EOFError:
                break
                
    def estimate_levels(self, nframe=10, sigma=5.0):
        """
        Estimate the standard deviation of the data.
        """
        
        raise NotImplementedError


class DRXFile(LDPFileBase):
    """
    Class to make it easy to interface with a DRX file.  DRX data consist of a
    time series of complex data a variable sample rate of up to 19.6 MHz from
    the beamformer.
    
    Methods defined for this class are:
      * get_info - Get information about the file's contents
      * get_remaining_frame_count - Get the number of frames remaining in the file
      * offset - Offset a specified number of seconds into the file
      * read_frame - Read and return a single `lsl.reader.drx.Frame` instance
      * read - Read a chunk of data in and return it as a numpy array
      * estimate_levels - Estimate the n-sigma level for the absolute value of the voltages 
    """
    
    def _ready_file(self):
        """
        Given an open file handle, find the start of valid DRX data.  This function:
        1) aligns on the first valid Mark 5C frame and
        2) skips over frames with a decimation of zero. 
        3) aligns the tuning/polarization timetags
        """
        
        # Align on the start of a Mark5C packet...
        while True:
            try:
                junkFrame = drx.read_frame(self.fh)
                try:
                    # ... this really is DRX data
                    assert(junkFrame.is_drx)
                    # ... that has a valid decimation
                    srate = junkFrame.sample_rate
                    # ... that it comes after 1980 (I don't know why this happens)
                    assert(junkFrame.payload.timetag > 61849368000000000)
                    break
                except (ZeroDivisionError, AssertionError):
                    pass
            except errors.SyncError:
                self.fh.seek(-drx.FRAME_SIZE+1, 1)
                
        self.fh.seek(-drx.FRAME_SIZE, 1)
        
        # Sometimes we get a file that has a large gap at the beginning for
        # reasons unknown.  If this seems to have happened jump over the gap
        # and issue a warning to the user.
        for checkpoint_size in (1, 2, 4, 8, 16):
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
                
            if filesize > checkpoint_size*1024**2:
                foffset = 0
                toffset = 0
                with FilePositionSaver(self.fh):
                    self.fh.seek(((checkpoint_size*1024**2-1)//drx.FRAME_SIZE+1)*drx.FRAME_SIZE, 1)
                    newFrame = drx.read_frame(self.fh)
                    toffset = newFrame.time - junkFrame.time
                    if toffset > 300:
                        foffset = self.fh.tell() - drx.FRAME_SIZE
                        
                if foffset > 0:
                    warnings.warn(colorfy("{{%%yellow Large (%.1f hr) gap at the beginning, skipping in %i B}}" % (toffset/3600, foffset)), RuntimeWarning)
                    self.fh.seek(foffset, 0)
                    break
                    
        # Line up the time tags for the various tunings/polarizations
        ids = []
        timetags = []
        for i in range(32):
            junkFrame = drx.read_frame(self.fh)
            b,t,p = junkFrame.id
            id = (t,p)
            if id not in ids:
                ids.append(id)
            timetags.append(junkFrame.payload.timetag)
        self.fh.seek(-32*drx.FRAME_SIZE, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the DRX file.
        """
        
        try:
            filesize = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            filesize = self.fh.size
        nFramesFile = (filesize - self.fh.tell()) // drx.FRAME_SIZE
        tunepols = drx.get_frames_per_obs(self.fh)
        tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepol
        bits = 4
        
        with FilePositionSaver(self.fh):
            beams = []
            tunes = []
            pols = []
            tuning1 = 0.0
            tuning2 = 0.0
            for i in range(32):
                try:
                    junkFrame0 = self.read_frame()
                    junkFrame1 = self.read_frame()
                except (errors.SyncError, errors.EOFError):
                    break
                for junkFrame in (junkFrame0, junkFrame1):
                    b,t,p = junkFrame.id
                    srate = junkFrame.sample_rate
                    if b not in beams:
                        beams.append(b)
                    if t not in tunes:
                        tunes.append(t)
                    if p not in pols:
                        pols.append(p)
                        
                    if t == 1:
                        tuning1 = junkFrame.central_freq
                    else:
                        tuning2 = junkFrame.central_freq
                        
                if i == 0:
                    start = junkFrame0.time
                    startRaw = junkFrame0.payload.timetag - junkFrame0.header.time_offset
                    
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': drx.FRAME_SIZE,
                            'nbeampol': beampols, 'beam': b, 
                            'sample_rate': srate, 'data_bits': bits, 
                            'start_time': start, 'start_time_samples': startRaw, 'freq1': tuning1, 'freq2': tuning2}
                        
        # Initialize the buffer as part of the description process
        self.buffer = DRXFrameBuffer(beams=beams, tunes=tunes, pols=pols, nsegments=LDP_CONFIG.get('drx_buffer_size'))
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open DRX file.  This function 
        returns the exact offset time.
        """
        
        # Figure out how far we need to offset inside the file
        junkFrame = drx.read_frame(self.fh)
        self.fh.seek(-drx.FRAME_SIZE, 1)
        
        # Get the initial time, sample rate, and beampols
        t0 = junkFrame.time
        if getattr(self, "_timetag", None) is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is not None:
                t0 = FrameTimestamp.from_dp_timetag(curr, junkFrame.header.time_offset)
        sample_rate = junkFrame.sample_rate
        beampols = drx.get_frames_per_obs(self.fh)
        beampols = sum(beampols)
        
        # Offset in frames for beampols beam/tuning/pol. sets
        ioffset = int(offset * sample_rate / 4096 * beampols)
        ioffset = int(1.0 * ioffset / beampols) * beampols
        self.fh.seek(ioffset*drx.FRAME_SIZE, 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.
        nattempt = 0
        diffs_used = deque([], 25)
        while True:
            junkFrame = drx.read_frame(self.fh)
            self.fh.seek(-drx.FRAME_SIZE, 1)
            
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            t1 = junkFrame.time
            sample_rate = junkFrame.sample_rate
            beampols = drx.get_frames_per_obs(self.fh)
            beampols = sum(beampols)
            
            ## See how far off the current frame is from the target
            tDiff = (t1 - t0) - offset
            diffs_used.append(tDiff)
            
            ## Half that to come up with a new seek parameter
            tCorr   = -tDiff / 2.0
            cOffset = int(tCorr * sample_rate / 4096 * beampols)
            cOffset = int(1.0 * cOffset / beampols) * beampols
            ioffset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset == 0:
                break
            try:
                self.fh.seek(cOffset*drx.FRAME_SIZE, 1)
                nattempt += 1
                assert(len(set(diffs_used)) > len(diffs_used)//4)
                assert(nattempt < 1000)
            except (IOError, AssertionError):
                warnings.warn(colorfy("{{%yellow Could not find the correct offset, giving up}}"), RuntimeWarning)

                break
                
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Zero out the time tag checker
        self._timetag = None
        
        return t1 - t0
        
    def read_frame(self, return_ci8=False):
        """
        Read and return a single `lsl.reader.drx.Frame` instance.  If
        `return_ci8` is True then the frame will contain `lsl.reader.base.CI8`
        data instead of numpy.complex64 data.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Zero out the time tag checker
        self._timetagSkip = None
        self._timetag = None
        
        drx_rf = drx.read_frame_ci8 if return_ci8 else drx.read_frame
        return drx_rf(self.fh)
        
    def read(self, duration, time_in_samples=False, return_ci8=False):
        """
        Given an open DRX file and an amount of data to read in in seconds, read 
        in the data and return a three-element tuple of the actual duration read 
        in, the time for the first sample, and the data as numpy 
        array.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.ndp.fS` if the 
        `time_in_samples' keyword is set.
        
        If `return_ci8` is True then the data are returned will contain 
        `lsl.reader.base.CI8` data instead of numpy.complex64.  The two
        dimensions are input by samples.
        
        ..note::
            This function always returns an array with the first dimension
            holding four elements.  These elements contain, in order:
             * Tuning 1, polarization X
             * Tuning 1, polarization Y
             * Tuning 2, polarization X
             * Tuning 2, polarization Y
        """
        
        # Make sure there is file left to read
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Covert the sample rate to an expected timetag skip
        if getattr(self, "_timetagSkip", None) is None:
            self._timetagSkip = int(4096 / self.description['sample_rate'] * fS)
            
        # Setup the read_frame version to use
        drx_rf = drx.read_frame_ci8 if return_ci8 else drx.read_frame
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = {0:0, 1:0, 2:0, 3:0}
            
        # Find out how many frames to read in
        frame_count = int(round(1.0 * duration * self.description['sample_rate'] / 4096))
        frame_count = frame_count if frame_count else 1
        
        # Setup the output arrays
        setTime = None
        if return_ci8:
            data = np.zeros((4,frame_count*4096), dtype=CI8)
            data_view = data.view(np.int16)
        else:
            data = np.zeros((4,frame_count*4096), dtype=np.complex64)
            data_view = data.view(np.float64)
            
        # Go!
        nFrameSets = 0
        eofFound = False
        count = {0:0, 1:0, 2:0, 3:0}
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            if not self.buffer.overfilled:
                cFrames = deque()
                for i in range(self.description['nbeampol']):
                    try:
                        cFrames.append( drx_rf(self.fh, verbose=False) )
                    except errors.EOFError:
                        eofFound = True
                        self.buffer.append(cFrames)
                        cFrames = []
                        break
                    except errors.SyncError:
                        continue
                self.buffer.append(cFrames)
                
            cTimetag = self.buffer.peek()
            if cTimetag is None:
                # Continue adding frames if nothing comes out.
                continue
            else:
                # Otherwise, make sure we are on track
                aStand = 0
                if self._timetag[aStand] == 0:
                    pass
                elif cTimetag != self._timetag[aStand]+self._timetagSkip:
                    missing = (cTimetag - self._timetag[aStand] - self._timetagSkip) / float(self._timetagSkip)
                    if int(missing) == missing and missing < LDP_CONFIG.get('drx_autofill_size'):
                        ## This is kind of black magic down here
                        for m in range(int(missing)):
                            m = self._timetag[aStand] + self._timetagSkip*(m+1)
                            try:
                                baseframe = copy.deepcopy(cFrames[0])
                            except NameError:
                                baseframe = copy.deepcopy(self.buffer.buffer[cTimetag][0])
                            baseframe.payload.timetag = m
                            baseframe.payload._data *= 0
                            self.buffer.append(baseframe)
            cFrames = self.buffer.get()
            
            # If something comes out, add it to the data array
            for cFrame in cFrames:
                b,t,p = cFrame.id
                aStand = 2*(t-1) + p
                cTimetag = cFrame.payload.timetag
                if self._timetag[aStand] == 0:
                    self._timetag[aStand] = cTimetag - self._timetagSkip
                if cTimetag != self._timetag[aStand]+self._timetagSkip:
                    actStep = cTimetag - self._timetag[aStand]
                    if self.ignore_timetag_errors:
                        warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i}}" % (self._timetagSkip, t, p, actStep)), RuntimeWarning)
                    else:
                        raise RuntimeError(f"Invalid timetag skip encountered, expected {self._timetagSkip} on tuning {t}, pol {p}, but found {actStep}")
                        
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag - cFrame.header.time_offset
                    else:
                        setTime = cFrame.time
                        
                data_view[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data.view(data_view.dtype)
                count[aStand] +=  1
                self._timetag[aStand] = cTimetag
            nFrameSets += 1
            
        # If we've hit the end of the file and haven't read in enough frames, 
        # flush the buffer
        if eofFound and nFrameSets != frame_count:
            for cFrames in self.buffer.flush():
                for cFrame in cFrames:
                    b,t,p = cFrame.id
                    aStand = 2*(t-1) + p
                    cTimetag = cFrame.payload.timetag
                    if self._timetag[aStand] == 0:
                        self._timetag[aStand] = cTimetag - self._timetagSkip
                    if cTimetag != self._timetag[aStand]+self._timetagSkip:
                        actStep = cTimetag - self._timetag[aStand]
                        if self.ignore_timetag_errors:
                            warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i}}" % (self._timetagSkip, t, p, actStep)), RuntimeWarning)
                        else:
                            raise RuntimeError(f"Invalid timetag skip encountered, expected {self._timetagSkip} on tuning {t}, pol {p}, but found {actStep}")
                            
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag - cFrame.header.time_offset
                        else:
                            setTime = cFrame.time
                            
                    data_view[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data.view(data_view.dtype)
                    count[aStand] +=  1
                    self._timetag[aStand] = cTimetag
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets * 4096 / self.description['sample_rate']
            
        return duration, setTime, data
        
    def estimate_levels(self, nframe=100, sigma=5.0):
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
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Sample the data
        with FilePositionSaver(self.fh):
            count = {0:0, 1:0, 2:0, 3:0}
            data = np.zeros((4, nframe*4096))
            for i in range(nframe):
                for j in range(self.description['nbeampol']):
                    # Read in the next frame and anticipate any problems that could occur
                    try:
                        cFrame = drx.read_frame(self.fh, verbose=False)
                    except errors.EOFError:
                        break
                    except errors.SyncError:
                        continue
                        
                    b,t,p = cFrame.id
                    aStand = 2*(t-1) + p
                    
                    data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = np.abs( cFrame.payload.data )
                    count[aStand] +=  1
                    
        # Statistics
        rv = norm()
        frac = rv.cdf(sigma) - rv.cdf(-sigma)
        index = int(round(data.shape[1]*frac))
        if index == data.shape[1]:
            index = data.shape[1] - 1
        
        levels = [0, 0, 0, 0]
        for i in range(4):
            data2 = sorted(data[i,:])
            levels[i] = data2[index]
            
        return levels


class DRX8File(LDPFileBase):
    """
    Class to make it easy to interface with a DRX8 file.  DRX8 data consist of a
    time series of complex data a variable sample rate of up to 19.6 MHz from
    the beamformer.
    
    Methods defined for this class are:
      * get_info - Get information about the file's contents
      * get_remaining_frame_count - Get the number of frames remaining in the file
      * offset - Offset a specified number of seconds into the file
      * read_frame - Read and return a single `lsl.reader.drx.Frame` instance
      * read - Read a chunk of data in and return it as a numpy array
      * estimate_levels - Estimate the n-sigma level for the absolute value of the voltages 
    """
    
    def _ready_file(self):
        """
        Given an open file handle, find the start of valid DRX8 data.  This function:
        1) aligns on the first valid Mark 5C frame and
        2) skips over frames with a decimation of zero. 
        3) aligns the tuning/polarization timetags
        """
        
        # Align on the start of a Mark5C packet...
        while True:
            try:
                junkFrame = drx8.read_frame(self.fh)
                try:
                    # ... this really is DRX8 data
                    assert(junkFrame.is_drx8)
                    # ... that has a valid decimation
                    srate = junkFrame.sample_rate
                    # ... that it comes after 1980 (I don't know why this happens)
                    assert(junkFrame.payload.timetag > 61849368000000000)
                    break
                except (ZeroDivisionError, AssertionError):
                    pass
            except errors.SyncError:
                self.fh.seek(-drx8.FRAME_SIZE+1, 1)
                
        self.fh.seek(-drx8.FRAME_SIZE, 1)
        
        # Sometimes we get a file that has a large gap at the beginning for
        # reasons unknown.  If this seems to have happened jump over the gap
        # and issue a warning to the user.
        for checkpoint_size in (1, 2, 4, 8, 16):
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
                
            if filesize > checkpoint_size*1024**2:
                foffset = 0
                toffset = 0
                with FilePositionSaver(self.fh):
                    self.fh.seek(((checkpoint_size*1024**2-1)//drx8.FRAME_SIZE+1)*drx8.FRAME_SIZE, 1)
                    newFrame = drx8.read_frame(self.fh)
                    toffset = newFrame.time - junkFrame.time
                    if toffset > 300:
                        foffset = self.fh.tell() - drx8.FRAME_SIZE
                        
                if foffset > 0:
                    warnings.warn(colorfy("{{%%yellow Large (%.1f hr) gap at the beginning, skipping in %i B}}" % (toffset/3600, foffset)), RuntimeWarning)
                    self.fh.seek(foffset, 0)
                    break
                    
        # Line up the time tags for the various tunings/polarizations
        ids = []
        timetags = []
        for i in range(32):
            junkFrame = drx8.read_frame(self.fh)
            b,t,p = junkFrame.id
            id = (t,p)
            if id not in ids:
                ids.append(id)
            timetags.append(junkFrame.payload.timetag)
        self.fh.seek(-32*drx8.FRAME_SIZE, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the DRX8 file.
        """
        
        try:
            filesize = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            filesize = self.fh.size
        nFramesFile = (filesize - self.fh.tell()) // drx8.FRAME_SIZE
        tunepols = drx8.get_frames_per_obs(self.fh)
        tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepol
        bits = 8
        
        with FilePositionSaver(self.fh):
            beams = []
            tunes = []
            pols = []
            tuning1 = 0.0
            tuning2 = 0.0
            for i in range(32):
                try:
                    junkFrame0 = self.read_frame()
                    junkFrame1 = self.read_frame()
                except (errors.SyncError, errors.EOFError):
                    break
                for junkFrame in (junkFrame0, junkFrame1):
                    b,t,p = junkFrame.id
                    srate = junkFrame.sample_rate
                    if b not in beams:
                        beams.append(b)
                    if t not in tunes:
                        tunes.append(t)
                    if p not in pols:
                        pols.append(p)
                        
                    if t == 1:
                        tuning1 = junkFrame.central_freq
                    else:
                        tuning2 = junkFrame.central_freq
                        
                if i == 0:
                    start = junkFrame0.time
                    startRaw = junkFrame0.payload.timetag - junkFrame0.header.time_offset
                    
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': drx8.FRAME_SIZE,
                            'nbeampol': beampols, 'beam': b, 
                            'sample_rate': srate, 'data_bits': bits, 
                            'start_time': start, 'start_time_samples': startRaw, 'freq1': tuning1, 'freq2': tuning2}
                        
        # Initialize the buffer as part of the description process
        self.buffer = DRX8FrameBuffer(beams=beams, tunes=tunes, pols=pols, nsegments=LDP_CONFIG.get('drx_buffer_size'))
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open DRX8 file.  This function 
        returns the exact offset time.
        """
        
        # Figure out how far we need to offset inside the file
        junkFrame = drx8.read_frame(self.fh)
        self.fh.seek(-drx8.FRAME_SIZE, 1)
        
        # Get the initial time, sample rate, and beampols
        t0 = junkFrame.time
        if getattr(self, "_timetag", None) is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is not None:
                t0 = FrameTimestamp.from_dp_timetag(curr, junkFrame.header.time_offset)
        sample_rate = junkFrame.sample_rate
        beampols = drx8.get_frames_per_obs(self.fh)
        beampols = sum(beampols)
        
        # Offset in frames for beampols beam/tuning/pol. sets
        ioffset = int(offset * sample_rate / 4096 * beampols)
        ioffset = int(1.0 * ioffset / beampols) * beampols
        self.fh.seek(ioffset*drx8.FRAME_SIZE, 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate. 
        nattempt = 0
        diffs_used = deque([], 25)
        while True:
            junkFrame = drx8.read_frame(self.fh)
            self.fh.seek(-drx8.FRAME_SIZE, 1)
            
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            t1 = junkFrame.time
            sample_rate = junkFrame.sample_rate
            beampols = drx8.get_frames_per_obs(self.fh)
            beampols = sum(beampols)
            
            ## See how far off the current frame is from the target
            tDiff = (t1 - t0) - offset
            diffs_used.append(tDiff)
            
            ## Half that to come up with a new seek parameter
            tCorr   = -tDiff / 2.0
            cOffset = int(tCorr * sample_rate / 4096 * beampols)
            cOffset = int(1.0 * cOffset / beampols) * beampols
            ioffset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset == 0:
                break
            try:
                self.fh.seek(cOffset*drx8.FRAME_SIZE, 1)
                nattempt += 1
                assert(len(set(diffs_used)) > len(diffs_used)//4)
                assert(nattempt < 1000)
            except (IOError, AssertionError):
                warnings.warn(colorfy("{{%yellow Could not find the correct offset, giving up}}"), RuntimeWarning)

                break
                
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Zero out the time tag checker
        self._timetag = None
        
        return t1 - t0
        
    def read_frame(self, return_ci8=False):
        """
        Read and return a single `lsl.reader.drx8.Frame` instance.  If
        `return_ci8` is True then the frame will contain `lsl.reader.base.CI8`
        data instead of numpy.complex64 data.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Zero out the time tag checker
        self._timetagSkip = None
        self._timetag = None
        
        drx8_rf = drx8.read_frame_ci8 if return_ci8 else drx8.read_frame
        return drx8_rf(self.fh)
        
    def read(self, duration, time_in_samples=False, return_ci8=False):
        """
        Given an open DRX8 file and an amount of data to read in in seconds, read 
        in the data and return a three-element tuple of the actual duration read 
        in, the time for the first sample, and the data as numpy 
        array.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.ndp.fS` if the 
        `time_in_samples' keyword is set.
        
        If `return_ci8` is True then the data are returned will contain 
        `lsl.reader.base.CI8` data instead of numpy.complex64.  The two
        dimensions are input by samples.
        
        ..note::
            This function always returns an array with the first dimension
            holding four elements.  These elements contain, in order:
             * Tuning 1, polarization X
             * Tuning 1, polarization Y
             * Tuning 2, polarization X
             * Tuning 2, polarization Y
        """
        
        # Make sure there is file left to read
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Covert the sample rate to an expected timetag skip
        if getattr(self, "_timetagSkip", None) is None:
            self._timetagSkip = int(4096 / self.description['sample_rate'] * fS)
            
        # Setup the read_frame version to use
        drx8_rf = drx8.read_frame_ci8 if return_ci8 else drx8.read_frame
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = {0:0, 1:0, 2:0, 3:0}
            
        # Find out how many frames to read in
        frame_count = int(round(1.0 * duration * self.description['sample_rate'] / 4096))
        frame_count = frame_count if frame_count else 1
        
        # Setup the output arrays
        setTime = None
        if return_ci8:
            data = np.zeros((4,frame_count*4096), dtype=CI8)
            data_view = data.view(np.int16)
        else:
            data = np.zeros((4,frame_count*4096), dtype=np.complex64)
            data_view = data.view(np.float64)
            
        # Go!
        nFrameSets = 0
        eofFound = False
        count = {0:0, 1:0, 2:0, 3:0}
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            if not self.buffer.overfilled:
                cFrames = deque()
                for i in range(self.description['nbeampol']):
                    try:
                        cFrames.append( drx8_rf(self.fh, verbose=False) )
                    except errors.EOFError:
                        eofFound = True
                        self.buffer.append(cFrames)
                        cFrames = []
                        break
                    except errors.SyncError:
                        continue
                self.buffer.append(cFrames)
                
            cTimetag = self.buffer.peek()
            if cTimetag is None:
                # Continue adding frames if nothing comes out.
                continue
            else:
                # Otherwise, make sure we are on track
                aStand = 0
                if self._timetag[aStand] == 0:
                    pass
                elif cTimetag != self._timetag[aStand]+self._timetagSkip:
                    missing = (cTimetag - self._timetag[aStand] - self._timetagSkip) / float(self._timetagSkip)
                    if int(missing) == missing and missing < LDP_CONFIG.get('drx_autofill_size'):
                        ## This is kind of black magic down here
                        for m in range(int(missing)):
                            m = self._timetag[aStand] + self._timetagSkip*(m+1)
                            try:
                                baseframe = copy.deepcopy(cFrames[0])
                            except NameError:
                                baseframe = copy.deepcopy(self.buffer.buffer[cTimetag][0])
                            baseframe.payload.timetag = m
                            baseframe.payload._data *= 0
                            self.buffer.append(baseframe)
            cFrames = self.buffer.get()
            
            # If something comes out, add it to the data array
            for cFrame in cFrames:
                b,t,p = cFrame.id
                aStand = 2*(t-1) + p
                cTimetag = cFrame.payload.timetag
                if self._timetag[aStand] == 0:
                    self._timetag[aStand] = cTimetag - self._timetagSkip
                if cTimetag != self._timetag[aStand]+self._timetagSkip:
                    actStep = cTimetag - self._timetag[aStand]
                    if self.ignore_timetag_errors:
                        warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i}}" % (self._timetagSkip, t, p, actStep)), RuntimeWarning)
                    else:
                        raise RuntimeError(f"Invalid timetag skip encountered, expected {self._timetagSkip} on tuning {t}, pol {p}, but found {actStep}")
                        
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag - cFrame.header.time_offset
                    else:
                        setTime = cFrame.time
                        
                data_view[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data.view(data_view.dtype)
                count[aStand] +=  1
                self._timetag[aStand] = cTimetag
            nFrameSets += 1
            
        # If we've hit the end of the file and haven't read in enough frames, 
        # flush the buffer
        if eofFound and nFrameSets != frame_count:
            for cFrames in self.buffer.flush():
                for cFrame in cFrames:
                    b,t,p = cFrame.id
                    aStand = 2*(t-1) + p
                    cTimetag = cFrame.payload.timetag
                    if self._timetag[aStand] == 0:
                        self._timetag[aStand] = cTimetag - self._timetagSkip
                    if cTimetag != self._timetag[aStand]+self._timetagSkip:
                        actStep = cTimetag - self._timetag[aStand]
                        if self.ignore_timetag_errors:
                            warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i}}" % (self._timetagSkip, t, p, actStep)), RuntimeWarning)
                        else:
                            raise RuntimeError(f"Invalid timetag skip encountered, expected {self._timetagSkip} on tuning {t}, pol {p}, but found {actStep}")
                            
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag - cFrame.header.time_offset
                        else:
                            setTime = cFrame.time
                            
                    data_view[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data.view(data_view.dtype)
                    count[aStand] +=  1
                    self._timetag[aStand] = cTimetag
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets * 4096 / self.description['sample_rate']
            
        return duration, setTime, data
        
    def estimate_levels(self, nframe=100, sigma=5.0):
        """
        Estimate the n-sigma level for the absolute value of the voltages.  
        Returns a list with indicies corresponding to:
         0)  Tuning 1, X pol.
         1)  Tuning 1, Y pol.
         2)  Tuning 2, X pol.
         3)  Tuning 2, Y pol.
        
        ..note::
            The returned list always has four items, regardless of whether 
            or not the input DRX8 file has one or two tunings.
        """
        
        # Make sure there is file left to read
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Sample the data
        with FilePositionSaver(self.fh):
            count = {0:0, 1:0, 2:0, 3:0}
            data = np.zeros((4, nframe*4096))
            for i in range(nframe):
                for j in range(self.description['nbeampol']):
                    # Read in the next frame and anticipate any problems that could occur
                    try:
                        cFrame = drx8.read_frame(self.fh, verbose=False)
                    except errors.EOFError:
                        break
                    except errors.SyncError:
                        continue
                        
                    b,t,p = cFrame.id
                    aStand = 2*(t-1) + p
                    
                    data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = np.abs( cFrame.payload.data )
                    count[aStand] +=  1
                    
        # Statistics
        rv = norm()
        frac = rv.cdf(sigma) - rv.cdf(-sigma)
        index = int(round(data.shape[1]*frac))
        if index == data.shape[1]:
            index = data.shape[1] - 1
        
        levels = [0, 0, 0, 0]
        for i in range(4):
            data2 = sorted(data[i,:])
            levels[i] = data2[index]
            
        return levels


class DRSpecFile(LDPFileBase):
    """
    Class to make it easy to interface with a DR Spectrometer file.  DR
    Spectrometer data contain DRX data that has been transformed to the Fourier
    domain, converted to power, and integrated on-the-fly by the data recorder. 
    These data can have various integration times, channel counts, and
    polarization products stored.
    
    Methods defined for this class are:
     * get_info - Get information about the file's contents
     * get_remaining_frame_count - Get the number of frames remaining in the file
     * offset - Offset a specified number of seconds into the file
     * read_frame - Read and return a single `lsl.reader.drspec.Frame` instance
     * read - Read a chunk of data in and return it as a numpy array
    """
    
    def _ready_file(self):
        """
        Ready the DRSpec file.
        """
        
        # Align on the start of a Mark5C packet...
        while True:
            try:
                junkFrame = drspec.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(1, 1)
        self.fh.seek(-drspec.get_frame_size(self.fh), 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the DRSpec file.
        """
        
        try:
            filesize = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            filesize = self.fh.size
        FRAME_SIZE = drspec.get_frame_size(self.fh)
        nFramesFile = filesize // FRAME_SIZE
        LFFT = drspec.get_transform_size(self.fh)
        with FilePositionSaver(self.fh):
            junkFrame = drspec.read_frame(self.fh)
            
        bits = 32
        beam = junkFrame.id
        beampols = 4
        srate = junkFrame.sample_rate
        nInt = junkFrame.header.nints
        tInt = nInt*LFFT/srate
        start = junkFrame.time
        startRaw = junkFrame.payload.timetag - junkFrame.header.time_offset
        tuning1, tuning2 = junkFrame.central_freq
        prod = junkFrame.data_products
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': FRAME_SIZE, 
                            'nbeampol': beampols, 'beam': beam, 
                            'sample_rate': srate, 'data_bits': bits, 
                            'start_time': start, 'start_time_samples': startRaw, 'freq1': tuning1, 'freq2': tuning2, 
                            'nint': nInt, 'tint': tInt, 'LFFT': LFFT, 
                            'nproduct': len(prod), 'data_products': prod}
                        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open DR spectrometer file.  This 
        function returns the exact offset time.
        """
        
        # Gather some basic information and read in the first frame
        junkFrame = drspec.read_frame(self.fh)
        self.fh.seek(-self.description['frame_size'], 1)
        
        # Get the initial time, sample rate, and integration time
        t0 = junkFrame.time
        
        # Offset in frames for beampols beam/tuning/pol. sets
        ioffset = int(round(offset / self.description['tint']))
        self.fh.seek(ioffset*self.description['frame_size'], 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.
        diffs_used = deque([], 25)
        while True:
            junkFrame = drspec.read_frame(self.fh)
            self.fh.seek(-self.description['frame_size'], 1)
            
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            t1 = junkFrame.time
            sample_rate = junkFrame.sample_rate
            LFFT = junkFrame.transform_size
            tInt = junkFrame.header.nints*LFFT/sample_rate
            
            ## See how far off the current frame is from the target
            tDiff = t1 - (t0 + offset)
            diffs_used.append(tDiff)
            
            ## Half that to come up with a new seek parameter
            tCorr   = -tDiff / 2.0
            cOffset = int(round(tCorr / tInt))
            ioffset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset == 0:
                break
            try:
                self.fh.seek(cOffset*self.description['frame_size'], 1)
                assert(len(set(diffs_used)) > len(diffs_used)//4)
            except (IOError, AssertionError):
                warnings.warn(colorfy("{{%yellow Could not find the correct offset, giving up}}"), RuntimeWarning)

                break
                
        # Update the file metadata
        self._describe_file()
        
        # Zero out the timetag checker
        self._timetag = None
        
        return t1 - t0
        
    def read_frame(self):
        """
        Read and return a single `lsl.reader.drspec.Frame` instance.
        """
        
        # Update the timetag checker
        self._timetag = None
        
        return drspec.read_frame(self.fh)
        
    def read(self, duration, time_in_samples=False):
        """
        Given an open DR spectrometer file and an amount of data read in in 
        seconds, read in the data and return a three-element tuple of the actual 
        duration read in, the times at the beginning of each stream, and the 
        data as numpy array.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.ndp.fS` if the 
        `time_in_samples' keyword is set.
        
        ..note::
            This function always returns a 3-D array with the first dimension
            indexing over data product, the second over time and the third over
            frequency channel.
        """
        
        # Make sure there is file left to read
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Covert the sample rate to an expected timetag skip
        timetagSkip = self.description['tint']
        
        # Setup the counter variables:  frame count and time tag count
        count = 0
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            junkFrame = drspec.read_frame(self.fh)
            self._timetag = junkFrame.time[0] + (junkFrame.time[1] - timetagSkip)
            self.fh.seek(-self.description['frame_size'], 1)
        
        # Find out how many frames to read in
        frame_count = int(round(1.0 * duration / self.description['tint']))
        frame_count = frame_count if frame_count else 1
        
        # Setup the output arrays
        data = np.zeros((2*self.description['nproduct'],frame_count,self.description['LFFT']), dtype=np.float32)
        
        # Go!
        nFrameSets = 0
        setTime = None
        for i in range(frame_count):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = drspec.read_frame(self.fh, verbose=False)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue
                
            cTimetag = cFrame.time
            if cTimetag > self._timetag + 1.001*timetagSkip:
                actStep = cTimetag - self._timetag
                if self.ignore_timetag_errors:
                    warnings.warn("Invalid timetag skip encountered, expected %i but found %i" % (timetagSkip, actStep), RuntimeWarning)
                else:
                    raise RuntimeError(f"Invalid timetag skip encountered, expected {timetagSkip} but found {actStep}")
                    
            if setTime is None:
                if time_in_samples:
                    setTime = cFrame.payload.timetag - cFrame.header.time_offset
                else:
                    setTime = cFrame.time
                    
            for j,p in enumerate(self.description['data_products']):
                data[j+0,                             count, :] = getattr(cFrame.payload, '%s0' % p, None)
                data[j+self.description['nproduct'], count, :] = getattr(cFrame.payload, '%s1' % p, None)
            count +=  1
            self._timetag = cTimetag
            nFrameSets += 1
            
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets * self.description['tint']
        
        return duration, setTime, data


class TBXFile(LDPFileBase):
    """
    Class to make it easy to interface with a TBT/TBS file.  TBF data are a complex
    frequency-domain product that contains blocks of channels from all antennas
    in the array.  Each channel has a bandwidth of f\ :sub:`C` (25 kHz) and
    there may be up to 3584 channels within a single recording.  The stand
    ordering is based on the input into the digital system rather than the stand
    number in the array.  
    
    Methods defined for this class are:
     * get_info - Get information about the file's contents
     * get_remaining_frame_count - Get the number of frames remaining in the file
     * offset - Offset a specified number of seconds into the file
     * read_frame - Read and return a single `lsl.reader.tbw.Frame` instance
     * read - Read in the capture and return it as a numpy array
    """
    
    def _ready_file(self):
        """
        Find the start of valid TBF data.  This function:
        1) Aligns on the first valid Mark 5C frame.
        """
        
        # Align on the start of a Mark5C packet
        while True:
            try:
                tbx.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(1, 1)
                
        # Skip over any DRX frames the start of the file
        i = 0
        while True:
            try:
                mark = self.fh.tell()
                tbx.read_frame(self.fh)
                frame_size = self.fh.tell() - mark
                break
            except errors.SyncError:
                i += 1
                self.fh.seek(1, 1)
        if i == 0:
            self.fh.seek(-frame_size, 1)
        self.fh.seek(-frame_size, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the TBF file.
        """
        
        with FilePositionSaver(self.fh):
            # Read in frame
            mark = self.fh.tell()
            junkFrame = tbx.read_frame(self.fh)
            frame_size = self.fh.tell() - mark
            self.fh.seek(-frame_size, 1)
            
            # Basic file information
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
            nFramesFile = (filesize - self.fh.tell()) // frame_size
            bits = 4
            nFramesPerObs = tbx.get_frames_per_obs(self.fh)
            nchan = tbx.get_channel_count(self.fh)
            firstFrameCount = tbx.get_first_frame_count(self.fh)
            
            # Pre-load the channel mapper
            self.mapper = tbx.get_first_channel(self.fh, all_frames=True)
            
            # Check for contiguous frequency coverage
            chan_steps = np.diff(self.mapper)
            channel_count = np.median(chan_steps)
            if not all(chan_steps == channel_count):
                bad_steps = np.where(chan_steps != channel_count)[0]
                warnings.warn(colorfy("{{%%yellow File appears to contain %i frequency gap(s) of size %s channels}}" % (len(bad_steps), ','.join([str(chan_steps[g]) for g in bad_steps]))), RuntimeWarning)
                
            # Find the "real" starttime
            while junkFrame.header.frame_count != firstFrameCount:
                junkFrame = tbx.read_frame(self.fh)
            start = junkFrame.time
            startRaw = junkFrame.payload.timetag
            
        # Calculate the frequencies
        freq = np.zeros(nchan)
        for i,c in enumerate(self.mapper):
            freq[i*channel_count:(i+1)*channel_count] = c + np.arange(channel_count)
        freq *= srate
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': frame_size,
                            'sample_rate': srate, 'data_bits': bits,
                            'nantenna': nstand*2, 'nchan': nchan, 'freq1': freq, 'start_time': start, 
                            'start_time_samples': startRaw, 'frame_channel_count': channel_count}
                        
        # Initialize the buffer as part of the description process
        self.buffer = TBXFrameBuffer(chans=self.mapper, nsegments=LDP_CONFIG.get('tbx_buffer_size'))
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open TBF file.  This function 
        returns the exact offset time.
        """
        
        frame_size = self.description['frame_size']
        
        # Find out where we really are taking into account the buffering
        buffer_offset = 0
        if getattr(self, "_timetag", None) is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is None:
                frame = tbx.read_frame(self.fh)
                self.fh.seek(-frame_size, 1)
                curr = frame.payload.time_tag
            buffer_offset = curr - self._timetag
            buffer_offset = buffer_offset / fS
            
        offset = offset - buffer_offset
        framesPerObs = self.description['nchan'] // self.description['frame_channel_count']
        frameOffset = int(offset * self.description['sample_rate'] * framesPerObs)
        frameOffset = int(1.0 * frameOffset / framesPerObs) * framesPerObs
        self.fh.seek(frameOffset*frame_size, 1)
        
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
        
        # Reset the timetag checker
        self._timetag = None
        
        return 1.0 * frameOffset / framesPerObs / self.description['sample_rate']
        
    def read_frame(self, return_ci8=False):
        """
        Read and return a single `lsl.reader.tbx.Frame` instance.  If
        `return_ci8` is True then the frame will contain `lsl.reader.base.CI8`
        data instead of numpy.complex64 data.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        tbx_rf = tbx.read_frame_ci8 if return_ci8 else tbx.read_frame
        return tbx_rf(self.fh)
        
    def read(self, duration, time_in_samples=False, return_ci8=False):
        """
        Given an amount of data to read in in seconds, read in the data from a
        TBT/TBS capture.  This function returns a three-element tuple with
        elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 3-D Numpy array of data.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common..ndp.fS` if the 
        `time_in_samples' keyword is set.
        
        The sorting order of the output data array is by 
        digitizer number - 1.
        """ 
        
        # Make sure there is file left to read
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Covert the sample rate to an expected timetag skip
        timetagSkip = int(1.0 / self.description['sample_rate'] * fS)
        
        # Setup the read_frame version to use
        tbx_rf = tbx.read_frame_ci8 if return_ci8 else tbx.read_frame
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            
        # Find out how many frames to read in
        framesPerObs = self.description['nchan'] // self.description['frame_channel_count']
        frame_count = int(round(1.0 * duration * self.description['sample_rate']))
        frame_count = frame_count if frame_count else 1
        duration = frame_count / self.description['sample_rate']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in range(framesPerObs)]
        if return_ci8:
            data = np.zeros((self.description['nantenna'], self.description['nchan'], frame_count), dtype=CI8)
            data_view = data.view(np.int16)
        else:
            data = np.zeros((self.description['nantenna'], self.description['nchan'], frame_count), dtype=np.complex64)
            data_view = data.view(np.float64)

        nSkip = 0
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            cFrames = deque()
            for i in range(framesPerObs):
                try:
                    cFrame = tbx_rf(self.fh, verbose=False)
                    if not cFrame.is_tbx:
                        nSkip += 1
                        continue
                    cFrames.append( cFrame )
                    nSkip = 0
                except errors.EOFError:
                    eofFound = True
                    self.buffer.append(cFrames)
                    cFrames = []
                    break
                except errors.SyncError:
                    nSkip += 1
                    if nSkip > 40000:
                        eofFound = True
                        self.buffer.append(cFrames)
                        cFrames = []
                        break
                    self.fh.seek(drx.FRAME_SIZE, 1)
                    continue
                    
            self.buffer.append(cFrames)
            cFrames = self.buffer.get()
            
            # Continue adding frames if nothing comes out.
            if cFrames is None:
                continue
                
            # If something comes out, add it to the data array
            cTimetag = cFrames[0].payload.timetag
            if self._timetag == 0:
                self._timetag = cTimetag - timetagSkip
            if cTimetag != self._timetag+timetagSkip:
                actStep = cTimetag - self._timetag
                if self.ignore_timetag_errors:
                    warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i, but found %i}}" % (timetagSkip, actStep)), RuntimeWarning)
                else:
                    raise RuntimeError(f"Invalid timetag skip encountered, expected {timetagSkip}, but found {actStep}")
            self._timetag = cFrames[0].payload.timetag
            
            for cFrame in cFrames:
                first_chan = cFrame.header.first_chan
                
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag
                    else:
                        setTime = cFrame.time
                        
                subData = cFrame.payload.data
                subData = subData.reshape(self.description['frame_channel_count'],-1)
                subData = subData.T
                
                aStand = self.mapper.index(first_chan)
                data_view[:,aStand*self.description['frame_channel_count']:(aStand+1)*self.description['frame_channel_count'],count[aStand]] = subData.view(data_view.dtype)
                count[aStand] += 1
            nFrameSets += 1
            
        # If we've hit the end of the file and haven't read in enough frames, 
        # flush the buffer
        if eofFound and nFrameSets != frame_count:
            for cFrames in self.buffer.flush():
                cTimetag = cFrames[0].payload.timetag
                if self._timetag == 0:
                    self._timetag = cTimetag - timetagSkip
                if cTimetag != self._timetag+timetagSkip:
                    actStep = cTimetag - self._timetag
                    if self.ignore_timetag_errors:
                        warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i, but found %i}}" % (timetagSkip, actStep)), RuntimeWarning)
                    else:
                        raise RuntimeError(f"Invalid timetag skip encountered, expected {timetagSkip}, but found {actStep}")
                self._timetag = cFrames[0].payload.timetag
                
                for cFrame in cFrames:
                    first_chan = cFrame.header.first_chan
                    
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag
                        else:
                            setTime = cFrame.time
                        
                    subData = cFrame.payload.data
                    subData = subData.reshape(self.description['frame_channel_count'],-1)
                    subData = subData.T
                    
                    aStand = self.mapper.index(first_chan)
                    data_view[:,aStand*self.description['frame_channel_count']:(aStand+1)*self.description['frame_channel_count'],count[aStand]] = subData.view(data_view.dtype)
                    count[aStand] += 1
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Sanity check at the end to see if we actually read anything.  
        if nFrameSets == 0 and duration > 0:
            raise errors.EOFError()
            
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets  / self.description['sample_rate']
        
        return duration, setTime, data


class CORFile(LDPFileBase):
    """
    Class to make it easy to interface with a COR file.  COR data contain full
    polarization complex visibility data from each baseline pair in the array.  
    Each channel has a bandwidth of f\ :sub:`C` (25 kHz) and there may be up to
    792 channels within a single recording.  The stand numbering for the
    baseline pair is based on the input into the digital system rather than the
    stand number in the array.
    
    Methods defined for this class are:
     * get_info - Get information about the file's contents
     * get_remaining_frame_count - Get the number of frames remaining in the file
     * offset - Offset a specified number of seconds into the file
     * read_frame - Read and return a single `lsl.reader.tbw.Frame` instance
     * read - Read in the capture and return it as a numpy array
    """
    
    def _ready_file(self):
        """
        Find the start of valid COR data.  This function:
        1) Aligns on the first valid Mark 5C frame.
        """
        
        # Align on the start of a Mark5C packet
        while True:
            try:
                cor.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(1, 1)
                
        # Skip over any DRX frames the start of the file
        i = 0
        while True:
            try:
                mark = self.fh.tell()
                cor.read_frame(self.fh)
                frame_size = self.fh.tell() - mark
                break
            except errors.SyncError:
                i += 1
                self.fh.seek(1, 1)
        if i == 0:
            self.fh.seek(-frame_size, 1)
        self.fh.seek(-frame_size, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the COR file.
        """
        
        # Read in frame
        with FilePositionSaver(self.fh):
            mark = self.fh.tell()
            junkFrame = cor.read_frame(self.fh)
            frame_size = self.fh.tell() - mark
            self.fh.seek(-frame_size, 1)
            
            # Basic file information
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
            nFramesFile = (filesize - self.fh.tell()) // frame_size
            adp_id = junkFrame.adp_id
            srate = fC
            if adp_id & 0x04:
                srate = fC
            bits = 32
            nFramesPerObs = cor.get_frames_per_obs(self.fh)
            nchan = cor.get_channel_count(self.fh)
            nBaseline = cor.get_baseline_count(self.fh)
            
            # Pre-load the baseline mapper
            # NOTE: This is done with a dictionary rather than a list since 
            #       the look-ups are much faster
            self.bmapperd = {}
            k = 0
            for i in range(1, 256+1):
                for j in range(i, 256+1):
                    self.bmapperd[(i,j)] = k
                    k += 1
                    
            # Pre-load the channel mapper
            self.cmapper = cor.get_first_channel(self.fh, all_frames=True)
            self.cmapper.sort()
            
        # Create a channel mapper dictionary
        self.cmapperd = {}
        for i,c in enumerate(self.cmapper):
            self.cmapperd[c] = i
            
        # Calculate the frequencies
        chan_steps = np.diff(self.mapper)
        channel_count = np.median(chan_steps) // 4
        freq = np.zeros(nchan)
        for i,c in enumerate(self.cmapper):
            freq[i*channel_count:(i+1)*channel_count] = c + np.arange(channel_count)
        freq *= srate
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': frame_size,
                            'sample_rate': srate, 'data_bits': bits, 
                            'nantenna': 512, 'nchan': nchan, 'freq1': freq, 'start_time': start, 
                            'start_time_samples': startRaw, 'nbaseline': nBaseline, 'frame_channel_count': channel_count, 'tint':cFrame.integration_time}
                        
        # Initialize the buffer as part of the description process
        self.buffer = CORFrameBuffer(chans=self.cmapper, reorder=False, nsegments=LDP_CONFIG.get('cor_buffer_size'))
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open COR file.  This function 
        returns the exact offset time.
        
        .. note::
            The offset provided by this function is relatively crude due to the
            structure of COR files.
            
        .. versionchanged:: 1.2.4
            Offsets are now relative to the current location in the file rather
            than to the start of the file
        """
        
        # Find out where we really are taking into account the buffering
        buffer_offset = 0
        if getattr(self, "_timetag", None) is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is None:
                frame = cor.read_frame(self.fh)
                self.fh.seek(-cor.FRAME_SIZE, 1)
                curr = frame.payload.time_tag
            buffer_offset = curr - self._timetag
            buffer_offset = buffer_offset / fS
            
        offset = offset - buffer_offset
        framesPerObs = self.description['nchan'] // self.description['frame_channel_count'] * self.description['nbaseline']
        frameOffset = int(offset / self.description['tint'] * framesPerObs)
        frameOffset = int(1.0 * frameOffset / framesPerObs) * framesPerObs
        self.fh.seek(frameOffset*cor.FRAME_SIZE, 1)
        
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
        
        # Reset the timetag checker
        self._timetag = None
        
        return 1.0 * frameOffset / framesPerObs * self.description['tint']
        
    def read_frame(self):
        """
        Read and return a single `lsl.reader.cor.Frame` instance.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        return cor.read_frame(self.fh)
        
    def read(self, duration=None, time_in_samples=False):
        """
        Read and return the entire COR capture.  This function returns 
        a three-element tuple with elements of:
        0) the actual duration of data read in, 
        1) the time tag for the first sample, and
        2) a 5-D Numpy array of data.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.ndp.fS` if the 
        `time_in_samples' keyword is set.
        
        The sorting order of the output data array is by 
        baseline.
        """ 
        
        # Make sure there is file left to read
        try:
            curr_size = os.fstat(self.fh.fileno()).st_size
        except AttributeError:
            curr_size = self.fh.size
        if self.fh.tell() == curr_size:
            try:
                if self.buffer.is_empty():
                    raise errors.EOFError()
            except AttributeError:
                raise errors.EOFError()
                
        # Covert the sample rate to an expected timetag skip
        timetagSkip = int(self.description['tint'] * fS)
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            
        # Find out how many frames to read in
        if duration is None:
            duration = self.description['nframe'] / framesPerObs * self.description['tint']
        framesPerObs = self.description['nchan'] // self.description['frame_channel_count'] * self.description['nbaseline']
        frame_count = int(round(1.0 * duration / self.description['tint']))
        frame_count = frame_count if frame_count else 1
        duration = frame_count * self.description['tint']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in range(framesPerObs)]
        data = np.zeros((self.description['nbaseline'], self.description['nchan'], 2, 2, frame_count), dtype=np.complex64)
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            cFrames = deque()
            for i in range(framesPerObs):
                try:
                    cFrame = cor.read_frame(self.fh, verbose=False)
                    if not cFrame.is_cor:
                        continue
                    cFrames.append( cFrame )
                except errors.EOFError:
                    eofFound = True
                    self.buffer.append(cFrames)
                    cFrames = []
                    break
                except errors.SyncError:
                    continue
                    
            self.buffer.append(cFrames)
            cFrames = self.buffer.get()
            
            # Continue adding frames if nothing comes out.
            if cFrames is None:
                continue
                
            # If something comes out, add it to the data array
            cTimetag = cFrames[0].payload.timetag
            if self._timetag == 0:
                self._timetag = cTimetag - timetagSkip
            if cTimetag != self._timetag+timetagSkip:
                actStep = cTimetag - self._timetag
                if self.ignore_timetag_errors:
                    warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i, but found %i}}" % (timetagSkip, actStep)), RuntimeWarning)
                else:
                    raise RuntimeError(f"Invalid timetag skip encountered, expected {timetagSkip}, but found {actStep}")
            self._timetag = cFrames[0].payload.timetag
            
            for cFrame in cFrames:
                first_chan = cFrame.header.first_chan
                
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag
                    else:
                        setTime = cFrame.time
                        
                aBase = self.bmapperd[cFrame.id]
                aChan = self.cmapperd[first_chan]
                aStand = aBase*len(self.cmapper) + aChan
                data[aBase,aChan*self.description['frame_channel_count']:(aChan+1)*self.description['frame_channel_count'],:,:,count[aStand]] = cFrame.payload.data
                count[aStand] += 1
            nFrameSets += 1
            
        # If we've hit the end of the file and haven't read in enough frames, 
        # flush the buffer
        if eofFound and nFrameSets != frame_count:
            for cFrames in self.buffer.flush():
                cTimetag = cFrames[0].payload.timetag
                if self._timetag == 0:
                    self._timetag = cTimetag - timetagSkip
                if cTimetag != self._timetag+timetagSkip:
                    actStep = cTimetag - self._timetag
                    if self.ignore_timetag_errors:
                        warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i, but found %i}}" % (timetagSkip, actStep)), RuntimeWarning)
                    else:
                        raise RuntimeError(f"Invalid timetag skip encountered, expected {timetagSkip}, but found {actStep}")
                self._timetag = cFrames[0].payload.timetag
                
                for cFrame in cFrames:
                    first_chan = cFrame.header.first_chan
                    
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag
                        else:
                            setTime = cFrame.time
                            
                    aBase = self.bmapperd[cFrame.id]
                    aChan = self.cmapperd[first_chan]
                    aStand = aBase*len(self.cmapper) + aChan
                    data[aBase,aChan*self.description['frame_channel_count']:(aChan+1)*self.description['frame_channel_count'],:,:,count[aStand]] = cFrame.payload.data
                    count[aStand] += 1
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Sanity check at the end to see if we actually read anything.  
        # This is needed because of how COR and DRX interact where COR
        # files can be padded at the end with DRX data
        if nFrameSets == 0 and duration > 0:
            raise errors.EOFError()
            
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets  * self.description['tint']
        
        return duration, setTime, data


def LWADataFile(filename=None, fh=None, ignore_timetag_errors=False, buffering=-1):
    """
    Wrapper around the various LWA data reader classes defined here that takes
    a file, determines the data type, and initializes and returns the 
    appropriate LDP class.
    """
    
    # Open the file as appropriate
    is_splitfile = False
    if fh is None:
        fh = open(filename, 'rb')
    else:
        filename = fh.name
        if not isinstance(fh, SplitFileWrapper):
            if fh.mode.find('b') == -1:
                fh.close()
                fh = open(filename, 'rb')
        else:
            is_splitfile = True
            
    # Read a bit of data to try to find the right type
    for mode in (drx, drx8, drspec, tbx, cor):
        ## Set if we find a valid frame marker
        foundMatch = False
        ## Set if we can read more than one valid successfully
        foundMode = False
        
        ## Sort out the frame size.  This is tricky because DR spectrometer files
        ## have frames of different sizes depending on the mode
        if mode in (drspec, tbx, cor):
            try:
                mfs = drspec.get_frame_size(fh)
            except:
                mfs = 0
        else:
            mfs = mode.FRAME_SIZE
            
        ## Loop over the frame size to try and find what looks like valid data.  If
        ## is is found, set 'foundMatch' to True.
        for i in range(mfs):
            try:
                junkFrame = mode.read_frame(fh)
                foundMatch = True
                break
            except errors.EOFError:
                break
            except errors.SyncError:
                fh.seek(-mfs+1, 1)
                
        ## Did we strike upon a valid frame?
        if foundMatch:
            ### Is so, we now need to try and read more frames to make sure we have 
            ### the correct type of file
            fh.seek(-mfs, 1)
            
            try:
                for i in range(2):
                    junkFrame = mode.read_frame(fh)
                foundMode = True
            except errors.EOFError:
                break
            except errors.SyncError:
                ### Reset for the next mode...
                fh.seek(0)
        else:
            ### Reset for the next mode...
            fh.seek(0)
            
        ## Did we read more than one valid frame?
        if foundMode:
            break
            
    fh.seek(0)
    if not is_splitfile:
        fh.close()
        fh = None
        
    # Raise an error if nothing is found
    if not foundMode:
        raise RuntimeError(f"File '{filename}' does not appear to be a valid LWA1 data file")
        
    # Otherwise, build and return the correct LDPFileBase sub-class
    if mode == drx:
        ldpInstance = DRXFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == drx8:
        ldpInstance = DRX8File(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == tbx:
        ldpInstance = TBXFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == cor:
        ldpInstance = CORFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    else:
        ldpInstance = DRSpecFile(filename=filename, fh=fh,
                                 ignore_timetag_errors=ignore_timetag_errors,
                                 buffering=buffering)
        
    # Done
    return ldpInstance


import atexit
atexit.register(_open_ldp_files.close_all)
