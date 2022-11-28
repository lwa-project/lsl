"""
LWA Development Primitives - A set of utilities that should make developing 
new analysis software easier.  These functions wrap the nitty gritty of the 
file reading and unpacking behind Python objects.

Data format objects included are:
  * TBWFile
  * TBNFile
  * DRXFile
  * DRSpecFile
  * TBFFile
  * CORFILE

Also included are the LWA1DataFile, LWASVDataFile, and LWADataFile functions 
that take a filename and try to determine the correct data format object to
use.

.. versionchanged:: 1.2.0
    Added support for LWA-SV ADP data
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import abc
import copy
import numpy
import warnings
from textwrap import fill as tw_fill
from scipy.stats import norm
from collections import deque, defaultdict

from lsl.common.dp import fS
from lsl.common.adp import fC
from lsl.reader import tbw, tbn, drx, drspec, tbf, cor, errors
from lsl.reader.buffer import TBNFrameBuffer, DRXFrameBuffer, TBFFrameBuffer, CORFrameBuffer
from lsl.reader.utils import *
from lsl.reader.base import FrameTimestamp, CI8
from lsl.common.color import colorfy

from lsl.config import LSL_CONFIG
LDP_CONFIG = LSL_CONFIG.view('ldp')

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.5'
__all__ = ['TBWFile', 'TBNFile', 'DRXFile', 'DRSpecFile', 'TBFFile', 'LWA1DataFile', 
           'LWASVDataFile', 'LWADataFile']


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
        return "%s @ %s" % (type(self).__name__, self.filename)
        
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
                raise ValueError("Unknown key '%s'" % key)
                
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


class TBWFile(LDPFileBase):
    """
    Class to make it easy to interface with a TBW file.  Methods defined for this class are:
     * get_info - Get information about the file's contents
     * get_remaining_frame_count - Get the number of frames remaining in the file
     * read_frame - Read and return a single `lsl.reader.tbw.Frame` instance
     * read - Read in the capture and return it as a numpy array
    """
    
    def _ready_file(self):
        """
        Find the start of valid TBW data.  This function:
         1) Aligns on the first valid Mark 5C frame and
         2) Skips over any TBN frames at the beginning of the file.
        """
        
        # Align on the start of a Mark5C packet
        while True:
            try:
                junkFrame = tbw.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(-tbw.FRAME_SIZE+1, 1)
                
        # Jump over any TBN data in the file
        while not junkFrame.header.is_tbw:
            try:
                junkFrame = tbw.read_frame(self.fh)
            except errors.SyncError:
                ## If we reached this then we are probably in an old TBW file that has
                ## a bunch of TBN frames at the beginning.  We need to seek backwards,
                ## realign on the sync word, and read forwards again.
                
                ## Jump back a TBW frame
                self.fh.seek(-tbw.FRAME_SIZE, 1)
                
                ## Find the sync word again
                while True:
                    try:
                        tbn.read_frame(self.fh)
                        break
                    except errors.SyncError:
                        self.fh.seek(-tbn.FRAME_SIZE+1, 1)
                        
                ## Find the end of the TBN data
                while True:
                    try:
                        tbn.read_frame(self.fh)
                    except errors.SyncError:
                        break
                self.fh.seek(-2*tbn.FRAME_SIZE, 1)
                junkFrame = tbw.read_frame(self.fh)
        self.fh.seek(-tbw.FRAME_SIZE, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the TBW file.
        """
        
        with FilePositionSaver(self.fh):
            junkFrame = self.read_frame()
            self.fh.seek(-tbw.FRAME_SIZE, 1)
            
            # Basic file information
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
            nFramesFile = (filesize - self.fh.tell()) // tbw.FRAME_SIZE
            srate = 196e6
            bits = junkFrame.data_bits
            
            # Trick to figure out how many antennas are in a file and the "real" 
            # start time.  For details of why this needs to be done, see the read()
            # function below.
            idsFound = []
            timesFound = []
            filePosRef = self.fh.tell()
            while True:
                try:
                    for i in range(26):
                        frame = tbw.read_frame(self.fh)
                        while not frame.header.is_tbw:
                            frame = tbw.read_frame(self.fh)
                        stand = frame.id
                        if stand not in idsFound:
                            idsFound.append(stand)
                        if frame.header.frame_count < 1000:
                            timesFound.append( (frame.header.frame_count-1, frame.payload.timetag) )
                    self.fh.seek(tbw.FRAME_SIZE*(30000-26), 1)
                except:
                    break
                    
        # What is that start time again?
        startTimeTag = None
        for fc,tt in timesFound:
            tt = tt - fc*(1200 if bits == 4 else 400)
            if startTimeTag is None or tt < startTimeTag:
                startTimeTag = tt
        start = startTimeTag / fS
        startRaw = startTimeTag
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': tbw.FRAME_SIZE,
                            'sample_rate': srate, 'data_bits': bits, 'nantenna': 2*len(idsFound), 
                            'start_time': start, 'start_time_samples': startRaw}
                        
    def read_frame(self):
        """
        Read and return a single `lsl.reader.tbw.Frame` instance.
        """
        
        frame = tbw.read_frame(self.fh)
        while not frame.header.is_tbw:
            frame = tbw.read_frame(self.fh)
            
        return frame
        
    def read(self, duration=None, time_in_samples=False):
        """
        Read and return the entire TBW capture.  This function returns 
        a three-element tuple with elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 2-D Numpy array of data.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.dp.fS` if the 
        `time_in_samples' keyword is set.
        
        The sorting order of the output data array is by 
        digitizer number - 1.
        
        .. note::
            Setting the 'duration' keyword has no effect on the read 
            process because the entire capture is always read in.
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
                
        # Get the data frame size
        dataSize = 400
        if self.description['data_bits'] == 4:
            dataSize = 1200
            
        # Find out how many frames to work with at a time
        nFrames = int(30000)
        
        # Initialize the time variables
        # Explination:
        #   This is needed to work out what the "real" start time is of the 
        #   capture due to buffering in the data recorder.  What can happen 
        #   is that the last ~4 MB of a previous capture could be stuck in 
        #   the data recoder's buffer and that the buffer won't get dumped 
        #   until the next capture is launch.  Thus, you can end up in a 
        #   situation where the first few valid TBW frames in a file are from 
        #   the previous capture.
        #   
        #   To get around this we use the frame count-correction time tag of 
        #   the lowest frame number found.  This skips over the trailing edge of 
        #   the previous capture (which should have a high frame count) while
        #   allowing the code to deal with files that may be missing the first
        #   frame from the first board to send a frame.
        setTime = None
        setTimeRef = 1000
        
        # Initialize the output data array
        data = numpy.zeros((self.description['nantenna'], nFrames*dataSize), dtype=numpy.int16)
        
        # Read in the next frame and anticipate any problems that could occur
        i = 0
        while i < ((self.description['nantenna']//2)*nFrames):
            try:
                cFrame = tbw.read_frame(self.fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue
                
            if not cFrame.header.is_tbw:
                continue
                
            stand = cFrame.header.id
            aStandX = 2*(stand-1) + 0
            aStandY = 2*(stand-1) + 1
            
            if cFrame.header.frame_count < setTimeRef:
                newSetTime = cFrame.payload.timetag - (cFrame.header.frame_count-1)*dataSize
                if setTime is None or cFrame.payload.timetag < setTime:
                    setTime = newSetTime
                    setTimeRef = cFrame.header.frame_count
                    
            try:
                cnt = cFrame.header.frame_count - 1
                data[aStandX, cnt*dataSize:(cnt+1)*dataSize] = cFrame.payload.data[0,:]
                data[aStandY, cnt*dataSize:(cnt+1)*dataSize] = cFrame.payload.data[1,:]
                
                i += 1
            except ValueError:
                pass
                
        # Deal with the time if we don't want it in samples
        if not time_in_samples:
            setTime = setTime / fS
            
        # Calculate the duration
        duration = data.shape[1]/self.description['sample_rate']
        
        return duration, setTime, data


class TBNFile(LDPFileBase):
    """
    Class to make it easy to interface with a TBN file.  Methods defined for this class are:
     * get_info - Get information about the file's contents
     * get_remaining_frame_count - Get the number of frames remaining in the file
     * offset - Offset a specified number of seconds into the file
     * read_frame - Read and return a single `lsl.reader.tbn.Frame` instance
     * read - Read a chunk of data in and return it as a numpy array
     * estimate_levels - Estimate the n-sigma level for the absolute value of the voltages 
    """
    
    def _ready_file(self):
        """
        Given an open file handle, find the start of valid TBN data.  This
        function:
         1) Aligns on the first valid Mark 5C frame and
         2) Skips over any TBW frames at the beginning of the file.
        """
        
        # Align on the start of a Mark5C packet
        while True:
            try:
                junkFrame = tbn.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(-tbn.FRAME_SIZE+1, 1)
                
        # Jump over any TBN data in the file
        while not junkFrame.header.is_tbn:
            junkFrame = tbn.read_frame(self.fh)
        self.fh.seek(-tbn.FRAME_SIZE, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the TBN file and initialize the frame circular buffer.
        """
        
        try:
            filesize = self.fh.size
        except AttributeError:
            filesize = os.fstat(self.fh.fileno()).st_size
        nFramesFile = (filesize - self.fh.tell()) // tbn.FRAME_SIZE
        framesPerObsX, framesPerObsY = tbn.get_frames_per_obs(self.fh)
        srate =  tbn.get_sample_rate(self.fh, nframe=((framesPerObsX+framesPerObsY)*3))
        bits = 8
        
        with FilePositionSaver(self.fh):
            junkFrame = self.read_frame()
        tuning1 = junkFrame.central_freq
        start = junkFrame.time
        startRaw = junkFrame.payload.timetag
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': tbn.FRAME_SIZE,
                            'nantenna': framesPerObsX+framesPerObsY, 
                            'sample_rate': srate, 'data_bits': bits, 
                            'start_time': start, 'start_time_samples': startRaw, 'freq1': tuning1}
                        
        # Initialize the buffer as part of the description process
        pols = []
        if framesPerObsX != 0:
            pols.append(0)
        if framesPerObsY != 0:
            pols.append(1)
        nAntenna = framesPerObsX + framesPerObsY
        
        self.buffer = TBNFrameBuffer(stands=range(1,nAntenna//len(pols)+1), pols=pols, nsegments=LDP_CONFIG.get('tbn_buffer_size'))
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open TBN file.  This function 
        returns the exact offset time.
        
        .. note::
            The offset provided by this function is relatively crude due to the
            structure of TBN files.
            
        .. versionchanged:: 1.2.4
            Offsets are now relative to the current location in the file rather
            than to the start of the file
        """
        
        # Figure out how far we need to offset inside the file
        junkFrame = tbn.read_frame(self.fh)
        self.fh.seek(-tbn.FRAME_SIZE, 1)
        
        # Get the initial time
        t0 = junkFrame.time
        if getattr(self, "_timetag", None) is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is not None:
                t0 = FrameTimestamp.from_dp_timetag(curr)
                
        # Offset in frames
        ioffset = int(offset * self.description['sample_rate'] / 512 * self.description['nantenna'])
        ioffset = int(1.0 * ioffset / self.description['nantenna']) * self.description['nantenna']
        self.fh.seek(ioffset*tbn.FRAME_SIZE, 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.  
        diffs_used = deque([], 25)
        while True:
            junkFrame = tbn.read_frame(self.fh)
            self.fh.seek(-tbn.FRAME_SIZE, 1)
            
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            t1 = junkFrame.time
            ## See how far off the current frame is from the target
            tDiff = (t1 - t0) - offset
            diffs_used.append(tDiff)
            
            ## Eighth that to come up with a new seek parameter
            tCorr   = -tDiff / 8.0
            cOffset = int(tCorr * self.description['sample_rate'] / 512 * self.description['nantenna'])
            cOffset = int(1.0 * cOffset / self.description['nantenna']) * self.description['nantenna']
            ioffset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset == 0:
                break
            try:
                self.fh.seek(cOffset*tbn.FRAME_SIZE, 1)
                assert(len(set(diffs_used)) > len(diffs_used)//4)
            except (IOError, AssertionError):
                warnings.warn(colorfy("{{%yellow Could not find the correct offset, giving up}}"), RuntimeWarning)

                break
                
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        return t1 - t0
        
    def read_frame(self, return_ci8=False):
        """
        Read and return a single `lsl.reader.tbn.Frame` instance.  If
        `return_ci8` is True then the frame will contain `lsl.reader.base.CI8`
        data instead of numpy.complex64 data.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        tbn_rf = tbn.read_frame_ci8 if return_ci8 else tbn.read_frame
        return tbn_rf(self.fh)
        
    def read(self, duration, time_in_samples=False, return_ci8=False):
        """
        Read in a chunk (in seconds) of TBN data.  This function returns 
        a three-element tuple with elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 2-D array of data (see below).
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.dp.fS` if the 
        `time_in_samples' keyword is set.
        
        If `return_ci8` is True then the data are returned will contain 
        `lsl.reader.base.CI8` data instead of numpy.complex64.  The two
        dimensions are input by samples.
        
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
        timetagSkip = int(512 / self.description['sample_rate'] * fS)
        
        # Setup the read_frame version to use
        tbn_rf = tbn.read_frame_ci8 if return_ci8 else tbn.read_frame
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            junkFrame = tbn_rf(self.fh)
            self._timetag = junkFrame.payload.timetag - timetagSkip
            self.fh.seek(-tbn.FRAME_SIZE, 1)
            
        # Find out how many frames to read in
        frame_count = int(round(1.0 * duration * self.description['sample_rate'] / 512))
        frame_count = frame_count if frame_count else 1
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in range(self.description['nantenna'])]
        if return_ci8:
            data = numpy.zeros((self.description['nantenna'], frame_count*512), dtype=CI8)
        else:
            data = numpy.zeros((self.description['nantenna'], frame_count*512), dtype=numpy.complex64)
        while True:
            if eofFound or nFrameSets == frame_count:
                break
            
            cFrames = deque()
            for i in range(self.description['nantenna']//2):
                try:
                    cFrames.append( tbn_rf(self.fh, verbose=False) )
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
            if cTimetag != self._timetag+timetagSkip:
                actStep = cTimetag - self._timetag
                if self.ignore_timetag_errors:
                    warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i, but found %i}}" % (timetagSkip, actStep)), RuntimeWarning)
                else:
                    raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
            self._timetag = cFrames[0].payload.timetag
            
            for cFrame in cFrames:
                stand,pol = cFrame.header.id
                aStand = 2*(stand-1)+pol
                
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag
                    else:
                        setTime = cFrame.time
                        
                data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                count[aStand] += 1
            nFrameSets += 1
            
        # If we've hit the end of the file and haven't read in enough frames, 
        # flush the buffer
        if eofFound and nFrameSets != frame_count:
            for cFrames in self.buffer.flush():
                cTimetag = cFrames[0].payload.timetag
                if cTimetag != self._timetag+timetagSkip:
                    actStep = cTimetag - self._timetag
                    if self.ignore_timetag_errors:
                        warnings.warn(colorfy("{{%%yellow Invalid timetag skip encountered, expected %i, but found %i}}" % (timetagSkip, actStep)), RuntimeWarning)
                    else:
                        raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
                self._timetag = cFrames[0].payload.timetag
                
                for cFrame in cFrames:
                    stand,pol = cFrame.header.id
                    aStand = 2*(stand-1)+pol
                    
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag
                        else:
                            setTime = cFrame.time
                        
                    data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                    count[aStand] += 1
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets * 512 / self.description['sample_rate']
        
        return duration, setTime, data
        
    def estimate_levels(self, nframe=100, sigma=5.0):
        """
        Estimate the n-sigma level for the absolute value of the voltages.  
        Returns a list with indicies that are the digitizer numbers minus one.
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
                
        # Go!
        with FilePositionSaver(self.fh):
            count = {}
            for i in range(self.description['nantenna']):
                count[i] = 0
            data = numpy.zeros((self.description['nantenna'], nframe*512))
            for i in range(nframe):
                for j in range(self.description['nantenna']):
                    # Read in the next frame and anticipate any problems that could occur
                    try:
                        cFrame = tbn.read_frame(self.fh, verbose=False)
                    except errors.EOFError:
                        break
                    except errors.SyncError:
                        continue
                        
                    s,p = cFrame.id
                    aStand = 2*(s-1) + p
                    
                    try:
                        data[aStand, count[aStand]*512:(count[aStand]+1)*512] = numpy.abs( cFrame.payload.data )
                        count[aStand] +=  1
                    except ValueError:
                        pass
                        
        # Statistics
        rv = norm()
        frac = rv.cdf(sigma) - rv.cdf(-sigma)
        index = int(round(data.shape[1]*frac))
        if index == data.shape[1]:
            index = data.shape[1] - 1
        
        levels = [0 for i in range(self.description['nantenna'])]
        for i in range(self.description['nantenna']):
            data2 = sorted(data[i,:])
            levels[i] = data2[index]
        
        return levels


class DRXFile(LDPFileBase):
    """
    Class to make it easy to interface with a DRX file.  Methods defined for this class are:
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
                assert(len(set(diffs_used)) > len(diffs_used)//4)
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
        tags can be returns as samples at `lsl.common.dp.fS` if the 
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
            data = numpy.zeros((4,frame_count*4096), dtype=CI8)
        else:
            data = numpy.zeros((4,frame_count*4096), dtype=numpy.complex64)
            
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
                            baseframe = copy.deepcopy(cFrames[0])
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
                        raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep))
                        
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag - cFrame.header.time_offset
                    else:
                        setTime = cFrame.time
                        
                data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data
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
                            raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep))
                            
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag - cFrame.header.time_offset
                        else:
                            setTime = cFrame.time
                            
                    data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.payload.data
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
            data = numpy.zeros((4, nframe*4096))
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
                    
                    data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = numpy.abs( cFrame.payload.data )
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
    Class to make it easy to interface with a DR Spectrometer file.  
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
        tags can be returns as samples at `lsl.common.dp.fS` if the 
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
        data = numpy.zeros((2*self.description['nproduct'],frame_count,self.description['LFFT']), dtype=numpy.float32)
        
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
                    raise RuntimeError("Invalid timetag skip encountered, expected %i but found %i" % (timetagSkip, actStep))
                    
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


def LWA1DataFile(filename=None, fh=None, ignore_timetag_errors=False, buffering=-1):
    """
    Wrapper around the various LWA1-related classes defined here that takes
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
    for mode in (drx, tbn, tbw, drspec):
        ## Set if we find a valid frame marker
        foundMatch = False
        ## Set if we can read more than one valid successfully
        foundMode = False
        
        ## Sort out the frame size.  This is tricky because DR spectrometer files
        ## have frames of different sizes depending on the mode
        if mode == drspec:
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
            
    # There is an ambiguity that can arise for TBW data such that it *looks* 
    # like TBN.  If the identified mode is TBN, skip halfway into the file and 
    # verify that it is still TBN.  We also need to catch the LWA-SV DRX vs.
    # TBF ambiguity since we could have been given an LWA-SV file by accident
    if mode in (tbn, drx):
        ## Sort out the frame size
        omfs = mode.FRAME_SIZE
        
        ## Seek half-way in
        if is_splitfile:
            nFrames = fh.size//omfs
        else:
            nFrames = os.path.getsize(filename)//omfs
        fh.seek(nFrames//2*omfs)
        
        ## Read a bit of data to try to find the right type
        for mode in (tbn, tbw, drx):
            ### Set if we find a valid frame marker
            foundMatch = False
            ### Set if we can read more than one valid successfully
            foundMode = False
            
            ### Sort out the frame size.
            mfs = mode.FRAME_SIZE
            
            ### Loop over the frame size to try and find what looks like valid data.  If
            ### is is found, set 'foundMatch' to True.
            for i in range(mfs):
                try:
                    junkFrame = mode.read_frame(fh)
                    foundMatch = True
                    break
                except errors.EOFError:
                    break
                except errors.SyncError:
                    fh.seek(-mfs+1, 1)
                    
            ### Did we strike upon a valid frame?
            if foundMatch:
                #### Is so, we now need to try and read more frames to make sure we have 
                #### the correct type of file
                fh.seek(-mfs, 1)
                
                try:
                    for i in range(4):
                        junkFrame = mode.read_frame(fh)
                    foundMode = True
                except errors.SyncError:
                    #### Reset for the next mode...
                    fh.seek(nFrames//2*omfs)
            else:
                #### Reset for the next mode...
                fh.seek(nFrames//2*omfs)
                
            ### Did we read more than one valid frame?
            if foundMode:
                break
                
    fh.seek(0)
    if not is_splitfile:
        fh.close()
        fh = None
        
    # Raise an error if nothing is found
    if not foundMode:
        raise RuntimeError("File '%s' does not appear to be a valid LWA1 data file" % filename)
        
    # Otherwise, build and return the correct LDPFileBase sub-class
    if mode == drx:
        ldpInstance = DRXFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == tbn:
        ldpInstance = TBNFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == tbw:
        ldpInstance = TBWFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    else:
        ldpInstance = DRSpecFile(filename=filename, fh=fh,
                                 ignore_timetag_errors=ignore_timetag_errors,
                                 buffering=buffering)
        
    # Done
    return ldpInstance


class TBFFile(LDPFileBase):
    """
    Class to make it easy to interface with a TBF file.  Methods defined for this class are:
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
                tbf.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(-tbf.FRAME_SIZE+1, 1)
                
        # Skip over any DRX frames the start of the file
        i = 0
        while True:
            try:
                tbf.read_frame(self.fh)
                break
            except errors.SyncError:
                i += 1
                self.fh.seek(-tbf.FRAME_SIZE+drx.FRAME_SIZE, 1)
        if i == 0:
            self.fh.seek(-tbf.FRAME_SIZE, 1)
        self.fh.seek(-tbf.FRAME_SIZE, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the TBF file.
        """
        
        with FilePositionSaver(self.fh):
            # Read in frame
            junkFrame = tbf.read_frame(self.fh)
            self.fh.seek(-tbf.FRAME_SIZE, 1)
            
            # Basic file information
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
            nFramesFile = (filesize - self.fh.tell()) // tbf.FRAME_SIZE
            srate = fC
            bits = 4
            nFramesPerObs = tbf.get_frames_per_obs(self.fh)
            nchan = tbf.get_channel_count(self.fh)
            firstFrameCount = tbf.get_first_frame_count(self.fh)
            firstChan = tbf.get_first_channel(self.fh)
            
            # Pre-load the channel mapper
            self.mapper = [firstChan+i*tbf.FRAME_CHANNEL_COUNT for i in range(nFramesPerObs)]
            
            # Find the "real" starttime
            while junkFrame.header.frame_count != firstFrameCount:
                junkFrame = tbf.read_frame(self.fh)
            start = junkFrame.time
            startRaw = junkFrame.payload.timetag
            
        # Calculate the frequencies
        freq = numpy.zeros(nchan)
        for i,c in enumerate(self.mapper):
            freq[i*tbf.FRAME_CHANNEL_COUNT:(i+1)*tbf.FRAME_CHANNEL_COUNT] = c + numpy.arange(tbf.FRAME_CHANNEL_COUNT)
        freq *= fC
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': tbf.FRAME_SIZE,
                            'sample_rate': srate, 'data_bits': bits, 
                            'nantenna': 512, 'nchan': nchan, 'freq1': freq, 'start_time': start, 
                            'start_time_samples': startRaw}
                        
        # Initialize the buffer as part of the description process
        self.buffer = TBFFrameBuffer(chans=self.mapper, nsegments=LDP_CONFIG.get('tbf_buffer_size'))
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open TBF file.  This function 
        returns the exact offset time.
        
        .. note::
            The offset provided by this function is relatively crude due to the
            structure of TBF files.
            
        .. versionchanged:: 1.2.4
            Offsets are now relative to the current location in the file rather
            than to the start of the file
        """
        
        # Find out where we really are taking into account the buffering
        buffer_offset = 0
        if getattr(self, "_timetag", None) is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is None:
                frame = tbf.read_frame(self.fh)
                self.fh.seek(-tbf.FRAME_SIZE, 1)
                curr = frame.payload.time_tag
            buffer_offset = curr - self._timetag
            buffer_offset = buffer_offset / fS
            
        offset = offset - buffer_offset
        framesPerObs = self.description['nchan'] // tbf.FRAME_CHANNEL_COUNT
        frameOffset = int(offset * self.description['sample_rate'] * framesPerObs)
        frameOffset = int(1.0 * frameOffset / framesPerObs) * framesPerObs
        self.fh.seek(frameOffset*tbf.FRAME_SIZE, 1)
        
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
        Read and return a single `lsl.reader.tbf.Frame` instance.  If
        `return_ci8` is True then the frame will contain `lsl.reader.base.CI8`
        data instead of numpy.complex64 data.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        tbf_rf = tbf.read_frame_ci8 if return_ci8 else tbf.read_frame
        return tbf_rf(self.fh)
        
    def read(self, duration=None, time_in_samples=False, return_ci8=False):
        """
        Read and return the entire TBF capture.  This function returns 
        a three-element tuple with elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 3-D Numpy array of data.
        
        The time tag is returned as seconds since the UNIX epoch as a 
        `lsl.reader.base.FrameTimestamp` instance by default.  However, the time 
        tags can be returns as samples at `lsl.common.dp.fS` if the 
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
        tbf_rf = tbf.read_frame_ci8 if return_ci8 else tbf.read_frame
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            
        # Find out how many frames to read in
        if duration is None:
            duration = self.description['nframe'] / framesPerObs / self.description['sample_rate']
        framesPerObs = self.description['nchan'] // tbf.FRAME_CHANNEL_COUNT
        frame_count = int(round(1.0 * duration * self.description['sample_rate']))
        frame_count = frame_count if frame_count else 1
        duration = frame_count / self.description['sample_rate']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in range(framesPerObs)]
        if return_ci8:
            data = numpy.zeros((self.description['nantenna'], self.description['nchan'], frame_count), dtype=CI8)
        else:
            data = numpy.zeros((self.description['nantenna'], self.description['nchan'], frame_count), dtype=numpy.complex64)
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            cFrames = deque()
            for i in range(framesPerObs):
                try:
                    cFrame = tbf_rf(self.fh, verbose=False)
                    if not cFrame.is_tbf:
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
                    raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
            self._timetag = cFrames[0].payload.timetag
            
            for cFrame in cFrames:
                first_chan = cFrame.header.first_chan
                
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag
                    else:
                        setTime = cFrame.time
                        
                subData = cFrame.payload.data
                subData.shape = (tbf.FRAME_CHANNEL_COUNT,512)
                subData = subData.T
                
                aStand = self.mapper.index(first_chan)
                data[:,aStand*tbf.FRAME_CHANNEL_COUNT:(aStand+1)*tbf.FRAME_CHANNEL_COUNT,count[aStand]] = subData
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
                        raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
                self._timetag = cFrames[0].payload.timetag
                
                for cFrame in cFrames:
                    first_chan = cFrame.header.first_chan
                    
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag
                        else:
                            setTime = cFrame.time
                        
                    subData = cFrame.payload.data
                    subData.shape = (tbf.FRAME_CHANNEL_COUNT,512)
                    subData = subData.T
                    
                    aStand = self.mapper.index(first_chan)
                    data[:,aStand*tbf.FRAME_CHANNEL_COUNT:(aStand+1)*tbf.FRAME_CHANNEL_COUNT,count[aStand]] = subData
                    count[aStand] += 1
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Sanity check at the end to see if we actually read anything.  
        # This is needed because of how TBF and DRX interact where TBF
        # files can be padded at the end with DRX data
        if nFrameSets == 0 and duration > 0:
            raise errors.EOFError()
            
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets  / self.description['sample_rate']
        
        return duration, setTime, data


class CORFile(LDPFileBase):
    """
    Class to make it easy to interface with a COR file.  Methods defined for this class are:
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
                self.fh.seek(-cor.FRAME_SIZE+1, 1)
                
        # Skip over any DRX frames the start of the file
        i = 0
        while True:
            try:
                cor.read_frame(self.fh)
                break
            except errors.SyncError:
                i += 1
                self.fh.seek(-cor.FRAME_SIZE+drx.FRAME_SIZE, 1)
        if i == 0:
            self.fh.seek(-cor.FRAME_SIZE, 1)
        self.fh.seek(-cor.FRAME_SIZE, 1)
        
        return True
        
    def _describe_file(self):
        """
        Describe the COR file.
        """
        
        # Read in frame
        with FilePositionSaver(self.fh):
            junkFrame = cor.read_frame(self.fh)
            self.fh.seek(-cor.FRAME_SIZE, 1)
            
            # Basic file information
            try:
                filesize = os.fstat(self.fh.fileno()).st_size
            except AttributeError:
                filesize = self.fh.size
            nFramesFile = (filesize - self.fh.tell()) // cor.FRAME_SIZE
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
            self.cmapper = []
            marker = self.fh.tell()
            firstFrameCount = 2**64-1
            while len(self.cmapper) < nchan/cor.FRAME_CHANNEL_COUNT:
                cFrame = cor.read_frame(self.fh)
                if cFrame.header.first_chan not in self.cmapper:
                    self.cmapper.append( cFrame.header.first_chan )
                if cFrame.header.frame_count < firstFrameCount:
                    firstFrameCount = cFrame.header.frame_count
                    start = junkFrame.time
                    startRaw = junkFrame.payload.timetag
            self.cmapper.sort()
            
        # Create a channel mapper dictionary
        self.cmapperd = {}
        for i,c in enumerate(self.cmapper):
            self.cmapperd[c] = i
            
        # Calculate the frequencies
        freq = numpy.zeros(nchan)
        for i,c in enumerate(self.cmapper):
            freq[i*cor.FRAME_CHANNEL_COUNT:(i+1)*cor.FRAME_CHANNEL_COUNT] = c + numpy.arange(cor.FRAME_CHANNEL_COUNT)
        freq *= fC
        
        self.description = {'size': filesize, 'nframe': nFramesFile, 'frame_size': cor.FRAME_SIZE,
                            'sample_rate': srate, 'data_bits': bits, 
                            'nantenna': 512, 'nchan': nchan, 'freq1': freq, 'start_time': start, 
                            'start_time_samples': startRaw, 'nbaseline': nBaseline, 'tint':cFrame.integration_time}
                        
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
        framesPerObs = self.description['nchan'] // cor.FRAME_CHANNEL_COUNT * self.description['nbaseline']
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
        tags can be returns as samples at `lsl.common.dp.fS` if the 
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
        framesPerObs = self.description['nchan'] // cor.FRAME_CHANNEL_COUNT * self.description['nbaseline']
        frame_count = int(round(1.0 * duration / self.description['tint']))
        frame_count = frame_count if frame_count else 1
        duration = frame_count * self.description['tint']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in range(framesPerObs)]
        data = numpy.zeros((self.description['nbaseline'], self.description['nchan'], 2, 2, frame_count), dtype=numpy.complex64)
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
                    raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
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
                data[aBase,aChan*cor.FRAME_CHANNEL_COUNT:(aChan+1)*cor.FRAME_CHANNEL_COUNT,:,:,count[aStand]] = cFrame.payload.data
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
                        raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
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
                    data[aBase,aChan*cor.FRAME_CHANNEL_COUNT:(aChan+1)*cor.FRAME_CHANNEL_COUNT,:,:,count[aStand]] = cFrame.payload.data
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


def LWASVDataFile(filename=None, fh=None, ignore_timetag_errors=False, buffering=-1):
    """
    Wrapper around the various LWA-SV-related classes defined here that takes
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
    for mode in (drx, tbn, tbf, cor, drspec):
        ## Set if we find a valid frame marker
        foundMatch = False
        ## Set if we can read more than one valid successfully
        foundMode = False
        
        ## Sort out the frame size.  This is tricky because DR spectrometer files
        ## have frames of different sizes depending on the mode
        if mode == drspec:
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
            
    # There is an ambiguity that can arise for TBF data such that it *looks* 
    # like DRX.  If the identified mode is DRX, skip halfway into the file and 
    # verify that it is still DRX.   We also need to catch the LWA1 TBN vs.
    # TBW ambiguity since we could have been given an LWA1 file by accident.
    if mode in (drx, tbn):
        ## Sort out the frame size
        omfs = mode.FRAME_SIZE
        
        ## Seek half-way in
        if is_splitfile:
            nFrames = fh.size//omfs
        else:
            nFrames = os.path.getsize(filename)//omfs
        fh.seek(nFrames//2*omfs)
        
        ## Read a bit of data to try to find the right type
        for mode in (tbn, drx, tbf):
            ### Set if we find a valid frame marker
            foundMatch = False
            ### Set if we can read more than one valid successfully
            foundMode = False
            
            ### Sort out the frame size.
            mfs = mode.FRAME_SIZE
            
            ### Loop over the frame size to try and find what looks like valid data.  If
            ### is is found, set 'foundMatch' to True.
            for i in range(mfs):
                try:
                    junkFrame = mode.read_frame(fh)
                    foundMatch = True
                    break
                except errors.EOFError:
                    break
                except errors.SyncError:
                    fh.seek(-mfs+1, 1)
                    
            ### Did we strike upon a valid frame?
            if foundMatch:
                #### Is so, we now need to try and read more frames to make sure we have 
                #### the correct type of file
                fh.seek(-mfs, 1)
                
                try:
                    for i in range(4):
                        junkFrame = mode.read_frame(fh)
                    foundMode = True
                except errors.SyncError:
                    #### Reset for the next mode...
                    fh.seek(nFrames//2*omfs)
            else:
                #### Reset for the next mode...
                fh.seek(nFrames//2*omfs)
                
            ### Did we read more than one valid frame?
            if foundMode:
                break
                
    fh.seek(0)
    if not is_splitfile:
        fh.close()
        fh = None
    
    # Raise an error if nothing is found
    if not foundMode:
        raise RuntimeError("File '%s' does not appear to be a valid LWA-SV data file" % filename)
        
    # Otherwise, build and return the correct LDPFileBase sub-class
    if mode == drx:
        ldpInstance = DRXFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == tbn:
        ldpInstance = TBNFile(filename=filename, fh=fh,
                              ignore_timetag_errors=ignore_timetag_errors,
                              buffering=buffering)
    elif mode == tbf:
        ldpInstance = TBFFile(filename=filename, fh=fh,
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


def LWADataFile(filename=None, fh=None, ignore_timetag_errors=False, buffering=-1):
    """
    Wrapper around the various classes defined here that takes a file, 
    determines the data type, and initializes and returns the appropriate
    LDP class.
    """
    
    found = False
    
    # LWA-1?
    if not found:
        try:
            ldpInstance = LWA1DataFile(filename=filename, fh=fh,
                                       ignore_timetag_errors=ignore_timetag_errors,
                                       buffering=buffering)
            found = True
        except RuntimeError:
            pass
            
    # LWA-SV?
    if not found:
        try:
            ldpInstance = LWASVDataFile(filename=filename, fh=fh,
                                       ignore_timetag_errors=ignore_timetag_errors,
                                       buffering=buffering)
            found = True
        except RuntimeError:
            pass
            
    # Failed?
    if not found:
        raise RuntimeError("File '%s' does not appear to be a valid LWA1 or LWA-SV data file" % filename)
        
    return ldpInstance


import atexit
atexit.register(_open_ldp_files.close_all)
