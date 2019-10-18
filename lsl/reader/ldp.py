# -*- coding: utf-8 -*-

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

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
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

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.4'
__revision__ = '$Rev$'
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
    
    def __init__(self, filename=None, fh=None, ignore_timetag_errors=False):
        # Make sure that we are given either a filename or an open file handle
        if filename is None and fh is None:
            raise RuntimeError("Must specify either a filename or open file instance")
            
        # Store a valid file handle and mark the object as ready
        if fh is None:
            self.filename = filename
            self.fh = open(filename, 'rb')
        else:
            self.filename = fh.name
            if not isinstance(fh, SplitFileWrapper):
                if fh.mode.find('b') == -1:
                    fh.close()
                    fh = open(self.filename, 'rb')
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
        return "%s @ %s" % (self.__name__, self.filename)
        
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
        
    def _ready_file(self):
        """
        Method for finding the start of valid data.  This will be over-
        ridden in the format-specific subclasses.
        """
        
        raise NotImplementedError
        
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
    def nframes_remaining(self):
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
        
    def read_frame(self):
        """
        Read a single frame from the data.
        """
        
        raise NotImplementedError
        
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
                
    def estimate_levels(self, nframes=10, sigma=5.0):
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
                        junkFrame = tbn.read_frame(self.fh)
                        break
                    except errors.SyncError:
                        self.fh.seek(-tbn.FRAME_SIZE+1, 1)
                        
                ## Find the end of the TBN data
                while True:
                    try:
                        junkFrame = tbn.read_frame(self.fh)
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
            start = junkFrame.time
            startRaw = junkFrame.payload.timetag
            
            # Trick to figure out how many antennas are in a file and the "real" 
            # start time.  For details of why this needs to be done, see the read()
            # function below.
            idsFound = []
            timesFound = []
            filePosRef = self.fh.tell()
            while True:
                try:
                    for i in xrange(26):
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
        
        self.description = {'size': filesize, 'nframes': nFramesFile, 'frame_size': tbw.FRAME_SIZE,
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
        
        The time tag is returned as seconds since the UNIX epoch by default.
        However, the time tags can be returns as samples at fS if the 
        time_in_samples keyword is set.
        
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
        srate =  tbn.get_sample_rate(self.fh, nframes=((framesPerObsX+framesPerObsY)*3))
        bits = 8
        
        with FilePositionSaver(self.fh):
            junkFrame = self.read_frame()
        tuning1 = junkFrame.central_freq
        start = junkFrame.time
        startRaw = junkFrame.payload.timetag
        
        self.description = {'size': filesize, 'nframes': nFramesFile, 'frame_size': tbn.FRAME_SIZE,
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
        
        self.buffer = TBNFrameBuffer(stands=range(1,nAntenna//len(pols)+1), pols=pols)
        
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
        
        # Find out where we really are taking into account the buffering
        buffer_offset = 0
        if self._timetag is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is None:
                frame = tbn.read_frame(self.fh)
                self.fh.seek(-tbn.FRAME_SIZE, 1)
                curr = frame.payload.time_time
            buffer_offset = curr - self._timetag
            buffer_offset = buffer_offset / fS
            
        offset = offset - buffer_offset
        frameOffset = int(offset * self.description['sample_rate'] / 512 * self.description['nantenna'])
        frameOffset = int(1.0 * frameOffset / self.description['nantenna']) * self.description['nantenna']
        self.fh.seek(frameOffset*tbn.FRAME_SIZE, 1)
        
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        return 1.0 * frameOffset / self.description['nantenna'] * 512 / self.description['sample_rate']
        
    def read_frame(self):
        """
        Read and return a single `lsl.reader.tbn.Frame` instance.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        return tbn.read_frame(self.fh)
        
    def read(self, duration, time_in_samples=False):
        """
        Read in a chunk (in seconds) of TBN data.  This function returns 
        a three-element tuple with elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 2-D Numpy array of data.
        
        The time tag is returned as seconds since the UNIX epoch by default.
        However, the time tags can be returns as samples at fS if the 
        time_in_samples keyword is set.
        
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
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            junkFrame = tbn.read_frame(self.fh)
            self._timetag = junkFrame.payload.timetag - timetagSkip
            self.fh.seek(-tbn.FRAME_SIZE, 1)
            
        # Find out how many frames to read in
        frame_count = int(round(1.0 * duration * self.description['sample_rate'] / 512))
        frame_count = frame_count if frame_count else 1
        duration = frame_count * 512 / self.description['sample_rate']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in xrange(self.description['nantenna'])]
        data = numpy.zeros((self.description['nantenna'], frame_count*512), dtype=numpy.complex64)
        while True:
            if eofFound or nFrameSets == frame_count:
                break
            
            cFrames = deque()
            for i in xrange(self.description['nantenna']//2):
                try:
                    cFrames.append( tbn.read_frame(self.fh, verbose=False) )
                except errors.EOFError:
                    eofFound = True
                    self.buffer.append(cFrames)
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
                    warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
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
                        setTime = sum(cFrame.time)
                        
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
                        warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
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
                            setTime = sum(cFrame.time)
                        
                    data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                    count[aStand] += 1
                nFrameSets += 1
                
                if nFrameSets == frame_count:
                    break
                    
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets * 512 / self.description['sample_rate']
        
        return duration, setTime, data
        
    def estimate_levels(self, nframes=100, sigma=5.0):
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
            for i in xrange(self.description['nantenna']):
                count[i] = 0
            data = numpy.zeros((self.description['nantenna'], nframes*512))
            for i in xrange(nframes):
                for j in xrange(self.description['nantenna']):
                    # Read in the next frame and anticipate any problems that could occur
                    try:
                        cFrame = tbn.read_frame(self.fh, verbose=False)
                    except errors.EOFError:
                        break
                    except errors.SyncError:
                        continue
                        
                    s,p = cFrame.id
                    aStand = 2*(s-1) + p
                    
                    data[aStand, count[aStand]*512:(count[aStand]+1)*512] = numpy.abs( cFrame.payload.data )
                    count[aStand] +=  1
                    
        # Statistics
        rv = norm()
        frac = rv.cdf(sigma) - rv.cdf(-sigma)
        index = int(round(data.shape[1]*frac))
        if index == data.shape[1]:
            index = data.shape[1] - 1
        
        levels = [0 for i in xrange(self.description['nantenna'])]
        for i in xrange(self.description['nantenna']):
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
                    break
                except ZeroDivisionError:
                    pass
            except errors.SyncError:
                self.fh.seek(-drx.FRAME_SIZE+1, 1)
                
        self.fh.seek(-drx.FRAME_SIZE, 1)
        
        # Line up the time tags for the various tunings/polarizations
        ids = []
        timetags = []
        for i in xrange(32):
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
        beams = drx.get_beam_count(self.fh)
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
            for i in xrange(32):
                junkFrame = self.read_frame()
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
                    start = junkFrame.time
                    startRaw = junkFrame.payload.timetag - junkFrame.header.time_offset
                    
        self.description = {'size': filesize, 'nframes': nFramesFile, 'frame_size': drx.FRAME_SIZE,
                            'beampols': beampols, 'beam': b, 
                            'sample_rate': srate, 'data_bits': bits, 
                            'start_time': start, 'start_time_samples': startRaw, 'freq1': tuning1, 'freq2': tuning2}
                        
        # Initialize the buffer as part of the description process
        self.buffer = DRXFrameBuffer(beams=beams, tunes=tunes, pols=pols)
        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open DRX file.  This function 
        returns the exact offset time.
        """
        
        # Figure out how far we need to offset inside the file
        junkFrame = drx.read_frame(self.fh)
        self.fh.seek(-drx.FRAME_SIZE, 1)
        
        # Get the initial time, sample rate, and beampols
        ti0, tf0 = junkFrame.time
        if self._timetag is not None:
            curr = self.buffer.peek(require_filled=False)
            if curr is not None:
                ti0 = (curr - junkFrame.header.timeOffset) // int(fS)
                tf0 = ((curr - junkFrame.header.timeOffset) % int(fS)) / fS
        sample_rate = junkFrame.sample_rate
        beampols = drx.get_frames_per_obs(self.fh)
        beampols = reduce(int.__add__, beampols)
        
        # Offset in frames for beampols beam/tuning/pol. sets
        ioffset = int(offset * sample_rate / 4096 * beampols)
        ioffset = int(1.0 * ioffset / beampols) * beampols
        self.fh.seek(ioffset*drx.FRAME_SIZE, 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.  
        while True:
            junkFrame = drx.read_frame(self.fh)
            self.fh.seek(-drx.FRAME_SIZE, 1)
            
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            ti1, tf1 = junkFrame.time
            sample_rate = junkFrame.sample_rate
            beampols = drx.get_frames_per_obs(self.fh)
            beampols = reduce(int.__add__, beampols)
            
            ## See how far off the current frame is from the target
            tDiff = ti1 - (ti0 + offset) + tf1 - tf0
            
            ## Half that to come up with a new seek parameter
            tCorr   = -tDiff / 2.0
            cOffset = int(tCorr * sample_rate / 4096 * beampols)
            cOffset = int(1.0 * cOffset / beampols) * beampols
            ioffset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset is 0:
                break
            try:
                self.fh.seek(cOffset*drx.FRAME_SIZE, 1)
            except IOError:
                warnings.warn("Could not find the correct offset, giving up", RuntimeWarning)
                break
                
        # Update the file metadata
        self._describe_file()
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Zero out the time tag checker
        self._timetag = None
        
        return t1 - t0
        
    def read_frame(self):
        """
        Read and return a single `lsl.reader.drx.Frame` instance.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Zero out the time tag checker
        self._timetagSkip = None
        self._timetag = None
        
        return drx.read_frame(self.fh)
        
    def read(self, duration, time_in_samples=False):
        """
        Given an open DRX file and an amount of data to read in in seconds, read 
        in the data and return a three-element tuple of the actual duration read 
        in, the time for the first sample, and the data as numpy 
        array.
        
        ..note::
            This function always returns a 2-D array with the first dimension
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
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = {0:0, 1:0, 2:0, 3:0}
            
        # Find out how many frames to read in
        frame_count = int(round(1.0 * duration * self.description['sample_rate'] / 4096))
        frame_count = frame_count if frame_count else 1
        duration = frame_count * 4096 / self.description['sample_rate']
        
        # Setup the output arrays
        setTime = None
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
                for i in xrange(self.description['beampols']):
                    try:
                        cFrames.append( drx.read_frame(self.fh, verbose=False) )
                    except errors.EOFError:
                        eofFound = True
                        self.buffer.append(cFrames)
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
                    if int(missing) == missing and missing < 50:
                        ## This is kind of black magic down here
                        for m in xrange(int(missing)):
                            m = self._timetag[aStand] + self._timetagSkip*(m+1)
                            baseframe = copy.deepcopy(cFrames[0])
                            baseframe.payload.timeTag = m
                            baseframe.payload.data *= 0
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
                        warnings.warn("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep), RuntimeWarning)
                    else:
                        raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep))
                        
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag - cFrame.header.time_offset
                    else:
                        setTime = sum(cFrame.time)
                        
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
                            warnings.warn("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep), RuntimeWarning)
                        else:
                            raise RuntimeError("Invalid timetag skip encountered, expected %i on tuning %i, pol %i, but found %i" % (self._timetagSkip, t, p, actStep))
                            
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag - cFrame.header.time_offset
                        else:
                            setTime = sum(cFrame.time)
                            
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
        
    def estimate_levels(self, nframes=100, sigma=5.0):
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
            data = numpy.zeros((4, nframes*4096))
            for i in xrange(nframes):
                for j in xrange(self.description['beampols']):
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
        for i in xrange(4):
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
            except errors.syncError:
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
        
        self.description = {'size': filesize, 'nframes': nFramesFile, 'frame_size': FRAME_SIZE, 
                            'beampols': beampols, 'beam': beam, 
                            'sample_rate': srate, 'data_bits': bits, 
                            'start_time': start, 'start_time_samples': startRaw, 'freq1': tuning1, 'freq2': tuning2, 
                            'nInt': nInt, 'tint': tInt, 'LFFT': LFFT, 
                            'nproducts': len(prod), 'data_products': prod}
                        
    def offset(self, offset):
        """
        Offset a specified number of seconds in an open DR spectrometer file.  This 
        function returns the exact offset time.
        """
        
        # Gather some basic information and read in the first frame
        junkFrame = drspec.read_frame(self.fh)
        self.fh.seek(-self.description['frame_size'], 1)
        
        # Get the initial time, sample rate, and integration time
        ti0, tf0 = junkFrame.time
        
        # Offset in frames for beampols beam/tuning/pol. sets
        ioffset = int(round(offset / self.description['tint']))
        self.fh.seek(ioffset*self.description['frame_size'], 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.  
        while True:
            junkFrame = drspec.read_frame(self.fh)
            self.fh.seek(-self.description['frame_size'], 1)
            
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            ti1, tf1 = junkFrame.time
            sample_rate = junkFrame.sample_rate
            LFFT = junkFrame.get_transform_size()
            tInt = junkFrame.header.nints*LFFT/sample_rate
            
            ## See how far off the current frame is from the target
            tDiff = ti1 - (ti0 + offset) + tf1 - tf0
            
            ## Half that to come up with a new seek parameter
            tCorr   = -tDiff / 2.0
            cOffset = int(round(tCorr / tInt))
            ioffset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset is 0:
                break
            self.fh.seek(cOffset*self.description['frame_size'], 1)
            
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
        duration = frame_count * self.description['tint']
        
        # Setup the output arrays
        data = numpy.zeros((2*self.description['nproducts'],frame_count,self.description['LFFT']), dtype=numpy.float32)
        
        # Go!
        nFrameSets = 0
        setTime = None
        for i in xrange(frame_count):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = drspec.read_frame(self.fh, verbose=False)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue
                
            cTimetag = sum(cFrame.time)
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
                    setTime = sum(cFrame.time)
                    
            for j,p in enumerate(self.description['data_products']):
                data[j+0,                             count, :] = getattr(cFrame.payload, '%s0' % p, None)
                data[j+self.description['nproducts'], count, :] = getattr(cFrame.payload, '%s1' % p, None)
            count +=  1
            self._timetag = cTimetag
            nFrameSets += 1
            
        # Adjust the duration to account for all of the things that could 
        # have gone wrong while reading the data
        duration = nFrameSets * self.description['tint']
        
        return duration, setTime, data


def LWA1DataFile(filename=None, fh=None, ignore_timetag_errors=False):
    """
    Wrapper around the various LWA1-related classes defined here that takes
    a file, determines the data type, and initializes and returns the 
    appropriate LDP class.
    """
    
    # Open the file as appropriate
    if fh is None:
        fh = open(filename, 'rb')
    else:
        filename = fh.name
        if fh.mode.find('b') == -1:
            fh.close()
            fh = open(filename, 'rb')
            
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
        for i in xrange(mfs):
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
                for i in xrange(2):
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
            for i in xrange(mfs):
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
                    for i in xrange(4):
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
                
    fh.close()
    
    # Raise an error if nothing is found
    if not foundMode:
        raise RuntimeError("File '%s' does not appear to be a valid LWA1 data file" % filename)
        
    # Otherwise, build and return the correct LDPFileBase sub-class
    if mode == drx:
        ldpInstance = DRXFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    elif mode == tbn:
        ldpInstance = TBNFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    elif mode == tbw:
        ldpInstance = TBWFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    else:
        ldpInstance = DRSpecFile(filename, ignore_timetag_errors=ignore_timetag_errors)
        
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
                junkFrame = tbf.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(-tbf.FRAME_SIZE+1, 1)
                
        # Skip over any DRX frames the start of the file
        i = 0
        while True:
            try:
                junkFrame = tbf.read_frame(self.fh)
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
            
            # Pre-load the channel mapper
            self.mapper = []
            firstFrameCount = 2**64-1
            while len(self.mapper) < nchan/tbf.FRAME_CHANNEL_COUNT:
                cFrame = tbf.read_frame(self.fh)
                if cFrame.header.first_chan not in self.mapper:
                    self.mapper.append( cFrame.header.first_chan )
                if cFrame.header.frame_count < firstFrameCount:
                    firstFrameCount = cFrame.header.frame_count
                    start = junkFrame.time
                    startRaw = junkFrame.payload.timetag
            self.mapper.sort()
            
        # Calculate the frequencies
        freq = numpy.zeros(nchan)
        for i,c in enumerate(self.mapper):
            freq[i*tbf.FRAME_CHANNEL_COUNT:(i+1)*tbf.FRAME_CHANNEL_COUNT] = c + numpy.arange(tbf.FRAME_CHANNEL_COUNT)
        freq *= fC
        
        self.description = {'size': filesize, 'nframes': nFramesFile, 'frame_size': tbf.FRAME_SIZE,
                            'sample_rate': srate, 'data_bits': bits, 
                            'nantenna': 512, 'nchan': nchan, 'freq1': freq, 'start_time': start, 
                            'start_time_samples': startRaw}
                        
        # Initialize the buffer as part of the description process
        self.buffer = TBFFrameBuffer(chans=self.mapper)
        
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
        if self._timetag is not None:
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
        
    def read_frame(self):
        """
        Read and return a single `lsl.reader.tbf.Frame` instance.
        """
        
        # Reset the buffer
        if hasattr(self, "buffer"):
            self.buffer.reset()
            
        # Reset the timetag checker
        self._timetag = None
        
        return tbf.read_frame(self.fh)
        
    def read(self, duration=None, time_in_samples=False):
        """
        Read and return the entire TBF capture.  This function returns 
        a three-element tuple with elements of:
         0) the actual duration of data read in, 
         1) the time tag for the first sample, and
         2) a 3-D Numpy array of data.
        
        The time tag is returned as seconds since the UNIX epoch by default.
        However, the time tags can be returns as samples at fS if the 
        time_in_samples keyword is set.
        
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
        
        # Setup the counter variables:  frame count and time tag count
        if getattr(self, "_timetag", None) is None:
            self._timetag = 0
            
        # Find out how many frames to read in
        if duration is None:
            duration = self.description['nframes'] / framesPerObs / self.description['sample_rate']
        framesPerObs = self.description['nchan'] // tbf.FRAME_CHANNEL_COUNT
        frame_count = int(round(1.0 * duration * self.description['sample_rate']))
        frame_count = frame_count if frame_count else 1
        duration = frame_count / self.description['sample_rate']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in xrange(framesPerObs)]
        data = numpy.zeros((self.description['nantenna'], self.description['nchan'], frame_count), dtype=numpy.complex64)
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            cFrames = deque()
            for i in xrange(framesPerObs):
                try:
                    cFrame = tbf.read_frame(self.fh, verbose=False)
                    if not cFrame.is_tbf:
                        continue
                    cFrames.append( cFrame )
                except errors.EOFError:
                    eofFound = True
                    self.buffer.append(cFrames)
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
                    warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
                else:
                    raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
            self._timetag = cFrames[0].payload.timetag
            
            for cFrame in cFrames:
                first_chan = cFrame.header.first_chan
                
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag
                    else:
                        setTime = sum(cFrame.time)
                        
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
                        warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
                    else:
                        raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
                self._timetag = cFrames[0].payload.timetag
                
                for cFrame in cFrames:
                    first_chan = cFrame.header.first_chan
                    
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag
                        else:
                            setTime = sum(cFrame.time)
                        
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
                junkFrame = cor.read_frame(self.fh)
                break
            except errors.SyncError:
                self.fh.seek(-cor.FRAME_SIZE+1, 1)
                
        # Skip over any DRX frames the start of the file
        i = 0
        while True:
            try:
                junkFrame = cor.read_frame(self.fh)
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
            for i in xrange(1, 256+1):
                for j in xrange(i, 256+1):
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
        
        self.description = {'size': filesize, 'nframes': nFramesFile, 'frame_size': cor.FRAME_SIZE,
                            'sample_rate': srate, 'data_bits': bits, 
                            'nantenna': 512, 'nchan': nchan, 'freq1': freq, 'start_time': start, 
                            'start_time_samples': startRaw, 'nbaseline': nBaseline, 'tint':cFrame.get_integration_time()}
                        
        # Initialize the buffer as part of the description process
        self.buffer = CORFrameBuffer(chans=self.cmapper, reorder=False)
        
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
        if self._timetag is not None:
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
        
        The time tag is returned as seconds since the UNIX epoch by default.
        However, the time tags can be returns as samples at fS if the 
        time_in_samples keyword is set.
        
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
            duration = self.description['nframes'] / framesPerObs * self.description['tint']
        framesPerObs = self.description['nchan'] // cor.FRAME_CHANNEL_COUNT * self.description['nbaseline']
        frame_count = int(round(1.0 * duration / self.description['tint']))
        frame_count = frame_count if frame_count else 1
        duration = frame_count * self.description['tint']
        
        nFrameSets = 0
        eofFound = False
        setTime = None
        count = [0 for i in xrange(framesPerObs)]
        data = numpy.zeros((self.description['nbaseline'], self.description['nchan'], 2, 2, frame_count), dtype=numpy.complex64)
        while True:
            if eofFound or nFrameSets == frame_count:
                break
                
            cFrames = deque()
            for i in xrange(framesPerObs):
                try:
                    cFrame = cor.read_frame(self.fh, verbose=False)
                    if not cFrame.isCOR():
                        continue
                    cFrames.append( cFrame )
                except errors.EOFError:
                    eofFound = True
                    self.buffer.append(cFrames)
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
                    warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
                else:
                    raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
            self._timetag = cFrames[0].payload.timetag
            
            for cFrame in cFrames:
                first_chan = cFrame.header.first_chan
                
                if setTime is None:
                    if time_in_samples:
                        setTime = cFrame.payload.timetag
                    else:
                        setTime = sum(cFrame.time)
                        
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
                        warnings.warn("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep), RuntimeWarning)
                    else:
                        raise RuntimeError("Invalid timetag skip encountered, expected %i, but found %i" % (timetagSkip, actStep))
                self._timetag = cFrames[0].payload.timetag
                
                for cFrame in cFrames:
                    first_chan = cFrame.header.first_chan
                    
                    if setTime is None:
                        if time_in_samples:
                            setTime = cFrame.payload.timetag
                        else:
                            setTime = sum(cFrame.time)
                            
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


def LWASVDataFile(filename=None, fh=None, ignore_timetag_errors=False):
    """
    Wrapper around the various LWA-SV-related classes defined here that takes
    a file, determines the data type, and initializes and returns the 
    appropriate LDP class.
    """
    
    # Open the file as appropriate
    if fh is None:
        fh = open(filename, 'rb')
    else:
        filename = fh.name
        if fh.mode.find('b') == -1:
            fh.close()
            fh = open(filename, 'rb')
            
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
        for i in xrange(mfs):
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
                for i in xrange(2):
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
            for i in xrange(mfs):
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
                    for i in xrange(4):
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
                
    fh.close()
    
    # Raise an error if nothing is found
    if not foundMode:
        raise RuntimeError("File '%s' does not appear to be a valid LWA-SV data file" % filename)
        
    # Otherwise, build and return the correct LDPFileBase sub-class
    if mode == drx:
        ldpInstance = DRXFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    elif mode == tbn:
        ldpInstance = TBNFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    elif mode == tbf:
        ldpInstance = TBFFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    elif mode == cor:
        ldpInstance = CORFile(filename, ignore_timetag_errors=ignore_timetag_errors)
    else:
        ldpInstance = DRSpecFile(filename, ignore_timetag_errors=ignore_timetag_errors)
        
    # Done
    return ldpInstance


def LWADataFile(filename=None, fh=None, ignore_timetag_errors=False):
    """
    Wrapper around the various classes defined here that takes a file, 
    determines the data type, and initializes and returns the appropriate
    LDP class.
    """
    
    found = False
    
    # LWA-1?
    if not found:
        try:
            ldpInstance = LWA1DataFile(filename=filename, fh=fh, ignore_timetag_errors=ignore_timetag_errors)
            found = True
        except RuntimeError:
            pass
            
    # LWA-SV?
    if not found:
        try:
            ldpInstance = LWASVDataFile(filename=filename, fh=fh, ignore_timetag_errors=ignore_timetag_errors)
            found = True
        except RuntimeError:
            pass
            
    # Failed?
    if not found:
        raise RuntimeError("File '%s' does not appear to be a valid LWA1 or LWA-SV data file" % filename)
        
    return ldpInstance


import atexit
atexit.register(_open_ldp_files.close_all)
