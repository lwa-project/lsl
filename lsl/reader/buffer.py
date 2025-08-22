"""
Buffer for dealing with out-of-order/missing frames.

.. versionchanged:: 0.5
    Removed support for DRX FrameBuffers since they shouldn't be needed.
    
.. versionchanged:: 0.6
    Removed support for TBW FrameBuffers since they didn't really work.
    
.. versionchanged:: 1.0.0
    Added back in the DRX FrameBuffers since they might be useful.
    Dropped TBW support from within the FrameBufferBase class.
    
.. versionchanged:: 1.0.2
    Added support for the LWA-SV ADP DRX8 mode
    
.. versionchanged:: 1.1.4
    Dropped support for the LWA-SV ADP DRX8 mode
    
.. versionchanged:: 1.2.1
    Added support for the LWA-SV TBF and COR modes
    
.. versionchanged:: 3.0.2
    Added support for the DRX8 mode (again).
"""

import abc
import copy
from functools import cmp_to_key
from collections import deque, OrderedDict

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '1.3'
__all__ = ['FrameBufferBase', 'TBNFrameBuffer', 'DRXFrameBuffer',
           'DRX8FrameBuffer', 'TBFFrameBuffer', 'CORFrameBuffer',
           'VDIFFrameBuffer']


def _cmp_frames(x, y):
    """
    Function to compare two frames and sort by stand/beam number (TBN, DRX)
    or by channel number (TBF).
    """
    
    # Parse if frame IDs to extract the stand/beam, tunning, and polarization
    # information (where appropriate)
    try:
        idsX = x.id
        idsY = y.id
    except AttributeError:
        idsX = x.header.first_chan
        idsY = y.header.first_chan
        
    # Do a try...except block to catch TBW vs. TBN/DRX
    try:
        len1 = len(idsX)
        if len1 == 2:
            sX = 2*(idsX[0]-1) + idsX[1]
            sY = 2*(idsY[0]-1) + idsY[1]
        else:
            sX = 16777216*(idsX[0]-1) + 4096*(idsX[1]-1) + idsX[2]
            sY = 16777216*(idsY[0]-1) + 4096*(idsY[1]-1) + idsY[2]
    except TypeError:
        sX = idsX
        sY = idsY
        
    # Do the comparison
    if sY > sX:
        return -1
    elif sX > sY:
        return 1
    else:
        return 0
    

class FrameBufferBase(object):
    """
    Frame buffer for re-ordering TBN and DRX frames in time order.  
    This class is filled with frames and a returns a frame list when 
    the 'nsegments' starts filling.  In that case, the oldest segment 
    is returned.

    The buffer also keeps track of what has already been read out so 
    that tardy frames are just dropped.  For buffers that are read out,
    missing frames are replaced with frames filled with zeros.
    
    .. note::
        Due to the nature of the buffer, it is possible that there will
        still be 'nsegments'-1 segements in the buffer that are either
        full or partially full.  This can be retrieved using the buffer's 
        'flush()' function.
    """
    
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, mode='TBN', stands=[], beams=[], tunes=[], pols=[], chans=[], threads=[], nsegments=6, reorder=False):
        """
        Initialize the buffer with a list of:
          * TBN
            * list of stands
            * list of pols
          * DRX
            * list of beams
            * list of tunings
            * list of pols.
          * TBF
            * list of start channels
          * COR
            * list of stands
            * list of start channels
          * VDIF
            * list of thread IDs
        By doing this, we should be able to keep up with when the buffer 
        is full and to help figure out which stands/beam/tunings/pols. are 
        missing.
        """
        
        # Input validation
        if mode.upper() not in ('TBN', 'DRX', 'DRX8', 'TBF', 'COR', 'VDIF'):
            raise RuntimeError(f"Invalid observing mode '{mode}'")
            
        if mode.upper() == 'TBN':
            for pol in pols:
                if pol not in (0, 1):
                    raise RuntimeError(f"Invalid polarization '{pol}'")
                    
        elif mode.upper() in ('DRX', 'DRX8'):
            for tune in tunes:
                if tune not in (1, 2):
                    raise RuntimeError(f"Invalid tuning '{tune}'")
            for pol in pols:
                if pol not in (0, 1):
                    raise RuntimeError(f"Invalid polarization '{pol}'")
                    
        elif mode.upper() == 'TBF':
            for chan in chans:
                if chan not in range(4096):
                    raise RuntimeError(f"Invalid start channel '{chan}'")
                    
        elif mode.upper() == 'COR':
            for chan in chans:
                if chan not in range(4096):
                    raise RuntimeError(f"Invalid start channel '{chan}'")
            
        else:
            for thread in threads:
                if thread not in range(1024):
                    raise RuntimeError(f"Invalid thread ID '{thread}'")
                    
        # The buffer itself
        self.nsegments = nsegments
        self.buffer = OrderedDict()
        self.done = deque([0,], maxlen=self.nsegments)
        
        # Buffer statistics
        self.full = 0		# Number of times a full buffer was emptied
        self.partial = 0	# Number of times a partial buffer was emptied
        self.dropped = 0	# Number of late or duplicate frames dropped
        self.missing = 0	# Number of missing frames
        
        # Information describing the observations
        self.mode = mode.upper()
        self.stands = stands
        self.beams = beams
        self.tunes = tunes
        self.pols = pols
        self.chans = chans
        self.threads = threads
        
        # If we should reorder the returned frames by stand/pol or not
        self.reorder = reorder
        
        # Figure out how many frames fill the buffer and the list of all
        # possible frames in the data set
        self.nFrames, self.possibleFrames = self.get_max_frames()
        self.possibleFrames = set(self.possibleFrames)
        
    @property
    def filled(self):
        """
        Indicates whether or not the ring buffer is full or not.
        
        .. versionadded:: 1.2.4
        """
        
        return False if len(self.buffer) < self.nsegments else True
        
    @property
    def overfilled(self):
        """
        Indicates whether or not the ring buffer has too many segements or not.
        
        .. versionadded:: 1.2.4
        """
        
        return True if len(self.buffer) > self.nsegments else False
        
    @abc.abstractmethod
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        
        This will be overridden by sub-classes of FrameBufferBase.
        """
        
        raise NotImplementedError
        
    @abc.abstractmethod
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for storing/sorting frames in the ring buffer.
        
        This will be overridden by sub-classes of FrameBufferBase.
        """
        
        raise NotImplementedError
        
    @abc.abstractmethod
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        
        This will be overridden by sub-classes of FrameBufferBase.
        """
        
        raise NotImplementedError
        
    def append(self, frames):
        """
        Append a new frame to the buffer with the appropriate time tag.  
        True is returned if the frame was added to the buffer and False if 
        the frame was dropped because it belongs to a buffer that has 
        already been returned.
        """
        
        # Convert input to a deque (if needed)
        typeName = type(frames).__name__
        if typeName == 'deque':
            pass
        elif typeName == 'list':
            frames = deque(frames)
        else:
            frames = deque([frames,])
            
        # Loop over frames
        for frame in frames:
            # Make sure that it is not in the `done' list.  If it is,
            # disgaurd the frame and make a note of it.
            fom = self.get_figure_of_merit(frame)
            if fom <= self.done[-1]:
                self.dropped += 1
                continue
                
            # If that time tag hasn't been done yet, add it to the 
            # buffer in the correct place.
            try:
                self.buffer[fom].append(frame)
            except KeyError:
                self.buffer[fom] = deque()
                self.buffer[fom].append(frame)
                
        return True
        
    def put(self, *args, **kwds):
        """
        Synonymous with 'append'.
        """
        
        self.append(*args, **kwds)
        
    def peek(self, require_filled=True):
        """
        Peek into the buffer to see what the next key to be retruned by a 
        call to get() will be.  Returns None if the buffer is not full and
+        `require_filled` is True.
        
        .. versionadded:: 1.2.4
        """
        
        # If the ring is not full, return nothing
        if require_filled and not self.filled:
            return None
            
        # Get the current status of the buffer
        keys = list(self.buffer.keys())
        try:
            key_to_return = min(keys)
        except ValueError:
            key_to_return = None
        return key_to_return
        
    def get(self, key_to_return=None):
        """
        Return a list of frames that consitute a 'full' buffer.  Afterwards, 
        delete that buffer and mark it as closed so that any missing frames that
        are recieved late are dropped.  If none of the buffers are ready to be 
        dumped, None is returned.
        """
        
        if key_to_return is None:
            key_to_return = self.peek()
            if key_to_return is None:
                return None
        returnCount = len(self.buffer[key_to_return])
        
        if returnCount == self.nFrames:
            ## Everything is good (Is it really???)
            self.full = self.full + 1
            
            output = self.buffer[key_to_return]
        elif returnCount < self.nFrames:
            ## There are too few frames
            self.partial = self.partial + 1
            self.missing = self.missing + (self.nFrames - returnCount)
            
            output = self.buffer[key_to_return]
            
            ## Fill in the missing frames
            output.extend( map(lambda x: self.create_fill(key_to_return, x), self._get_missing_list(key_to_return)) )
        else:
            ## There are too many frames
            self.full = self.full + 1
            
            output = []
            frameIDs = []
            for frame in self.buffer[key_to_return]:
                try:
                    newID = frame.id
                except AttributeError:
                    newID = frame.header.first_chan
                if newID not in frameIDs:
                    output.append(frame)
                    frameIDs.append(newID)
                else:
                    self.dropped += 1
                    
        del(self.buffer[key_to_return])
        self.done.append(key_to_return)
        
        # Sort and return
        if self.reorder:
            output = list(output)
            output.sort(key=cmp_to_key(_cmp_frames))
        return output
        
    def flush(self):
        """
        Generator to return all of the remaining frames in the buffer 
        from buffers that are considered 'full'.  Afterwards, delete all 
        buffers.  This is useful for emptying the buffer after reading 
        in all of the data.
        
        .. note::
            It is possible for this function to return list of packets that
            are mostly invalid.
            
        .. versionchanged::1.2.1
            Converted to a generator to make sure that we can get all of 
            the LWA-SV COR frames out of the buffer
        """
        
        for key in list(self.buffer):
            yield self.get(key_to_return=key)
            
    def reset(self):
        """
        Emtpy the contents of the buffer and reset it to a clean state.
        
        .. versionadded::1.2.1
        """
        
        self.buffer.clear()
        self.done.clear()
        self.done.append(0)
        
    def is_empty(self):
        """
        Determine if there is anything in the buffer or not.  Returns False 
        if there is, True if there is not.
        
        .. versionadded::1.2.1
        """
        
        return False if len(self.buffer) else True
        
    def _get_missing_list(self, key):
        """
        Create a list of tuples of missing frame information.
        """
        
        # Find out what frames we have
        if hasattr(self, 'frameID'):
            fnc = lambda x: self.frameID(x)
        else:
            if self.mode == 'VDIF':
                fnc = lambda x: x.id[1]
            elif self.mode == 'TBF':
                fnc = lambda x: x.header.first_chan
            elif self.mode == 'COR':
                fnc = lambda x: x.id+(x.header.first_chan,)
            else:
                fnc = lambda x: x.id
        frameList = set(map(fnc, self.buffer[key]))
        
        # Compare the existing list with the possible list stored in the 
        # FrameBufferBase object to build a list of missing frames.
        missingList = self.possibleFrames.difference(frameList)
        
        return missingList
        
    def status(self):
        """
        Print out the status of the buffer.  This contains information about:
          1. The current buffer fill level
          2. The numer of full and partial buffer dumps preformed
          3. The number of missing frames that fill packets needed to be created
             for
          4. The number of frames that arrived too late to be incorporated into 
             one of the ring buffers
        """
        
        nf = 0
        for key in self.buffer:
            nf = nf + len(self.buffer[key])
            
        outString = ''
        outString = '\n'.join([outString, "Current buffer level:  %i frames" % nf])
        outString = '\n'.join([outString, "Buffer dumps:  %i full / %i partial" % (self.full, self.partial)])
        outString = '\n'.join([outString, "--"])
        outString = '\n'.join([outString, "Missing frames:  %i" % self.missing])
        outString = '\n'.join([outString, "Dropped frames:  %i" % self.dropped])
        
        print(outString)


class TBNFrameBuffer(FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with TBN frames.
    See :class:`lsl.reader.buffer.FrameBufferBase` for a description of how the 
    buffering is implemented.
    
    Keywords:
      stands
        list of stand to expect packets for
    
      pols
        list of polarizations to expect packets for
    
      nsegments
        number of ring segments to use for the buffer (default is 20)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        stand/polarization (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+------------------------------------------------+
    |          |                TBN Filter Code                 |
    | Segments +------+------+------+------+------+------+------+
    |          |  1   |  2   |  3   |  4   |  5   |  6   |  7   |
    +----------+------+------+------+------+------+------+------+
    |    10    |  5.1 |  1.6 |  0.8 |  0.4 |  0.2 |  0.1 | 0.05 |
    +----------+------+------+------+------+------+------+------+
    |    20    | 10.2 |  3.3 |  1.6 |  0.8 |  0.4 |  0.2 | 0.10 |
    +----------+------+------+------+------+------+------+------+
    |    30    | 15.4 |  4.9 |  2.5 |  1.2 |  0.6 |  0.3 | 0.15 |
    +----------+------+------+------+------+------+------+------+
    |    40    | 20.5 |  6.6 |  3.3 |  1.6 |  0.8 |  0.4 | 0.20 |
    +----------+------+------+------+------+------+------+------+
    |    50    | 25.6 |  8.2 |  4.1 |  2.0 |  1.0 |  0.5 | 0.26 |
    +----------+------+------+------+------+------+------+------+
    |   100    | 51.2 | 16.4 |  8.2 |  4.1 |  2.0 |  1.0 | 0.51 |
    +----------+------+------+------+------+------+------+------+
    
    """
    
    def __init__(self, stands=[], pols=[0, 1], nsegments=20, reorder=False):
        FrameBufferBase.__init__(self, mode='TBN', stands=stands, pols=pols, nsegments=nsegments, reorder=reorder)
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.stands)*len(self.pols)
        for stand in self.stands:
            for pol in self.pols:
                frameList.append((stand,pol))
                
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For TBN it is:
            <frame timetag in ticks>
        """
        
        return frame.payload.timetag
        
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        stand, pol = frameParameters
        fillFrame.header.tbn_id = 2*(stand-1) + pol + 1
        
        # Zero the data for the fill packet
        try:
            fillFrame.payload._data *= 0
        except TypeError:
            fillFrame.payload._data[...] = 0
            
        # Invalidate the frame
        fillFrame.valid = False
        
        return fillFrame


class DRXFrameBuffer(FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with DRX frames.
    See :class:`lsl.reader.buffer.FrameBufferBase` for a description of how the 
    buffering is implemented.
    
    Keywords:
      beams
        list of beam to expect packets for
    
      tunes
        list of tunings to expect packets for
    
      pols
        list of polarizations to expect packets for
    
      nsegments
        number of ring segments to use for the buffer (default is 20)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        stand/polarization (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------------------------------------------------+
    |          |                 DRX Filter Code                  |
    | Segments +------+------+------+------+------+-------+-------+
    |          |  1   |  2   |  3   |  4   |  5   |  6    |  7    |
    +----------+------+------+------+------+------+-------+-------+
    |    10    | 0.16 | 0.08 | 0.04 | 0.02 | 0.01 | 0.004 | 0.002 |
    +----------+------+------+------+------+------+-------+-------+
    |    20    | 0.33 | 0.16 | 0.08 | 0.04 | 0.02 | 0.008 | 0.004 |
    +----------+------+------+------+------+------+-------+-------+
    |    30    | 0.49 | 0.25 | 0.12 | 0.06 | 0.03 | 0.013 | 0.006 |
    +----------+------+------+------+------+------+-------+-------+
    |    40    | 0.66 | 0.33 | 0.16 | 0.08 | 0.03 | 0.017 | 0.008 |
    +----------+------+------+------+------+------+-------+-------+
    |    50    | 0.82 | 0.41 | 0.20 | 0.10 | 0.04 | 0.021 | 0.010 |
    +----------+------+------+------+------+------+-------+-------+
    |   100    | 1.64 | 0.82 | 0.41 | 0.20 | 0.08 | 0.042 | 0.021 |
    +----------+------+------+------+------+------+-------+-------+
    
    """
    
    def __init__(self, beams=[], tunes=[1,2], pols=[0, 1], nsegments=20, reorder=False):
        FrameBufferBase.__init__(self, mode='DRX', beams=beams, tunes=tunes, pols=pols, nsegments=nsegments, reorder=reorder)
        
        # Make sure that we ignore packets before 1980
        # (although I don't really know why this happens sometimes)
        self.done[0] = 61849368000000000
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.beams)*len(self.tunes)*len(self.pols)
        for beam in self.beams:
            for tune in self.tunes:
                for pol in self.pols:
                    frameList.append((beam,tune,pol))
                    
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For DRX it is:
            <frame timetag in ticks>
        """
        
        return frame.payload.timetag
        
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        beam, tune, pol = frameParameters
        fillFrame.header.drx_id = (beam & 7) | ((tune & 7) << 3) | ((pol & 1) << 7)
        
        # Zero the data for the fill packet
        try:
            fillFrame.payload._data *= 0
        except TypeError:
            fillFrame.payload._data[...] = 0
        
        # Invalidate the frame
        fillFrame.valid = False
        
        return fillFrame


class DRX8FrameBuffer(FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with DRX8 frames.
    See :class:`lsl.reader.buffer.FrameBufferBase` for a description of how the 
    buffering is implemented.
    
    Keywords:
      beams
        list of beam to expect packets for
    
      tunes
        list of tunings to expect packets for
    
      pols
        list of polarizations to expect packets for
    
      nsegments
        number of ring segments to use for the buffer (default is 20)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        stand/polarization (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------------------------------------------------+
    |          |                DRX8 Filter Code                  |
    | Segments +------+------+------+------+------+-------+-------+
    |          |  1   |  2   |  3   |  4   |  5   |  6    |  7    |
    +----------+------+------+------+------+------+-------+-------+
    |    10    | 0.16 | 0.08 | 0.04 | 0.02 | 0.01 | 0.004 | 0.002 |
    +----------+------+------+------+------+------+-------+-------+
    |    20    | 0.33 | 0.16 | 0.08 | 0.04 | 0.02 | 0.008 | 0.004 |
    +----------+------+------+------+------+------+-------+-------+
    |    30    | 0.49 | 0.25 | 0.12 | 0.06 | 0.03 | 0.013 | 0.006 |
    +----------+------+------+------+------+------+-------+-------+
    |    40    | 0.66 | 0.33 | 0.16 | 0.08 | 0.03 | 0.017 | 0.008 |
    +----------+------+------+------+------+------+-------+-------+
    |    50    | 0.82 | 0.41 | 0.20 | 0.10 | 0.04 | 0.021 | 0.010 |
    +----------+------+------+------+------+------+-------+-------+
    |   100    | 1.64 | 0.82 | 0.41 | 0.20 | 0.08 | 0.042 | 0.021 |
    +----------+------+------+------+------+------+-------+-------+
    
    """
    
    def __init__(self, beams=[], tunes=[1,2], pols=[0, 1], nsegments=20, reorder=False):
        FrameBufferBase.__init__(self, mode='DRX8', beams=beams, tunes=tunes, pols=pols, nsegments=nsegments, reorder=reorder)
        
        # Make sure that we ignore packets before 1980
        # (although I don't really know why this happens sometimes)
        self.done[0] = 61849368000000000
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.beams)*len(self.tunes)*len(self.pols)
        for beam in self.beams:
            for tune in self.tunes:
                for pol in self.pols:
                    frameList.append((beam,tune,pol))
                    
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For DRX it is:
            <frame timetag in ticks>
        """
        
        return frame.payload.timetag
        
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        beam, tune, pol = frameParameters
        fillFrame.header.drx_id = (beam & 7) | ((tune & 7) << 3) | 0x40 | ((pol & 1) << 7)
        
        # Zero the data for the fill packet
        try:
            fillFrame.payload._data *= 0
        except TypeError:
            fillFrame.payload._data[...] = 0
        
        # Invalidate the frame
        fillFrame.valid = False
        
        return fillFrame


class TBFFrameBuffer(FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with TBF frames.
    See :class:`lsl.reader.buffer.FrameBufferBase` for a description of how the 
    buffering is implemented.
    
    Keywords:
      chans
        list of start channel numbers to expect data for
    
      nsegments
        number of ring segments to use for the buffer (default is 25)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        start channel (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------+
    | Segments |  Time  |
    +----------+--------+
    |    10    | 0.0004 |
    +----------+--------+
    |    25    | 0.001  |
    +----------+--------+
    |    50    | 0.002  |
    +----------+--------+
    |   100    | 0.004  |
    +----------+--------+
    |   200    | 0.008  |
    +----------+--------+
    
    """
    
    def __init__(self, chans, nsegments=25, reorder=False):
        FrameBufferBase.__init__(self, mode='TBF', chans=chans, nsegments=nsegments, reorder=reorder)
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.chans)
        for chans in self.chans:
            frameList.append(chans)
            
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For TBF this is:
        frame.payload.timetag
        """
        
        return frame.payload.timetag
        
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        chan = frameParameters
        fillFrame.header.first_chan = chan
        
        # Zero the data for the fill packet
        try:
            fillFrame.payload._data *= 0
        except TypeError:
            fillFrame.payload._data[...] = 0
            
        # Invalidate the frame
        fillFrame.valid = False
        
        return fillFrame


class CORFrameBuffer(FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with COR frames.
    See :class:`lsl.reader.buffer.FrameBufferBase` for a description of how the 
    buffering is implemented.
    
    Keywords:
      chans
        list of start channel numbers to expect data for
    
      nsegments
        number of ring segments to use for the buffer (default is 5)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        start channel (default is False)
    
    The number of segements in the ring can be converted to a buffer time in 
    seconds:
    
    +----------+--------+
    | Segments |  Time  |
    +----------+--------+
    |     1    |   10   |
    +----------+--------+
    |     2    |   20   |
    +----------+--------+
    |     5    |   50   |
    +----------+--------+
    
    """
    
    def __init__(self, chans, nsegments=5, reorder=False):
        FrameBufferBase.__init__(self, mode='COR', stands=list(range(1,256+1)), chans=chans, nsegments=nsegments, reorder=reorder)
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.stands)*(len(self.stands)+1)//2 * len(self.chans)
        for stand0 in self.stands:
            for stand1 in self.stands:
                if stand1 < stand0:
                    continue
                for chan in self.chans:
                    frameList.append((stand0,stand1,chan))
                    
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For TBF this is:
        frame.payload.timetag
        """
        
        return frame.payload.timetag
        
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        stand0, stand1, chan = frameParameters
        fillFrame.header.first_chan = chan
        fillFrame.payload.stand0 = stand0
        fillFrame.payload.stand1 = stand1
        
        # Zero the data for the fill packet
        fillFrame.payload._data *= 0
        
        # Invalidate the frame
        fillFrame.valid = False
        
        return fillFrame


class VDIFFrameBuffer(FrameBufferBase):
    """
    A sub-type of FrameBufferBase specifically for dealing with VDIF frames.
    See :class:`lsl.reader.buffer.FrameBufferBase` for a description of how the 
    buffering is implemented.
    
    Keywords:
      threads
        list of thread IDs to expect data for
    
      nsegments
        number of ring segments to use for the buffer (default is 10)
    
      reorder
        whether or not to reorder frames returned by get() or flush() by 
        stand/polarization (default is False)
    """
    
    def __init__(self, threads=[0,1], nsegments=10, reorder=False):
        FrameBufferBase.__init__(self, mode='VDIF', threads=threads, nsegments=nsegments, reorder=reorder)
        
    def get_max_frames(self):
        """
        Calculate the maximum number of frames that we expect from 
        the setup of the observations and a list of tuples that describes
        all of the possible stand/pol combination.
        """
        
        nFrames = 0
        frameList = []
        
        nFrames = len(self.threads)
        for thread in self.threads:
            frameList.append(thread)
            
        return (nFrames, frameList)
        
    def get_figure_of_merit(self, frame):
        """
        Figure of merit for sorting frames.  For VIDF this is:
        seconds_from_epoch * 1000000 + frame_in_second
        """
        
        return frame.header.seconds_from_epoch*1000000 + frame.header.frame_in_second
        
    def create_fill(self, key, frameParameters):
        """
        Create a 'fill' frame of zeros using an existing good
        packet as a template.
        """

        # Get a template based on the first frame for the current buffer
        fillFrame = copy.deepcopy(self.buffer[key][0])
        
        # Get out the frame parameters and fix-up the header
        thread = frameParameters
        fillFrame.header.thread_id = thread
        
        # Zero the data for the fill packet
        fillFrame.payload._data *= 0
        
        # Invalidate the frame
        fillFrame.valid = False
        
        return fillFrame
