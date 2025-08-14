"""
Unit test for the lsl.reader.buffer module.
"""

import os
import unittest

from lsl.reader import drx
from lsl.reader import errors
from lsl.reader import buffer


__version__  = "0.2"
__author__    = "Jayce Dowell"


drxFile = os.path.join(os.path.dirname(__file__), 'data', 'drx-test.dat')


class buffer_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader.buffer
    module."""
    
    #
    # DRX
    #
    
    def test_drx_basic(self):
        """Test the DRX ring buffer."""
        
        fh = open(drxFile, 'rb')
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
        fh.seek(-drx.FRAME_SIZE, 1)
        
        # Create the FrameBuffer instance
        frameBuffer = buffer.DRXFrameBuffer(beams=[b,], tunes=[1, 2], pols=[0, 1], nsegments=2)
        
        # Go
        dumped = []
        while True:
            try:
                cFrame = drx.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue

            frameBuffer.append(cFrame)
            cFrames = frameBuffer.get()

            if cFrames is None:
                continue
                
            dumped.append( cFrames[0].payload.timetag )
            
        fh.close()
        
        # Make sure we have the right number of frames in the buffer
        nFrames = 0
        for key in frameBuffer.buffer.keys():
            nFrames = nFrames + len(frameBuffer.buffer[key])
        self.assertEqual(nFrames, 1)
        self.assertEqual(nFrames+len(dumped)*4, 32+1)
        
        # Make sure nothing has happened that shouldn't have
        self.assertEqual(frameBuffer.full,    7)
        self.assertEqual(frameBuffer.partial, 1)
        self.assertEqual(frameBuffer.missing, 1)
        self.assertEqual(frameBuffer.dropped, 0)
        
        # Make sure we have the right keys
        for key in dumped:
            self.assertTrue(key in (257355782095018376, 257355782095059336, 257355782095100296, 257355782095141256, 257355782095182216, 257355782095223176, 257355782095264136, 257355782095305096))
            
        for key in frameBuffer.buffer.keys():
            self.assertTrue(key in (257355782095346056,))
        
        # Make sure the buffer keys have the right sizes
        self.assertEqual(len(frameBuffer.buffer[257355782095346056]), 1)
        
    def test_drx_reorder(self):
        """Test the reorder function of the TBN ring buffer."""
        
        fh = open(drxFile, 'rb')
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
        fh.seek(-drx.FRAME_SIZE, 1)
        
        # Create the FrameBuffer instance
        frameBuffer = buffer.DRXFrameBuffer(beams=[b,], tunes=[1, 2], pols=[0, 1], nsegments=2, reorder=True)
        
        # Go
        while True:
            try:
                cFrame = drx.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue

            frameBuffer.append(cFrame)
            cFrames = frameBuffer.get()

            if cFrames is None:
                continue
                
            # Make sure it has the right number of frames
            self.assertEqual(len(cFrames), 4)
            
            # Check the order
            for i in range(1, len(cFrames)):
                pB, pT, pP = cFrames[i-1].id
                cB, cT, cP = cFrames[i].id
                
                pID = 4*pB + 2*(pT-1) + pP
                cID = 4*cB + 2*(cT-1) + cP
                self.assertTrue(cID > pID)
                
        fh.close()
        
    def test_drx_buffer_flush(self):
        """Test the DRX ring buffer's flush() function."""
        
        fh = open(drxFile, 'rb')
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
        fh.seek(-drx.FRAME_SIZE, 1)
        
        # Create the FrameBuffer instance
        frameBuffer = buffer.DRXFrameBuffer(beams=[b,], tunes=[1, 2], pols=[0, 1], nsegments=2)
        
        # Go
        while True:
            try:
                cFrame = drx.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue

            frameBuffer.append(cFrame)
            cFrames = frameBuffer.get()

            if cFrames is None:
                continue
                
        fh.close()
        
        # Flush the buffer
        for cFrames in frameBuffer.flush():
            # Make sure the dump has one of the expected time tags
            self.assertTrue(cFrames[0].payload.timetag in (257355782095346056,))
            
            # Make sure it has the right number of frames
            self.assertEqual(len(cFrames), 4)
            
    def test_drx_buffer_reorder_flush(self):
        """Test the DRX ring buffer's flush() function with reordering."""
        
        fh = open(drxFile, 'rb')
        junkFrame = drx.read_frame(fh)
        b,t,p = junkFrame.id
        fh.seek(-drx.FRAME_SIZE, 1)
        
        # Create the FrameBuffer instance
        frameBuffer = buffer.DRXFrameBuffer(beams=[b,], tunes=[1, 2], pols=[0, 1], nsegments=2, reorder=True)
        
        # Go
        while True:
            try:
                cFrame = drx.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                continue

            frameBuffer.append(cFrame)
            cFrames = frameBuffer.get()

            if cFrames is None:
                continue
                
            # Make sure it has the right number of frames
            self.assertEqual(len(cFrames), 4)
            
            # Check the order
            for i in range(1, len(cFrames)):
                pB, pT, pP = cFrames[i-1].id
                cB, cT, cP = cFrames[i].id
                
                pID = 4*pB + 2*(pT-1) + pP
                cID = 4*cB + 2*(cT-1) + cP
                self.assertTrue(cID > pID)
                
        fh.close()
        
        # Flush the buffer
        for cFrames in frameBuffer.flush():
            # Make sure the dump has one of the expected time tags
            self.assertTrue(cFrames[0].payload.timetag in (257355782095346056,))
            
            # Make sure it has the right number of frames
            self.assertEqual(len(cFrames), 4)
            
            # Check the order
            for i in range(1, len(cFrames)):
                pB, pT, pP = cFrames[i-1].id
                cB, cT, cP = cFrames[i].id
                
                pID = 4*pB + 2*(pT-1) + pP
                cID = 4*cB + 2*(cT-1) + cP
                self.assertTrue(cID > pID)
                

class buffer_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader.buffer 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(buffer_tests)) 


if __name__ == '__main__':
    unittest.main()
