"""
Unit tests for future formats in the lsl.reader modules.
"""

import os
import numpy as np
import unittest
from datetime import timedelta

from lsl.reader import drx8, errors
from lsl.reader.base import FrameTimestamp


__version__  = "0.1"
__author__    = "Jayce Dowell"


drx8File = os.path.join(os.path.dirname(__file__), 'data', 'drx8-sim-test.dat')


class reader_future_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for future formats in the
    lsl.reader modules."""
    
    ### DRX8 ###
    
    def test_drx8_read(self):
        """Test reading in a frame from a DRX8 file."""
        
        fh = open(drx8File, 'rb')
        # First frame is really DRX and stores the IDs
        frame1 = drx8.read_frame(fh)
        beam, tune, pol = frame1.id
        self.assertEqual(beam, 4)
        self.assertEqual(tune, 1)
        self.assertEqual(pol,  1)
        # Second frame
        frame2 = drx8.read_frame(fh)
        beam, tune, pol = frame2.id
        self.assertEqual(beam, 4)
        self.assertEqual(tune, 2)
        self.assertEqual(pol,  0)
        fh.close()
        
    def test_drx8_read_ci8(self):
        """Test reading in a frame from a DRX file, ci8 style."""
        
        fh = open(drx8File, 'rb')
        frame1 = drx8.read_frame(fh)
        frame2 = drx8.read_frame(fh)
        fh.close()
        
        fh = open(drx8File, 'rb')
        frame3 = drx8.read_frame_ci8(fh)
        frame4 = drx8.read_frame_ci8(fh)
        fh.close()
        
        # Compare
        data1 = frame3.payload.data['re'] + 1j*frame3.payload.data['im']
        data2 = frame4.payload.data['re'] + 1j*frame4.payload.data['im']
        for i in range(512):
            self.assertAlmostEqual(frame1.payload.data[i], data1[i], 1e-6)
            self.assertAlmostEqual(frame2.payload.data[i], data2[i], 1e-6)
            
    def test_drx_errors(self):
        """Test reading in all frames from a truncated DRX file."""
        
        fh = open(drx8File, 'rb')
        # Frames 1 through 32
        for i in range(1,33):
            frame = drx8.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, drx8.read_frame, fh)
        fh.close()
        
        # If we offset in the file by 1 byte, we should be a 
        # sync error (errors.SyncError).
        fh = open(drx8File, 'rb')
        fh.seek(1)
        self.assertRaises(errors.SyncError, drx8.read_frame, fh)
        fh.close()
        
    def test_drx_beam(self):
        """Test finding out how many beams are present in a DRX file."""
        
        fh = open(drx8File, 'rb')
        nBeam = drx8.get_beam_count(fh)
        self.assertEqual(nBeam, 1)
        fh.close()
        
    def test_drx_block(self):
        """Test finding out how many tunings/pols. per beam are in a DRX file."""
        
        fh = open(drx8File, 'rb')
        b1, b2, b3, b4 = drx8.get_frames_per_obs(fh)
        self.assertEqual(b1, 0)
        self.assertEqual(b2, 0)
        self.assertEqual(b3, 0)
        self.assertEqual(b4, 4)
        fh.close()
        
    def test_drx_rate(self):
        """Test finding out the DRX sample rate."""
        
        fh = open(drx8File, 'rb')
        cFrame = drx8.read_frame(fh)
        fh.seek(0)
        
        # Sample rate
        self.assertEqual(cFrame.sample_rate, drx8.get_sample_rate(fh))
        
        # Filter code
        self.assertEqual(cFrame.filter_code, drx8.get_sample_rate(fh, filter_code=True))
        fh.close()
        
    def test_drx_comps(self):
        """Test the DRX frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(drx8File, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(drx8.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertFalse(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertFalse(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_drx_sort(self):
        """Test sorting DRX frames by time tags."""
        
        fh = open(drx8File, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(drx8.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in range(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_drx_math(self):
        """Test mathematical operations on DRX frame data via frames."""
        
        fh = open(drx8File, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(drx8.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        np.testing.assert_allclose(frameT.payload.data, 2*frames[0].payload.data, atol=1e-6)
        frameT *= 2.0
        np.testing.assert_allclose(frameT.payload.data, 4*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] * frames[1]
        np.testing.assert_allclose(frameT.payload.data, frames[0].payload.data*frames[1].payload.data, atol=1e-6)
        
        # Division
        frameT = frames[0] / 2.0
        np.testing.assert_allclose(frameT.payload.data, 0.5*frames[0].payload.data, atol=1e-6)
        frameT /= 2.0
        np.testing.assert_allclose(frameT.payload.data, 0.25*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] / frames[1]
        np.testing.assert_allclose(frameT.payload.data, frames[0].payload.data/frames[1].payload.data, atol=1e-6)
        
        # Addition
        frameA = frames[0] + 2.0
        np.testing.assert_allclose(frameA.payload.data, 2+frames[0].payload.data, atol=1e-6)
        frameA += 2.0
        np.testing.assert_allclose(frameA.payload.data, 4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] + frames[1]
        np.testing.assert_allclose(frameA.payload.data, frames[0].payload.data+frames[1].payload.data, atol=1e-6)
        
        # Subtraction
        frameA = frames[0] - 2.0
        np.testing.assert_allclose(frameA.payload.data, -2+frames[0].payload.data, atol=1e-6)
        frameA -= 2.0
        np.testing.assert_allclose(frameA.payload.data, -4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] - frames[1]
        np.testing.assert_allclose(frameA.payload.data, frames[0].payload.data-frames[1].payload.data, atol=1e-6)


class reader_future_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all future formats lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(reader_future_tests)) 


if __name__ == '__main__':
    unittest.main()
