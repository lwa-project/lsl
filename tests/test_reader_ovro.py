"""
Unit tests for the lsl.reader modules that are for OVRO-LWA.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import numpy
import unittest

from lsl.common.paths import DATA_BUILD
from lsl.reader import ovro
from lsl.reader import errors


__version__  = "0.1"
__author__    = "Jayce Dowell"


ovroFile = os.path.join(DATA_BUILD, 'tests', 'ovro-test.dat')


class reader_ovro_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### Triggered voltage buffer dump ###
    
    def test_ovro_read(self):
        """Test reading in a frame from a OVRO-LWA triggered voltage dump file."""
        
        fh = open(ovroFile, 'rb')
        # First frame stores the first channel
        frame1 = ovro.read_frame(fh)
        self.assertEqual(frame1.header.first_chan, 2288)
        # Second frame
        frame2 = ovro.read_frame(fh)
        self.assertEqual(frame2.header.first_chan, 2288)
        fh.close()
        
    def test_ovro_errors(self):
        """Test OVRO-LWA triggered voltage dump reading errors."""
        
        fh = open(ovroFile, 'rb')
        # Frames 1 through 10
        for i in range(1,11):
            frame = ovro.read_frame(fh)
            
        # Last frame should be an error (errors.EOFError)
        self.assertRaises(errors.EOFError, ovro.read_frame, fh)
        fh.close()
        
    def test_ovro_comps(self):
        """Test the OVRO-LWA triggered voltage dump frame comparison operators (>, <, etc.) for time tags."""
        
        fh = open(ovroFile, 'rb')
        # Frames 1 through 10
        frames = []
        for i in range(1,11):
            frames.append(ovro.read_frame(fh))
        fh.close()
        
        self.assertTrue(0 < frames[0])
        self.assertFalse(0 > frames[0])
        self.assertTrue(frames[-1] >= frames[0])
        self.assertTrue(frames[-1] <= frames[0])
        self.assertTrue(frames[0] == frames[0])
        self.assertTrue(frames[0] == frames[-1])
        self.assertFalse(frames[0] != frames[0])
        
    def test_ovro_sort(self):
        """Test sorting OVRO-LWA triggered voltage dump frames by time tags."""
        
        fh = open(ovroFile, 'rb')
        # Frames 1 through 9
        frames = []
        for i in range(1,10):
            frames.append(ovro.read_frame(fh))
        fh.close()
        
        frames.sort()
        frames = frames[::-1]
        
        for i in range(1,len(frames)):
            self.assertTrue( frames[i-1] >= frames[i] )
            
    def test_ovro_math(self):
        """Test mathematical operations on OVRO-LWA triggered voltage dump frame data via frames."""
        
        fh = open(ovroFile, 'rb')
        # Frames 1 through 9
        frames = []
        for i in range(1,10):
            frames.append(ovro.read_frame(fh))
        fh.close()
        
        # Multiplication
        frameT = frames[0] * 2.0
        numpy.testing.assert_allclose(frameT.payload.data, 2*frames[0].payload.data, atol=1e-6)
        frameT *= 2.0
        numpy.testing.assert_allclose(frameT.payload.data, 4*frames[0].payload.data, atol=1e-6)
        frameT = frames[0] * frames[1]
        numpy.testing.assert_allclose(frameT.payload.data, frames[0].payload.data*frames[1].payload.data, atol=1e-6)
        
        # Addition
        frameA = frames[0] + 2.0
        numpy.testing.assert_allclose(frameA.payload.data, 2+frames[0].payload.data, atol=1e-6)
        frameA += 2.0
        numpy.testing.assert_allclose(frameA.payload.data, 4+frames[0].payload.data, atol=1e-6)
        frameA = frames[0] + frames[1]
        numpy.testing.assert_allclose(frameA.payload.data, frames[0].payload.data+frames[1].payload.data, atol=1e-6)


class reader_ovro_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(reader_ovro_tests)) 


if __name__ == '__main__':
    unittest.main()
