"""
Unit test for the lsl.sim.tbn/drx modules.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import copy
import unittest
import tempfile
import shutil
import numpy

from lsl.common.paths import DATA_BUILD
from lsl.reader import tbn as tbnReader
from lsl.sim import tbn as tbnWriter
from lsl.reader import drx as drxReader
from lsl.sim import drx as drxWriter


__version__  = "0.4"
__author__    = "Jayce Dowell"


tbnFile = os.path.join(DATA_BUILD, 'tests', 'tbn-test.dat')
drxFile = os.path.join(DATA_BUILD, 'tests', 'drx-test.dat')


class fake_TBN_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.tbn
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""
        
        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-fakedata-', suffix='.tmp')
        
    def test_sim_frame(self):
        """Test the tbn.SimFrame class."""
        
        # Read in a TBN frame from the test file
        fh = open(tbnFile, 'rb')
        origFrame = tbnReader.read_frame(fh)
        fh.close()
        
        fakeFrame = tbnWriter.SimFrame()
        fakeFrame.load_frame(origFrame)
        # Test the validity of the SimFrame
        self.assertTrue(fakeFrame.is_valid())
        
    def test_sim_frame_ci8(self):
        """Test the tbn.SimFrame class, ci8 style."""
        
        # Read in a TBN frame from the test file
        fh = open(tbnFile, 'rb')
        origFrame = tbnReader.read_frame_ci8(fh)
        fh.close()
        
        fakeFrame = tbnWriter.SimFrame()
        fakeFrame.load_frame(origFrame)
        # Test the validity of the SimFrame
        self.assertTrue(fakeFrame.is_valid())
        
    def test_write_frame(self):
        """Test that the TBN data writer works."""
        
        testFile = os.path.join(self.testPath, 'tbn-test-W.dat')
        
        nFrames = os.path.getsize(tbnFile) // tbnReader.FRAME_SIZE
        
        # Read in a TBN frame from the test file
        fh = open(tbnFile, 'rb')
        origFrames = []
        for i in range(nFrames):
            origFrames.append( tbnReader.read_frame(fh) )
        fh.close()
        
        # Write the data to a TBN test frame
        fh = open(testFile, 'wb')
        for origFrame in origFrames:
            rawFrame = tbnWriter.frame_to_frame(origFrame)
            rawFrame.tofile(fh)
        fh.close()
        
        # Read in the 
        fh = open(testFile, 'rb')
        fakeFrames = []
        for i in range(nFrames):
            fakeFrames.append( tbnReader.read_frame(fh) )
        fh.close()
        
        for fakeFrame,origFrame in zip(fakeFrames, origFrames):
            # Test values returned by info functions
            self.assertEqual(fakeFrame.id[0], origFrame.id[0])
            self.assertEqual(fakeFrame.id[1], origFrame.id[1])
            
            # Test raw header values
            self.assertTrue(fakeFrame.header.is_tbn)
            self.assertEqual(fakeFrame.header.frame_count, origFrame.header.frame_count)
            self.assertEqual(fakeFrame.header.tuning_word, origFrame.header.tuning_word)
            
            # Test raw data values
            self.assertEqual(fakeFrame.payload.timetag, origFrame.payload.timetag)
            for i in range(512):
                self.assertEqual(fakeFrame.payload.data[i].real, origFrame.payload.data[i].real)
                self.assertEqual(fakeFrame.payload.data[i].imag, origFrame.payload.data[i].imag)
                
    def test_write_frame_ci8(self):
        """Test that the TBN data writer works, ci8 style."""
        
        testFile = os.path.join(self.testPath, 'tbn-test-W.dat')
        
        nFrames = os.path.getsize(tbnFile) // tbnReader.FRAME_SIZE
        
        # Read in a TBN frame from the test file
        fh = open(tbnFile, 'rb')
        origFrames = []
        for i in range(nFrames):
            origFrames.append( tbnReader.read_frame_ci8(fh) )
        fh.close()
        
        # Write the data to a TBN test frame
        fh = open(testFile, 'wb')
        for origFrame in origFrames:
            rawFrame = tbnWriter.frame_to_frame(origFrame)
            rawFrame.tofile(fh)
        fh.close()
        
        # Read in the 
        fh = open(testFile, 'rb')
        fakeFrames = []
        for i in range(nFrames):
            fakeFrames.append( tbnReader.read_frame_ci8(fh) )
        fh.close()
        
        for fakeFrame,origFrame in zip(fakeFrames, origFrames):
            # Test values returned by info functions
            self.assertEqual(fakeFrame.id[0], origFrame.id[0])
            self.assertEqual(fakeFrame.id[1], origFrame.id[1])
            
            # Test raw header values
            self.assertTrue(fakeFrame.header.is_tbn)
            self.assertEqual(fakeFrame.header.frame_count, origFrame.header.frame_count)
            self.assertEqual(fakeFrame.header.tuning_word, origFrame.header.tuning_word)
            
            # Test raw data values
            self.assertEqual(fakeFrame.payload.timetag, origFrame.payload.timetag)
            for i in range(512):
                self.assertEqual(fakeFrame.payload.data[i]['re'], origFrame.payload.data[i]['re'])
                self.assertEqual(fakeFrame.payload.data[i]['im'], origFrame.payload.data[i]['im'])
                
    def test_frame_data_errors(self):
        """Test the data error scenarios when validating a TBN SimFrame."""
        
        # Read in a TBN frame from the test file
        fh = open(tbnFile, 'rb')
        origFrame = tbnReader.read_frame(fh)
        fh.close()
        
        # Try to validate frame with the wrong data type
        fakeFrame = tbnWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.data = fakeFrame.payload.data.real
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
        # Try to validate frame with the wrong data size
        fakeFrame = tbnWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.data = None
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        fakeFrame = tbnWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.data = fakeFrame.payload.data[0:50]
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
    def test_frame_header_errors(self):
        """Test the header error scenarios when validating a TBN SimFrame."""
        
        # Read in a TBN frame from the test file
        fh = open(tbnFile, 'rb')
        origFrame = tbnReader.read_frame(fh)
        fh.close()
        
        # Try to validate frame with the wrong stand number
        fakeFrame = tbnWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.stand = 300
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
    def tearDown(self):
        """Remove the test path directory and its contents"""
        
        shutil.rmtree(self.testPath, ignore_errors=True)


class fake_DRX_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.drx
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""
        
        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-fakedata-', suffix='.tmp')
        
    def test_sim_frame(self):
        """Test the drx.SimFrame class."""
        
        # Read in a DRX frame from the test file
        fh = open(drxFile, 'rb')
        origFrame = drxReader.read_frame(fh)
        fh.close()
        
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(origFrame)
        # Test the validity of the SimFrame
        self.assertTrue(fakeFrame.is_valid())
        
    def test_sim_frame_ci8(self):
        """Test the drx.SimFrame class, ci8 style."""
        
        # Read in a DRX frame from the test file
        fh = open(drxFile, 'rb')
        origFrame = drxReader.read_frame_ci8(fh)
        fh.close()
        
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(origFrame)
        # Test the validity of the SimFrame
        self.assertTrue(fakeFrame.is_valid())
        
    def test_write_frame(self):
        """Test that the DRX data writer works."""
        
        testFile = os.path.join(self.testPath, 'drx-test-W.dat')
        
        nFrames = os.path.getsize(drxFile) // drxReader.FRAME_SIZE
        
        # Read in a TBN frame from the test file
        fh = open(drxFile, 'rb')
        origFrames = []
        for i in range(nFrames):
            origFrames.append( drxReader.read_frame(fh) )
        fh.close()
        
        # Write the data to a TBN test frame
        fh = open(testFile, 'wb')
        for origFrame in origFrames:
            rawFrame = drxWriter.frame_to_frame(origFrame)
            rawFrame.tofile(fh)
        fh.close()
        
        # Read in the 
        fh = open(testFile, 'rb')
        fakeFrames = []
        for i in range(nFrames):
            fakeFrames.append( drxReader.read_frame(fh) )
        fh.close()
        
        for fakeFrame,origFrame in zip(fakeFrames, origFrames):
            # Test values returned by info functions
            self.assertEqual(fakeFrame.id[0], origFrame.id[0])
            self.assertEqual(fakeFrame.id[1], origFrame.id[1])
            self.assertEqual(fakeFrame.id[2], origFrame.id[2])
            self.assertAlmostEqual(fakeFrame.sample_rate, origFrame.sample_rate, 4)
            
            # Test raw header values
            self.assertEqual(fakeFrame.header.second_count, origFrame.header.second_count)
            self.assertEqual(fakeFrame.header.decimation, origFrame.header.decimation)
            self.assertEqual(fakeFrame.header.time_offset, origFrame.header.time_offset)
            
            # Test raw data values
            self.assertEqual(fakeFrame.payload.timetag, origFrame.payload.timetag)
            self.assertEqual(fakeFrame.payload.flags, origFrame.payload.flags)
            for i in range(4096):
                self.assertEqual(fakeFrame.payload.data[i].real, origFrame.payload.data[i].real)
                self.assertEqual(fakeFrame.payload.data[i].imag, origFrame.payload.data[i].imag)
                
    def test_write_frame_ci8(self):
        """Test that the DRX data writer works, ci8 style."""
        
        testFile = os.path.join(self.testPath, 'drx-test-W.dat')
        
        nFrames = os.path.getsize(drxFile) // drxReader.FRAME_SIZE
        
        # Read in a TBN frame from the test file
        fh = open(drxFile, 'rb')
        origFrames = []
        for i in range(nFrames):
            origFrames.append( drxReader.read_frame_ci8(fh) )
        fh.close()
        
        # Write the data to a TBN test frame
        fh = open(testFile, 'wb')
        for origFrame in origFrames:
            rawFrame = drxWriter.frame_to_frame(origFrame)
            rawFrame.tofile(fh)
        fh.close()
        
        # Read in the 
        fh = open(testFile, 'rb')
        fakeFrames = []
        for i in range(nFrames):
            fakeFrames.append( drxReader.read_frame_ci8(fh) )
        fh.close()
        
        for fakeFrame,origFrame in zip(fakeFrames, origFrames):
            # Test values returned by info functions
            self.assertEqual(fakeFrame.id[0], origFrame.id[0])
            self.assertEqual(fakeFrame.id[1], origFrame.id[1])
            self.assertEqual(fakeFrame.id[2], origFrame.id[2])
            self.assertAlmostEqual(fakeFrame.sample_rate, origFrame.sample_rate, 4)
            
            # Test raw header values
            self.assertEqual(fakeFrame.header.second_count, origFrame.header.second_count)
            self.assertEqual(fakeFrame.header.decimation, origFrame.header.decimation)
            self.assertEqual(fakeFrame.header.time_offset, origFrame.header.time_offset)
            
            # Test raw data values
            self.assertEqual(fakeFrame.payload.timetag, origFrame.payload.timetag)
            self.assertEqual(fakeFrame.payload.flags, origFrame.payload.flags)
            for i in range(4096):
                self.assertEqual(fakeFrame.payload.data[i]['re'], origFrame.payload.data[i]['re'])
                self.assertEqual(fakeFrame.payload.data[i]['im'], origFrame.payload.data[i]['im'])
                
    def test_frame_data_errors(self):
        """Test the data error scenarios when validating a DRX SimFrame."""
        
        # Read in a DRX frame from the test file
        fh = open(drxFile, 'rb')
        origFrame = drxReader.read_frame(fh)
        fh.close()
        
        # Try to validate frame with the wrong data type
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.data = fakeFrame.payload.data.real
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
        # Try to validate frame with the wrong data size
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.data = None
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.data = fakeFrame.payload.data[0:50]
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
    def test_frame_header_errors(self):
        """Test the header error scenarios when validating a DRX SimFrame."""
        
        # Read in a DRX frame from the test file
        fh = open(drxFile, 'rb')
        origFrame = drxReader.read_frame(fh)
        fh.close()
        
        # Try to validate frame with the wrong beam number
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.beam = 5
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
        # Try to validate frame with the wrong tuning number
        fakeFrame = drxWriter.SimFrame()
        fakeFrame.load_frame(copy.deepcopy(origFrame))
        fakeFrame.tune = 3
        self.assertRaises(ValueError, fakeFrame.is_valid, raise_errors=True)
        
    def tearDown(self):
        """Remove the test path directory and its contents"""
        
        shutil.rmtree(self.testPath, ignore_errors=True)


class fakedata_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(fake_TBN_tests))
        self.addTests(loader.loadTestsFromTestCase(fake_DRX_tests))


if __name__ == '__main__':
    unittest.main()
