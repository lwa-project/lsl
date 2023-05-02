"""
Unit test for the lsl.sim.dp module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import unittest
import numpy
import tempfile
import shutil

import aipy

from lsl.sim import dp
from lsl.reader import tbn
from lsl.reader import drx
from lsl.common import dp as dp_common
from lsl.common import stations as lwa_common


__version__  = "0.3"
__author__    = "Jayce Dowell"


class simdp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.dp
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-simdp-', suffix='.tmp')
        self.src = aipy.src.get_catalog(['cyg',])
        
    def test_basic_tbn(self):
        """Test building a basic TBN signal"""

        testFile = os.path.join(self.testPath, 'tbn.dat')

        station = lwa_common.lwa1
        antennas = station.antennas
        
        fh = open(testFile, 'wb')
        dp.basic_signal(fh, antennas[:8], 2000, station=station, mode='TBN', filter=7, start_time=1000)
        fh.close()

        # Check the file size
        fileSize = os.path.getsize(testFile)
        nSamples = fileSize // tbn.FRAME_SIZE
        self.assertEqual(nSamples, 2000*8)

        # Check the time of the first frame
        fh = open(testFile, 'rb')
        frame = tbn.read_frame(fh)
        fh.close()
        self.assertEqual(frame.payload.timetag, 1000*dp_common.fS)
        
    def test_point_tbn(self):
        """Test building a point source TBN signal"""

        testFile = os.path.join(self.testPath, 'tbn.dat')
        
        station = lwa_common.lwa1
        antennas = station.antennas

        fh = open(testFile, 'wb')
        dp.point_source(fh, antennas[:8], self.src, 4, station=station, mode='TBN', filter=7, start_time=1000)
        fh.close()

        # Check the file size
        fileSize = os.path.getsize(testFile)
        nSamples = fileSize // tbn.FRAME_SIZE
        self.assertEqual(nSamples, 4*8)

        # Check the time of the first frame
        fh = open(testFile, 'rb')
        frame = tbn.read_frame(fh)
        fh.close()
        self.assertEqual(frame.payload.timetag, 1000*dp_common.fS)

    def test_basic_drx(self):
        """Test building a basic DRX signal"""

        testFile = os.path.join(self.testPath, 'drx.dat')

        fh = open(testFile, 'wb')
        dp.basic_signal(fh, numpy.array([1,2,3,4]), 10, mode='DRX', filter=6, ntuning=2, start_time=1000)
        fh.close()

        # Check the file size
        fileSize = os.path.getsize(testFile)
        nSamples = fileSize // drx.FRAME_SIZE
        self.assertEqual(nSamples, 10*4*2*2)

        # Check the file size
        fh = open(testFile, 'rb')
        frame = drx.read_frame(fh)
        fh.close()
        self.assertEqual(frame.payload.timetag, 1000*dp_common.fS)
        self.assertEqual(frame.header.frame_count, 0)
        self.assertEqual(frame.header.second_count, 0)

    def tearDown(self):
        """Remove the test path directory and its contents"""

        shutil.rmtree(self.testPath, ignore_errors=True)


class  simdp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.sim.vis units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(simdp_tests)) 


if __name__ == '__main__':
    unittest.main()
