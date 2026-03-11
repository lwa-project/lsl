"""
Unit test for the lsl.sim.ndp module.
"""

import os
import unittest
import numpy as np
import tempfile
import shutil

from lsl.sim import ndp
from lsl.reader import drx
from lsl.common import ndp as ndp_common
from lsl.common import stations as lwa_common


__version__  = "0.1"
__author__    = "Jayce Dowell"


class simdp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.ndp
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""
        
        np.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-simndp-', suffix='.tmp')
        
    def test_basic_drx(self):
        """Test building a basic DRX signal"""
        
        testFile = os.path.join(self.testPath, 'drx.dat')
        
        fh = open(testFile, 'wb')
        ndp.basic_signal(fh, np.array([1,2,3,4]), 10, mode='DRX', filter=6, ntuning=2, start_time=1000)
        fh.close()
        
        # Check the file size
        fileSize = os.path.getsize(testFile)
        nSamples = fileSize // drx.FRAME_SIZE
        self.assertEqual(nSamples, 10*4*2*2)
        
        # Check the file size
        fh = open(testFile, 'rb')
        frame = drx.read_frame(fh)
        fh.close()
        self.assertEqual(frame.payload.timetag, 1000*ndp_common.fS)
        self.assertEqual(frame.header.frame_count, 0)
        self.assertEqual(frame.header.second_count, 0)
        
    def tearDown(self):
        """Remove the test path directory and its contents"""
        
        shutil.rmtree(self.testPath, ignore_errors=True)


class  simndp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.sim.ndp units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(simndp_tests)) 


if __name__ == '__main__':
    unittest.main()
