"""
Unit test for the lsl.writer.vdif module.
"""

import os
import time
import unittest
import tempfile
import shutil
import numpy as np

from lsl.reader import base, vdif as vrdr, errors
from lsl.writer import vdif


__version__  = "0.3"
__author__    = "Jayce Dowell"


class vdif_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.vdif
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""
        
        np.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-vdif-', suffix='.tmp')
        
    def _get_real(self):
        """Private function to create a set of real data."""
        
        srate = 196000000
        fsize = 1000
        
        frames = []
        tt0 = int(time.time()) * srate
        for i in range(32):
            tt = tt0 + i*fsize*srate
            data = np.random.randn(fsize)*128
            data = np.clip(data, -128, 127).astype(np.int8)
            frames.append({'timetag': base.FrameTimestamp(tt//srate, (tt%srate)/srate),
                           'sample_rate': srate,
                           'bits': 8,
                           'data': data
                          })
            
        return frames
        
    def _get_complex(self):
        """Private function to create a set of complex data."""
        
        srate = 19600000
        fsize = 4000
        
        frames = []
        tt0 = int(time.time()) * srate
        for i in range(32):
            tt = tt0 + i*fsize*srate
            data = np.random.randn(fsize)*8 + 1j*np.random.randn(fsize)*8
            data = np.round(data)
            data.real = np.clip(data.real, -8, 7)
            data.imag = np.clip(data.imag, -8, 7)
            frames.append({'timetag': base.FrameTimestamp(tt//srate, (tt%srate)/srate),
                           'sample_rate': srate,
                           'bits': 4,
                           'data': data
                          })
            
        return frames
        
    def test_vdif_real(self):
        """Test writing real data to VDIF format."""
        
        # Setup the file names
        testFile = os.path.join(self.testPath, 'real-test-W.fits')
                                 
        # Get some data
        frames = self._get_real()

        # Write the data
        with open(testFile, 'wb') as fh:
            for frame in frames:
                vFrame = vdif.Frame(1, frame['timetag'], bits=frame['bits'],
                                    data=frame['data'], sample_rate=frame['sample_rate'])
                vFrame.write_raw_frame(fh)
                
        # Read it back in
        with open(testFile, 'rb') as fh:
            for tFrame in frames:
                vFrame = vrdr.read_frame(fh, sample_rate=frames[0]['sample_rate'])
                self.assertAlmostEqual(vFrame.time, tFrame['timetag'], 6)
                np.testing.assert_allclose((vFrame.payload.data*256-1)/2, tFrame['data'], atol=1e-6)
                
    def test_vdif_complex(self):
        """Test writing complex data to VDIF format."""
        
        # Setup the file names
        testFile = os.path.join(self.testPath, 'complex-test-W.fits')
        
        # Get some data
        frames = self._get_complex()
        
        # Write the data
        with open(testFile, 'wb') as fh:
            for frame in frames:
                ## We need an integer number of frame per second so adjust the sample rate
                vFrame = vdif.Frame(2, frame['timetag'], bits=frame['bits'],
                                    data=frame['data'], sample_rate=frame['sample_rate'])
                vFrame.write_raw_frame(fh)
                
        # Read it back in
        with open(testFile, 'rb') as fh:
            for tFrame in frames:
                vFrame = vrdr.read_frame(fh, sample_rate=frames[0]['sample_rate'])
                self.assertAlmostEqual(vFrame.time, tFrame['timetag'], 6)
                vData = vFrame.payload.data
                vData.real = vData.real*2.95
                vData.imag = vData.imag*2.95
                np.testing.assert_allclose(vData, tFrame['data'], atol=1e-6)
                
    def tearDown(self):
        """Remove the test path directory and its contents"""

        shutil.rmtree(self.testPath, ignore_errors=True)


class  vdif_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.writer.vdif units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(vdif_tests)) 


if __name__ == '__main__':
    unittest.main()
