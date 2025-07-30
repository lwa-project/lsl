"""
Unit test for the lsl.writer.vdif module.
"""

import os
import time
import unittest
import tempfile
import shutil
import numpy as np

from lsl.reader import tbw, tbn, vdif as vrdr, errors
from lsl.writer import vdif


__version__  = "0.2"
__author__    = "Jayce Dowell"

tbwFile = os.path.join(os.path.dirname(__file__), 'data', 'tbw-test.dat')
tbnFile = os.path.join(os.path.dirname(__file__), 'data', 'tbn-test.dat')


class vdif_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.vdif
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        np.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-vdif-', suffix='.tmp')

    def _get_tbw(self):
        """Private function to load in the test TBW data and get the frames."""

        with open(tbwFile, 'rb') as fh:
            # Frames 1 through 8
            frames = []
            for i in range(1,9):
                frames.append(tbw.read_frame(fh))
                
        return frames

    def _get_tbn(self):
        """Private function to load in the test TBN data and get the frames.  If 
        the keyword 'vanilla' is set to True, gain, sample rate, and frequency meta-
        data are not added to the frames."""

        with open(tbnFile, 'rb') as fh:
            # Frames 1 through 8
            frames = []
            for i in range(1,9):
                frames.append(tbn.read_frame(fh))
                # WAR for bad TBN timetags
                frames[-1].payload.timetag += 1673741823*196000000
                frames[-1].payload.timetag //= 19600000
                frames[-1].payload.timetag *= 19600000
                
        return frames

    def test_vdif_real(self):
        """Test writing real data to VDIF format."""

        # Setup the file names
        testFile = os.path.join(self.testPath, 'tbw-test-W.fits')

        # Get some data
        frames = self._get_tbw()

        # Write the data
        with open(testFile, 'wb') as fh:
            for frame in frames:
                vFrame = vdif.Frame(frame.id, frame.time, bits=8, data=frame.payload.data[0,:].astype(np.int8), sample_rate=196e6)
                vFrame.write_raw_frame(fh)
                
        # Read it back in
        with open(testFile, 'rb') as fh:
            for tFrame in frames:
                vFrame = vrdr.read_frame(fh, sample_rate=196e6)
                self.assertAlmostEqual(vFrame.time, tFrame.time, 6)
                np.testing.assert_allclose((vFrame.payload.data*256-1)/2, tFrame.payload.data[0,:].astype(np.int8), atol=1e-6)
                
    def test_vdif_complex(self):
        """Test writing complex data to VIDF format."""

        # Setup the file names
        testFile = os.path.join(self.testPath, 'tbn-test-W.fits')

        # Get some data
        frames = self._get_tbn()
        
        # Write the data
        with open(testFile, 'wb') as fh:
            for frame in frames:
                stand, pol = frame.id
                if pol == 1:
                    continue
                ## We need an integer number of frame per second so adjust the sample rate
                vFrame = vdif.Frame(stand, frame.time, bits=8, data=frame.payload.data, sample_rate=102.4e3)
                vFrame.write_raw_frame(fh)
                
        # Read it back in
        with open(testFile, 'rb') as fh:
            for tFrame in frames:
                stand, pol = tFrame.id
                if pol == 1:
                    continue
                vFrame = vrdr.read_frame(fh, sample_rate=102.4e3)
                self.assertAlmostEqual(vFrame.time, tFrame.time, 6)
                vData = vFrame.payload.data
                vData.real = (vData.real*256-1).astype(np.int8)//2
                vData.imag = (vData.imag*256-1).astype(np.int8)//2
                np.testing.assert_allclose(vData, tFrame.payload.data, atol=1e-6)
                
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
