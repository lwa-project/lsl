# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import division

"""Unit test for the lsl.writer.tsfits modules."""

import os
import time
import unittest
import tempfile
import numpy

from lsl.common.paths import DATA_BUILD
from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import errors
from lsl.writer import vdif


__revision__ = "$Rev$"
__version__  = "0.2"
__author__    = "Jayce Dowell"

tbwFile = os.path.join(DATA_BUILD, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(DATA_BUILD, 'tests', 'tbn-test.dat')


class vdif_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.vdif
    module."""

    testPath = None

    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-vdif-', suffix='.tmp')

    def __getTBW(self):
        """Private function to load in the test TBW data and get the frames."""

        fh = open(tbwFile, 'rb')

        # Frames 1 through 8
        frames = []
        for i in range(1,9):
            frames.append(tbw.read_frame(fh))

        fh.close()
        return frames

    def __getTBN(self, vanilla=False):
        """Private function to load in the test TBN data and get the frames.  If 
        the keyword 'vanilla' is set to True, gain, sample rate, and frequency meta-
        data are not added to the frames."""

        fh = open(tbnFile, 'rb')

        # Frames 1 through 8
        frames = []
        for i in range(1,9):
            frames.append(tbn.read_frame(fh))

        if not vanilla:
            # Set some values for the other meta-data
            for frame in frames:
                frame.sample_rate = 100000

        fh.close()
        return frames

    def test_vdif_real(self):
        """Test writing real data to VDIF format."""

        # Setup the file names
        testFile = os.path.join(self.testPath, 'tbw-test-W.fits')

        # Get some data
        frames = self.__getTBW()

        # Write the data
        fh = open(testFile, 'wb')
        for frame in frames:
            vFrame = vdif.Frame(frame.id, frame.get_time(), bits=12, data=frame.data.xy[0,:])
            vFrame.write_raw_frame(fh)
        fh.close()

        # Read it back in
        fh = open(testFile, 'rb')
        junk = numpy.fromfile(fh, dtype=numpy.uint8, count=16)
        frameSize = (junk[10]<<16) | (junk[9]<<8) | (junk[8])
        dataBits = (junk[15]>>2) & 15
        dataSize = (frameSize - 4) * 8
        fh.seek(0)
        for tFrame in frames:
            vFrame = numpy.fromfile(fh, dtype=numpy.uint8, count=frameSize*8)
            data = numpy.zeros(400, dtype=numpy.int16)
            j = 0
            for i in range(32,vFrame.shape[0],4):
                word = vFrame[i] | (vFrame[i+1]<<8) | (vFrame[i+2]<<16) | (vFrame[i+3]<<24)
                for k in range(32//dataBits):
                    data[j] = ((word>>(dataBits*k)) & (2**dataBits-1)) - 2**(dataBits-1)
                    j = j + 1
            
            for v,t in zip(data, tFrame.data.xy[0,:]):
                self.assertEqual(v, t)

        fh.close()

    def test_vdif_complex(self):
        """Test writing complex data to VIDF format."""

        # Setup the file names
        testFile = os.path.join(self.testPath, 'tbn-test-W.fits')

        # Get some data
        frames = self.__getTBN()
        
        # Write the data
        fh = open(testFile, 'wb')
        for frame in frames:
            stand, pol = frame.id
            if pol == 1:
                continue
            vFrame = vdif.Frame(stand, frame.get_time(), bits=8, data=frame.data.iq)
            vFrame.write_raw_frame(fh)
        fh.close()

        # Read it back in
        fh = open(testFile, 'rb')
        junk = numpy.fromfile(fh, dtype=numpy.uint8, count=16)
        frameSize = (junk[10]<<16) | (junk[9]<<8) | (junk[8])
        dataBits = (junk[15]>>2) & 15
        dataSize = (frameSize - 4) * 8
        fh.seek(0)
        
        for tFrame in frames[::2]:
            vFrame = numpy.fromfile(fh, dtype=numpy.uint8, count=frameSize*8)
            data = numpy.zeros(512, dtype=numpy.singlecomplex)
            j = 0
            for i in range(32,vFrame.shape[0],4):
                word = vFrame[i] | (vFrame[i+1]<<8) | (vFrame[i+2]<<16) | (vFrame[i+3]<<24)
                for k in range(32//dataBits):
                    if k % 2 == 0:
                        data.real[j] = ((word>>(dataBits*k)) & (2**dataBits-1)) - 2**(dataBits-1)
                    else:
                        data.imag[j] = ((word>>(dataBits*k)) & (2**dataBits-1)) - 2**(dataBits-1)
                        j = j + 1

            for v,t in zip(data, tFrame.data.iq):
                self.assertAlmostEqual(v, t, 6)
        
        fh.close()

    def tearDown(self):
        """Remove the test path directory and its contents"""

        tempFiles = os.listdir(self.testPath)
        for tempFile in tempFiles:
            os.unlink(os.path.join(self.testPath, tempFile))
        os.rmdir(self.testPath)
        self.testPath = None


class  vdif_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.writer.vdif units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(vdif_tests)) 


if __name__ == '__main__':
    unittest.main()
