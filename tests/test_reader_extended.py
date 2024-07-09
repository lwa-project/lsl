"""
Extended unit test for the lsl.reader modules
"""

import os
import glob
import unittest
import subprocess

from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import drx
from lsl.reader import vdif
from lsl.reader import drspec
from lsl.reader import errors
from lsl.reader.base import FrameTimestamp


__version__  = "0.1"
__author__    = "Jayce Dowell"


run_extended_tests = False
if os.getenv('GITHUB_ACTIONS', None) is not None:
    run_extended_tests = True


_VDIF_URL = 'https://fornax.phys.unm.edu/lwa/data/eLWA_test_small_raw.tar.gz'


@unittest.skipUnless(run_extended_tests, "requires appropriate environment variable to be set")
class extended_reader_tests(unittest.TestCase):
    """An extended unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    def setUp(self):
        """Create the temporary file directory."""
        
        if not os.path.exists('eLWA_test_raw.tar.gz'):
            subprocess.check_call(['curl', _VDIF_URL,
                                   '--silent',
                                   '-o', 'eLWA_test_raw.tar.gz'])
            subprocess.check_call(['tar', 'xzf', 'eLWA_test_raw.tar.gz'])
            
    def test_vdif_guppi_header(self):
        """Test reading a VDIF file with a GUPPI header."""
        
        filename = glob.glob('*.vdif')[0]
        
        # Open the file
        fh = open(filename, 'rb')
        
        # Make sure it is there
        self.assertTrue(vdif.has_guppi_header(fh))
        
        # Read the GUPPI header
        header = vdif.read_guppi_header(fh)
        self.assertEqual(header['NBITS'], 4)
        self.assertEqual(header['PKTSIZE'], 5000)
        self.assertAlmostEqual(header['OBSBW'], 8.0e6, 6)
        
        # And now read some frames to make sure that everything is ok
        frame0 = vdif.read_frame(fh)
        frame1 = vdif.read_frame(fh)
        self.assertEqual(frame0.payload.data.size, header['PKTSIZE']*8//header['NBITS'])
        
        fh.close()
        
    def test_vdif_frames_per_second(self):
        """Test reading a VDIF file with a GUPPI header."""
        
        filename = glob.glob('*.vdif')[0]
        
        # Open the file
        fh = open(filename, 'rb')
        
        # Read the GUPPI header and a frame
        header = vdif.read_guppi_header(fh)
        frame0 = vdif.read_frame(fh)
        
        # Find the number of frames per second
        fps = vdif.get_frames_per_second(fh)
        
        # Does it make sense?
        bw = header['OBSBW']
        bw *= 1 if frame0.header.is_complex else 2
        calculated_fps = int(bw / frame0.payload.data.size)
        self.assertEqual(fps, calculated_fps)
        
        # Find the sample rate
        sample_rate = vdif.get_sample_rate(fh)
        self.assertAlmostEqual(sample_rate, bw, 6)
        
        fh.close()


class reader_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(reader_tests)) 


if __name__ == '__main__':
    unittest.main()
