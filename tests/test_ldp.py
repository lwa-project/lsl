"""
Unit test for the lsl.reader.ldp module.
"""

import os
import unittest
import tempfile
import shutil
import numpy as np

from lsl.reader import ldp
from lsl.reader import errors
from lsl.reader.utils import SplitFileWrapper


__version__  = "0.3"
__author__    = "Jayce Dowell"


drxFile = os.path.join(os.path.dirname(__file__), 'data', 'drx-test.dat')
drspecFile = os.path.join(os.path.dirname(__file__), 'data', 'drspec-test.dat')

tbfFile = os.path.join(os.path.dirname(__file__), 'data', 'tbf-test.dat')


class ldp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    def setUp(self):
        """Create the temporary file directory."""

        self.testPath = tempfile.mkdtemp(prefix='test-ldp-', suffix='.tmp')
        
    ### DRX ###
    
    def test_ldp_drx(self):
        """Test the LDP interface for a DRX file."""
        
        f = ldp.DRXFile(drxFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 32)
        self.assertEqual(f.get_info('nbeampol'), 4)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 32)
        self.assertEqual(f.nbeampol, 4)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(0.005)
        
        # Reset
        f.reset()
        
        # Offset and read a chunk - short
        tSkip = f.offset(0.0001)
        tInt, tStart, data = f.read(0.005)
        
        # Reset
        f.reset()
        
        # Read a chunk - long
        tInt, tStart, data = f.read(1.00)
        
        # Reset
        f.reset()
        
        # Estimate levels
        f.estimate_levels()
        
        # Close it out
        f.close()
        
        # 'with' statement support
        with ldp.DRXFile(drxFile) as f:
            ## File info
            self.assertEqual(f.get_info("sample_rate"), 19.6e6)
            self.assertEqual(f.get_info("data_bits"), 4)
            self.assertEqual(f.get_info('nframe'), 32)
            self.assertEqual(f.get_info('nbeampol'), 4)
            
            self.assertEqual(f.sample_rate, 19.6e6)
            self.assertEqual(f.data_bits, 4)
            self.assertEqual(f.nframe, 32)
            self.assertEqual(f.nbeampol, 4)
            
            ## Read a frame
            frame = f.read_frame()
            
            ## Get the remaining frame count
            self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
            self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
            
        # generator support
        f = ldp.DRXFile(drxFile)
        i = 0
        for (tInt2, tStart2, data2) in f.read_sequence(1.0):
            self.assertEqual(tInt, tInt2)
            self.assertEqual(tStart, tStart2)
            self.assertEqual(data.shape, data2.shape)
            i += 1
        self.assertEqual(i, 1)
        f.close()
        
        # both at the same time
        with ldp.DRXFile(drxFile) as f:
            i = 0
            for (tInt2, tStart2, data2) in f.read_sequence(1.0):
                self.assertEqual(tInt, tInt2)
                self.assertEqual(tStart, tStart2)
                self.assertEqual(data.shape, data2.shape)
                i += 1
            self.assertEqual(i, 1)
            
    def test_ldp_drx_nocheck(self):
        """Test the LDP interface for a DRX file."""
        
        f = ldp.DRXFile(drxFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 32)
        self.assertEqual(f.get_info('nbeampol'), 4)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 32)
        self.assertEqual(f.nbeampol, 4)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(0.005)
        
        # Reset
        f.reset()
        
        # Read a chunk - long
        tInt, tStart, data = f.read(1.00)
        
        # Close it out
        f.close()
        
    def test_ldp_drx_ci8(self):
        """Test the LDP interface for a DRX file, ci8 style."""
        
        f = ldp.DRXFile(drxFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 32)
        self.assertEqual(f.get_info('nbeampol'), 4)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 32)
        self.assertEqual(f.nbeampol, 4)
        
        # Read a frame
        frame = f.read_frame(return_ci8=True)
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(0.005, return_ci8=True)
        self.assertEqual(len(data.shape), 2)
        
        # Reset
        f.reset()
        
        # Offset and read a chunk - short
        tSkip = f.offset(0.0001)
        tInt, tStart, data = f.read(0.005, return_ci8=True)
        self.assertEqual(len(data.shape), 2)
        
        # Reset
        f.reset()
        
        # Read a chunk - long
        tInt, tStart, data = f.read(1.00, return_ci8=True)
        self.assertEqual(len(data.shape), 2)
        
        # Go back and try it again without ci8 support
        f.reset()
        _, _, data2 = f.read(1.00, return_ci8=False)
        data = data['re'] + 1j*data['im']
        np.testing.assert_equal(data, data2)
        
        # Close it out
        f.close()
        
    ### DR Spectrometer ###
    
    def test_ldp_drspec(self):
        """Test the LDP interface for a DR Spectrometer file."""
        
        f = ldp.DRSpecFile(drspecFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 32)
        self.assertEqual(f.get_info('nframe'), 7)
        self.assertEqual(f.get_info('nbeampol'), 4)
        self.assertEqual(f.get_info('nproduct'), 2)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 32)
        self.assertEqual(f.nframe, 7)
        self.assertEqual(f.nbeampol, 4)
        self.assertEqual(f.nproduct, 2)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(0.01)
        
        # Reset
        f.reset()
        
        # Offset and read a chunk - short
        skip = f.offset(0.01)
        tInt, tStart, data = f.read(0.01)

        # Reset
        f.reset()
        
        # Read a chunk - long
        tInt, tStart, data = f.read(5.00)
        
        # Close it out
        f.close()
        
        # 'with' statement support
        with ldp.DRSpecFile(drspecFile) as f:
            ## File info
            self.assertEqual(f.get_info("sample_rate"), 19.6e6)
            self.assertEqual(f.get_info("data_bits"), 32)
            self.assertEqual(f.get_info('nframe'), 7)
            self.assertEqual(f.get_info('nbeampol'), 4)
            self.assertEqual(f.get_info('nproduct'), 2)
            
            self.assertEqual(f.sample_rate, 19.6e6)
            self.assertEqual(f.data_bits, 32)
            self.assertEqual(f.nframe, 7)
            self.assertEqual(f.nbeampol, 4)
            self.assertEqual(f.nproduct, 2)
            
            ## Read a frame
            frame = f.read_frame()
            
            ## Get the remaining frame count
            self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
            self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
            
        # generator support
        f = ldp.DRSpecFile(drspecFile)
        i = 0
        for (tInt2, tStart2, data2) in f.read_sequence(5.0):
            self.assertEqual(tInt, tInt2)
            self.assertEqual(tStart, tStart2)
            self.assertEqual(data.shape, data2.shape)
            i += 1
        self.assertEqual(i, 1)
        f.close()
        
        # both at the same time
        with ldp.DRSpecFile(drspecFile) as f:
            i = 0
            for (tInt2, tStart2, data2) in f.read_sequence(5.0):
                self.assertEqual(tInt, tInt2)
                self.assertEqual(tStart, tStart2)
                self.assertEqual(data.shape, data2.shape)
                i += 1
            self.assertEqual(i, 1)
            
    def test_ldp_drspec_nocheck(self):
        """Test the LDP interface for a DR Spectrometer file."""
        
        f = ldp.DRSpecFile(drspecFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 32)
        self.assertEqual(f.get_info('nframe'), 7)
        self.assertEqual(f.get_info('nbeampol'), 4)
        self.assertEqual(f.get_info('nproduct'), 2)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 32)
        self.assertEqual(f.nframe, 7)
        self.assertEqual(f.nbeampol, 4)
        self.assertEqual(f.nproduct, 2)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(0.01)
        
        # Reset
        f.reset()
        
        # Read a chunk - long
        tInt, tStart, data = f.read(5.00)
        
        # Close it out
        f.close()
        
    ### File Type Discovery ###
    
    def test_ldp_discover_drx(self):
        """Test the LDP LWA1DataFile function of DRX."""
        # DRX
        f = ldp.LWADataFile(drxFile)
        self.assertEqual(type(f), ldp.DRXFile)
        
    def test_ldp_discover_drspec(self):
        """Test the LDP LWA1DataFile function of DR Spectrometer."""
        # DR Spectrometer
        f = ldp.LWADataFile(drspecFile)
        self.assertEqual(type(f), ldp.DRSpecFile)
        
    def tearDown(self):
        """Cleanup"""
        for handler in list(ldp._open_ldp_files.handlers):
            handler.close()
        shutil.rmtree(self.testPath, ignore_errors=True)


class ldp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader.ldp 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(ldp_tests)) 


if __name__ == '__main__':
    unittest.main()
