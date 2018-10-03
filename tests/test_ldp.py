# -*- coding: utf-8 -*-

"""Unit test for lsl.reader.ldp module"""

import os
import unittest

from lsl.common.paths import DATA_BUILD
from lsl.reader import ldp
from lsl.reader import errors


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(DATA_BUILD, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(DATA_BUILD, 'tests', 'tbn-test.dat')
drxFile = os.path.join(DATA_BUILD, 'tests', 'drx-test.dat')
drspecFile = os.path.join(DATA_BUILD, 'tests', 'drspec-test.dat')

tbfFile = os.path.join(DATA_BUILD, 'tests', 'tbf-test.dat')


class ldp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### TBW ###
    
    def test_ldp_tbw(self):
        """Test the LDP interface for a TBW file."""
        
        f = ldp.TBWFile(tbwFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 196e6)
        self.assertEqual(f.get_info("data_bits"), 12)
        self.assertEqual(f.get_info("nframes"), 8)
        
        self.assertEqual(f.sample_rate, 196e6)
        self.assertEqual(f.data_bits, 12)
        self.assertEqual(f.nframes, 8)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
        # Reset
        f.reset()
        
        # Read more
        tInt, tStart, data = f.read()
        
        # Close it out
        f.close()
        
    def test_ldp_tbw_nocheck(self):
        """Test the LDP interface for a TBW file."""
        
        f = ldp.TBWFile(tbwFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 196e6)
        self.assertEqual(f.get_info("data_bits"), 12)
        self.assertEqual(f.get_info("nframes"), 8)
        
        self.assertEqual(f.sample_rate, 196e6)
        self.assertEqual(f.data_bits, 12)
        self.assertEqual(f.nframes, 8)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
        # Reset
        f.reset()
        
        # Read more
        tInt, tStart, data = f.read()
        
        # Close it out
        f.close()
        
    ### TBN ###
    
    def test_ldp_tbn(self):
        """Test the LDP interface for a TBN file."""
        
        f = ldp.TBNFile(tbnFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 100e3)
        self.assertEqual(f.get_info("data_bits"), 8)
        self.assertEqual(f.get_info("nframes"), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframes, 29)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
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
        
    def test_ldp_tbn_nocheck(self):
        """Test the LDP interface for a TBN file."""
        
        f = ldp.TBNFile(tbnFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 100e3)
        self.assertEqual(f.get_info("data_bits"), 8)
        self.assertEqual(f.get_info("nframes"), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframes, 29)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
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
        
    ### DRX ###
    
    def test_ldp_drx(self):
        """Test the LDP interface for a DRX file."""
        
        f = ldp.DRXFile(drxFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info("nframes"), 32)
        self.assertEqual(f.get_info("beampols"), 4)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframes, 32)
        self.assertEqual(f.beampols, 4)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
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
        
    def test_ldp_drx_nocheck(self):
        """Test the LDP interface for a DRX file."""
        
        f = ldp.DRXFile(drxFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info("nframes"), 32)
        self.assertEqual(f.get_info("beampols"), 4)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframes, 32)
        self.assertEqual(f.beampols, 4)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
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
        
    ### DR Spectrometer ###
    
    def test_ldp_drspec(self):
        """Test the LDP interface for a DR Spectrometer file."""
        
        f = ldp.DRSpecFile(drspecFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 32)
        self.assertEqual(f.get_info("nframes"), 7)
        self.assertEqual(f.get_info("beampols"), 4)
        self.assertEqual(f.get_info("nproducts"), 2)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 32)
        self.assertEqual(f.nframes, 7)
        self.assertEqual(f.beampols, 4)
        self.assertEqual(f.nproducts, 2)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
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
        
    def test_ldp_drspec_nocheck(self):
        """Test the LDP interface for a DR Spectrometer file."""
        
        f = ldp.DRSpecFile(drspecFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 19.6e6)
        self.assertEqual(f.get_info("data_bits"), 32)
        self.assertEqual(f.get_info("nframes"), 7)
        self.assertEqual(f.get_info("beampols"), 4)
        self.assertEqual(f.get_info("nproducts"), 2)
        
        self.assertEqual(f.sample_rate, 19.6e6)
        self.assertEqual(f.data_bits, 32)
        self.assertEqual(f.nframes, 7)
        self.assertEqual(f.beampols, 4)
        self.assertEqual(f.nproducts, 2)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info("nframes")-1)
        self.assertEqual(f.nframes_remaining, f.get_info("nframes")-1)
        
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
    
    def test_ldp_discover_tbw(self):
        """Test the LDP LWA1DataFile function of TBW."""
        # TBW
        f = ldp.LWA1DataFile(tbwFile)
        self.assertEqual(type(f), ldp.TBWFile)
        
    def test_ldp_discover_tbn(self):
        """Test the LDP LWA1DataFile function of TBN."""
        # TBN
        f = ldp.LWA1DataFile(tbnFile)
        self.assertEqual(type(f), ldp.TBNFile)
        
    def test_ldp_discover_drx(self):
        """Test the LDP LWA1DataFile function of DRX."""
        # DRX
        f = ldp.LWA1DataFile(drxFile)
        self.assertEqual(type(f), ldp.DRXFile)
        
    def test_ldp_discover_drspec(self):
        """Test the LDP LWA1DataFile function of DR Spectrometer."""
        # DR Spectrometer
        f = ldp.LWA1DataFile(drspecFile)
        self.assertEqual(type(f), ldp.DRSpecFile)
        
    def test_ldp_discover_tbf(self):
        """Test the LDP LWA1DataFile function of TBF."""
        # TBF
        self.assertRaises(RuntimeError, ldp.LWA1DataFile, tbfFile)
        
    def tearDown(self):
        """Cleanup"""
        for handler in list(ldp._open_ldp_files.handlers):
            handler.close()


class ldp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader.ldp 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(ldp_tests)) 


if __name__ == '__main__':
    unittest.main()
