"""
Unit test for the OVRO-LWA portion of the lsl.reader.ldp module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import unittest

from lsl.common.paths import DATA_BUILD
from lsl.reader import ldp
from lsl.reader import errors
from lsl.reader.utils import SplitFileWrapper


__version__  = "0.1"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(DATA_BUILD, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(DATA_BUILD, 'tests', 'tbn-test.dat')
drxFile = os.path.join(DATA_BUILD, 'tests', 'drx-test.dat')
drspecFile = os.path.join(DATA_BUILD, 'tests', 'drspec-test.dat')

tbfFile = os.path.join(DATA_BUILD, 'tests', 'tbf-test.dat')
corFile = os.path.join(DATA_BUILD, 'tests', 'cor-test.dat')

ovroFile = os.path.join(DATA_BUILD, 'tests', 'ovro-test.dat')


class ldp_ovro_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### Triggered voltage buffer dump ###
    
    def test_ldp_ovro(self):
        """Test the LDP interface for a triggered voltage buffer dump."""
        
        f = ldp.OVROFile(ovroFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 8192/196e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 10)
        
        self.assertEqual(f.sample_rate, 8192/196e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 10)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(0.1)
        
        # Close it out
        f.close()
        
    def test_ldp_ovro_nocheck(self):
        """Test the LDP interface for a triggered voltage buffer dump."""
        
        f = ldp.OVROFile(ovroFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 8192/196e6)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 10)
        
        self.assertEqual(f.sample_rate, 8192/196e6)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 10)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Close it out
        f.close()


class ldp_ovro_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader.ldp 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(ldp_ovro_tests)) 


if __name__ == '__main__':
    unittest.main()
