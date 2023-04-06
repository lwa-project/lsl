"""
Unit test for the ADP portion of the lsl.reader.ldp module.
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


__version__  = "0.2"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(DATA_BUILD, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(DATA_BUILD, 'tests', 'tbn-test.dat')
drxFile = os.path.join(DATA_BUILD, 'tests', 'drx-test.dat')
drspecFile = os.path.join(DATA_BUILD, 'tests', 'drspec-test.dat')

tbfFile = os.path.join(DATA_BUILD, 'tests', 'tbf-test.dat')
corFile = os.path.join(DATA_BUILD, 'tests', 'cor-test.dat')

ovroFile = os.path.join(DATA_BUILD, 'tests', 'ovro-test.dat')


class ldp_adp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    ### TBF ###
    
    def test_ldp_tbf(self):
        """Test the LDP interface for a TBF file."""
        
        f = ldp.TBFFile(tbfFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 25e3)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 5)
        
        self.assertEqual(f.sample_rate, 25e3)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 5)
        
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
        
    def test_ldp_tbf_nocheck(self):
        """Test the LDP interface for a TBF file."""
        
        f = ldp.TBFFile(tbfFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 25e3)
        self.assertEqual(f.get_info("data_bits"), 4)
        self.assertEqual(f.get_info('nframe'), 5)
        
        self.assertEqual(f.sample_rate, 25e3)
        self.assertEqual(f.data_bits, 4)
        self.assertEqual(f.nframe, 5)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Close it out
        f.close()
        
    ### COR ###
    
    def test_ldp_cor(self):
        """Test the LDP interface for a COR file."""
        
        f = ldp.CORFile(corFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info('nchan'), 72)
        self.assertEqual(f.get_info('nbaseline'), 32896)
        self.assertEqual(f.get_info('nframe'), 65)
        
        self.assertEqual(f.nchan, 72)
        self.assertEqual(f.nbaseline, 32896)
        self.assertEqual(f.nframe, 65)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(5)
        
        # Close it out
        f.close()
        
    def test_ldp_cor_nocheck(self):
        """Test the LDP interface for a COR file."""
        
        f = ldp.CORFile(corFile)
        
        # File info
        self.assertEqual(f.get_info('nchan'), 72)
        self.assertEqual(f.get_info('nbaseline'), 32896)
        self.assertEqual(f.get_info('nframe'), 65)
        
        self.assertEqual(f.nchan, 72)
        self.assertEqual(f.nbaseline, 32896)
        self.assertEqual(f.nframe, 65)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
        # Reset
        f.reset()
        
        # Read a chunk - short
        tInt, tStart, data = f.read(5)
        
        # Close it out
        f.close()
        
    ### File Type Discovery ###
    
    def test_ldp_discover_tbw(self):
        """Test the LDP LWA1DataFile function of TBW."""
        # TBW
        self.assertRaises(RuntimeError, ldp.LWASVDataFile, tbwFile)
        
    def test_ldp_discover_tbn(self):
        """Test the LDP LWASVDataFile function of TBN."""
        # TBN
        f = ldp.LWASVDataFile(tbnFile)
        self.assertEqual(type(f), ldp.TBNFile)
        
    def test_ldp_discover_drx(self):
        """Test the LDP LWASVDataFile function of DRX."""
        # DRX
        f = ldp.LWASVDataFile(drxFile)
        self.assertEqual(type(f), ldp.DRXFile)
        
    def test_ldp_discover_drspec(self):
        """Test the LDP LWASVDataFile function of DR Spectrometer."""
        # DR Spectrometer
        f = ldp.LWASVDataFile(drspecFile)
        self.assertEqual(type(f), ldp.DRSpecFile)
        
    def test_ldp_discover_tbf(self):
        """Test the LDP LWASVDataFile function of TBF."""
        # TBF
        f = ldp.LWASVDataFile(tbfFile)
        self.assertEqual(type(f), ldp.TBFFile)
        
    def test_ldp_discover_cor(self):
        """Test the LDP LWASVDataFile function of COR."""
        # TBF
        f = ldp.LWASVDataFile(corFile)
        self.assertEqual(type(f), ldp.CORFile)
        
    def test_ldp_discover_ovro(self):
        """Test the LDP LWASVDataFile function of OVRO-LWA triggered voltage buffer dump."""
        # TBF
        self.assertRaises(RuntimeError, ldp.LWASVDataFile, ovroFile)
        
    def test_ldp_discover_all_tbw(self):
        """Test the LDP LWADataFile function of TBW."""
        # TBW
        f = ldp.LWADataFile(tbwFile)
        self.assertEqual(type(f), ldp.TBWFile)
        
    def test_ldp_discover_all_tbn(self):
        """Test the LDP LWADataFile function of TBN."""
        # TBN
        f = ldp.LWADataFile(tbnFile)
        self.assertEqual(type(f), ldp.TBNFile)
        
    def test_ldp_discover_all_drx(self):
        """Test the LDP LWADataFile function of DRX."""
        # DRX
        f = ldp.LWADataFile(drxFile)
        self.assertEqual(type(f), ldp.DRXFile)
        
    def test_ldp_discover_all_drspec(self):
        """Test the LDP LWADataFile function of DR Spectrometer."""
        # DR Spectrometer
        f = ldp.LWADataFile(drspecFile)
        self.assertEqual(type(f), ldp.DRSpecFile)
        
    def test_ldp_discover_all_tbf(self):
        """Test the LDP LWADataFile function of TBF."""
        # TBF
        f = ldp.LWADataFile(tbfFile)
        self.assertEqual(type(f), ldp.TBFFile)
        
    def test_ldp_discover_all_cor(self):
        """Test the LDP LWADataFile function of COR."""
        # TBF
        f = ldp.LWADataFile(corFile)
        self.assertEqual(type(f), ldp.CORFile)
        
    
    ### SplitFileWrapper ###
    
    def test_ldp_splitfilewrapper_discover(self):
        """Test the LDP interface for type discover with a SplitFileWrapper."""
        
        w = SplitFileWrapper([tbnFile,])
        f = ldp.LWASVDataFile(fh=w)
        
        # File info
        self.assertTrue(isinstance(f, ldp.TBNFile))
        self.assertEqual(f.get_info("sample_rate"), 100e3)
        self.assertEqual(f.get_info("data_bits"), 8)
        self.assertEqual(f.get_info('nframe'), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframe, 29)
        
        # Read a frame
        frame = f.read_frame()
        
        f.close()
        w.close()
        
    def test_ldp_splitfilewrapper_mixed(self):
        """Test the LDP interface for a SplitFileWrapper in a contorted way."""
        
        w = SplitFileWrapper([drxFile, tbfFile], sort=False)
        f = ldp.LWASVDataFile(fh=w)
        
        # File info
        self.assertTrue(isinstance(f, ldp.DRXFile))
        
        # Read some
        frames = []
        while True:
            try:
                frame = f.read_frame()
                frames.append(frame)
            except errors.SyncError:
                continue
            except errors.EOFError:
                break
        self.assertEqual(len(frames), 32+1)  # There is data at the end of the "file" now      
        
    def test_ldp_splitfilewrapper_mixed2(self):
        """Test the LDP interface for a SplitFileWrapper in a contorted way."""
        
        w = SplitFileWrapper([drxFile, tbfFile], sort=False)
        f = ldp.TBFFile(fh=w)
        
        # Read some
        frames = []
        while True:
            try:
                frame = f.read_frame()
                frames.append(frame)
            except errors.SyncError:
                continue
            except errors.EOFError:
                break
        self.assertEqual(len(frames), 5-1)  # We loose part of the first frame to DRX   
        
    def tearDown(self):
        """Cleanup"""
        for handler in list(ldp._open_ldp_files.handlers):
            handler.close()


class ldp_adp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader.ldp 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(ldp_adp_tests)) 


if __name__ == '__main__':
    unittest.main()
