"""
Unit test for the lsl.reader.ldp module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import unittest
import tempfile
import shutil

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


class ldp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.reader
    modules."""
    
    def setUp(self):
        """Create the temporary file directory."""

        self.testPath = tempfile.mkdtemp(prefix='test-ldp-', suffix='.tmp')
        
    ### TBW ###
    
    def test_ldp_tbw(self):
        """Test the LDP interface for a TBW file."""
        
        f = ldp.TBWFile(tbwFile)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 196e6)
        self.assertEqual(f.get_info("data_bits"), 12)
        self.assertEqual(f.get_info('nframe'), 8)
        
        self.assertEqual(f.sample_rate, 196e6)
        self.assertEqual(f.data_bits, 12)
        self.assertEqual(f.nframe, 8)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
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
        self.assertEqual(f.get_info('nframe'), 8)
        
        self.assertEqual(f.sample_rate, 196e6)
        self.assertEqual(f.data_bits, 12)
        self.assertEqual(f.nframe, 8)
        
        # Read a frame
        frame = f.read_frame()
        
        # Get the remaining frame count
        self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
        self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
        
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
        self.assertEqual(f.get_info('nframe'), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframe, 29)
        
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
        tSkip = f.offset(0.005)
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
        with ldp.TBNFile(tbnFile) as f:
            ## File info
            self.assertEqual(f.get_info("sample_rate"), 100e3)
            self.assertEqual(f.get_info("data_bits"), 8)
            self.assertEqual(f.get_info('nframe'), 29)
            
            self.assertEqual(f.sample_rate, 100e3)
            self.assertEqual(f.data_bits, 8)
            self.assertEqual(f.nframe, 29)
            
            ## Read a frame
            frame = f.read_frame()
            
            ## Get the remaining frame count
            self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
            self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
            
        # generator support
        f = ldp.TBNFile(tbnFile)
        i = 0
        for (tInt2, tStart2, data2) in f.read_sequence(1.0):
            self.assertEqual(tInt, tInt2)
            self.assertEqual(tStart, tStart2)
            self.assertEqual(data.shape, data2.shape)
            i += 1
        self.assertEqual(i, 1)
        f.close()
        
        # both at the same time
        with ldp.TBNFile(tbnFile) as f:
            i = 0
            for (tInt2, tStart2, data2) in f.read_sequence(1.0):
                self.assertEqual(tInt, tInt2)
                self.assertEqual(tStart, tStart2)
                self.assertEqual(data.shape, data2.shape)
                i += 1
            self.assertEqual(i, 1)
            
    def test_ldp_tbn_nocheck(self):
        """Test the LDP interface for a TBN file."""
        
        f = ldp.TBNFile(tbnFile, ignore_timetag_errors=True)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 100e3)
        self.assertEqual(f.get_info("data_bits"), 8)
        self.assertEqual(f.get_info('nframe'), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframe, 29)
        
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
        
    ### SplitFileWrapper ###
    
    def testldp_splitfilewrapper(self):
        """Test the LDP interface for a SplitFileWrapper."""
        
        # Split up the TBN file into many many parts
        size = os.path.getsize(tbnFile)
        block = size // 16
        
        splitFiles = []
        f = open(tbnFile, 'rb')
        while size > 0:
            splitFiles.append(os.path.join(self.testPath, 'filepart%03i' % len(splitFiles)))
            with open(splitFiles[-1], 'wb') as w:
                w.write(f.read(block))
                size -= block
        f.close()
        
        w = SplitFileWrapper(splitFiles)
        f = ldp.TBNFile(fh=w)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 100e3)
        self.assertEqual(f.get_info("data_bits"), 8)
        self.assertEqual(f.get_info('nframe'), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframe, 29)
        
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
        w.close()
        
        # 'with' statement support
        with SplitFileWrapper(splitFiles) as w:
            with ldp.TBNFile(fh=w) as f:
                ## File info
                self.assertEqual(f.get_info("sample_rate"), 100e3)
                self.assertEqual(f.get_info("data_bits"), 8)
                self.assertEqual(f.get_info('nframe'), 29)
        
                self.assertEqual(f.sample_rate, 100e3)
                self.assertEqual(f.data_bits, 8)
                self.assertEqual(f.nframe, 29)
        
                ## Read a frame
                frame = f.read_frame()
        
                ## Get the remaining frame count
                self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
                self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
                
        # generator support
        w = SplitFileWrapper(splitFiles)
        f = ldp.TBNFile(fh=w)
        i = 0
        for (tInt2, tStart2, data2) in f.read_sequence(1.0):
            self.assertEqual(tInt, tInt2)
            self.assertEqual(tStart, tStart2)
            self.assertEqual(data.shape, data2.shape)
            i += 1
        self.assertEqual(i, 1)
        f.close()
        w.close()
        
        # both at the same time
        with SplitFileWrapper(splitFiles) as w:
            with ldp.TBNFile(fh=w) as f:
                i = 0
                for (tInt2, tStart2, data2) in f.read_sequence(1.0):
                    self.assertEqual(tInt, tInt2)
                    self.assertEqual(tStart, tStart2)
                    self.assertEqual(data.shape, data2.shape)
                    i += 1
                self.assertEqual(i, 1)
                
    def test_ldp_splitfilewrapper_single(self):
        """Test the LDP interface for a SplitFileWrapper with a single file."""
        
        w = SplitFileWrapper([tbnFile,])
        f = ldp.TBNFile(fh=w)
        
        # File info
        self.assertEqual(f.get_info("sample_rate"), 100e3)
        self.assertEqual(f.get_info("data_bits"), 8)
        self.assertEqual(f.get_info('nframe'), 29)
        
        self.assertEqual(f.sample_rate, 100e3)
        self.assertEqual(f.data_bits, 8)
        self.assertEqual(f.nframe, 29)
        
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
        w.close()
        
        # 'with' statement support
        with SplitFileWrapper([tbnFile,]) as w:
            with ldp.TBNFile(fh=w) as f:
                ## File info
                self.assertEqual(f.get_info("sample_rate"), 100e3)
                self.assertEqual(f.get_info("data_bits"), 8)
                self.assertEqual(f.get_info('nframe'), 29)
        
                self.assertEqual(f.sample_rate, 100e3)
                self.assertEqual(f.data_bits, 8)
                self.assertEqual(f.nframe, 29)
        
                ## Read a frame
                frame = f.read_frame()
        
                ## Get the remaining frame count
                self.assertEqual(f.get_remaining_frame_count(), f.get_info('nframe')-1)
                self.assertEqual(f.nframe_remaining, f.get_info('nframe')-1)
                
        # generator support
        w = SplitFileWrapper([tbnFile,])
        f = ldp.TBNFile(fh=w)
        i = 0
        for (tInt2, tStart2, data2) in f.read_sequence(1.0):
            self.assertEqual(tInt, tInt2)
            self.assertEqual(tStart, tStart2)
            self.assertEqual(data.shape, data2.shape)
            i += 1
        self.assertEqual(i, 1)
        f.close()
        w.close()
        
        # both at the same time
        with SplitFileWrapper([tbnFile,]) as w:
            with ldp.TBNFile(fh=w) as f:
                i = 0
                for (tInt2, tStart2, data2) in f.read_sequence(1.0):
                    self.assertEqual(tInt, tInt2)
                    self.assertEqual(tStart, tStart2)
                    self.assertEqual(data.shape, data2.shape)
                    i += 1
                self.assertEqual(i, 1)
                
    def test_ldp_splitfilewrapper_discover(self):
        """Test the LDP interface for type discover with a SplitFileWrapper."""
        
        w = SplitFileWrapper([tbnFile,])
        f = ldp.LWA1DataFile(fh=w)
        
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
        
        w = SplitFileWrapper([tbwFile, tbnFile], sort=False)
        f = ldp.LWA1DataFile(fh=w)
        
        # File info
        self.assertTrue(isinstance(f, ldp.TBWFile))
        
        # Read some
        frames = []
        while True:
            try:
                frame = f.read_frame()
                if frame.is_tbw:
                    frames.append(frame)
            except errors.SyncError:
                continue
            except errors.EOFError:
                break
        self.assertEqual(len(frames), 8+1)  # There is data at the end of the "file" now      
        
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
