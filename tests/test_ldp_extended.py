"""
Extended unit test for the lsl.reader.ldp module.
"""

import os
import unittest
import tempfile
import shutil
import subprocess

from lsl.reader import ldp
from lsl.reader import errors
from lsl.reader.utils import SplitFileWrapper


__version__  = "0.1"
__author__    = "Jayce Dowell"


run_extended_tests = False
if os.getenv('GITHUB_ACTIONS', None) is not None:
    run_extended_tests = True


_TBN_URL = 'https://lda10g.alliance.unm.edu/tutorial/Meteors/056761_000099453'
_DRX_URL = 'https://lda10g.alliance.unm.edu/tutorial/UnknownPulsar/056227_000024985_DRX.dat'
_SPC_URL = 'https://lda10g.alliance.unm.edu/tutorial/B0329+54/056770_000044687'

tbnFile = os.path.join(os.path.dirname(__file__), 'data', 'tbn-extended.dat')
drxFile = os.path.join(os.path.dirname(__file__), 'data', 'drx-extended.dat')
drspecFile = os.path.join(os.path.dirname(__file__), 'data', 'drspec-extended.dat')


@unittest.skipUnless(run_extended_tests, "requires appropriate environment variable to be set")
class extended_ldp_tests(unittest.TestCase):
    """An extended unittest.TestCase collection of unit tests for the lsl.reader.ldp
    module."""
    
    def setUp(self):
        """Download the files."""
        
        for filename,url in zip((tbnFile, drxFile, drspecFile), (_TBN_URL,_DRX_URL,_SPC_URL)):
            if not os.path.exists(filename):
                subprocess.check_call(['curl', url,
                                   '--silent',
                                   '--range', '0-%i' % (250*1024*1024),
                                   '-o', filename])
                
    def test_tbn_estimate(self):
        """Test estimating power levels in a TBN file."""
        
        idf = ldp.TBNFile(tbnFile)
        offset = idf.offset(0.1)       # This file needs a skip at the beginning
        levels = idf.estimate_levels()
        self.assertEqual(len(levels), 520)
        
    def test_tbn_read(self):
        """Test more involved reading from a TBN file."""
        
        idf = ldp.TBNFile(tbnFile)
        offset = idf.offset(0.1)       # This file needs a skip at the beginning
        
        # Read some
        for i in range(21):
            tInt, tStart, data = idf.read(0.1)
            
        idf.close()
        
    def test_tbn_offset(self):
        """Test offsetting inside a TBN file."""
        
        idf = ldp.TBNFile(tbnFile)
        offset = idf.offset(0.1)       # This file needs a skip at the beginning
        
        # Jump forwards
        fileStart = idf.start_time
        offset = idf.offset(0.1)
        
        # Read
        tInt, tStart, data = idf.read(0.11)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        # Jump forwards
        fileStart = tStart + tInt
        offset = idf.offset(0.1)
        
        # Read
        tInt, tStart, data = idf.read(0.12)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        # Jump backwards
        fileStart = tStart + tInt
        offset = idf.offset(-0.15)
        
        # Read
        tInt, tStart, data = idf.read(0.1)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        idf.close()
        
    def test_drx_estimate(self):
        """Test estimating power levels in a DRX file."""
        
        idf = ldp.DRXFile(drxFile)
        levels = idf.estimate_levels()
        self.assertEqual(len(levels), 4)
        
    def test_drx_read(self):
        """Test more involved reading from a DRX file."""
        
        idf = ldp.DRXFile(drxFile)
        
        # Read some
        for i in range(21):
            tInt, tStart, data = idf.read(0.1)
            
        idf.close()
        
    def test_drx_offset(self):
        """Test offsetting inside a DRX file."""
        
        idf = ldp.DRXFile(drxFile)
        
        # Jump forwards
        fileStart = idf.start_time
        offset = idf.offset(0.1)
        
        # Read
        tInt, tStart, data = idf.read(0.09)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        # Jump forwards
        fileStart = tStart + tInt
        offset = idf.offset(0.1)
        
        # Read
        tInt, tStart, data = idf.read(0.11)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        # Jump backwards
        fileStart = tStart + tInt
        offset = idf.offset(-0.15)
        
        # Read
        tInt, tStart, data = idf.read(0.12)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        idf.close()
        
    def test_drspec_read(self):
        """Test more involved reading from a DR Spectrometer file."""
        
        idf = ldp.DRSpecFile(drspecFile)
        
        # Read some
        for i in range(21):
            tInt, tStart, data = idf.read(0.1)
            
        idf.close()
        
    def test_drspec_offset(self):
        """Test offsetting inside a DR Spectrometer file."""
        
        idf = ldp.DRSpecFile(drspecFile)
        
        # Jump forwards
        fileStart = idf.start_time
        offset = idf.offset(0.1)
        
        # Read
        tInt, tStart, data = idf.read(0.1)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        # Jump forwards
        fileStart = tStart + tInt
        offset = idf.offset(0.1)
        
        # Read
        tInt, tStart, data = idf.read(0.1)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        # Jump backwards
        fileStart = tStart + tInt
        offset = idf.offset(-0.15)
        
        # Read
        tInt, tStart, data = idf.read(0.1)
        self.assertAlmostEqual(tStart, fileStart+offset, 9)
        
        idf.close()
        
    def testDown(self):
        """Cleanup"""
        
        for handler in list(ldp._open_ldp_files.handlers):
            handler.close()


class extended_ldp_test_suite(unittest.TestSuite):
    """An extended unittest.TestSuite class which contains all of the lsl.reader.ldp 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(extended_ldp_tests)) 


if __name__ == '__main__':
    unittest.main()
