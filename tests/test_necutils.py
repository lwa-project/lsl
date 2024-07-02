"""
Unit test for the lsl.sim.necutils module.
"""

import unittest
import tempfile
import shutil
import os
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from lsl.sim import necutils
from lsl.common.data_access import DataAccess


__version__   = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class necutils_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.necutils
    module."""
    
    def setUp(self):
        """Setup unit tests."""
        
        # get a reference to the input test file
        self.nec_name = os.path.join(os.path.dirname(__file__), 'data', 'bigblade_imp.out')
    
    def test_NECImpedance_init(self):
        """Test necutils.NECImpedance constructor method."""
        
        imp = necutils.NECImpedance(self.nec_name)
        
    def test_open_and_get_nec_freq(self):
        """Test necutils.open_and_get_nec_freq() function."""
        
        (fh, freq) = necutils.open_and_get_nec_freq(self.nec_name)   
        fh.close()
        
    def test_change_nec_freq(self):
        """Test the necutils.change_nec_freq() function."""
        
        testPath = tempfile.mkdtemp(prefix='test-necutils-', suffix='.tmp')
        for freq in (25.6, 38.7, 54.6, 75.02):
            with DataAccess.open('antenna/nec/lwa1_xep_1.nec', 'r') as fh:
                with open(os.path.join(testPath, 'lwa1_xep_1.nec'), 'w') as oh:
                    oh.write(fh.read())
                    
            filename = os.path.join(testPath, 'lwa1_xep_1.nec')
            necutils.change_nec_freq(filename, freq)
        
            with open(filename, 'r') as fh:
                found = False
                for line in fh:
                    if line[:2] == 'FR':
                        found = True
                        self.assertTrue(line.find('%.2f' % freq) != -1)
                self.assertTrue(found)
            os.unlink(filename)
        shutil.rmtree(testPath)
        
    def test_calculate_ime(self):
        """Test necutils.calculate_ime() function."""
        
        (freqs, ime) = necutils.calculate_ime(self.nec_name)
    
    def test_NECPattern_init(self):
        """Test necutils.NECPattern constructor method."""
        
        pat = necutils.NECPattern(self.nec_name, 5.0)
        self.assertRaises(ValueError, necutils.NECPattern, self.nec_name, 0.0, False)

    
class necutils_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lwa_user.necutils
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(necutils_tests))        
        
        
if __name__ == '__main__':
    unittest.main()
