"""
Unit test for the lsl.sim.necutil module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import unittest
import os
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from lsl.sim import necutil
from lsl.common.paths import DATA_BUILD


__version__   = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class necutil_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.necutil
    module."""
    
    def setUp(self):
        """Setup unit tests."""
        
        # get a reference to the input test file
        self.nec_name = os.path.join(DATA_BUILD, 'tests', 'bigblade_imp.out')
    
    def test_NECImpedance_init(self):
        """Test necutil.NECImpedance constructor method."""
        
        imp = necutil.NECImpedance(self.nec_name)
        
    def test_open_and_get_nec_freq(self):
        """Test necutil.open_and_get_nec_freq() function."""
        
        (fh, freq) = necutil.open_and_get_nec_freq(self.nec_name)   
        fh.close()
        
    def test_calculate_ime(self):
        """Test necutil.calculate_ime() function."""
        
        (freqs, ime) = necutil.calculate_ime(self.nec_name)
    
    def test_NECPattern_init(self):
        """Test necutil.NECPattern constructor method."""
        
        pat = necutil.NECPattern(self.nec_name, 5.0)
        self.assertRaises(ValueError, necutil.NECPattern, self.nec_name, 0.0, False)

    
class necutil_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lwa_user.necutil
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(necutil_tests))        
        
        
if __name__ == '__main__':
    unittest.main()

