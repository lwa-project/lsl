"""
Unit test for regressions in the lsl.common.busy module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import sys
import time
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
from lsl.common import busy

__version__  = "0.1"
__author__    = "Jayce Dowell"


class busy_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the regressions in LSL."""
    
    def test_default(self):
        """Test the busy indicator, default options."""
        
        bi = busy.BusyIndicator()
        
        sys.stdout = StringIO()
        
        bi.start()
        time.sleep(1)
        bi.stop()
        
        sys.stdout = sys.__stdout__


class busy_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.busy unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(busy_tests)) 


if __name__ == '__main__':
    unittest.main()
