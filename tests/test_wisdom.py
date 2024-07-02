"""
Unit test for the lsl.misc.wisdom module.
"""

import os
import time
import warnings
import unittest

from lsl.misc import wisdom
import lsl.testing


__version__  = "0.1"
__author__    = "Jayce Dowell"

class wisdom_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.statistics.robust
    module."""
    
    def test_show(self):
        """Test wisdom.show()"""
        
        with lsl.testing.SilentVerbose():
            wisdom.show()


class wisdom_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.wisdom 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(wisdom_tests)) 


if __name__ == '__main__':
    unittest.main()
