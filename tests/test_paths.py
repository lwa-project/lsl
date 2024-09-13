"""
Unit test for the lsl.common.paths module.
"""

import os
import unittest

from lsl.common import paths


__version__  = "0.1"
__author__    = "Jayce Dowell"


class paths_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.paths tests
    module."""

    def test_module_path(self):
        """Test the paths.module variable."""
        
        modPath = paths.MODULE_BUILD

        astroFile = os.path.join(modPath, 'astro.py')
        self.assertTrue(os.path.exists(astroFile))

        tbwFile = os.path.join(modPath, 'reader', 'tbw.py')
        self.assertTrue(os.path.exists(tbwFile))

    def test_data_path(self):
        """Test the paths.data variable."""
        
        DATA_PATH = paths.DATA_BUILD

        ssmif = os.path.join(DATA_PATH, 'lwa1-ssmif.txt')
        self.assertTrue(os.path.exists(ssmif))

        timeFile = os.path.join(DATA_PATH, 'astro', 'Leap_Second.dat')
        self.assertTrue(os.path.exists(timeFile))


class paths_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.paths
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(paths_tests)) 


if __name__ == '__main__':
    unittest.main()
