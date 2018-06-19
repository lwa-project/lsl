# -*- coding: utf-8 -*-

"""Unit test for lsl.common.paths module."""

import os
import unittest

from lsl.common import paths


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class paths_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.paths tests
    module."""

    def test_module_path(self):
        """Test the paths.module variable."""
        
        modPath = paths.moduleBuild

        astroFile = os.path.join(modPath, 'astro.py')
        self.assertTrue(os.path.exists(astroFile))

        tbwFile = os.path.join(modPath, 'reader', 'tbw.py')
        self.assertTrue(os.path.exists(tbwFile))

    def test_data_path(self):
        """Test the paths.data variable."""
        
        dataPath = paths.dataBuild

        ssmif = os.path.join(dataPath, 'lwa1-ssmif.txt')
        self.assertTrue(os.path.exists(ssmif))

        timeFile = os.path.join(dataPath, 'astro', 'tai-utc.dat')
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
