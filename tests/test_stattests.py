# -*- coding: utf-8 -*-

"""Unit test for the lsl.statistics.tests module."""

import os
import time
import unittest
import numpy

from lsl.statistics import stattests


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class stattests_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.statistics.tests
    module."""

    def test_wald_wolfowitz(self):
        """Test the Wald-Wolfowitz (runs) test"""

        data = numpy.random.randn(1024)
        pValue = stattests.wald_wolfowitz(data)


class stattests_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.statistics.tests 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(stattests_tests)) 


if __name__ == '__main__':
    unittest.main()
