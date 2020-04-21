"""
Unit test for the lsl.misc.dedispersion module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import time
import warnings
import unittest
import numpy

from lsl.misc import dedispersion


__version__  = "0.1"
__author__    = "Jayce Dowell"

class dedispersion_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.dedispersion
    module."""
    
    def test_delay(self):
        """Test calculating the dispersion delay"""
        
        d0 = dedispersion.delay(74e6, 12.455)
        d1 = dedispersion.delay([74e6,], 12.455)
        d2 = dedispersion.delay([74e6, 80e6], 12.455)
        d3 = dedispersion.delay(numpy.array([74e6, 80e6]), 12.455)
        self.assertAlmostEqual(d0, d1, 6)
        self.assertAlmostEqual(d2[0], d3[0], 6)
        
    def test_incoherent(self):
        """Test incoherent dedispersion"""
        
        t = numpy.arange(3000)*0.1
        freq = numpy.linspace(-9.8e6, 9.8e6, 1024) + 74e6
        data = numpy.random.randn(t.size, freq.size)
        data[100,:] += 10
        data = dedispersion.incoherent(freq, data, t[1]-t[0], -12.455)
        data2 = dedispersion.incoherent(freq, data, t[1]-t[0],  2.455)
        data = dedispersion.incoherent(freq, data, t[1]-t[0],  12.455)
        self.assertAlmostEqual(data[100,:].mean(), 10.0, 1)
        self.assertTrue(data2[100,:].mean() < 10-1)


class dedispersion_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.dedispersion 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(dedispersion_tests)) 


if __name__ == '__main__':
    unittest.main()
