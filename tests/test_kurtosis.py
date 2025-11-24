"""
Unit test for the lsl.statistics.kurtosis module.
"""

import os
import time
import warnings
import unittest
import numpy as np

from lsl.statistics import kurtosis


__version__  = "0.1"
__author__    = "Jayce Dowell"

class kurtosis_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.statistics.kurtosis
    module."""
    
    def test_limits(self):
        """Test the limits returned by get_limits()"""
        
        N = 10
        M = 300
        sigma = 3
        
        lower, upper = kurtosis.get_limits(sigma, M, N=N)
        # Table 1 of Nita & Gary (2010, MNRAS 406, L60)
        self.assertAlmostEqual(lower, 0.76648, 5)
        self.assertAlmostEqual(upper, 1.28313, 5)
        
    def test_kurtosis(self):
        """Test that spectral kurtosis runs"""
        
        data = np.random.randn(300) + np.random.randn(300)*1j
        sk = kurtosis.spectral_fft(data)
        lower, upper = kurtosis.get_limits(6, data.size, N=1)
        self.assertTrue(sk > lower)
        self.assertTrue(sk < upper)
        
    def test_kurtosis_axis(self):
        """Test that spectral kurtosis runs along an axis"""
        
        data = np.random.randn(300, 128) + np.random.randn(300, 128)*1j
        sk0 = np.zeros(data.shape[1])
        for i in range(data.shape[1]):
            sk0[i] = kurtosis.spectral_fft(data[:,i])
        sk1 = kurtosis.spectral_fft(data, axis=0)
        lower, upper = kurtosis.get_limits(6, data.shape[0], N=1)
        
        self.assertEqual(sk0.size, sk1.size)
        for s0,s1 in zip(sk0, sk1):
            self.assertAlmostEqual(s0, s1, 6)
            self.assertTrue(s0 > lower)
            self.assertTrue(s1 < upper)


class kurtosis_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.statistics.kurtosis 
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(kurtosis_tests)) 


if __name__ == '__main__':
    unittest.main()
