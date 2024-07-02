"""
Unit test for the lsl.misc.scattering module.
"""

import os
import time
import warnings
import unittest
import numpy as np

from lsl.misc import scattering


__version__  = "0.1"
__author__    = "Jayce Dowell"

class scattering_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.scattering
    module."""
    
    def setUp(self):
        self.t = np.linspace(0.0, 1.333, 512)
        self.intrinsic = np.zeros(self.t.shape, dtype=np.float32)
        self.intrinsic[128] += 10
        self.intrinsic[128+32] += 5
        screen = scattering.thin(self.t, 0.015)
        self.scattered = np.abs(np.fft.ifft(np.fft.fft(self.intrinsic)*np.fft.fft(screen[::-1]).conj()))
        
    def test_thin(self):
        """Test thin screen descattering"""
        
        tau, merit, profile = scattering.unscatter(self.t, self.scattered, 0.005, 0.030, 0.005, 
                                                   screen=scattering.thin, verbose=False)
        self.assertAlmostEqual(tau, 0.015, 3)
        
    def test_thick(self):
        """Test thick screen descattering"""
        
        tau, merit, profile = scattering.unscatter(self.t, self.scattered, 0.005, 0.030, 0.005,
                                                   screen=scattering.thick, verbose=False)
        
    def test_uniform(self):
        """Test uniforms screen descattering"""
        
        tau, merit, profile = scattering.unscatter(self.t, self.scattered, 0.005, 0.030, 0.005,
                                                   screen=scattering.uniform, verbose=False)


class scattering_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.scattering 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(scattering_tests)) 


if __name__ == '__main__':
    unittest.main()
