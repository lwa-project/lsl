"""
Unit test for the lsl.misc.rfutils module.
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

from lsl.misc import rfutils


__version__  = "0.1"
__author__    = "Jayce Dowell"

class rfutils_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.rfutils
    module."""
    
    def test_dBd(self):
        """Test the dBd conversions"""
        
        # https://www.g4urh.co.uk/amateur_radio/watt_converter.php
        self.assertAlmostEqual(rfutils.dBi_to_dBd(10), 7.851561519523021, 5)
        self.assertAlmostEqual(rfutils.dBd_to_dBi(7.851561519523021), 10, 5)
        
    def test_area(self):
        """Test the effective area calculator"""
        
        # https://www.atnf.csiro.au/people/Tasso.Tzioumis/sms2014/presentations/Clegg(RF_Engineering).pptx.pdf
        self.assertAlmostEqual(rfutils.calculate_effective_area(0.36), 994, 0)
        
    def test_gain(self):
        """Test the gain calculator"""
        
        # https://www.atnf.csiro.au/people/Tasso.Tzioumis/sms2014/presentations/Clegg(RF_Engineering).pptx.pdf
        self.assertAlmostEqual(rfutils.gain_to_dBi(2.5e-5, 1.8e9), 15, 0)
        
    def test_sefd(self):
        """Test the SEFD calculator"""
        
        # https://web.njit.edu/~gary/728/Lecture5.html
        sefd = rfutils.calculate_sefd(100.0, gain=0.04, area=numpy.pi*10**2)
        self.assertAlmostEqual(sefd, 100/0.04, 6)
        
    def test_Jy(self):
        """Test the Jy conversions"""
        
        # https://www.atnf.csiro.au/people/Tasso.Tzioumis/sms2014/presentations/Clegg(RF_Engineering).pptx.pdf
        gain = rfutils.dBi_to_gain(0.0, 1.8e9)
        self.assertAlmostEqual(rfutils.Jy_to_dBm(1.0, 10e6, gain), -187, 0)
        self.assertAlmostEqual(rfutils.dBm_to_Jy(-187.0, 10e6, gain), 1.0, 0)


class rfutils_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.rfutils 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(rfutils_tests)) 


if __name__ == '__main__':
    unittest.main()
