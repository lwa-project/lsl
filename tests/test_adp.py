"""
Unit test for the lsl.common.adp module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import warnings
import unittest
import numpy

from lsl.common import adp
from lsl.common import stations


__version__  = "0.1"
__author__    = "Jayce Dowell"


class adp_bandpass_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.adp
    module function for the bandpass."""
    
    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
    
    def test_tbn_bandpass(self):
        """Test that the TBN bandpass generator actually runs."""
        
        fnc = adp.tbn_filter(sample_rate=1e5, npts=256)
        junk = fnc(1e3)

    def test_drx_bandpass(self):
        """Test that the DRX bandpass generator actually runs."""
        
        fnc = adp.drx_filter(sample_rate=19.6e6, npts=256)
        junk = fnc(1e3)


class adp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.adp
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(adp_bandpass_tests)) 


if __name__ == '__main__':
    unittest.main()
