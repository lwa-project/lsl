"""
Unit test for the lsl.correlator.filterbank module.
"""

import os
import numpy as np
import unittest

from lsl.correlator import filterbank
from lsl.correlator.fx import null_window
import lsl.testing


__version__  = "0.3"
__author__    = "Jayce Dowell"


class filterbank_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.filterbank
    module."""
    
    def run_filterbank_test(self, dtype, ntap=2, nchan=256, window=null_window):
        data = np.random.rand(nchan*ntap*4)
        data = data.astype(dtype)
        fnc = "fft%i" % ntap
        if hasattr(filterbank, fnc):
            getattr(filterbank, fnc)(data, nchan, window=window)
        else:
            filterbank.fft(data, nchan, P=ntap, window=window)
            
    def test_filterbank_real(self):
        """Test that the filterbank works on real-valued data."""
        
        for ntap in (1, 2, 4, 8, 16):
            for dtype in (np.float32, np.float64):
                with self.subTest(ntap=ntap, dtype=dtype):
                    self.run_filterbank_test(dtype, ntap=ntap)
                    
    def test_filterbank_complex(self):
        """Test that the filterbank works on real-valued data."""
        
        for ntap in (1, 2, 4, 8, 16):
            for dtype in (np.complex64, np.complex128):
                with self.subTest(ntap=ntap, dtype=dtype):
                    self.run_filterbank_test(dtype, ntap=ntap)
                    
    def test_filterbank_window(self):
        """Test that window functions can be passed to the filterbank."""
        
        #
        # Real
        #
        
        for ntap in (1, 2, 4, 8, 16):
            for dtype in (np.float32, np.float64):
                with self.subTest(ntap=ntap, dtype=dtype):
                    self.run_filterbank_test(dtype, ntap=ntap, window=np.blackman)
                    
        #
        # Complex
        #
        
        for ntap in (1, 2, 4, 8, 16):
            for dtype in (np.float32, np.float64):
                with self.subTest(ntap=ntap, dtype=dtype):
                    self.run_filterbank_test(dtype, ntap=ntap, window=np.hamming)


class filterbank_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.correlator.filterbank
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(filterbank_tests)) 


if __name__ == '__main__':
    unittest.main()
