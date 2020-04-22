"""
Unit test for regressions in the lsl.common.progress module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import numpy
import unittest

from lsl.common import progress

__version__  = "021"
__author__    = "Jayce Dowell"


class progress_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the regressions in LSL."""
    
    def test_bar_default(self):
        """Test the progress bar, default options."""
        
        pbar = progress.ProgressBar()
        for i in range(101):
            pbar.inc(2)
            pbar.dec(1)
        
    def test_bar_attributes(self):
        """Test the progress bar's attributes."""
        
        pbar2 = progress.ProgressBar()
        for i in range(101):
            pbar2 += 2
            pbar2 -= 1
            
            pbar2 = pbar2 + 1
            pbar2 = pbar2 - 1
            
            self.assertEqual(pbar2.amount, i+1)
            
    def test_bar_show(self):
        """Test the progress bar's rendering."""
        
        # With percentage
        pbar = progress.ProgressBar()
        for i in range(101):
            pbar.inc(1)
            pbar.show()
            
        # Without percentage
        pbar = progress.ProgressBar(print_percent=False)
        for i in range(101):
            pbar.inc(1)
            pbar.show()
            
    def test_bar_plus_default(self):
        """Test the progress bar, default options."""
        
        pbar = progress.ProgressBarPlus()
        for i in range(101):
            pbar.inc(2)
            pbar.dec(1)
        
    def test_bar_plus_attributes(self):
        """Test the progress bar's attributes."""
        
        pbar2 = progress.ProgressBarPlus()
        for i in range(101):
            pbar2 += 2
            pbar2 -= 1
            
            pbar2 = pbar2 + 1
            pbar2 = pbar2 - 1
            
            self.assertEqual(pbar2.amount, i+1)
            
    def test_bar_plus_show(self):
        """Test the progress bar's rendering."""
        
        # With percentage
        pbar = progress.ProgressBarPlus()
        for i in range(101):
            pbar.inc(1)
            pbar.show()
            
        # Without percentage
        pbar = progress.ProgressBarPlus(print_percent=False)
        for i in range(101):
            pbar.inc(1)
            pbar.show()


class progress_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.progress unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(progress_tests)) 


if __name__ == '__main__':
    unittest.main()
