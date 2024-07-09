"""
Unit test for regressions in the lsl.common.progress module.
"""

import os
import unittest

from lsl.common import progress
import lsl.testing

__version__  = "0.3"
__author__    = "Jayce Dowell"


class progress_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the regressions in LSL."""
    
    def run_bar_default_test(self, bartype, **kwds):
        pbar = bartype(**kwds)
        for i in range(101):
            pbar.inc(2)
            pbar.dec(1)
    
    def test_bar_default(self):
        """Test the progress bar."""
        
        for bartype in (progress.ProgressBar, progress.ProgressBarPlus, progress.DownloadBar):
            with self.subTest(bartype=bartype):
                self.run_bar_default_test(bartype)
                self.run_bar_default_test(bartype, color='green')
                
    def run_bar_attibutes_test(self, bartype):
        pbar2 = bartype()
        for i in range(101):
            pbar2 += 2
            pbar2 -= 1
            
            pbar2 = pbar2 + 1
            pbar2 = pbar2 - 1
            
            self.assertEqual(pbar2.amount, i+1)
            
    def test_bar_attributes(self):
        """Test the progress bar's attributes."""
        
        for bartype in (progress.ProgressBar, progress.ProgressBarPlus, progress.DownloadBar):
            with self.subTest(bartype=bartype):
                self.run_bar_attibutes_test(bartype)
                
    def run_bar_show_test(self, bartype, **kwds):
        # With percentage
        pbar = bartype(**kwds)
        for i in range(101):
            pbar.inc(1)
            pbar.show()
            
        # Without percentage
        pbar = bartype(print_percent=False, **kwds)
        for i in range(101):
            pbar.inc(1)
            pbar.show()
            
    def test_bar_show(self):
        """Test the progress bar's rendering."""
        
        for bartype in (progress.ProgressBar, progress.ProgressBarPlus, progress.DownloadBar):
            with self.subTest(bartype=bartype):
                self.run_bar_show_test(bartype)
                self.run_bar_show_test(bartype, color='green')


class progress_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.progress unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(progress_tests)) 


if __name__ == '__main__':
    unittest.main()
