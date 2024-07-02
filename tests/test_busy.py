"""
Unit test for regressions in the lsl.common.busy module.
"""

import sys
import time
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
from lsl.common import busy
import lsl.testing

__version__  = "0.3"
__author__    = "Jayce Dowell"


class busy_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the regressions in LSL."""
    
    def test_default(self):
        """Test the standard busy indicator, default options."""
        
        with lsl.testing.SilentVerbose():
            bi = busy.BusyIndicator()
            
            bi.start()
            time.sleep(1)
            bi.stop()
            
    def test_color(self):
        """Test the standard busy indicator, default options and color."""
        
        with lsl.testing.SilentVerbose():
            bi = busy.BusyIndicator(color='green')
            
            bi.start()
            time.sleep(1)
            bi.stop()
            
    def test_context(self):
        """Test the standard busy indicator as a context manager."""
        
        with lsl.testing.SilentVerbose():
            with busy.BusyIndicator() as bi:
                time.sleep(1)
                
    def test_styles(self):
        """Test the various styles of the busy indicator plus."""
        
        for style in ('boomerang', 'pingpong', 'flow'):
            for color in (None, 'red'):
                with self.subTest(style=style, color=color):
                    with lsl.testing.SilentVerbose():
                        with busy.BusyIndicatorPlus(style=style, color=color) as bi:
                            time.sleep(1)
                            
        with self.assertRaises(ValueError):
            busy.BusyIndicatorPlus(style='testing?')


class busy_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.busy unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(busy_tests)) 


if __name__ == '__main__':
    unittest.main()
