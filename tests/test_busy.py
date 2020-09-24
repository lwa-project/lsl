"""
Unit test for regressions in the lsl.common.busy module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import sys
import time
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
from lsl.common import busy

__version__  = "0.2"
__author__    = "Jayce Dowell"


class busy_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the regressions in LSL."""
    
    def test_default(self):
        """Test the standard busy indicator, default options."""
        
        bi = busy.BusyIndicator()
        
        sys.stdout = StringIO()
        
        bi.start()
        time.sleep(1)
        bi.stop()
        
        sys.stdout = sys.__stdout__
        
    def test_color(self):
        """Test the standard busy indicator, default options and color."""
        
        bi = busy.BusyIndicator(color='green')
        
        sys.stdout = StringIO()
        
        bi.start()
        time.sleep(1)
        bi.stop()
        
        sys.stdout = sys.__stdout__
        
    def test_context(self):
        """Test the standard busy indicator as a context manager."""
        
        sys.stdout = StringIO()
        
        with busy.BusyIndicator() as bi:
            time.sleep(1)
            
        sys.stdout = sys.__stdout__
        
    def test_styles(self):
        """Test the various styles of the busy indicator plus."""
        
        for style in ('boomerang', 'pingpong', 'flow'):
            for color in (None, 'red'):
                sys.stdout = StringIO()
                
                with busy.BusyIndicatorPlus(style=style, color=color) as bi:
                    time.sleep(1)
                    
                sys.stdout = sys.__stdout__
                
         self.assertRaises(ValueError, busy.BusyIndicatorPlus, (), {'style':'testing?'})            


class busy_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.busy unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(busy_tests)) 


if __name__ == '__main__':
    unittest.main()
