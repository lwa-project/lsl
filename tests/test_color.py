"""
Unit tests for the lsl.common.color module.
"""

import unittest

from lsl.common import color

__version__  = "0.1"
__author__    = "Jayce Dowell"


class color_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.color
    module."""
    
    def test_colorfy(self):
        """Test colorfy"""
        
        color.colorfy("This is a plain string")
        color.colorfy("This is a {{%red string}}")
        color.colorfy("This is an {{%underline string}}")
        color.colorfy("This is an {{%yellow {{%underline string}}}}")


class color_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.color units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(color_tests)) 


if __name__ == '__main__':
    unittest.main()
    
