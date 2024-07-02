"""
Unit test for the lsl.misc.telemetry module.
"""

import os
import time
import warnings
import unittest

from lsl.misc import telemetry


__version__  = "0.1"
__author__    = "Jayce Dowell"


@telemetry.track_function
def _trivial(x):
    return x + 5


@telemetry.track_function_timed
def _trivial2(x):
    return x + 6


class _TrivialClass(object):
    @telemetry.track_method
    def one(self, x):
        return x + 5
        
    @telemetry.track_method_timed
    def two(self, x):
        return x + 6


class telemetry_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.statistics.robust
    module."""
    
    def test_script(self):
        """Test telemetry tracking of scripts"""
        
        telemetry.track_script()
        
    def test_function(self):
        """Test telemetry tracking of functions"""
        
        for i in (0, 5, 10, 15, 20):
            _trivial(i)
            
    def test_function_timed(self):
        """Test telemetry tracking of functions with timing information"""
            
        for i in (0, 5, 10, 15, 20):
            _trivial2(i)
            
    def test_method(self):
        """Test telemetry tracking of methods"""
        
        c = _TrivialClass()
        for i in (0, 5, 10, 15, 20):
            c.one(i)
            
    def test_method_timed(self):
        """Test telemetry tracking of methods with timing information"""
        
        c = _TrivialClass()
        for i in (0, 5, 10, 15, 20):
            c.two(i)
    
    def tearDown(self):
        telemetry._telemetry_client.clear()


class telemetry_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.telemetry 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(telemetry_tests)) 


if __name__ == '__main__':
    unittest.main()
