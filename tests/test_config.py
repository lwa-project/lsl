"""
Unit test for the lsl.config module.
"""

import os
import unittest
import tempfile
import shutil
import importlib
import inspect

from lsl.config import LSL_CONFIG


__version__  = "0.1"
__author__    = "Jayce Dowell"


class config_tests(unittest.TestCase):
    def test_save_load(self):
        """Test the LSL configuration save and load functions."""
        
        LSL_CONFIG._save_config()
        
        # Force a reload
        mod = inspect.getmodule(LSL_CONFIG)
        importlib.reload(mod)
        
    def test_str(self):
        """Test converting the LSL configuration to a string."""
        str(LSL_CONFIG)
        
    def test_repr(self):
        """Test building a text representation of the LSL configuration."""
        repr(LSL_CONFIG)
        
    def test_set(self):
        """Test the set function of the LSL configuration."""
        
        ref_value = LSL_CONFIG.get('download.timeout')*1
        
        LSL_CONFIG.set('download.timeout', 300)
        self.assertTrue(LSL_CONFIG.get('download.timeout') != ref_value)
        
        LSL_CONFIG.set('download.timeout', ref_value)
        self.assertEqual(LSL_CONFIG.get('download.timeout'), ref_value)
        
    def test_set_temp(self):
        """Test the set_temp function of the LSL configuration."""
        
        ref_value = LSL_CONFIG.get('download.timeout')*1
        
        with LSL_CONFIG.set_temp('download.timeout', 1000000):
            self.assertTrue(LSL_CONFIG.get('download.timeout') != ref_value)
            
        self.assertEqual(LSL_CONFIG.get('download.timeout'), ref_value)
        
    def test_sub_str(self):
        """Test converting a section of the LSL configuration to a string."""
        
        dconfig = LSL_CONFIG.view('download')
        str(dconfig)
        
    def test_sub_repr(self):
        """Test building a text representation for a section of the LSL configuration."""
        
        dconfig = LSL_CONFIG.view('download')
        repr(dconfig)
        
    def test_sub_set(self):
        """Test the set function of the LSL configuration."""
        
        dconfig = LSL_CONFIG.view('download')
        ref_value = dconfig.get('timeout')*1
        
        dconfig.set('timeout', 300)
        self.assertTrue(dconfig.get('timeout') != ref_value)
        
        dconfig.set('timeout', ref_value)
        self.assertEqual(dconfig.get('timeout'), ref_value)
        
    def test_sub_set_temp(self):
        """Test the set_temp function for a section of the LSL configuration."""
        
        dconfig = LSL_CONFIG.view('download')
        ref_value = dconfig.get('timeout')*1
        
        with dconfig.set_temp('timeout', 1000000):
            self.assertTrue(dconfig.get('timeout') != ref_value)
            
        self.assertEqual(dconfig.get('timeout'), ref_value)


class config_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.config units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(config_tests)) 


if __name__ == '__main__':
    unittest.main()
