"""
Unit test for regressions in the lsl.sim.beam module.
"""

import os
import unittest
import numpy as np

from lsl.sim import beam as simbeam

__version__  = "0.1"
__author__    = "Jayce Dowell"


class sim_beam_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the regressions in LSL."""
    
    def test_model_list(self):
        """Test listing all available dipole response models."""
        
        names = simbeam.get_avaliable_models()
        self.assertEqual(len(names), 7)
        
    def test_response(self):
        """Test that the beam_response function actually runs."""
        
        az0 = 0.0
        alt0 = 45.0
        
        az = np.arange(360)
        alt = np.arange(360)*0 + 45
        
        az2 = np.arange(360)
        alt2 = np.arange(90)
        az2, alt2 = np.meshgrid(az2, alt2)
        
        for model in simbeam.get_avaliable_models():
            with self.subTest(model=model, dim=0):
                pattern = simbeam.beam_response(model, 'XX', az0, alt0, frequency=70e6)
                self.assertTrue(isinstance(pattern, float))
                self.assertTrue(np.isfinite(pattern))
                
            with self.subTest(model=model, dim=1):
                pattern = simbeam.beam_response(model, 'XX', az, alt, frequency=70e6)
                self.assertEqual(pattern.shape, az.shape)
                self.assertTrue(np.isfinite(pattern[0]))
                
            with self.subTest(model=model, dim=2):
                pattern = simbeam.beam_response(model, 'XX', az2, alt2, frequency=70e6)
                self.assertEqual(pattern.shape, az2.shape)
                self.assertTrue(np.isfinite(pattern[0,0]))
                
        for bad_model in ('bad'):
            with self.subTest(bad_model=bad_model):
                self.assertRaises(ValueError, simbeam.beam_response, bad_model, 'XX', az, alt, frequency=70e6)
                
    def test_response_pols(self):
        """
        Test that the beam_response function works with across polarizations."""
        
        az = np.arange(360)
        alt = np.arange(360)*0 + 45
        
        az2 = np.arange(360)
        alt2 = np.arange(90)
        az2, alt2 = np.meshgrid(az2, alt2)
        
        for model in simbeam.get_avaliable_models():
            for pol in ('XX', 'XY', 'YX', 'YY', 'I', 'Q', 'U', 'V'):
                with self.subTest(pol=pol, dim=1):
                    pattern = simbeam.beam_response(model, pol, az, alt, frequency=70e6)
                    self.assertEqual(pattern.shape, az.shape)
                    
                with self.subTest(pol=pol, dim=2):
                    pattern = simbeam.beam_response(model, pol, az2, alt2, frequency=70e6)
                    self.assertEqual(pattern.shape, az2.shape)
                    
            for bad_pol in ('X', 'Z'): 
                with self.subTest(bad_pol=bad_pol):
                    self.assertRaises(ValueError, simbeam.beam_response, model, bad_pol, az, alt, frequency=70e6)
                    
            break
            
    def test_mueller(self):
        """
        Test that the mueller_matrix function actually runs.
        """
        
        az = np.arange(360)
        alt = np.arange(360)*0 + 45
        
        az2 = np.arange(360)
        alt2 = np.arange(90)
        az2, alt2 = np.meshgrid(az2, alt2)
        
        for model in simbeam.get_avaliable_models():
            with self.subTest(model=model, dim=1):
                pattern = simbeam.mueller_matrix(model, az, alt, frequency=70e6)
                self.assertEqual(pattern.shape, (4,4)+az.shape)
                
            with self.subTest(model=model, dim=2):
                pattern = simbeam.mueller_matrix(model, az2, alt2, frequency=70e6)
                self.assertEqual(pattern.shape, (4,4)+az2.shape)


class sim_beam_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.sim.beam unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(sim_beam_tests)) 


if __name__ == '__main__':
    unittest.main()
