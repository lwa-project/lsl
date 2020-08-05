"""
Unit test for the lsl.misc.beamformer module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import unittest
import numpy

from lsl.misc import beamformer
from lsl.common import stations


__version__  = "0.2"
__author__    = "Jayce Dowell"


class beamformer_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.dp
    module."""
    
    def test_calc_delay(self):
        """Check that the beamformer.calc_delay function actually runs"""
        
        station = stations.lwa1
        antennas = station.antennas
        
        out = beamformer.calc_delay(antennas[:3])
        self.assertEqual(len(out), 3)
        
        out = beamformer.calc_delay(antennas[:3], freq=49.0e6)
        self.assertEqual(len(out), 3)
        
        out = beamformer.calc_delay(antennas[:3], freq=49.0e6, azimuth=45, elevation=30)
        self.assertEqual(len(out), 3)
        
    def test_pointing_limits(self):
        """Test that beamformer.calc_delay respects the pointing limits"""
        
        station = stations.lwa1
        antennas = station.antennas
        
        # Azimuth  checks
        self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, -5, 30)
        self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, 365, 30)
        
        # Elevation checks
        self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, 45, -5)
        self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, 45, 95)
        
    def test_int_delay_and_sum(self):
        """Check that the beamformer.int_delay_and_sum function actually runs"""
        
        station = stations.lwa1
        antennas = station.antennas
        data = numpy.random.rand(3, 1000)
        
        beam = beamformer.int_delay_and_sum(antennas[:3], data, azimuth=45.0, elevation=30.0)
        
    def test_int_beam_shape(self):
        """Check that the beamformer.int_beam_shape function actually runs"""
        
        station = stations.lwa1
        antennas = station.antennas
        
        # Get a list of valid antennas
        ants = []
        for ant in antennas:
            if ant.pol == 1:
                continue
            if ant.combined_status != 33:
                continue
                
            ants.append(ant)
            if len(ants) == 16:
                break
                
        # Multithreaded test for accuracy
        out = beamformer.int_beam_shape(ants, azimuth=135.0, elevation=60.0)
        
        i = out.argmax()
        azDiff = numpy.abs(135.0 - i / 90)
        elDiff = numpy.abs(60.0 - i % 90)
        self.assertTrue(azDiff <= 1)
        self.assertTrue(elDiff <= 1)
        
        # Single threaded test for coverage
        out = beamformer.int_beam_shape(ants[:1], azimuth=135.0, elevation=60.0, disable_pool=True)
            
    def test_phase_and_sum(self):
        """Check that the beamformer.phase_and_sum function actually runs"""
        
        station = stations.lwa1
        antennas = station.antennas
        data = numpy.random.rand(3, 1000)
        
        beam = beamformer.phase_and_sum(antennas[:3], data, azimuth=45.0, elevation=30.0)
        
    def test_phase_beam_shape(self):
        """Check that the beamformer.phase_beam_shape function actually runs"""
        
        station = stations.lwa1
        antennas = station.antennas
        
        # Get a list of valid antennas
        ants = []
        for ant in antennas:
            if ant.pol == 1:
                continue
            if ant.combined_status != 33:
                continue
                
            ants.append(ant)
            if len(ants) == 16:
                break
                
        # Test for accuracy
        out = beamformer.phase_beam_shape(ants, azimuth=135.0, elevation=60.0)
        
        i = out.argmax()
        azDiff = numpy.abs(135.0 - i / 90)
        elDiff = numpy.abs(60.0 - i % 90)
        self.assertTrue(azDiff <= 1)
        self.assertTrue(elDiff <= 1)


class  beamformer_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.sim.vis units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(beamformer_tests)) 


if __name__ == '__main__':
    unittest.main()
