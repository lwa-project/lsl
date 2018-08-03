# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for the lsl.misc.beamformer module."""

import os
import unittest
import numpy

from lsl.misc import beamformer
from lsl.common import stations


__revision__ = "$Rev$"
__version__  = "0.2"
__author__    = "Jayce Dowell"


class beamformer_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.dp
    module."""
    
    def test_calcDelay(self):
        """Check that the beamformer.calcDelay function actually runs"""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        
        out = beamformer.calcDelay(antennas[:3])
        self.assertEqual(len(out), 3)
        
        out = beamformer.calcDelay(antennas[:3], freq=49.0e6)
        self.assertEqual(len(out), 3)
        
        out = beamformer.calcDelay(antennas[:3], freq=49.0e6, azimuth=45, elevation=30)
        self.assertEqual(len(out), 3)
        
    def test_pointing_limits(self):
        """Test that beamformer.calcDelay respects the pointing limits"""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        
        # Azimuth  checks
        self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, antennas[:3], 49.0e6, -5, 30)
        self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, antennas[:3], 49.0e6, 365, 30)
        
        # Elevation checks
        self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, antennas[:3], 49.0e6, 45, -5)
        self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, antennas[:3], 49.0e6, 45, 95)
        
    def test_intDelayAndSum(self):
        """Check that the beamformer.intDelayAndSum function actually runs"""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        data = numpy.random.rand(3, 1000)
        
        beam = beamformer.intDelayAndSum(antennas[:3], data, azimuth=45.0, elevation=30.0)
        
    def test_intBeamShape(self):
        """Check that the beamformer.intBeamShape function actually runs"""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        
        # Get a list of valid antennas
        ants = []
        for ant in antennas:
            if ant.pol == 1:
                continue
            if ant.getStatus() != 33:
                continue
                
            ants.append(ant)
            if len(ants) == 15:
                break
                
        # Multithreaded test for accuracy
        out = beamformer.intBeamShape(ants, azimuth=135.0, elevation=60.0)
        
        i = out.argmax()
        azDiff = numpy.abs(135.0 - i / 90)
        elDiff = numpy.abs(60.0 - i % 90)
        self.assertTrue(azDiff <= 1)
        self.assertTrue(elDiff <= 1)
        
        # Single threaded test for coverage
        out = beamformer.intBeamShape(ants[:1], azimuth=135.0, elevation=60.0, DisablePool=True)
            
    def test_phaseAndSum(self):
        """Check that the beamformer.phaseAndSum function actually runs"""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        data = numpy.random.rand(3, 1000)
        
        beam = beamformer.phaseAndSum(antennas[:3], data, azimuth=45.0, elevation=30.0)
        
    def test_phaseBeamShape(self):
        """Check that the beamformer.phaseBeamShape function actually runs"""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        
        # Get a list of valid antennas
        ants = []
        for ant in antennas:
            if ant.pol == 1:
                continue
            if ant.getStatus() != 33:
                continue
                
            ants.append(ant)
            if len(ants) == 15:
                break
                
        # Test for accuracy
        out = beamformer.phaseBeamShape(ants, azimuth=135.0, elevation=60.0)
        
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
