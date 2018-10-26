# -*- coding: utf-8 -*-

"""Unit test for the lsl.sim.vis module."""

import os
import unittest
import numpy

from lsl.sim import vis
from lsl.imaging.data import VisibilityData
from lsl.common import stations as lwa_common


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class simvis_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.vis
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings."""
        
        numpy.seterr(all='ignore')
        
    def test_build_aa_flat(self):
        """Test building a antenna array object with uniform sky response."""
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        
        aa = vis.build_sim_array(lwa1, antennas, freqs, force_flat=True)
        
    def test_build_aa(self):
        """Test building a antenna array object with realistic sky response."""
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        # Check the number of stands
        self.assertEqual(len(aa.ants), len(antennas))
        
        # Check the frequencies comming out
        for fo, fi in zip(aa.get_afreqs(), freqs):
            self.assertAlmostEqual(fo, fi/1e9, 6)
            
    def test_build_data(self):
        """Test building simulated visibility data"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Do a check of frequencies
        for fa, fq in zip(out.freq, freqs):
            self.assertAlmostEqual(fa, fq, 6)
            
        # Do a check to make sure that the polarizations
        self.assertEqual(out.npol, 4)
        self.assertTrue('XX' in out.pols)
        self.assertTrue('XY' in out.pols)
        self.assertTrue('YX' in out.pols)
        self.assertTrue('YY' in out.pols)
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Do a check of frequencies
        for fa, fq in zip(out.freq, freqs):
            self.assertAlmostEqual(fa, fq, 6)
            
        # Do a check to make sure that the polarizations
        self.assertEqual(out.npol, 4)
        self.assertTrue('XX' in out.pols)
        self.assertTrue('XY' in out.pols)
        self.assertTrue('YX' in out.pols)
        self.assertTrue('YY' in out.pols)
        
    def test_build_data_res(self):
        """Test building simulated visibility data with resolved sources"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, resolve_src=True)
        
        # Do a check of keys
        # Do a check of frequencies
        for fa, fq in zip(out.freq, freqs):
            self.assertAlmostEqual(fa, fq, 6)
            
        # Do a check to make sure that the polarizations
        self.assertEqual(out.npol, 4)
        self.assertTrue('XX' in out.pols)
        self.assertTrue('XY' in out.pols)
        self.assertTrue('YX' in out.pols)
        self.assertTrue('YY' in out.pols)
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, resolve_src=True)
        
        # Do a check of frequencies
        for fa, fq in zip(out.freq, freqs):
            self.assertAlmostEqual(fa, fq, 6)
            
        # Do a check to make sure that the polarizations
        self.assertEqual(out.npol, 4)
        self.assertTrue('XX' in out.pols)
        self.assertTrue('XY' in out.pols)
        self.assertTrue('YX' in out.pols)
        self.assertTrue('YY' in out.pols)
        
    def test_scale_data(self):
        """Test that we can scale a data dictionary without error"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Scale
        amp = vis.scale_data(out, numpy.ones(len(antennas))*2, numpy.zeros(len(antennas)))
        # Delay
        phs = vis.scale_data(out, numpy.ones(len(antennas)), numpy.ones(len(antennas)))
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Scale
        amp = vis.scale_data(out, numpy.ones(len(antennas))*2, numpy.zeros(len(antennas)))
        # Delay
        phs = vis.scale_data(out, numpy.ones(len(antennas)), numpy.ones(len(antennas)))
        
        #
        # VisibilityData test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        out2 = VisibilityData(out)
        
        # Scale
        amp2 = vis.scale_data(out2, numpy.ones(len(antennas))*2, numpy.zeros(len(antennas)))
        # Delay
        phs2 = vis.scale_data(out2, numpy.ones(len(antennas)), numpy.ones(len(antennas)))
        
        
    def test_shift_data(self):
        """Test that we can shift the uvw coordinates of a data dictionary 
        without error"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Shift
        sft = vis.shift_data(out, aa)
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Shift
        sft = vis.shift_data(out, aa)
        
        #
        # VisibilityData test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        out2 = VisibilityData()
        out2.append( out )
        
        # Shift
        sft2 = vis.shift_data(out2, aa)
        
    def test_add_noise(self):
        """Test that we can add baseline noise to a data dictionary without
        error"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Add in the noise
        na = vis.add_baseline_noise(out, 15e3, 0.061)
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Add in the noise
        na = vis.add_baseline_noise(out, 15e3, 0.061, bandwidth=1e6)
        
        #
        # VisibilityData test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = numpy.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        out2 = VisibilityData()
        out2.append( out )
        
        # Add in the noise
        na2 = vis.add_baseline_noise(out2, 15e3, 0.061)


class  simvis_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.sim.vis units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(simvis_tests)) 


if __name__ == '__main__':
    unittest.main()
