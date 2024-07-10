"""
Unit test for the lsl.sim.vis module.
"""

import os
import unittest
import ephem
import numpy as np

import aipy

from lsl.sim import vis
from lsl.imaging.data import VisibilityData
from lsl.common import stations as lwa_common


__version__  = "0.1"
__author__    = "Jayce Dowell"


class simvis_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.sim.vis
    module."""
    
    def setUp(self):
        """Turn off all numpy warnings."""
        
        np.seterr(all='ignore')
        
    def test_earth_satellite(self):
        tle = ["ISS (ZARYA)  ",
               "1 25544U 98067A   20118.40744172 -.00000643  00000-0 -34717-5 0  9993",
               "2 25544  51.6437 239.7966 0001371 201.5257 265.1663 15.49320404224148"]
        iss = vis.RadioEarthSatellite(tle, 0.150, tpower=1.0, tbw=100e5)
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)

        aa = vis.build_sim_array(lwa1, antennas, freqs)
        aa.set_unixtime(1588026422.0)
        iss.compute(aa)
        
        iss.get_crds('eq')
        iss.get_crds('top')
        
    def test_build_aa_flat(self):
        """Test building a antenna array object with uniform sky response."""
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 10e6)
        
        aa = vis.build_sim_array(lwa1, antennas, freqs, force_flat=True)
        aa[0].bm_response((0,0,1), pol='x')
        aa[0].bm_response((0,0,1), pol='y')
        bm1 = aa[0].get_beam_shape()
        
        az = np.zeros((360,90))
        for i in range(360):
            az[i,:] = i*np.pi/180.0
        alt = np.zeros((360,90))
        for i in range(90):
            alt[:,i] = i*np.pi/180.0
        xyz = aipy.coord.azalt2top(np.concatenate([[az],[alt]]))
        bm2 = aa[0].bm_response(xyz, pol='x')
        np.testing.assert_allclose(bm1.transpose(2,0,1), bm2)
        
    def test_build_aa_gaussian(self):
        """Test building a antenna array object with Gaussian sky response."""
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 10e6)

        aa = vis.build_sim_array(lwa1, antennas, freqs, force_gaussian=True)
        aa[0].bm_response((0,0,1), pol='x')
        aa[0].bm_response((0,0,1), pol='y')
        bm1 = aa[0].get_beam_shape()
        
        az = np.zeros((360,90))
        for i in range(360):
            az[i,:] = i*np.pi/180.0
        alt = np.zeros((360,90))
        for i in range(90):
            alt[:,i] = i*np.pi/180.0
        xyz = aipy.coord.azalt2top(np.concatenate([[az],[alt]]))
        bm2 = aa[0].bm_response(xyz, pol='x')
        np.testing.assert_allclose(bm1.transpose(2,0,1), bm2)
        
    def test_build_aa(self):
        """Test building a antenna array object with realistic sky response."""
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        aa.set_asp_filter('split')
        aa[0].bm_response((0,0,1))
        # Check the number of stands
        self.assertEqual(len(aa.ants), len(antennas))
        
        # Check the frequencies comming out
        np.testing.assert_allclose(aa.get_afreqs(), freqs/1e9)
        
        # Check that other methods even run
        aa.get_baseline_fast(0, 1)
        aa.gen_uvw_fast(0, 1)
        aa.gen_phs_fast('z', 0, 1)
        
        aa.set_unixtime(1588026422.0)
        vis.SOURCES['crab'].compute(aa)
        aa.gen_phs_fast(vis.SOURCES['crab'], 0, 1)
            
    def test_build_data(self):
        """Test building simulated visibility data"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Do a check of frequencies
        np.testing.assert_allclose(out.freq, freqs)
        
        # Do a check to make sure that the polarizations
        self.assertEqual(out.npol, 4)
        self.assertTrue('XX' in out.pols)
        self.assertTrue('XY' in out.pols)
        self.assertTrue('YX' in out.pols)
        self.assertTrue('YY' in out.pols)
        
        # Try a simulation on a single baselines
        aa.sim(0, 1)
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Do a check of frequencies
        np.testing.assert_allclose(out.freq, freqs)
        
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
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, resolve_src=True)
        
        # Do a check of keys
        # Do a check of frequencies
        np.testing.assert_allclose(out.freq, freqs)
        
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
        freqs = np.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, resolve_src=True)
        
        # Do a check of frequencies
        np.testing.assert_allclose(out.freq, freqs)
        
        # Do a check to make sure that the polarizations
        self.assertEqual(out.npol, 4)
        self.assertTrue('XX' in out.pols)
        self.assertTrue('XY' in out.pols)
        self.assertTrue('YX' in out.pols)
        self.assertTrue('YY' in out.pols)
        
    def test_build_data_phase_center(self):
        """Test building simulated visibility data with different phase center types."""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary - a few times
        with self.subTest(type='str'):
            out = vis.build_sim_data(aa, vis.SOURCES, phase_center='z')
            self.assertRaises(ValueError, vis.build_sim_data, aa, vis.SOURCES, phase_center='notgoingtowork')
            
        with self.subTest(type='ephem.Body'):
            bdy = ephem.FixedBody()
            bdy._ra = '1:02:03'
            bdy._dec = '+89:00:00'
            bdy._epoch = ephem.J2000
            
            out = vis.build_sim_data(aa, vis.SOURCES, phase_center=bdy)
            self.assertAlmostEqual(out.phase_center._ra, bdy._ra, 6)
            self.assertAlmostEqual(out.phase_center._dec, bdy._dec, 6)
            self.assertAlmostEqual(out.phase_center._epoch, bdy._epoch, 6)
            
    def test_scale_data(self):
        """Test that we can scale a data dictionary without error"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Scale
        amp = vis.scale_data(out, np.ones(len(antennas))*2, np.zeros(len(antennas)))
        # Delay
        phs = vis.scale_data(out, np.ones(len(antennas)), np.ones(len(antennas)))
        
        #
        # Single-channel test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.array([30e6,])
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        
        # Scale
        amp = vis.scale_data(out, np.ones(len(antennas))*2, np.zeros(len(antennas)))
        # Delay
        phs = vis.scale_data(out, np.ones(len(antennas)), np.ones(len(antennas)))
        
        #
        # VisibilityData test
        #
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES)
        out2 = VisibilityData(out)
        
        # Scale
        amp2 = vis.scale_data(out2, np.ones(len(antennas))*2, np.zeros(len(antennas)))
        # Delay
        phs2 = vis.scale_data(out2, np.ones(len(antennas)), np.ones(len(antennas)))
        
        
    def test_shift_data(self):
        """Test that we can shift the uvw coordinates of a data dictionary 
        without error"""
        
        # Setup
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
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
        freqs = np.array([30e6,])
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
        freqs = np.arange(30e6, 50e6, 1e6)
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
        freqs = np.arange(30e6, 50e6, 1e6)
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
        freqs = np.array([30e6,])
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
        freqs = np.arange(30e6, 50e6, 1e6)
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
