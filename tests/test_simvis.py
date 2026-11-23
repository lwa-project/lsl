"""
Unit test for the lsl.sim.vis module.
"""

import os
import copy
import unittest
import ephem
import numpy as np

import aipy

from astropy.coordinates import SkyCoord

from lsl.sim import vis
from lsl.imaging.data import VisibilityData, VisibilityDataSet, PolarizationDataSet
from lsl.correlator import fx
from lsl.common import stations as lwa_common
import lsl.testing


__version__  = "0.1"
__author__    = "Jayce Dowell"


@lsl.testing.quiet_lsl_logging
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
            
        with self.subTest(type='astropy.coordinates.SkyCoord'):
            sc = SkyCoord('1h02m03s', '+89d00m00s', frame='fk5', equinox='J2000')
            
            out = vis.build_sim_data(aa, vis.SOURCES, phase_center=sc)
            self.assertAlmostEqual(out.phase_center._ra, sc.ra.rad, 6)
            self.assertAlmostEqual(out.phase_center._dec, sc.dec.rad, 6)
            
    def _expected_scale(self, freq, baselines, amps, delays, phase_offsets):
        """Per-baseline complex gain a_i*a_j*exp(2j*pi*f*(d_i-d_j) + 1j*(p_i-p_j))."""
        
        fq = freq / 1e9
        return np.array([amps[i]*amps[j]*np.exp(2j*np.pi*fq*(delays[i]-delays[j]) + 1j*(phase_offsets[i]-phase_offsets[j])) \
                         for i,j in baselines])
        
    def test_scale_data(self):
        """Test that scaling a data dictionary applies the expected per-baseline gains"""
        
        lwa1 = lwa_common.lwa1
        antennas = lwa1.antennas[0:20]
        N = len(antennas)
        rng = np.random.default_rng(42)
        amps = rng.uniform(0.5, 2.0, N)
        delays = rng.uniform(-5.0, 5.0, N)
        phases = rng.uniform(-np.pi, np.pi, N)
        
        for label,freqs in (('multi-channel', np.arange(30e6, 50e6, 1e6)), ('single-channel', np.array([30e6,]))):
            with self.subTest(case=label):
                aa = vis.build_sim_array(lwa1, antennas, freqs)
                out = vis.build_sim_data(aa, vis.SOURCES)
                orig = {pds.polarization: pds.data.copy() for pds in out}
                
                scl = vis.scale_data(out, amps, delays, phases)
                g = self._expected_scale(out.freq, out.baselines, amps, delays, phases)
                for pds in scl:
                    with self.subTest(pol=pds.polarization):
                        np.testing.assert_allclose(pds.data, orig[pds.polarization]*g, rtol=1e-5, atol=1e-5*np.abs(orig[pds.polarization]).max())
                        
                # The input should not be modified
                for pds in out:
                    np.testing.assert_array_equal(pds.data, orig[pds.polarization])
                    
                # Without phase offsets
                scl = vis.scale_data(out, amps, delays)
                g = self._expected_scale(out.freq, out.baselines, amps, delays, np.zeros(N))
                np.testing.assert_allclose(scl.XX.data, orig['XX']*g, rtol=1e-5, atol=1e-5*np.abs(orig['XX']).max())
                
                # Undoing the scaling should recover the original data
                rnd = vis.scale_data(vis.scale_data(out, amps, delays, phases), 1/amps, -delays, -phases)
                for pds in rnd:
                    np.testing.assert_allclose(pds.data, orig[pds.polarization], rtol=1e-5, atol=1e-5*np.abs(orig[pds.polarization]).max())
                    
        with self.subTest(case='VisibilityData'):
            freqs = np.arange(30e6, 50e6, 1e6)
            aa = vis.build_sim_array(lwa1, antennas, freqs)
            out = vis.build_sim_data(aa, vis.SOURCES)
            out_later = out.copy(include_pols=True)
            out_later.jd += 10/86400.0
            out2 = VisibilityData([out, out_later])
            
            scl2 = vis.scale_data(out2, amps, delays, phases)
            self.assertEqual(len(scl2), 2)
            g = self._expected_scale(out.freq, out.baselines, amps, delays, phases)
            for data_set in scl2:
                np.testing.assert_allclose(data_set.XX.data, out.XX.data*g, rtol=1e-5, atol=1e-5*np.abs(out.XX.data).max())
                
        with self.subTest(case='bad input'):
            self.assertRaises(TypeError, vis.scale_data, out.XX.data, amps, delays)
            
    def test_scale_data_convention(self):
        """Test that a positive delay in scale_data removes a signal lag of the same size"""
        
        # Three copies of the same antenna (same position and cable) so that
        # the correlator applies no relative geometric or cable delays
        ant = lwa_common.lwa1.antennas[0]
        ants = [copy.deepcopy(ant) for i in range(3)]
        for i,a in enumerate(ants):
            a.id = i + 1
            a.stand.id = 1000 + i
            
        # Antenna 1 sees the signal D samples later than antennas 0 and 2.  For
        # a delay before the mixer the baseband signal also picks up a phase
        # of -2*pi*fc*tau.
        fs, fc, D, nsamp = 19.6e6, 74e6, 3, 2**16
        tau = D/fs
        rng = np.random.default_rng(7)
        s = (rng.normal(size=nsamp+D) + 1j*rng.normal(size=nsamp+D)).astype(np.complex64)
        signals = np.array([s[D:], s[:-D]*np.exp(-2j*np.pi*fc*tau), s[D:]]).astype(np.complex64)
        baselines, freq, cross = fx.FXMaster(signals, ants, LFFT=64, sample_rate=fs, central_freq=fc, return_baselines=True)
        
        ds = VisibilityDataSet(2460000.5, freq, [(0,1),(0,2),(1,2)], np.zeros((3,3,freq.size)))
        ds.append(PolarizationDataSet('XX', data=cross))
        lag = np.array([0.0, tau*1e9, 0.0])
        cal = vis.scale_data(ds, np.ones(3), lag)
        
        # After the correction all baselines should have a flat, zero phase
        for b in range(3):
            np.testing.assert_allclose(np.angle(cal.XX.data[b]), 0.0, atol=0.1)
            
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
