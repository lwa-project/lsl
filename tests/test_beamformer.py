"""
Unit test for the lsl.misc.beamformer module.
"""

import os
import unittest
import ephem
import numpy as np

from astropy.coordinates import Angle

from lsl.misc import beamformer, _beamformer
from lsl.common import stations
import lsl.testing


__version__  = "0.3"
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
        
        out = beamformer.calc_delay(antennas[:3], freq=49.0e6, azimuth=45, altitude=30)
        self.assertEqual(len(out), 3)
        
        # Angle support
        az = 45
        alt = 60
        out0 = beamformer.calc_delay(antennas[:3], freq=49.0e6, azimuth=az, altitude=alt)
        
        az = Angle('45deg')
        alt = Angle('60d00m00s')
        out1 = beamformer.calc_delay(antennas[:3], freq=49.0e6, azimuth=az, altitude=alt)
        np.testing.assert_allclose(out0, out1)
        
        az = ephem.degrees('45')
        alt = ephem.degrees('60')
        out1 = beamformer.calc_delay(antennas[:3], freq=49.0e6, azimuth=az, altitude=alt)
        np.testing.assert_allclose(out0, out1)
        
        # Multi-frequency support
        out0 = beamformer.calc_delay(antennas[:3], freq=49.0e6, azimuth=az, altitude=alt)
        out1 = beamformer.calc_delay(antennas[:3], freq=[49.0e6], azimuth=az, altitude=alt)
        np.testing.assert_allclose(out0, out1[:,0])
        
        out0 = beamformer.calc_delay(antennas[:3], freq=[49.0e6, 49.0e6], azimuth=az, altitude=alt)
        np.testing.assert_allclose(out0[:,0], out0[:,1])
        
    def test_pointing_limits(self):
        """Test that beamformer.calc_delay respects the pointing limits"""
        
        station = stations.lwa1
        antennas = station.antennas
        
        # Azimuth  checks
        with self.subTest(type='azimuth'):
            self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, -5, 30)
            self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, 365, 30)
            
        # Altitude checks
        with self.subTest(type='altitude'):
            self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, 45, -5)
            self.assertRaises(ValueError, beamformer.calc_delay, antennas[:3], 49.0e6, 45, 95)
            
    def test_phase_and_sum(self):
        """Check that the beamformer.phase_and_sum function actually runs"""
        
        station = stations.lwa1
        antennas = station.antennas
        data = np.random.rand(3, 1000)
        
        with self.subTest(type='time domain'):
            beam = beamformer.phase_and_sum(antennas[:3], data, azimuth=45.0, altitude=30.0)
            
            
        _freqs = [38e6,56e6,64e6,74e6]
        for n in (1,4):
            freqs = [_freqs[i] for i in range(n)]
            data = np.random.rand(3, len(freqs), 1000)
            with self.subTest(type='frequency domain', nchan=n):
                beam = beamformer.phase_and_sum(antennas[:3], data, central_freq=freqs,
                                                azimuth=45.0, altitude=30.0)
                
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
        out = beamformer.phase_beam_shape(ants, azimuth=135.0, altitude=60.0)
        
        i = out.argmax()
        azDiff = np.abs(135.0 - i / 90)
        elDiff = np.abs(60.0 - i % 90)
        self.assertTrue(azDiff <= 1)
        self.assertTrue(elDiff <= 1)
        
    def test_accuracy(self):
        """Check that the extension and numpy.einsum give the same results."""
        
        ninput = 64
        nchan = 16
        ntime = 516
        
        coeffs = np.random.randn(ninput,nchan) + 1j*np.random.randn(ninput,nchan)
        data = np.random.rand(ninput,nchan,ntime) + 1j*np.random.rand(ninput,nchan,ntime) + 5
        
        # Special holders for the ci8 data tests - numpy first then the LSL ci8 type
        data_r = np.clip(data.real*127, -128, 127).astype(np.int8)
        data_i = np.clip(data.imag*127, -128, 127).astype(np.int8)
        n_data = data_r + 1j * data_i
        n_data = n_data.astype(np.complex64)
        
        c_data = np.zeros((ninput,nchan,ntime,2), dtype=np.int8)
        c_data[...,0] = data_r
        c_data[...,1] = data_i
        
        for ctype in (np.complex64, np.complex128):
            t_coeffs = coeffs.astype(ctype)
            
            for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64, np.complex64, np.complex128):
                if np.dtype(dtype).kind == 'c':
                    t_data = data.astype(dtype)
                else:
                    t_data = data.real.astype(dtype)
                    
                with self.subTest(coeffs=ctype, data=dtype):
                    valid0 = np.zeros(ninput, dtype=np.uint8)
                    valid1 = np.zeros(ninput, dtype=np.uint8)
                    valid0[0::2] = 1
                    valid1[1::2] = 1
                    
                    ex_output = _beamformer.BEngine(t_coeffs, t_data, valid0, valid1)
                        
                    np_output = np.zeros((2,nchan,ntime), dtype=np.complex64)
                    np_output[0,...] = np.einsum('sc,sct->ct', t_coeffs[0::2,:], t_data[0::2,:,:], optimize=True) / (ninput // 2)
                    np_output[1,...] = np.einsum('sc,sct->ct', t_coeffs[1::2,:], t_data[1::2,:,:], optimize=True) / (ninput // 2)
                    
                    np.testing.assert_allclose(ex_output, np_output, rtol=1e-5, atol=1e-5)
                    
            with self.subTest(coeffs=ctype, data='ci8'):
                ex_output = _beamformer.BEngine(t_coeffs, c_data, valid0, valid1)
                
                np_output = np.zeros((2,nchan,ntime), dtype=np.complex64)
                np_output[0,...] = np.einsum('sc,sct->ct', t_coeffs[0::2,:], n_data[0::2,:,:], optimize=True) / (ninput // 2)
                np_output[1,...] = np.einsum('sc,sct->ct', t_coeffs[1::2,:], n_data[1::2,:,:], optimize=True) / (ninput // 2)
                
                np.testing.assert_allclose(ex_output, np_output, rtol=1e-5, atol=1e-5)


class  beamformer_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.sim.vis unit
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(beamformer_tests)) 


if __name__ == '__main__':
    unittest.main()
