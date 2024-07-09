"""
Unit test for the lsl.correlator.uvutils module.
"""

import warnings
import unittest
import ephem
import numpy as np

from astropy.coordinates import Angle

from lsl.correlator import uvutils
from lsl.common import stations
import lsl.testing


__version__  = "0.6"
__author__    = "Jayce Dowell"


class uvutils_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.uvutils
    module."""
    
    def setUp(self):
        """Turn off all numpy and python warnings."""

        np.seterr(all='ignore')
        warnings.simplefilter('ignore')

    def test_baseline_gen(self):
        """Test that the generated baselines contain the correct numbers of elements."""

        standList = np.array([100, 101, 102, 103])

        bl = uvutils.get_baselines(standList, include_auto=False)
        self.assertEqual(len(bl), 6)
        bl = uvutils.get_baselines(standList, include_auto=True)
        self.assertEqual(len(bl), 10)
        
    def test_antenna_lookup(self):
        """Test baseline number to antenna lookup function."""
        
        standList = np.array([100, 101, 102, 103])

        bl = uvutils.get_baselines(standList, include_auto=False)
        ind = uvutils.baseline_to_antennas(0, standList)
        self.assertEqual(ind[0], 100)
        self.assertEqual(ind[1], 101)
        
        ind = uvutils.baseline_to_antennas(1, standList, baseline_list=bl)
        self.assertEqual(ind[0], 100)
        self.assertEqual(ind[1], 102)
        
    def test_baseline_lookup(self):
        """Test antennas to baseline lookup function."""
        
        standList = np.array([100, 101, 102, 103])
        bl = uvutils.get_baselines(standList, include_auto=False)
        
        ind = uvutils.antennas_to_baseline(100, 101, standList, include_auto=False)
        self.assertEqual(ind, 0)
        
        ind = uvutils.antennas_to_baseline(100, 102, standList, baseline_list=bl)
        self.assertEqual(ind, 1)
        
        ind = uvutils.antennas_to_baseline(100, 103, standList, include_auto=False)
        self.assertEqual(ind, 2)
        
    def run_compute_uvw_test(self, antennas, freq, HA=0.0, dec=34.0):
        out = uvutils.compute_uvw(antennas, HA=HA, dec=dec, freq=freq)
        
        nbl = len(antennas)*(len(antennas)-1)//2
        try:
            expected_shape = freq.shape
        except AttributeError:
            try:
                expected_shape = (len(freq),)
            except TypeError:
                expected_shape = (1,)
        expected_shape = (nbl,3)+expected_shape
        self.assertEqual(out.shape, expected_shape)
        return out
        
    def test_compute_uvw(self):
        """Test the the compute_uvw function runs."""
        
        station = stations.lwa1
        antennas = station.antennas
        
        # Frequency is a scalar
        with self.subTest(mode='scalar'):
            freq = 45e6
            out = self.run_compute_uvw_test(antennas[0:60:2], freq)
            
        # Frequency is a list
        with self.subTest(mode='list'):
            freq = [45e6, 60e6]
            out = self.run_compute_uvw_test(antennas[0:60:2], freq)

        # Frequency is an array
        ## 1-D
        with self.subTest(mode='1D'):
            freq = np.linspace(45e6, 60e6, 1024)
            out0 = self.run_compute_uvw_test(antennas[0:60:2], freq)
            
        ## 2-D
        with self.subTest(mode='2D'):
            freq.shape = (512, 2)
            out1 = self.run_compute_uvw_test(antennas[0:60:2], freq)
        
        ## 3-D
        with self.subTest(mode='3D'):
            freq.shape = (128, 4, 2)
            out2 = self.run_compute_uvw_test(antennas[0:60:2], freq)
        
        # Make sure the values are the same
        out1.shape = out0.shape
        out2.shape = out0.shape
        np.testing.assert_allclose(out0, out1)
        np.testing.assert_allclose(out0, out2)
        
        # Pass in a list of baselines
        bl = uvutils.get_baselines(antennas[0:30:2], include_auto=False)
        out3 = uvutils.compute_uvw(bl, freq=freq)
        self.assertEqual(len(out3), len(bl))
        
        # Angle support
        freq = [45e6, 60e6]
        HA = 1.0
        dec = 34.0
        out0 = self.run_compute_uvw_test(antennas[0:60:2], freq, HA=HA, dec=dec)
        
        HA = Angle('1h0m0s')
        dec = Angle('34deg')
        out1 = self.run_compute_uvw_test(antennas[0:60:2], freq, HA=HA, dec=dec)
        np.testing.assert_allclose(out0, out1)
        
        HA = ephem.hours('1:00:00')
        dec = ephem.degrees('34:00:00')
        out1 = self.run_compute_uvw_test(antennas[0:60:2], freq, HA=HA, dec=dec)
        np.testing.assert_allclose(out0, out1)
        
    def test_compute_uv_track(self):
        """Test that the compute_uv_track function runs."""
        
        station = stations.lwa1
        antennas = station.antennas
        
        out = uvutils.compute_uv_track(antennas[0:60:2])
        
        # Make sure we have the right dimensions
        self.assertEqual(out.shape, (435,2,512))
        
        # Pass in a list of baselines
        bl = uvutils.get_baselines(antennas[0:30:2], include_auto=False)
        out0 = uvutils.compute_uv_track(bl)
        self.assertEqual(len(out0), len(bl))
        
        # Angle support
        dec = 34.0
        out0 = uvutils.compute_uv_track(antennas[0:60:2], dec=dec)
        
        dec = Angle('34deg')
        out1 = uvutils.compute_uv_track(antennas[0:60:2], dec=dec)
        np.testing.assert_allclose(out0, out1)
        
        dec = ephem.degrees('34:00:00')
        out1 = uvutils.compute_uv_track(antennas[0:60:2], dec=dec)
        np.testing.assert_allclose(out0, out1)

        
class uvutils_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(uvutils_tests)) 


if __name__ == '__main__':
    unittest.main()
