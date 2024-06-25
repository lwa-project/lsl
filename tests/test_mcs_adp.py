"""
Unit test for the lsl.common.mcsADP module.
"""

import os
import unittest
from datetime import datetime

from lsl.common import mcsADP


__version__  = "0.1"
__author__    = "Jayce Dowell"


class mcs_adp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.mcs
    module."""
    
    def test_delay_conversion(self):
        """Test the MCS delay conversion"""
        
        delay_in = 5.0
        delay_out = mcsADP.mcsd_to_delay(mcsADP.delay_to_mcsd(delay_in))
        self.assertTrue(abs(delay_in-delay_out) < 1e9/196e6/16) # Within 1/16 of a sample
        
    def test_gain_conversion(self):
        """Test the MCS gain conversion"""
        
        gain_in = 0.5
        gain_out = mcsADP.mcsg_to_gain(mcsADP.gain_to_mcsg(0.5))
        self.assertTrue(abs(gain_in-gain_out) < 1/2.**15)
        
    def test_datetime(self):
        """Test the datetime to MJD, MPM conversion"""
        
        dt = datetime.strptime("2012-06-15 06:34:09", "%Y-%m-%d %H:%M:%S")
        mjd, mpm = mcsADP.datetime_to_mjdmpm(dt)
        
        self.assertEqual(mjd, 56093)
        self.assertEqual(mpm, 23649000)
    
    def test_mjdmpm(self):
        """Test the MJD, MPM to datetime conversion"""
        
        mjd, mpm = 56093, 23649000
        dt = mcsADP.mjdmpm_to_datetime(mjd, mpm)
        
        self.assertEqual(dt.strftime("%Y-%m-%d %H:%M:%S"), "2012-06-15 06:34:09")
        
    def test_summary_limits(self):
        """Test valid summary values"""
        
        for i in range(0, 6+1):
            mcsADP.summary_to_string(i)
        self.assertRaises(ValueError, mcsADP.summary_to_string, 7)
        
    def test_sid_limits(self):
        """Test valid subsystem ID values"""
        
        for i in range(1, 19+1):
            mcsADP.sid_to_string(i)
        self.assertRaises(ValueError, mcsADP.sid_to_string, 0)
        self.assertRaises(ValueError, mcsADP.sid_to_string, 20)
        
    def test_cid_limits(self):
        """Test valid command ID values"""
        
        for i in range(0, 41+1):
            mcsADP.cid_to_string(i)
        self.assertRaises(ValueError, mcsADP.cid_to_string, 42)
        
    def test_mode_limits(self):
        """Test valid observing mode values"""
        
        for i in range(1, 9+1):
            mcsADP.mode_to_string(i)
        self.assertRaises(ValueError, mcsADP.mode_to_string, 0)
        self.assertRaises(ValueError, mcsADP.mode_to_string, 10)
        
    def test_pointing_correction(self):
        """Test the pointing correction function"""
        
        az = 63.4
        el = 34.2
        
        # No rotation
        theta = 0.0
        phi = 0.0
        psi = 0.0
        azP, elP = mcsADP.apply_pointing_correction(az, el, theta, phi, psi, degrees=True)
        self.assertAlmostEqual(azP, az, 1)
        self.assertAlmostEqual(elP, el, 1)
        
        # Azimuth only
        theta = 0.0
        phi = 0.0
        psi = 1.0
        azP, elP = mcsADP.apply_pointing_correction(az, el, theta, phi, psi, degrees=True)
        self.assertAlmostEqual(azP, az-1.0, 1)
        self.assertAlmostEqual(elP, el, 1)
        
        # Something random
        theta = 23.0
        phi = 10.0
        psi = 1.5
        azP, elP = mcsADP.apply_pointing_correction(az, el, theta, phi, psi, degrees=True)
        self.assertAlmostEqual(azP, 62.40, 1)
        self.assertAlmostEqual(elP, 34.37, 1)

    
class mcs_adp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.mcsADP
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(mcs_adp_tests))        
        
        
if __name__ == '__main__':
    unittest.main()
