# -*- coding: utf-8 -*-

"""Unit test for lsl.common.mcs"""

import os
import unittest
from datetime import datetime

from lsl.common import mcs


__revision__ = "$Rev: 839 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class mcs_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.mcs
	module."""
	
	def test_datetime(self):
		"""Test the datetime to MJD, MPM conversion"""
		
		dt = datetime.strptime("2012-06-15 06:34:09", "%Y-%m-%d %H:%M:%S")
		mjd, mpm = mcs.datetime2mjdmpm(dt)
		
		self.assertEqual(mjd, 56093)
		self.assertEqual(mpm, 23649000)
	
	def test_mjdmpm(self):
		"""Test the MJD, MPM to datetime conversion"""
		
		mjd, mpm = 56093, 23649000
		dt = mcs.mjdmpm2datetime(mjd, mpm)
		
		self.assertEqual(dt.strftime("%Y-%m-%d %H:%M:%S"), "2012-06-15 06:34:09")
		
	def test_pointing_correction(self):
		"""Test the pointing correction function"""
		
		az = 63.4
		el = 34.2
		
		# No rotation
		theta = 0.0
		phi = 0.0
		psi = 0.0
		azP, elP = mcs.applyPointingCorrection(az, el, theta, phi, psi, degrees=True)
		self.assertAlmostEqual(azP, az, 1)
		self.assertAlmostEqual(elP, el, 1)
		
		# Azimuth only
		theta = 0.0
		phi = 0.0
		psi = 1.0
		azP, elP = mcs.applyPointingCorrection(az, el, theta, phi, psi, degrees=True)
		self.assertAlmostEqual(azP, az-1.0, 1)
		self.assertAlmostEqual(elP, el, 1)
		
		# Something random
		theta = 23.0
		phi = 10.0
		psi = 1.5
		azP, elP = mcs.applyPointingCorrection(az, el, theta, phi, psi, degrees=True)
		self.assertAlmostEqual(azP, 62.40, 1)
		self.assertAlmostEqual(elP, 34.37, 1)

    
class mcs_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.mcs
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(mcs_tests))        
        
        
if __name__ == '__main__':
	unittest.main()
