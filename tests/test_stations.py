# -*- coding: utf-8 -*-

"""Unit test for lsl.common.stations module."""

import unittest

from lsl.common import stations


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class stations_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.stations
	module."""

	def test_station(self):
		"""Test retrieving a stations from the stations module."""

		lwa1 = stations.lwa1()
		self.assertTrue(isinstance(lwa1, stations.LWAStation))

	def test_ecef_conversion(self):
		"""Test the stations.geo2ecef() function."""

		lat = 0.0
		lng = 0.0
		elev = 0.0
		x, y, z = stations.geo2ecef(lat, lng, elev)
		self.assertAlmostEqual(x, 6378137.0)
		self.assertAlmostEqual(y, 0.0)
		self.assertAlmostEqual(z, 0.0)


class stations_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.stations
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(stations_tests)) 


if __name__ == '__main__':
	unittest.main()
