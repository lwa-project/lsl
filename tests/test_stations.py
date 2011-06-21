# -*- coding: utf-8 -*-

"""Unit test for lsl.common.stations module."""

import os
import unittest
from datetime import datetime

from lsl.common.paths import dataBuild as dataPath
from lsl.common import stations


__revision__ = "$ Revision: 2 $"
__version__  = "0.2"
__author__    = "Jayce Dowell"


class stations_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.stations
	module."""

	def test_station(self):
		"""Test retrieving a stations from the stations module."""

		lwa1 = stations.lwa1
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
		
	def test_prototype(self):
		"""Test retrieving a PrototypeStation from the stations module."""
		
		proto = stations.prototypeSystem
		self.assertTrue(isinstance(proto, stations.PrototypeStation))
		
	def test_prototype_ants(self):
		"""Test retrieving antennas from a prototype system."""
		
		proto = stations.prototypeSystem
		
		# Check that we get the right number of antennas for the system
		ants = proto.getAntennas(datetime(2011, 4, 4, 0, 0, 0))
		self.assertEqual(len(ants), 20)
		
		# Again
		ants = proto.getAntennas(datetime(2011, 1, 1, 0, 0, 0))
		self.assertEqual(len(ants), 20)
		
		# And check that we actually get out what we need in the right order
		antExpected = [(206,0), (183,0), (153,0), (174,0), (38,0), (34,0), (67,0), (181,0), (80,0), (14,0), 
					(254,0), (118,0), (246,0), (9,0), (69,0), (168,0), (258,0), (4,0), (158,0), (205,0)]
		for i in xrange(len(ants)):
			pair = (ants[i].stand.id, ants[i].pol)
			self.assertEqual(pair, antExpected[i])
	
	def test_ssmif_text(self):
		"""Test the text SSMIF parser."""
		
		ssmifFile = os.path.join(dataPath, 'lwa1-ssmif.txt')
		out = stations.parseSSMIF(ssmifFile)
	
	#def test_ssmif_binary(self):
		#"""Test the binary SSMIF parser."""
		
		## This is pretty fake down here.  We need a *real* binary SSMIF file to do this
		#ssmifFile = os.path.join(dataPath, 'lwa1-ssmif.txt')
		#out = stations.parseSSMIF(ssmifFile)


class stations_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.stations
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(stations_tests)) 


if __name__ == '__main__':
	unittest.main()
