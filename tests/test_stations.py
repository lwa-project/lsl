# -*- coding: utf-8 -*-

"""Unit test for lsl.common.stations module."""

import os
import ephem
import pickle
import unittest
from datetime import datetime

from lsl.common.paths import dataBuild as dataPath
from lsl.common import stations


__revision__ = "$Rev$"
__version__  = "0.3"
__author__    = "Jayce Dowell"


class stations_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.stations
	module."""

	def test_station(self):
		"""Test retrieving a stations from the stations module."""

		lwa1 = stations.lwa1
		self.assertTrue(isinstance(lwa1, stations.LWAStation))
		
	def test_observer(self):
		"""Test the ephem.Observer portion of an LWAStation."""
		
		lwa1 = stations.lwa1
		jov = ephem.Jupiter()
		
		lwa1.date = '2013/7/10 22:07:07'
		lwa1.compute(jov)
		
		# RA/Dec
		self.assertAlmostEqual(jov.ra,  ephem.hours('6:14:41.01'))
		self.assertAlmostEqual(jov.dec, ephem.degrees('23:11:49.1'))
		
		#Az/Alt
		self.assertAlmostEqual(jov.az,  ephem.degrees('274:40:25.3'))
		self.assertAlmostEqual(jov.alt, ephem.degrees('37:24:09.8'))
		
	def test_pickle(self):
		"""Test pickling of LWAStation instances."""
		
		lwa1 = stations.lwa1
		
		# Pickle and re-load
		out  = pickle.dumps(lwa1)
		lwa1Prime = pickle.loads(out)
		
		# Test similarity
		self.assertAlmostEqual(lwa1.lat, lwa1Prime.lat)
		self.assertAlmostEqual(lwa1.long, lwa1Prime.long)
		self.assertAlmostEqual(lwa1.elev, lwa1Prime.elev)
		for i in xrange(520):
			self.assertEqual(lwa1.antennas[i].id, lwa1Prime.antennas[i].id)
			self.assertEqual(lwa1.antennas[i].stand.id, lwa1Prime.antennas[i].stand.id)
			self.assertEqual(lwa1.antennas[i].digitizer, lwa1Prime.antennas[i].digitizer)
			
		# Check independence
		lwa1Prime.antennas[100].stand.id = 888
		self.assertTrue(lwa1.antennas[100].stand.id != lwa1Prime.antennas[100].stand.id)
		
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
	
	def test_ssmif_binary(self):
		"""Test the binary SSMIF parser."""
		
		ssmifFile = os.path.join(dataPath, 'tests', 'ssmif.dat')
		out = stations.parseSSMIF(ssmifFile)


class stations_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.stations
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(stations_tests)) 


if __name__ == '__main__':
	unittest.main()
