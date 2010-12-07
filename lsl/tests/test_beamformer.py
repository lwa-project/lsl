# -*- coding: utf-8 -*-

"""Unit test for the lsl.misc.beamformer module."""

import os
import unittest
import numpy

from lsl.misc import beamformer
from lsl.common import stations as lwa_common


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class beamformer_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.dp
	module."""

	def test_calcDelay(self):
		"""Check that the beamformer.calcDelay function actually runs"""

		out = beamformer.calcDelay(numpy.array([1,2,3]))
		self.assertEqual(len(out), 3)
		
		out = beamformer.calcDelay(numpy.array([1,2,3]), freq=49.0e6)
		self.assertEqual(len(out), 3)

		out = beamformer.calcDelay(numpy.array([1,2,3]), freq=49.0e6, azimuth=45, elevation=30)
		self.assertEqual(len(out), 3)

	def test_pointing_limits(self):
		"""Test that beamformer.calcDelay respects the pointing limits"""

		# Azimuth  checks
		self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, numpy.array([1,2,3]), 49.0e6, -5, 30)
		self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, numpy.array([1,2,3]), 49.0e6, 365, 30)

		# Elevation checks
		self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, numpy.array([1,2,3]), 49.0e6, 45, -5)
		self.assertRaises(beamformer.BeamformingError, beamformer.calcDelay, numpy.array([1,2,3]), 49.0e6, 45, 95)

	def test_intDelayAndSum(self):
		"""Check that the beamformer.intDelayAndSum function actually runs"""
		
		stands = numpy.array([1,2,3])
		data = numpy.random.rand(3, 1000)
		
		beam = beamformer.intDelayAndSum(stands, data, azimuth=45.0, elevation=30.0)

	def test_intBeamShape(self):
		"""Check that the beamformer.intBeamShape function actually runs"""

		lwa1 = lwa_common.lwa1()
		stands = lwa1.getStands()

		out = beamformer.intBeamShape(stands, azimuth=135.0, elevation=60.0)
		i = out.argmax()
		az = i / 90
		el = i % 90
		self.assertAlmostEqual(az, 135.0, 1)
		self.assertAlmostEqual(el, 60.0, 1)
		

class  beamformer_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.sim.vis units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(beamformer_tests)) 


if __name__ == '__main__':
	unittest.main()
