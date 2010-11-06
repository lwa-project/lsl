# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.uvUtils module."""

import unittest
import numpy

from lsl.correlator import uvUtils


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class uvUtils_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.uvUtils
	module."""

	def test_position(self):
		"""Test positions directly and via the position cache object."""

		# One stand
		xyz = uvUtils.getXYZ(1)
		self.assertAlmostEqual(xyz[0,0], -0.67)
		self.assertAlmostEqual(xyz[0,1], -54.63)
		self.assertAlmostEqual(xyz[0,2], 1.63)

		cache = uvUtils.PositionCache()
		xyz = cache.getXYZ(1)
		self.assertAlmostEqual(xyz[0,0], -0.67)
		self.assertAlmostEqual(xyz[0,1], -54.63)
		self.assertAlmostEqual(xyz[0,2], 1.63)

	def test_relative(self):
		"""Test baseline positions directly and via the position cache object."""

		x, y, z = uvUtils.getRelativeXYZ(1, 1)
		self.assertAlmostEqual(x, 0)
		self.assertAlmostEqual(y, 0)
		self.assertAlmostEqual(z, 0)

		x, y, z = uvUtils.getRelativeXYZ(1, 2)
		self.assertAlmostEqual(x, 1.21)
		self.assertAlmostEqual(y, 5.30)
		self.assertAlmostEqual(z, -0.13)

		cache = uvUtils.PositionCache()
		x, y, z = cache.getRelativeXYZ(1, 1)
		self.assertAlmostEqual(x, 0)
		self.assertAlmostEqual(y, 0)
		self.assertAlmostEqual(z, 0)

		x, y, z = cache.getRelativeXYZ(1, 2)
		self.assertAlmostEqual(x, 1.21)
		self.assertAlmostEqual(y, 5.30)
		self.assertAlmostEqual(z, -0.13)

	def test_position_bounds(self):
		"""Test that stand ID number bounds are respected by the position routines."""

		self.assertRaises(uvUtils.uvUtilsError, uvUtils.getXYZ, -10)
		self.assertRaises(uvUtils.uvUtilsError, uvUtils.getXYZ, 0)
		self.assertRaises(uvUtils.uvUtilsError, uvUtils.getXYZ, 500)

		cache = uvUtils.PositionCache()
		self.assertRaises(uvUtils.uvUtilsError, cache.getXYZ, -10)
		self.assertRaises(uvUtils.uvUtilsError, cache.getXYZ, 0)
		self.assertRaises(uvUtils.uvUtilsError, cache.getXYZ, 500)

	def test_cable(self):
		"""Test cable delays directly and via the cable cahce object."""

		d = uvUtils.cableDelay(1, 10e6)
		self.assertAlmostEqual(d, 557.913e-9)

		cache = uvUtils.CableCache(10e6)
		d = cache.cableDelay(1)
		self.assertAlmostEqual(d, 557.913e-9)

	def test_cable_bounds(self):
		"""Test that stand ID number bound are respected by the cable routines."""

		self.assertRaises(uvUtils.uvUtilsError, uvUtils.cableDelay, -10, 10e6)
		self.assertRaises(uvUtils.uvUtilsError, uvUtils.cableDelay, 0, 10e6)
		self.assertRaises(uvUtils.uvUtilsError, uvUtils.cableDelay, 500, 10e6)

		cache = uvUtils.CableCache(10e6)
		self.assertRaises(uvUtils.uvUtilsError, cache.cableDelay, -10)
		self.assertRaises(uvUtils.uvUtilsError, cache.cableDelay, 0)
		self.assertRaises(uvUtils.uvUtilsError, cache.cableDelay, 500)

	def test_attenuation(self):
		"""Test that the cable attenuation is working."""

		a = uvUtils.cableAttenuation(1)
		self.assertAlmostEqual(a, 4.4524597279)

	def test_attenuation_bounds(self):
		"""Test that stand ID number bound are respected by the attenuation routine."""

		self.assertRaises(uvUtils.uvUtilsError, uvUtils.cableAttenuation, -10)
		self.assertRaises(uvUtils.uvUtilsError, uvUtils.cableAttenuation, 0)
		self.assertRaises(uvUtils.uvUtilsError, uvUtils.cableAttenuation, 500)

	def test_baseline_gen(self):
		"""Test that the generated baselines contian the correct numbers of elements."""

		standList = numpy.array([100, 101, 102, 103])

		bl = uvUtils.getBaselines(standList, IncludeAuto=False, Indicies=False)
		self.assertEqual(len(bl), 6)
		bl = uvUtils.getBaselines(standList, IncludeAuto=True, Indicies=False)
		self.assertEqual(len(bl), 10)

	def test_baseline_ind(self):
		"""Test that the baselines generated with Incidies do return indicies and vice
		versa."""

		standList = numpy.array([100, 101, 102, 103])

		bl = uvUtils.getBaselines(standList, IncludeAuto=False, Indicies=False)
		bl = numpy.array(bl)
		self.assertTrue(bl.min() == 100)
		bl = uvUtils.getBaselines(standList, IncludeAuto=True, Indicies=False)
		bl = numpy.array(bl)
		self.assertTrue(bl.min() == 100)

		bl = uvUtils.getBaselines(standList, IncludeAuto=False, Indicies=True)
		bl = numpy.array(bl)
		self.assertTrue(bl.max() < 100)
		bl = uvUtils.getBaselines(standList, IncludeAuto=True, Indicies=True)
		bl = numpy.array(bl)
		self.assertTrue(bl.max() < 100)
		
class uvUtils_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(uvUtils_tests)) 


if __name__ == '__main__':
	unittest.main()
