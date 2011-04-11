# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.uvUtils module."""

import unittest
import numpy

from lsl.correlator import uvUtils


__version__  = "0.4"
__revision__ = "$ Revision: 2 $"
__author__    = "Jayce Dowell"


class uvUtils_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.uvUtils
	module."""

	def test_baseline_gen(self):
		"""Test that the generated baselines contain the correct numbers of elements."""

		standList = numpy.array([100, 101, 102, 103])

		bl = uvUtils.getBaselines(standList, IncludeAuto=False, Indicies=False)
		self.assertEqual(len(bl), 6)
		bl = uvUtils.getBaselines(standList, IncludeAuto=True, Indicies=False)
		self.assertEqual(len(bl), 10)

	def test_baseline_ind(self):
		"""Test that the baselines generated with Indicies do return indicies and vice
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
