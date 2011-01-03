# -*- coding: utf-8 -*-

"""Unit test for the lsl.statistics.robust module."""

import os
import time
import unittest
import numpy

from lsl.statistics import robust


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"

class robust_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.statistics.robust
	module."""

	def test_mean(self):
		"""Test the outlier-resistant mean function."""

		# Make sure that it can do simple averages
		a = numpy.random.rand(512)
		self.assertAlmostEqual(a.mean(), robust.mean(a, Cut=100), 6)

		b = numpy.random.randn(512)
		self.assertAlmostEqual(b.mean(), robust.mean(b, Cut=100), 6)

		# Make sure it can reject obvious points
		b = 1.0*a
		b[10] = 1e6
		self.assertTrue(robust.mean(b) < b.mean())

	def test_std(self):
		"""Test the outlier-resistant standard deviation function."""

		# Make sure it can reject obvious points
		b = numpy.random.randn(512)
		b[10] = 1e6
		self.assertTrue(robust.std(b) < b.std())

	def test_biweight(self):
		"""Test the outlier-resistant biweighted mean function."""

		# Make sure that it can do simple averages
                a = numpy.random.rand(512)
                self.assertAlmostEqual(a.mean(), robust.biweightMean(a, Cut=100), 6)

                b = numpy.random.randn(512)
                self.assertAlmostEqual(b.mean(), robust.biweightMean(b, Cut=100), 6)

                # Make sure it can reject obvious points
                b = 1.0*a
                b[10] = 1e6
                self.assertTrue(robust.biweightMean(b) < b.mean())


class robust_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.statistics.robust 
	units tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(robust_tests)) 


if __name__ == '__main__':
	unittest.main()
