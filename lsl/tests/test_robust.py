# -*- coding: utf-8 -*-

"""Unit test for the lsl.statistics.robust module."""

import os
import time
import unittest
import numpy

from lsl.statistics import robust


__revision__ = "$Revision: 4 $"
__version__  = "0.2"
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

	def test_linefit(self):
		"""Test the outlier-resistant line fitter."""

		x = numpy.linspace(-10, 10, 256)
		y = 1.34*x + 0.56

		cc = robust.linefit(x, y)
		self.assertAlmostEqual(cc[0], 1.34, 2)
		self.assertAlmostEqual(cc[1], 0.56, 2)

		cc = robust.linefit(x, y, Bisector=True)
		self.assertAlmostEqual(cc[0], 1.34, 2)
		self.assertAlmostEqual(cc[1], 0.56, 2)

		x = numpy.linspace(-10, 10, 2048)
		y = 2.86*x - 0.56

		cc = robust.linefit(x, y)
		self.assertAlmostEqual(cc[0], 2.86,  2)
		self.assertAlmostEqual(cc[1], -0.56, 2)

		cc = robust.linefit(x, y, Bisector=True)
		self.assertAlmostEqual(cc[0], 2.86,  2)
		self.assertAlmostEqual(cc[1], -0.56, 2)

	def test_polyfit(self):
		"""Test the outlier-resistant polynomial fitter."""

		x = numpy.linspace(-10, 10, 256)
		y = 0.012*x**3 + 0.34*x**2 + 1.34*x + 0.56

		cc = robust.polyfit(x, y, 3)
		self.assertAlmostEqual(cc[0], 0.012, 3)
		self.assertAlmostEqual(cc[1], 0.340, 3)
		self.assertAlmostEqual(cc[2], 1.340, 3)
		self.assertAlmostEqual(cc[3], 0.560, 3)

		x = numpy.linspace(-10, 10, 2048)
		y = 0.003*x**4 + 0.012*x**3 + 0.34*x**2 + 1.34*x + 0.56

		cc = robust.polyfit(x, y, 4)
		self.assertAlmostEqual(cc[0], 0.003, 3)
		self.assertAlmostEqual(cc[1], 0.012, 3)
		self.assertAlmostEqual(cc[2], 0.340, 3)
		self.assertAlmostEqual(cc[3], 1.340, 3)
		self.assertAlmostEqual(cc[4], 0.560, 3)


class robust_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.statistics.robust 
	units tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(robust_tests)) 


if __name__ == '__main__':
	unittest.main()
