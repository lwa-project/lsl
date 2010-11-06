# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.visUtils module."""

import unittest
import numpy

from lsl.correlator import visUtils


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class visUtils_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.visUtils
	module."""
	
	def test_argument(self):
		"""Test the argument function."""

		cValue = 5*numpy.exp(-2j*numpy.pi*0.455)
		arg = visUtils.argument(cValue)
		self.assertAlmostEqual(arg, -2*numpy.pi*0.455)

	def test_unmask_data(self):
		"""Test teh unmaskCompressedData function."""

		data = numpy.random.rand(100)
		mask = numpy.zeros_like(data)
		for i in range(0,100,17):
			mask[i] = 0
		maskedData = 

class visUtils_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.correlator.visUtils units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(visUtils_tests)) 


if __name__ == '__main__':
	unittest.main()
