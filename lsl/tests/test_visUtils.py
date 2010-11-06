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
		"""Test the unmaskCompressedData function."""

		data = numpy.random.rand(100)
		mask1 = numpy.zeros_like(data)
		mask2 = numpy.zeros_like(data)
		for i in range(0,100,17):
			mask1[i] = 0
		maskedData1 = numpy.ma.array(data, mask=mask1)	# every 17th value masked
		maskedData1 = maskedData1.compressed()
		maskedData2 = numpy.ma.array(data, mask=mask2)	# none masked
		maskedData2 = maskedData2.compressed()
		
		unmaskedData1 = visUtils.unmaskCompressedData(maskedData1, mask1)
		unmaskedData2 = visUtils.unmaskCompressedData(maskedData2, mask2)
		
		# Make sure both arrays are unmasked to the correct size
		self.assertEqual(len(data), len(unmaskedData1))
		self.assertEqual(len(data), len(unmaskedData2))
		
	def test_unwrap_phase(self):
		"""Test the unwrap function."""
		
		# Range from 0 to 10 pi
		x = numpy.linspace(0, 10*numpy.pi, num=50)
		# Wrap
		x = x % (2*numpy.pi)
		# Unwrap
		y = visUtils.unwrap(x)
		
		self.assertAlmostEqual(y[-1], 10*numpy.pi)
		
		
class visUtils_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.correlator.visUtils units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(visUtils_tests)) 


if __name__ == '__main__':
	unittest.main()
