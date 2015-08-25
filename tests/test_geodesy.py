# -*- coding: utf-8 -*-

"""Unit test for lsl.misc.geodesy module."""

import os
import unittest

from lsl.misc import geodesy


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class geodesy_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.misc.geodesy module."""
	
	def test_read_MAIA(self):
		"""Test that a line from one of the MAIA files can be read."""
		
		ACCURACY = 6
		
		line = '92 1 1 48622.00 I  0.182969 0.000672  0.168817 0.000345  I-0.1251669 0.0000207  1.8335 0.0201  I   -10.437     .507     -.917     .165   .182400   .167900  -.1253000    -9.900    -1.700\n'
		newEOP = geodesy.EOP()
		newEOP.fromMAIA(line)
		
		self.assertAlmostEqual(newEOP.x, 0.182969, ACCURACY)
		self.assertAlmostEqual(newEOP.y, 0.168817, ACCURACY)
		self.assertAlmostEqual(newEOP.utDiff, -0.1251669, ACCURACY)
		self.assertEqual(newEOP.type, 'final')
		
	def test_read_mjd(self):
		"""Test reading a specific MJD via getEOP."""
		
		ACCURACY = 6
		
		# Read in the data for January 2, 1973
		eop = geodesy.getEOP(41684.0)
		
		# Check
		self.assertAlmostEqual(eop.x, 0.120679, ACCURACY)
		self.assertAlmostEqual(eop.y, 0.137008, ACCURACY)
		self.assertAlmostEqual(eop.utDiff,  0.8084136, ACCURACY)
		
		# Read in the data for January 2, 1993
		eops = geodesy.getEOP(48989.0)
		eop = eops
		
		# Check
		self.assertAlmostEqual(eop.x, 0.208317, ACCURACY)
		self.assertAlmostEqual(eop.y, 0.356582, ACCURACY)
		self.assertAlmostEqual(eop.utDiff,  0.0594726, ACCURACY)
		
		# Read in three values and check the output
		mjdList = [41684.0, 41685.0, 41686.0]
		eops = geodesy.getEOP(mjdList)
		self.assertEqual(len(eops), 3)
		for i in range(3):
			self.assertTrue(mjdList[i] in eops)


class geodesy_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.misc.geodesy
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(geodesy_tests)) 


if __name__ == '__main__':
	unittest.main()

