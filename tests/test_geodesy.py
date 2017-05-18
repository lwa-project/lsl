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
		
		## http://toshi.nofs.navy.mil/ser7/finals2000A.all
		line = '92 1 1 48622.00 I  0.182985 0.000672  0.168782 0.000345  I-0.1251614 0.0000207  1.8335 0.0201  I    -0.086    0.202     0.130    0.165   .182400   .167900  -.1253000     0.129    -0.653\n'
		newEOP = geodesy.EOP()
		newEOP.fromMAIA(line)
		
		self.assertAlmostEqual(newEOP.x, 0.182985, ACCURACY)
		self.assertAlmostEqual(newEOP.y, 0.168782, ACCURACY)
		self.assertAlmostEqual(newEOP.utDiff, -0.1251614, ACCURACY)
		self.assertEqual(newEOP.type, 'final')
		
	def test_read_mjd(self):
		"""Test reading a specific MJD via getEOP."""
		
		ACCURACY = 6
		
		# Read in the data for January 2, 1973
		eop = geodesy.getEOP(41684.0)
		
		# Check
		## http://toshi.nofs.navy.mil/ser7/finals2000A.all
		## 73 1 2 41684.00 I  0.120727 0.009786  0.136963 0.015902  I 0.8084318 0.0002710  0.0000 0.1916  P    -0.766    0.199    -0.720    0.300   .143000   .137000   .8075000   -18.637    -3.667 
		self.assertAlmostEqual(eop.x, 0.120727, ACCURACY)
		self.assertAlmostEqual(eop.y, 0.136963, ACCURACY)
		self.assertAlmostEqual(eop.utDiff,  0.8084318, ACCURACY)
		
		# Read in the data for January 2, 1993
		eops = geodesy.getEOP(48989.0)
		eop = eops
		
		# Check
		## http://toshi.nofs.navy.mil/ser7/finals2000A.all
		## 93 1 2 48989.00 I  0.208331 0.000166  0.356585 0.000241  I 0.0594794 0.0000108  2.6575 0.0094  I     0.187    0.233     0.020    0.180   .208700   .355900   .0593100    -0.020    -0.087
		self.assertAlmostEqual(eop.x, 0.208331, ACCURACY)
		self.assertAlmostEqual(eop.y, 0.356585, ACCURACY)
		self.assertAlmostEqual(eop.utDiff,  0.0594794, ACCURACY)
		
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

