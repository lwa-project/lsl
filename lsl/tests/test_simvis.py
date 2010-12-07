# -*- coding: utf-8 -*-

"""Unit test for the lsl.sim.vis module."""

import os
import unittest
import numpy

from lsl.sim import vis
from lsl.common import stations as lwa_common


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class simvis_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.vis
	module."""

	def test_build_aa_flat(self):
		"""Test building a antenna array object with uniform sky response."""

		lwa1 = lwa_common.lwa1()
		stands = lwa1.getStands()
		freqs = numpy.arange(30e6, 50e6, 1e6)

		aa = vis.buildSimArray(lwa1, stands, freqs, ForceFlat=True)

	def test_build_aa(self):
		"""Test building a antenna array object with realistic sky response."""

		lwa1 = lwa_common.lwa1()
		stands = lwa1.getStands()
		freqs = numpy.arange(30e6, 50e6, 1e6)

		aa = vis.buildSimArray(lwa1, stands, freqs)
		# Check the number of stands
		self.assertEqual(len(aa.ants), len(stands))

		# Check the frequencies comming out
		for fo, fi in zip(aa.get_afreqs(), freqs):
			self.assertAlmostEqual(fo, fi/1e9, 6)
		

class  simvis_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.sim.vis units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(simvis_tests)) 


if __name__ == '__main__':
	unittest.main()
