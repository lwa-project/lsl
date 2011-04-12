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

	def setUp(self):
		"""Turn off all numpy warnings."""

		numpy.seterr(all='ignore')

	def test_build_aa_flat(self):
		"""Test building a antenna array object with uniform sky response."""

		lwa1 = lwa_common.lwa1
		antennas = lwa1.getAntennas()[0:20]
		freqs = numpy.arange(30e6, 50e6, 1e6)

		aa = vis.buildSimArray(lwa1, antennas, freqs, ForceFlat=True)

	def test_build_aa(self):
		"""Test building a antenna array object with realistic sky response."""

		lwa1 = lwa_common.lwa1
		antennas = lwa1.getAntennas()[0:20]
		freqs = numpy.arange(30e6, 50e6, 1e6)

		aa = vis.buildSimArray(lwa1, antennas, freqs)
		# Check the number of stands
		self.assertEqual(len(aa.ants), len(antennas))

		# Check the frequencies comming out
		for fo, fi in zip(aa.get_afreqs(), freqs):
			self.assertAlmostEqual(fo, fi/1e9, 6)

	def test_build_data(self):
		"""Test building simulated visibility data"""

		# Setup
		lwa1 = lwa_common.lwa1
		antennas = lwa1.getAntennas()[0:20]
		freqs = numpy.arange(30e6, 50e6, 1e6)
		aa = vis.buildSimArray(lwa1, antennas, freqs)

		# Build the data dictionary
		out = vis.buildSimData(aa, vis.srcs)

		# Do a check of keys
		keyList = out.keys()
		for key in ['freq', 'isMasked', 'bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			self.assertTrue(key in keyList)

		# Do a check of frequencies
		for fa, fq in zip(out['freq'], freqs):
			self.assertAlmostEqual(fa, fq, 6)

		# Do a check to make sure that the entries with secondary keys have them
		for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			secondaryKeyList = out[key].keys()
			for key2 in ['xx', 'yy', 'xy', 'yx']:
				self.assertTrue(key2 in secondaryKeyList)

		# Do a check to make sure that the entries with secondary keys also 
		# have lists in them
		for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			secondaryKeyList = out[key].keys()
			for key2 in ['xx', 'yy', 'xy', 'yx']:
				self.assertTrue(type(out[key][key2]).__name__ == 'list')
		
	def test_scale_data(self):
		"""Test that we can scale a data dictionary without error"""

		# Setup
		lwa1 = lwa_common.lwa1
		antennas = lwa1.getAntennas()[0:20]
		freqs = numpy.arange(30e6, 50e6, 1e6)
		aa = vis.buildSimArray(lwa1, antennas, freqs)

		# Build the data dictionary
		out = vis.buildSimData(aa, vis.srcs)

		# Scale
		amp = vis.scaleData(out, numpy.ones(len(antennas))*2, numpy.zeros(len(antennas)))
		# Delay
		phs = vis.scaleData(out, numpy.ones(len(antennas)), numpy.ones(len(antennas)))

	def test_shift_data(self):
		"""Test that we can shift the uvw coordinates of a data dictionary 
		without error"""

		# Setup
		lwa1 = lwa_common.lwa1
		antennas = lwa1.getAntennas()[0:20]
		freqs = numpy.arange(30e6, 50e6, 1e6)
		aa = vis.buildSimArray(lwa1, antennas, freqs)

		# Build the data dictionary
		out = vis.buildSimData(aa, vis.srcs)

		# Shift
		sft = vis.shiftData(out, aa)


class  simvis_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.sim.vis units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)

		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(simvis_tests)) 


if __name__ == '__main__':
	unittest.main()
