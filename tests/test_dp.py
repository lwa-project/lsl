# -*- coding: utf-8 -*-

"""Unit test for lsl.common.dp module."""

import warnings
import unittest
import numpy

from lsl.common import dp


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class dp_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.dp
	module."""
	
	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')
	
	def test_tbn_bandpass(self):
		"""Test that the TBN bandpass generator actually runs."""
		
		fnc = dp.tbnFilter(sampleRate=1e5)
		junk = fnc(1e3)

	def test_drx_bandpass(self):
		"""Test that the DRX bandpass generator actually runs."""
		
		fnc = dp.drxFilter(sampleRate=19.6e6)
		junk = fnc(1e3)


class dp_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.dp
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(dp_tests)) 


if __name__ == '__main__':
	unittest.main()
