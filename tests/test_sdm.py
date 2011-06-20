# -*- coding: utf-8 -*-

"""
Unit test for the lsl.common.sdf module.
"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.common import sdm


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


# This is pretty fake down here
sdmFile = os.path.join(dataPath, 'tests', 'bigblade_imp.out')


class sdm_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.sdm
	module."""

	def test_sdm_parse(self):
		"""Test that the SDM parser runs."""
		
		junk = sdm.parseSDM(sdmFile)


class sdm_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.sdm units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(sdm_tests)) 


if __name__ == '__main__':
	unittest.main()