# -*- coding: utf-8 -*-

"""Unit test for lsl.sim.nec_util module."""

import unittest
import logging
import os
import StringIO

from lsl.sim import nec_util
from lsl.common.paths import dataBuild as dataPath


__revision__  = "$Revision$"
__version__   = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class nec_util_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.nec_util
	module."""
	
	def setUp(self):
		"""Setup unit tests."""
		
		# disable logger since we intentionally generate errors 
		logging.basicConfig(stream = StringIO.StringIO())
		
		# get a reference to the input test file
		self.nec_name = os.path.join(dataPath, 'tests', 'bigblade_imp.out')
	
	def test_NECImpedance_init(self):
		"""Test nec_util.NECImpedance constructor method."""
		
		imp = nec_util.NECImpedance(self.nec_name)
		
	def test_open_and_get_nec_freq(self):
		"""Test nec_util.open_and_get_nec_freq() function."""
		
		(f, freq) = nec_util.open_and_get_nec_freq(self.nec_name)   
		
	def test_calcIME(self):
		"""Test nec_util.calcIME() function."""
		
		(freqs, ime) = nec_util.calcIME(self.nec_name)
	
	def test_NECPattern_init(self):
		"""Test nec_util.NECPattern constructor method."""
		
		pat = nec_util.NECPattern(self.nec_name, 5.0)
		self.assertRaises(ValueError, nec_util.NECPattern, self.nec_name, 0.0, False)
   
    
class nec_util_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lwa_user.nec_util
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(nec_util_tests))        
        
        
if __name__ == '__main__':
	unittest.main()
 
