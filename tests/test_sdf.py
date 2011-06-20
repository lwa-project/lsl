# -*- coding: utf-8 -*-

"""
Unit test for the lsl.common.sdf module.
"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.common import sdf


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(dataPath, 'tests', 'tbw-sdf.txt')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-sdf.txt')
drxFile = os.path.join(dataPath, 'tests', 'drx-sdf.txt')


class sdf_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.sdf
	module."""

	### TBW ###

	def test_tbw_parse(self):
		"""Test reading in a TBW SDF file."""
		
		project = sdf.parseSDF(tbwFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode, 'TBW')
		self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[0].mpm,      0)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TBW')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  70000)
	
	def test_tbw_write(self):
		"""Test writing a TBW SDF file."""
		
		project = sdf.parseSDF(tbwFile)
		out = project.render()
	
	def test_tbw_errors(self):
		"""Test various TBW SDF errors."""
		
		project = sdf.parseSDF(tbwFile)
		
		# Bad number of TBW bits
		project.sessions[0].observations[0].bits = 6
		self.assertFalse(project.validate())
		
		# Bad number of TBW samples
		project.sessions[0].observations[0].bits = 4
		project.sessions[0].observations[0].samples = 72000000
		self.assertFalse(project.validate())
		
		project.sessions[0].observations[0].bits = 12
		project.sessions[0].observations[0].samples = 72000000
		self.assertFalse(project.validate())
	
	### TBN ###
	
	def test_tbn_parse(self):
		"""Test reading in a TBN SDF file."""
		
		project = sdf.parseSDF(tbnFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode, 'TBN')
		self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[0].mpm,      0)
		self.assertEqual(project.sessions[0].observations[0].dur,  10000)
		self.assertEqual(project.sessions[0].observations[0].freq1, 438261968)
		self.assertEqual(project.sessions[0].observations[0].filter, 7)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TBN')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
		self.assertEqual(project.sessions[0].observations[1].dur,  10000)
		self.assertEqual(project.sessions[0].observations[1].freq1, 832697741)
		self.assertEqual(project.sessions[0].observations[1].filter, 7)
	
	def test_tbn_write(self):
		"""Test writing a TBN SDF file."""
		
		project = sdf.parseSDF(tbnFile)
		out = project.render()
	
	def test_tbn_errors(self):
		"""Test various TBN SDF errors."""
		
		project = sdf.parseSDF(tbnFile)
		
		# Bad filter
		project.sessions[0].observations[0].filter = 10
		self.assertFalse(project.validate())
		
		# Bad frequency
		project.sessions[0].observations[0].filter = 7
		project.sessions[0].observations[0].frequency1 = 90.0e6
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		# Bad duration
		project.sessions[0].observations[0].frequency1 = 38.0e6
		project.sessions[0].observations[0].duration = '96:00:00.000'
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
	
	### DRX ###
	
	def test_drx_parse(self):
		"""Test reading in a DRX SDF file."""
		
		project = sdf.parseSDF(drxFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode, 'TRK_RADEC')
		self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[0].mpm,      0)
		self.assertEqual(project.sessions[0].observations[0].dur,  10000)
		self.assertEqual(project.sessions[0].observations[0].freq1,  438261968)
		self.assertEqual(project.sessions[0].observations[0].freq2, 1928352663)
		self.assertEqual(project.sessions[0].observations[0].filter, 7)
		self.assertAlmostEqual(project.sessions[0].observations[0].ra, 5.6, 6)
		self.assertAlmostEqual(project.sessions[0].observations[0].dec, 22.0, 6)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_RADEC')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
		self.assertEqual(project.sessions[0].observations[1].dur,  10000)
		self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
		self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
		self.assertEqual(project.sessions[0].observations[1].filter, 7)
		self.assertAlmostEqual(project.sessions[0].observations[1].ra, 5.6, 6)
		self.assertAlmostEqual(project.sessions[0].observations[1].dec, 22.0, 6)
		
	def test_drx_write(self):
		"""Test writing a DRX SDF file."""
		
		project = sdf.parseSDF(drxFile)
		out = project.render()
		
	def test_drx_errors(self):
		"""Test various DRX SDF errors."""
		
		project = sdf.parseSDF(drxFile)
		
		# Bad filter
		project.sessions[0].observations[0].filter = 10
		self.assertFalse(project.validate())
		
		# Bad frequency
		project.sessions[0].observations[0].filter = 7
		project.sessions[0].observations[0].frequency1 = 90.0e6
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		project.sessions[0].observations[0].frequency1 = 38.0e6
		project.sessions[0].observations[0].frequency2 = 90.0e6
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		# Bad duration
		project.sessions[0].observations[0].frequency2 = 38.0e6
		project.sessions[0].observations[0].duration = '96:00:00.000'
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		# Bad pointing
		project.sessions[0].observations[0].duration = '00:00:01.000'
		project.sessions[0].observations[0].dec = -72.0
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())


class sdf_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.sdf units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(sdf_tests)) 


if __name__ == '__main__':
	unittest.main()
