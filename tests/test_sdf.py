# -*- coding: utf-8 -*-

"""
Unit test for the lsl.common.sdf module.
"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.common import sdf
from lsl.common.stations import lwa1


__revision__ = "$Rev$"
__version__  = "0.2"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(dataPath, 'tests', 'tbw-sdf.txt')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-sdf.txt')
drxFile = os.path.join(dataPath, 'tests', 'drx-sdf.txt')
solFile = os.path.join(dataPath, 'tests', 'sol-sdf.txt')
jovFile = os.path.join(dataPath, 'tests', 'jov-sdf.txt')
stpFile = os.path.join(dataPath, 'tests', 'stp-sdf.txt')


class sdf_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.sdf
	module."""
	
	### General ###
	def test_time(self):
		"""Test the sdf.parseTimeString() function."""
		
		# Different realizations of the same thing
		s1 = "EST 2011-01-01 12:13:14.567"
		s2 = "EST 2011 01 01 12:13:14.567"
		s3 = "EST 2011 Jan 01 12:13:14.567"
		
		self.assertEqual(sdf.parseTimeString(s1), sdf.parseTimeString(s2))
		self.assertEqual(sdf.parseTimeString(s1), sdf.parseTimeString(s3))
		self.assertEqual(sdf.parseTimeString(s2), sdf.parseTimeString(s3))
		
		# Month name and month number agreement
		for n,m in enumerate(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
			s1 = "UTC 2011-%s-14 12:13:14.000" % m
			s2 = "UTC 2011-%02i-14 12:13:14.000" % (n+1)
			self.assertEqual(sdf.parseTimeString(s1), sdf.parseTimeString(s2))
			
		# Time zone agreement - UTC
		s1 = "2011-01-01 12:13:14.567"
		s2 = "2011 01 01 12:13:14.567"
		self.assertEqual(sdf.parseTimeString(s1), sdf.parseTimeString(s2))
		
		# Time zone agreement - local
		for o,z in enumerate(['EST', 'CST', 'MST', 'PST']):
			h = 12
			o = -5 - o
			s1 = "%s 2011 01 01 %02i:13:14.567" % ('UTC', h)
			s2 = "%s 2011 01 01 %02i:13:14.567" % (z, h+o)
			self.assertEqual(sdf.parseTimeString(s1), sdf.parseTimeString(s2))
			
		# Details
		s1 = "2011 01 02 03:04:05.678"
		out = sdf.parseTimeString(s1)
		## Date
		self.assertEqual(out.year, 2011)
		self.assertEqual(out.month, 1)
		self.assertEqual(out.day, 2)
		## Time
		self.assertEqual(out.hour, 3)
		self.assertEqual(out.minute, 4)
		self.assertEqual(out.second, 5)
		self.assertEqual(out.microsecond, 678000)
		
		# LST at LWA1
		s1 = "LST 2013-01-08 19:42:00.000"
		s2 = "UTC 2013-01-08 19:38:26.639"
		self.assertEqual(sdf.parseTimeString(s1, site=lwa1), sdf.parseTimeString(s2))

	### TBW ###

	def test_tbw_parse(self):
		"""Test reading in a TBW SDF file."""
		
		project = sdf.parseSDF(tbwFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode,  'TBW')
		self.assertEqual(project.sessions[0].observations[0].mjd,   55616)
		self.assertEqual(project.sessions[0].observations[0].mpm,       0)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode,  'TBW')
		self.assertEqual(project.sessions[0].observations[1].mjd,   55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  700000)
	
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
		self.assertEqual(project.sessions[0].observations[0].filter,   7)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TBN')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
		self.assertEqual(project.sessions[0].observations[1].dur,  10000)
		self.assertEqual(project.sessions[0].observations[1].freq1, 832697741)
		self.assertEqual(project.sessions[0].observations[1].filter,   7)
	
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
		project.sessions[0].observations[0].frequency1 = 95.0e6
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		# Bad duration
		project.sessions[0].observations[0].frequency1 = 38.0e6
		project.sessions[0].observations[0].duration = '96:00:00.000'
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
	
	### DRX - TRK_RADEC ###
	
	def test_drx_parse(self):
		"""Test reading in a TRK_RADEC SDF file."""
		
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
		self.assertEqual(project.sessions[0].observations[0].filter,   7)
		self.assertAlmostEqual(project.sessions[0].observations[0].ra, 5.6, 6)
		self.assertAlmostEqual(project.sessions[0].observations[0].dec, 22.0, 6)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_RADEC')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
		self.assertEqual(project.sessions[0].observations[1].dur,  10000)
		self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
		self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
		self.assertEqual(project.sessions[0].observations[1].filter,   7)
		self.assertAlmostEqual(project.sessions[0].observations[1].ra, 5.6, 6)
		self.assertAlmostEqual(project.sessions[0].observations[1].dec, 22.0, 6)
		
	def test_drx_write(self):
		"""Test writing a TRK_RADEC SDF file."""
		
		project = sdf.parseSDF(drxFile)
		out = project.render()
		
	def test_drx_errors(self):
		"""Test various TRK_RADEC SDF errors."""
		
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
		
	### DRX - TRK_SOL ###
	
	def test_sol_parse(self):
		"""Test reading in a TRK_SOL SDF file."""
		
		project = sdf.parseSDF(solFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode, 'TRK_SOL')
		self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[0].mpm,      0)
		self.assertEqual(project.sessions[0].observations[0].dur,  10000)
		self.assertEqual(project.sessions[0].observations[0].freq1,  438261968)
		self.assertEqual(project.sessions[0].observations[0].freq2, 1928352663)
		self.assertEqual(project.sessions[0].observations[0].filter,   7)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_SOL')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
		self.assertEqual(project.sessions[0].observations[1].dur,  10000)
		self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
		self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
		self.assertEqual(project.sessions[0].observations[1].filter,   7)
		
	def test_sol_write(self):
		"""Test writing a TRK_SOL SDF file."""
		
		project = sdf.parseSDF(solFile)
		out = project.render()
		
	def test_sol_errors(self):
		"""Test various TRK_SOL SDF errors."""
		
		project = sdf.parseSDF(solFile)
		
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
		
	### DRX - TRK_JOV ###
	
	def test_jov_parse(self):
		"""Test reading in a TRK_JOV SDF file."""
		
		project = sdf.parseSDF(jovFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode, 'TRK_JOV')
		self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[0].mpm,      0)
		self.assertEqual(project.sessions[0].observations[0].dur,  10000)
		self.assertEqual(project.sessions[0].observations[0].freq1,  438261968)
		self.assertEqual(project.sessions[0].observations[0].freq2, 1928352663)
		self.assertEqual(project.sessions[0].observations[0].filter,   7)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_JOV')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
		self.assertEqual(project.sessions[0].observations[1].dur,  10000)
		self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
		self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
		self.assertEqual(project.sessions[0].observations[1].filter,   7)
		
	def test_jov_write(self):
		"""Test writing a TRK_JOV SDF file."""
		
		project = sdf.parseSDF(jovFile)
		out = project.render()
		
	def test_jov_errors(self):
		"""Test various TRK_JOV SDF errors."""
		
		project = sdf.parseSDF(jovFile)
		
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
		
	### DRX - STEPPED ###
	
	def test_stp_parse(self):
		"""Test reading in a STEPPED SDF file."""
		
		project = sdf.parseSDF(stpFile)
		
		# Basic file structure
		self.assertEqual(len(project.sessions), 1)
		self.assertEqual(len(project.sessions[0].observations), 2)
		
		# Observational setup - 1
		self.assertEqual(project.sessions[0].observations[0].mode, 'STEPPED')
		self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[0].mpm, 440000)
		self.assertEqual(project.sessions[0].observations[0].dur, 300000)
		self.assertEqual(project.sessions[0].observations[0].filter,   7)
		self.assertEqual(project.sessions[0].observations[0].obsFEE[0], [1,1])
		self.assertEqual(project.sessions[0].observations[0].aspFlt[0], 2)
		self.assertEqual(project.sessions[0].observations[0].aspAT1[0], 10)
		self.assertEqual(project.sessions[0].observations[0].aspAT2[0], 12)
		self.assertEqual(project.sessions[0].observations[0].aspATS[0], 14)
		
		# Steps - 1
		self.assertEqual(len(project.sessions[0].observations[0].steps), 4)
		for i in xrange(4):
			self.assertEqual(project.sessions[0].observations[0].steps[i].RADec, project.sessions[0].observations[0].RADec)
			self.assertEqual(project.sessions[0].observations[0].steps[i].freq1,  832697741)
			self.assertEqual(project.sessions[0].observations[0].steps[i].freq2, 1621569285)
		self.assertAlmostEqual(project.sessions[0].observations[0].steps[0].c1, 90.0, 6)
		self.assertAlmostEqual(project.sessions[0].observations[0].steps[0].c2, 45.0, 6)
		self.assertEqual(project.sessions[0].observations[0].steps[0].dur, 60000)
		self.assertAlmostEqual(project.sessions[0].observations[0].steps[-1].c1, 0.0, 6)
		self.assertAlmostEqual(project.sessions[0].observations[0].steps[-1].c2, 1.0, 6)
		self.assertEqual(project.sessions[0].observations[0].steps[-1].dur, 120000)
		
		# Observational setup - 2
		self.assertEqual(project.sessions[0].observations[1].mode, 'STEPPED')
		self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
		self.assertEqual(project.sessions[0].observations[1].mpm, 800000)
		self.assertEqual(project.sessions[0].observations[1].dur, 180000)
		self.assertEqual(project.sessions[0].observations[1].filter,   7)
		self.assertEqual(project.sessions[0].observations[1].obsFEE[0], [1,0])
		self.assertEqual(project.sessions[0].observations[1].aspFlt[0], 1)
		self.assertEqual(project.sessions[0].observations[1].aspAT1[0], 11)
		self.assertEqual(project.sessions[0].observations[1].aspAT2[0], 13)
		self.assertEqual(project.sessions[0].observations[1].aspATS[0], 15)
		
		# Steps - 2
		self.assertEqual(len(project.sessions[0].observations[1].steps), 2)
		for i in xrange(2):
			self.assertEqual(project.sessions[0].observations[1].steps[i].RADec, project.sessions[0].observations[1].RADec)
			self.assertEqual(project.sessions[0].observations[1].steps[i].freq1,  832697741)
			self.assertEqual(project.sessions[0].observations[1].steps[i].freq2, 1621569285)
		self.assertAlmostEqual(project.sessions[0].observations[1].steps[0].c1, 0.0, 6)
		self.assertAlmostEqual(project.sessions[0].observations[1].steps[0].c2, 90.0, 6)
		self.assertEqual(project.sessions[0].observations[1].steps[0].dur, 60000)
		self.assertAlmostEqual(project.sessions[0].observations[1].steps[-1].c1, 12.0, 6)
		self.assertAlmostEqual(project.sessions[0].observations[1].steps[-1].c2, 80.0, 6)
		self.assertEqual(project.sessions[0].observations[1].steps[-1].dur, 120000)
		
	def test_stp_write(self):
		"""Test writing a STEPPED SDF file."""
		
		project = sdf.parseSDF(stpFile)
		out = project.render()
		
	def test_stp_errors(self):
		"""Test various STEPPED SDF errors."""
		
		project = sdf.parseSDF(stpFile)
		
		# Bad filter
		project.sessions[0].observations[0].filter = 10
		self.assertFalse(project.validate())
		
		# Bad frequency
		project.sessions[0].observations[0].filter = 7
		project.sessions[0].observations[0].steps[0].frequency1 = 90.0e6
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		project.sessions[0].observations[0].steps[0].frequency1 = 38.0e6
		project.sessions[0].observations[0].steps[1].frequency2 = 90.0e6
		project.sessions[0].observations[0].update()
		self.assertFalse(project.validate())
		
		# Bad duration
		project.sessions[0].observations[0].steps[1].frequency2 = 38.0e6
		project.sessions[0].observations[0].steps[2].duration = '96:00:00.000'
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
