# -*- coding: utf-8 -*-

"""Unit test for lsl.reader.ldp module"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import ldp
from lsl.reader import errors


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')
drxFile = os.path.join(dataPath, 'tests', 'drx-test.dat')
drspecFile = os.path.join(dataPath, 'tests', 'drspec-test.dat')

tbfFile = os.path.join(dataPath, 'tests', 'tbf-test.dat')


class ldp_adp_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.reader
	modules."""
	
	### TBF ###
	
	def test_ldp_tbf(self):
		"""Test the LDP interface for a TBF file."""
		
		f = ldp.TBFFile(tbfFile)
		
		# File info
		self.assertEqual(f.getInfo("sampleRate"), 25e3)
		self.assertEqual(f.getInfo("dataBits"), 4)
		self.assertEqual(f.getInfo("nFrames"), 5)
		
		# Read a frame
		frame = f.readFrame()
		
		# Get the remaining frame count
		self.assertEqual(f.getRemainingFrameCount(), f.getInfo("nFrames")-1)
		
		# Reset
		f.reset()
		
		# Close it out
		f.close()
		
	def test_ldp_tbf_nocheck(self):
		"""Test the LDP interface for a TBF file."""
		
		f = ldp.TBFFile(tbfFile, ignoreTimeTagErrors=True)
		
		# File info
		self.assertEqual(f.getInfo("sampleRate"), 25e3)
		self.assertEqual(f.getInfo("dataBits"), 4)
		self.assertEqual(f.getInfo("nFrames"), 5)
		
		# Read a frame
		frame = f.readFrame()
		
		# Get the remaining frame count
		self.assertEqual(f.getRemainingFrameCount(), f.getInfo("nFrames")-1)
		
		# Reset
		f.reset()
		
		# Close it out
		f.close()
		
	### File Type Discovery ###
	
	def test_ldp_discover_tbw(self):
		"""Test the LDP LWA1DataFile function of TBW."""
		# TBW
		self.assertRaises(RuntimeError, ldp.LWASVDataFile, tbwFile)
		
	def test_ldp_discover_tbn(self):
		"""Test the LDP LWASVDataFile function of TBN."""
		# TBN
		f = ldp.LWASVDataFile(tbnFile)
		self.assertEqual(type(f), ldp.TBNFile)
		
	def test_ldp_discover_drx(self):
		"""Test the LDP LWASVDataFile function of DRX."""
		# DRX
		f = ldp.LWASVDataFile(drxFile)
		self.assertEqual(type(f), ldp.DRXFile)
		
	def test_ldp_discover_drspec(self):
		"""Test the LDP LWASVDataFile function of DR Spectrometer."""
		# DR Spectrometer
		f = ldp.LWASVDataFile(drspecFile)
		self.assertEqual(type(f), ldp.DRSpecFile)
		
	def test_ldp_discover_tbf(self):
		"""Test the LDP LWASVDataFile function of TBF."""
		# TBF
		f = ldp.LWASVDataFile(tbfFile)
		self.assertEqual(type(f), ldp.TBFFile)


class ldp_adp_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader.ldp 
	unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(ldp_adp_tests)) 


if __name__ == '__main__':
	unittest.main()
