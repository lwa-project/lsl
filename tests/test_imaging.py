# -*- coding: utf-8 -*-

"""Unit test for lsl.imaging modules"""

import os
import copy
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.imaging import utils, selfCal


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


idiFile = os.path.join(dataPath, 'tests', 'idi-test.fits')


class imaging_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.imaging
	modules."""
	
	def test_CorrelatedData(self):
		"""Test the utils.CorrelatedData class."""
		
		# Open the FITS IDI file
		idi = utils.CorrelatedData(idiFile)
		
		# Dates
		self.assertEqual(idi.dateObs.strftime("%Y-%m-%dT%H:%M:%S"), "2013-01-14T23:30:37")

		# Stand and baseline counts
		self.assertEqual(len(idi.stands), 5)
		self.assertEqual(idi.totalBaselineCount, 5*(5+1)/2)

		# Basic functions (just to see that they run)
		junk = idi.getAntennaArray()
		junk = idi.getObserver()
		junk = idi.getDataSet(1)

		# Error checking
		self.assertRaises(RuntimeError, idi.getDataSet, 2)

	def test_sort(self):
		"""Test the utils.sortDataDict function."""
		
		# Open the FITS IDI file
		idi = utils.CorrelatedData(idiFile)
		
		# Get some data to sort
		ds = idi.getDataSet(1, sort=False)
		
		# Sort
		dss = copy.deepcopy(ds)
		utils.sortDataDict(dss)
		for prop in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			for pol in ds['bls'].keys():
				self.assertEqual(len(dss[prop][pol]), len(ds[prop][pol]))
		
	def test_prune(self):
		"""Test the utils.pruneBaselineRange function."""
		
		# Open the FITS IDI file
		idi = utils.CorrelatedData(idiFile)
		
		# Get some data to sort
		ds = idi.getDataSet(1)
		
		# Prune
		dsp1 = utils.pruneBaselineRange(ds, uvMin=10)
		for prop in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			for pol in ds['bls'].keys():
				self.assertTrue(len(dsp1[prop][pol]) < len(ds[prop][pol]))
				
		# Auto-prune
		dsp2 = idi.getDataSet(1, uvMin=10)
		for prop in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			for pol in ds['bls'].keys():
				self.assertEqual(len(dsp1[prop][pol]), len(dsp2[prop][pol]))

		# Auto-prune that should result in no baselines
		dsp3 = idi.getDataSet(1, uvMin=100)
		for prop in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			for pol in ds['bls'].keys():
				self.assertEqual(len(dsp3[prop][pol]), 0)
				
	def test_gridding(self):
		"""Test building a image from a visibility data set."""
		
		# Open the FITS IDI file
		idi = utils.CorrelatedData(idiFile)
		
		# Build the image
		ds = idi.getDataSet(1)
		junk = utils.buildGriddedImage(ds)

		# Error checking
		self.assertRaises(RuntimeError, utils.buildGriddedImage, ds, pol='xy')
		
	def test_selfcal(self):
		"""Test running a simple self calibration."""
		
		# Open the FITS IDI file
		idi = utils.CorrelatedData(idiFile)
		
		# Go for it!
		aa = idi.getAntennaArray()
		ds = idi.getDataSet(1)
		junk = selfCal.selfCal(aa, ds, ds, 173, 'xx')
		
		# Error checking
		self.assertRaises(RuntimeError, selfCal.selfCal, aa, ds, ds, 173, 'yx', refAnt=0  )
		self.assertRaises(RuntimeError, selfCal.selfCal, aa, ds, ds, 173, 'yx', refAnt=564)


class imaging_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.imaging units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(imaging_tests)) 


if __name__ == '__main__':
	unittest.main()
