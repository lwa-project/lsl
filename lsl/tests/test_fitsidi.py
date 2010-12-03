# -*- coding: utf-8 -*-

"""Unit test for the lsl.writer.fitsidi modules."""

import os
import time
import unittest
import tempfile
import numpy
import pyfits

from lsl.common import stations as lwa_common
from lsl.correlator import uvUtils
from lsl.writer import fitsidi


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class fitsidi_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.writer.fitsidi
	module."""

	def __initData(self):
		"""Private function to generate a random set of data for writing a FITS
		IDI file.  The data is returned as a dictionary with keys:
		* freq - frequency array in Hz
		* site - lwa.common.stations object
		* stands - array of stand numbers
		* bl - list of baseline pairs in real stand numbers
		* vis - array of visibility data in baseline x freq format
		"""

		# Frequency range
		freq = numpy.arange(0,512)*20e6/512 + 40e6
		# Site and stands
		site = lwa_common.lwa1()
		stands = site.getStands('2010/11/25 00:00:00')
		# Set baselines and data
		blList = uvUtils.getBaselines(stands, IncludeAuto=True, Indicies=False)
		visData = numpy.random.rand(len(blList), len(freq))
		visData = visData.astype(numpy.complex64)

		return {'freq': freq, 'site': site, 'stands': stands, 'bl': blList, 'vis': visData}

	def test_write_tables(self):
		"""Test if the FITS IDI writer writes all of the tables."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		# Check that all of the extensions are there
		extNames = [hdu.name for hdu in hdulist]
		for ext in ['ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'BANDPASS', 'SOURCE', 'UV_DATA']:
			self.assertTrue(ext in extNames)
		hdulist.close()

	def test_array_geometry(self):
		"""Test the ARRAY_GEOMETRY table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		ag = hdulist['ARRAY_GEOMETRY'].data
		# Correct number of stands
		self.assertEqual(len(data['stands']), len(ag.field('NOSTA')))

		# Correct stand names
		names = ['LWA%03i' % stand for stand in data['stands']]
		for name, anname in zip(names, ag.field('ANNAME')):
			self.assertEqual(name, anname)

		hdulist.close()

	def test_frequency(self):
		"""Test the FREQUENCY table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		fq = hdulist['FREQUENCY'].data
		# Correct number of FREQIDs
		self.assertEqual(len(fq.field('FREQID')), 1)

		# Correct channel width
		self.assertAlmostEqual(fq.field('CH_WIDTH')[0], numpy.abs(data['freq'][1]-data['freq'][0]), 4)

		# Correct bandwidth
		self.assertAlmostEqual(fq.field('TOTAL_BANDWIDTH')[0], numpy.abs(data['freq'][-1]-data['freq'][0]).astype(numpy.float32), 4)

		# Correct sideband
		self.assertEqual(fq.field('SIDEBAND')[0], 1)

		hdulist.close()

	def test_antenna(self):
		"""Test the ANTENNA table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		an = hdulist['ANTENNA'].data
		# Correct number of stands
		self.assertEqual(len(data['stands']), len(an.field('ANTENNA_NO')))

		# Correct FREQIDs
		for freqid in an.field('FREQID'):
			self.assertEqual(freqid, 1)

		hdulist.close()

	def test_bandpass(self):
		"""Test the BANDPASS table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		bp = hdulist['BANDPASS'].data
		# Correct number of entries
		self.assertEqual(len(data['stands']), len(bp.field('ANTENNA_NO')))

		# Correct Source ID number
		for src in bp.field('SOURCE_ID'):
			self.assertEqual(src, 0)

		# Correct FREQIDs
		for freqid in bp.field('FREQID'):
			self.assertEqual(freqid, 1)

		hdulist.close()

	def test_source(self):
		"""Test the SOURCE table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		so = hdulist['SOURCE'].data
		# Correct number of entries
		self.assertEqual(len(so.field('SOURCE_ID')), 1)

		# Correct Source ID number
		self.assertEqual(so.field('SOURCE_ID'), 1)

		hdulist.close()

	def test_uvdata(self):
		"""Test the UV_DATA table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		uv = hdulist['UV_DATA'].data

		# Load the mapper
		try:
			mp = hdulist['NOSTA_MAPPER'].data
			nosta = mp.field('NOSTA')
			noact = mp.field('NOACT')
		except KeyError:
			ag = hdulist['ARRAY_GEOMETRY'].data
			nosta = ag.field('NOSTA')
			noact = ag.field('NOACT')
		mapper = {}
		for s,a in zip(nosta, noact):
			mapper[s] = a

		# Correct number of visibilities
		self.assertEqual(len(uv.field('FLUX')), data['vis'].shape[0])
		
		# Correct number of frequencies
		for vis in uv.field('FLUX'):
			self.assertEqual(len(vis), 2*len(data['freq']))

		# Correct values
		for bl, vis in zip(uv.field('BASELINE'), uv.field('FLUX')):
			# Convert mapped stands to real stands
			stand1 = mapper[(bl >> 8) & 255]
			stand2 = mapper[bl & 255]

			# Find out which visibility set in the random data corresponds to the 
			# current visibility
			i = 0
			for s1,s2 in data['bl']:
				if s1 == stand1 and s2 == stand2:
					break
				else:
					i = i + 1
			
			# Extract the data and run the comparison
			visData = numpy.zeros(len(data['freq']), dtype=numpy.complex64)
			visData.real = vis[0::2]
			visData.imag = vis[1::2]
			for vd, sd in zip(visData, data['vis'][i,:]):
				self.assertAlmostEqual(vd, sd, 8)
			i = i + 1
		
		hdulist.close()

	def test_mapper(self):
		"""Test the NOSTA_MAPPER table."""

		testTime = time.time()
		testPath = tempfile.mkdtemp(prefix='fitsidi-', suffix='.tmp')
		testFile = os.path.join(testPath, 'idi-test.fits')
		
		# Get some data
		data = self.__initData()
		
		# Start the file
		fits = fitsidi.IDI(testFile, refTime=testTime)
		fits.setStokes(['xx'])
		fits.setFrequency(data['freq'])
		fits.setGeometry(data['site'], data['stands'])
		fits.addDataSet(testTime, 6.0, data['bl'], data['vis'])
		fits.write()

		# Open the file and examine
		hdulist = pyfits.open(testFile)
		extNames = [hdu.name for hdu in hdulist]
		if data['stands'].max() > 255:
			self.assertTrue('NOSTA_MAPPER' in extNames)

			# Make sure the mapper makes sense
			mp = hdulist['NOSTA_MAPPER'].data
			ag = hdulist['ARRAY_GEOMETRY'].data
			mNoSta = mp.field('NOSTA')
			aNoSta = ag.field('NOSTA')
			mNoAct = mp.field('NOACT')
			aAnNam = ag.field('ANNAME')
			for msta, mact, asta, anam in zip(mNoSta, mNoAct, aNoSta, aAnNam):
				self.assertEqual(msta, asta)
				self.assertEqual(mact, int(anam[3:]))


class fitsidi_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(fitsidi_tests))


if __name__ == '__main__':
	unittest.main()