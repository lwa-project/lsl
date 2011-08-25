# -*- coding: utf-8 -*-

"""Unit test for the lsl.writer.sdfits module."""

import os
import unittest
import numpy
import tempfile
import pyfits

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import errors
from lsl.writer import sdfits


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"

tbwFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')


class sdfits_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.writer.sdfits
	module."""

	testPath = None

	def setUp(self):
		"""Turn off all numpy warnings and create the temporary file directory."""

		numpy.seterr(all='ignore')
		self.testPath = tempfile.mkdtemp(prefix='test-sdfits-', suffix='.tmp')

	def __getTBW(self):
		"""Private function to load in the test TBW data and get the frames."""

		fh = open(tbwFile, 'rb')

		# Frames 1 through 8
		frames = []
		for i in range(1,9):
			frames.append(tbw.readFrame(fh))

		return frames

	def __getTBN(self, vanilla=False):
		"""Private function to load in the test TBN data and get the frames.  If 
		the keyword 'vanilla' is set to True, gain, sample rate, and frequency meta-
		data are not added to the frames."""

		fh = open(tbnFile, 'rb')

		# Frames 1 through 8
		frames = []
		for i in range(1,9):
			frames.append(tbn.readFrame(fh))

		if not vanilla:
			# Set some values for the other meta-data
			for frame in frames:
				frame.setSampleRate(100000)

		return frames

	#def test_write_tbw_single(self):
		#"""Test that TBW data can be written to a SD FITS file one at a time."""

		## Setup the file names
		#testFile = os.path.join(self.testPath, 'tbw-test-WS.fits')

		## Get some data and organize it into a nice numpy array
		#frames = self.__getTBW()
		#junk = [frame.parseID() for frame in frames]
		#stands = []
		#for j in junk:
			#if j not in stands:
				#stands.append(j)
		#count = {1: 0, 2: 0}
		#frameData = numpy.zeros((len(stands), 4, 2, 400), dtype=numpy.int16)
		#for frame in frames:
			#frameData[frame.parseID()-1, count[frame.parseID()], :, :] = frame.data.xy
			#count[frame.parseID()] = count[frame.parseID()] + 1

		## Load it in
		#fits = sdfits.TBW(testFile, LFFT=128, UseQueue=False)
		#for frame in frames:
			#fits.addStandData(frame)
		#fits.close()
	
		## Read it back in
		#hdulist = pyfits.open(testFile)
		## Check if we have the correct number of extensions
		#self.assertEqual(len(stands)+1, len(hdulist))
		
		#for hdu in hdulist[1:]:
			#stand = hdu.data
			#sNumb = hdu.header['STAND']

			## Check the pol. values
			#for pol in stand.field('pol'):
				#self.assertTrue(pol in [0, 1])

			## Check the data length
			#x, y = stand.field('data').shape
			#self.assertEqual(x, 8, "Shape mis-match on dimension 1: %i != 8" % x)
			#self.assertEqual(y, 127, "Shape mis-match on dimension 2: %i != 127" % y)

			#for data in stand.field('data'):
				#self.assertEqual(data.dtype.kind, 'f')
				#self.assertEqual(len(data), 127)
			
			### Check the data itself
			##count = {0: 0, 1: 0}
			##for pol, data in zip(stand.field('pol'), stand.field('data')):
				##for ts, fs in zip(data, frameData[sNumb-1,count[pol],pol,:]):
					##self.assertAlmostEqual(ts, fs, 6)
				##count[pol] = count[pol] + 1
		
		#hdulist.close

	def test_write_tbw_queue(self):
		"""Test that TBW data can be written to a SD FITS file with a queue."""

		# Setup the file names
		testFile = os.path.join(self.testPath, 'tbw-test-WQ.fits')

		# Get some data and organize it into a nice numpy array
		frames = self.__getTBW()
		junk = [frame.parseID() for frame in frames]
		stands = []
		for j in junk:
			if j not in stands:
				stands.append(j)
		count = {1: 0, 2: 0}
		frameData = numpy.zeros((len(stands), 4, 2, 400), dtype=numpy.int16)
		for frame in frames:
			frameData[frame.parseID()-1, count[frame.parseID()], :, :] = frame.data.xy
			count[frame.parseID()] = count[frame.parseID()] + 1

		# Load it in
		fits = sdfits.TBW(testFile, LFFT=128, UseQueue=True)
		for frame in frames:
			fits.addStandData(frame)
		fits.close()

		# Read it back in
		hdulist = pyfits.open(testFile)
		# Check if we have the correct number of extensions
		self.assertEqual(len(stands)+1, len(hdulist))
		
		for hdu in hdulist[1:]:
			stand = hdu.data
			sNumb = hdu.header['STAND']

			# Check the pol. values
			for pol in stand.field('pol'):
				self.assertTrue(pol in [0, 1])

			# Check the data length
			x, y = stand.field('data').shape
			self.assertEqual(x, 8, "Shape mis-match on dimension 1: %i != 8" % x)
			self.assertEqual(y, 127, "Shape mis-match on dimension 2: %i != 127" % y)

			for data in stand.field('data'):
				self.assertEqual(data.dtype.kind, 'f')
				self.assertEqual(len(data), 127)
			
			## Check the data itself
			#count = {0: 0, 1: 0}
			#for pol, data in zip(stand.field('pol'), stand.field('data')):
				#for ts, fs in zip(data, frameData[sNumb-1,count[pol],pol,:]):
					#self.assertAlmostEqual(ts, fs, 6)
				#count[pol] = count[pol] + 1
		
		hdulist.close

	def test_write_tbn_single(self):
		"""Test that TBN data can be written to a SD FITS file one at a time."""

		# Setup the file names
		testFile = os.path.join(self.testPath, 'tbn-test-WS.fits')

		# Get some data and organize it into a nice numpy array
		frames = self.__getTBN()
		junk = [frame.parseID()[0] for frame in frames]
		stands = []
		for j in junk:
			if j not in stands:
				stands.append(j)
		count = {1: 0, 2: 0, 3: 0, 4: 0}
		frameData = numpy.zeros((len(stands), 4, 2, 512), dtype=numpy.singlecomplex)
		for frame in frames:
			stand, pol = frame.parseID()
			frameData[stand-1, count[stand], pol, :] = frame.data.iq
			if pol == 1:
				count[stand] = count[stand] + 1

		# Load it in
		fits = sdfits.TBN(testFile, LFFT=256, UseQueue=False)
		fits.setSampleRate(100000)
		for frame in frames:
			fits.addStandData(frame)
		fits.close()

		# Read it back in
		hdulist = pyfits.open(testFile)
		# Check if we have the correct number of extensions
		self.assertEqual(len(stands)+1, len(hdulist))
		
		for hdu in hdulist[1:]:
			stand = hdu.data
			sNumb = hdu.header['STAND']

			# Check the pol. values
			for pol in stand.field('pol'):
				self.assertTrue(pol in [0, 1])

			# Check the data length
			x, y = stand.field('data').shape
			self.assertEqual(x, 2, "Shape mis-match on dimension 1: %i != 2" % x)
			self.assertEqual(y, 255, "Shape mis-match on dimension 2: %i != 255" % y)

			for data in stand.field('data'):
				self.assertEqual(data.dtype.kind, 'f')
				self.assertEqual(len(data), 255)
			
			## Check the data itself
			#count = {0: 0, 1: 0}
			#for pol, data in zip(stand.field('pol'), stand.field('data')):
				#for ts, fs in zip(data, frameData[sNumb-1,count[pol],pol,:]):
					#self.assertAlmostEqual(ts, fs, 6)
				#count[pol] = count[pol] + 1
		
		hdulist.close

	def test_write_tbn_queue(self):
		"""Test that TBN data can be written to a SD FITS file with a queue"""
		
		# Setup the file names
		testFile = os.path.join(self.testPath, 'tbn-test-WQ.fits')

		# Get some data
		frames = self.__getTBN()

		# Load it in
		junk = [frame.parseID()[0] for frame in frames]
		stands = []
		for j in junk:
			if j not in stands:
				stands.append(j)
		count = {1: 0, 2: 0, 3: 0, 4: 0}
		frameData = numpy.zeros((len(stands), 4, 2, 512), dtype=numpy.singlecomplex)
		for frame in frames:
			stand, pol = frame.parseID()
			frameData[stand-1, count[stand], pol, :] = frame.data.iq
			if pol == 1:
				count[stand] = count[stand] + 1

		# Load it in
		fits = sdfits.TBN(testFile, LFFT=256, UseQueue=False)
		fits.setSampleRate(100000)
		for frame in frames:
			fits.addStandData(frame)
		fits.close()

		# Read it back in
		hdulist = pyfits.open(testFile)
		# Check if we have the correct number of extensions
		self.assertEqual(len(stands)+1, len(hdulist))
		
		for hdu in hdulist[1:]:
			stand = hdu.data
			sNumb = hdu.header['STAND']

			# Check the pol. values
			for pol in stand.field('pol'):
				self.assertTrue(pol in [0, 1])

			# Check the data length
			x, y = stand.field('data').shape
			self.assertEqual(x, 2, "Shape mis-match on dimension 1: %i != 2" % x)
			self.assertEqual(y, 255, "Shape mis-match on dimension 2: %i != 255" % y)

			for data in stand.field('data'):
				self.assertEqual(data.dtype.kind, 'f')
				self.assertEqual(len(data), 255)
			
			## Check the data itself
			#count = {0: 0, 1: 0}
			#for pol, data in zip(stand.field('pol'), stand.field('data')):
				#for ts, fs in zip(data, frameData[sNumb-1,count[pol],pol,:]):
					#self.assertAlmostEqual(ts, fs, 6)
				#count[pol] = count[pol] + 1
		
		hdulist.close

	def tearDown(self):
		"""Remove the test path directory and its contents"""

		tempFiles = os.listdir(self.testPath)
		for tempFile in tempFiles:
			os.unlink(os.path.join(self.testPath, tempFile))
		os.rmdir(self.testPath)
		self.testPath = None


class  sdfits_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.sim.vis units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)

		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(sdfits_tests)) 


if __name__ == '__main__':
	unittest.main()
