# -*- coding: utf-8 -*-

"""Unit test for the lsl.sim.s60/tbw/tbn/drx modules."""

import os
import unittest
import tempfile
import numpy

from lsl.sim import errors


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class fake_S60_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.s60
	module."""

	from lsl.reader import s60 as reader
	from lsl.sim import s60 as writer

	testPath = tempfile.mkdtemp(prefix='fakeS60-', suffix='.tmp')
	testFile = os.path.join(testPath, 's60-test.dat')
	dataFile = os.path.join(dataPath, 'tests', 's60-test.dat')

	def test_read_write_frame(self):
		"""Test that the S60 data writer works for frames."""

		# Read in a S60 frame from the test file
		fh = open(dataFile, 'rb')
		origData = reader.readFrame(fh)
		fh.close()
		# Write the data to a S60 test frame
		fh = open(testFile, 'wb')
		rawFrame = writer.frame2frame(origData)
		rawFrame.tofile(fh)
		fh.close()
		# Read in the 
		fh = open(testFile, 'rb')
		fakeData = reader.readFrame(fh)
		fh.close()

		for i in range(25):
			self.assertAlmostEqual(fakeData[i].real, origData[i].real, 4)
			self.assertAlmostEqual(fakeData[i].imag, origData[i].imag, 4)

	def test_read_write_chunk(self):
		"""Test that the S60 data writer works for chunks."""

		# Read in a S60 frame from the test file
		fh = open(dataFile, 'rb')
		origData = reader.readChunk(fh, Chunk=2202)
		fh.close()
		# Write the data to a S60 test frame
		fh = open(testFile, 'wb')
		rawFrame = writer.chunk2frame(origData)
		rawFrame.tofile(fh)
		fh.close()
		# Read in the 
		fh = open(testFile, 'rb')
		fakeData = reader.readChunk(fh, Chunk=2202)
		fh.close()

		self.assertEqual(len(fakeData), len(origData))
		for i in range(25):
			self.assertAlmostEqual(fakeData[i].real, origData[i].real, 4)
			self.assertAlmostEqual(fakeData[i].imag, origData[i].imag, 4)

	def test_frame_errors(self):
		


class fake_TBW_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.tbw
	module."""

	from lsl.reader import tbw as reader
	from lsl.sim import tbw as writer

	testPath = tempfile.mkdtemp(prefix='fakeTBW-', suffix='.tmp')
	testFile = os.path.join(testPath, 'tbw-test.dat')
	dataFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')

	def test_read_write_frame(self):
		"""Test that the TBW data writer works."""

		# Read in a TBW frame from the test file
		fh = open(dataFile, 'rb')
		origFrame = reader.readFrame(fh)
		fh.close()
		# Write the data to a TBW test frame
		fh = open(testFile, 'wb')
		rawFrame = writer.frame2frame(origData)
		rawFrame.tofile(fh)
		fh.close()
		# Read in the 
		fh = open(testFile, 'rb')
		fakeFrame = reader.readFrame(fh)
		fh.close()

		# Test values returned by info functions
		self.assertEqual(fakeFrame.parseID(), origFrame.parseID())
		self.assertEqual(fakeFrame.getDataBits(), origFrame.getDataBits())
		# Test raw header values
		self.assertTrue(fakeFrame.header.isTBW())
		self.assertEqual(fakeFrame.header.frameCount, origFrame.header.frameCount)
		self.assertEqual(fakeFrame.header.secondsCount, origFrame.header.secondsCount)
		# Test raw data values
		self.assertEqual(fakeFrame.data.timeTag, origFrame.data.timeTag)
		for i in range(25):
			self.assertAlmostEqual(fakeFrame.data.xy[0,i], origFrame.data.xy[0,i], 4)
			self.assertAlmostEqual(fakeFrame.data.xy[1,i], origFrame.data.xy[1,i], 4)


class fake_TBN_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.tbn
	module."""

	from lsl.reader import tbn as reader
	from lsl.sim import tbn as writer

	testPath = tempfile.mkdtemp(prefix='fakeTBN-', suffix='.tmp')
	testFile = os.path.join(testPath, 'tbn-test.dat')
	dataFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')

	def test_read_write(self):
		"""Test that the TBN data writer works."""

		# Read in a TBN frame from the test file
		fh = open(dataFile, 'rb')
		origData = reader.readFrame(fh)
		fh.close()
		# Write the data to a TBN test frame
		fh = open(testFile, 'wb')
		rawFrame = writer.frame2frame(origData)
		rawFrame.tofile(fh)
		fh.close()
		# Read in the 
		fh = open(testFile, 'rb')
		fakeData = reader.readFrame(fh)
		fh.close()

		# Test values returned by info functions
		self.assertEqual(fakeFrame.parseID()[0], origFrame.parseID()[0])
		self.assertEqual(fakeFrame.parseID()[1], origFrame.parseID()[1])
		# Test raw header values
		self.assertTrue(fakeFrame.header.isTBN())
		self.assertEqual(fakeFrame.header.frameCount, origFrame.header.frameCount)
		self.assertEqual(fakeFrame.header.secondsCount, origFrame.header.secondsCount)
		# Test raw data values
		self.assertEqual(fakeFrame.data.timeTag, origFrame.data.timeTag)
		for i in range(25):
			self.assertAlmostEqual(fakeData[i].real, origData[0].real, 4)
			self.assertAlmostEqual(fakeData[i].imag, origData[i].imag, 4)


class fakedata_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(fake_S60_tests))
		self.addTests(loader.loadTestsFromTestCase(fake_TBW_tests))
		self.addTests(loader.loadTestsFromTestCase(fake_TBN_tests))


if __name__ == '__main__':
	unittest.main()
