# -*- coding: utf-8 -*-

"""Unit test for the lsl.sim.s60/tbw/tbn/drx modules."""

import os
import copy
import unittest
import tempfile
import numpy

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import s60 as s60Reader
from lsl.sim import s60 as s60Writer
from lsl.reader import tbw as tbwReader
from lsl.sim import tbw as tbwWriter
from lsl.reader import tbn as tbnReader
from lsl.sim import tbn as tbnWriter
from lsl.reader import drx as drxReader
from lsl.sim import drx as drxWriter
from lsl.sim import errors


__revision__ = "$ Revision: 4 $"
__version__  = "0.2"
__author__    = "Jayce Dowell"


s60File = os.path.join(dataPath, 'tests', 's60-test.dat')
tbwFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')
drxFile = os.path.join(dataPath, 'tests', 'drx-test.dat')


class fake_S60_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.s60
	module."""

	testPath = None

	def setUp(self):
		"""Turn off all numpy warnings and create the temporary file directory."""

		numpy.seterr(all='ignore')
		self.testPath = tempfile.mkdtemp(prefix='test-fakedata-', suffix='.tmp')

	def test_write_frame(self):
		"""Test that the S60 data writer works for frames."""

		testFile = os.path.join(self.testPath, 's60-test-W.dat')

		# Read in a S60 frame from the test file
		fh = open(s60File, 'rb')
		origData = s60Reader.readFrame(fh)
		fh.close()
		# Write the data to a S60 test frame
		fh = open(testFile, 'wb')
		rawFrame = s60Writer.frame2frame(origData)
		rawFrame.tofile(fh)
		fh.close()
		# Read in the 
		fh = open(testFile, 'rb')
		fakeData = s60Reader.readFrame(fh)
		fh.close()

		# Test the array size and values
		self.assertEqual(len(fakeData), s60Reader.FrameSize/2)
		for i in range(s60Reader.FrameSize/2):
			self.assertAlmostEqual(fakeData[i].real, origData[i].real, 4)
			self.assertAlmostEqual(fakeData[i].imag, origData[i].imag, 4)

	def test_write_chunk(self):
		"""Test that the S60 data writer works for chunks."""

		testFile = os.path.join(self.testPath, 's60-test-C.dat')

		# Read in a S60 frame from the test file
		fh = open(s60File, 'rb')
		origData = s60Reader.readChunk(fh, Chunk=2202)
		fh.close()
		
		# Write the data to a S60 test frame
		fh = open(testFile, 'wb')
		rawFrame = s60Writer.chunk2frame(origData)
		rawFrame.tofile(fh)
		fh.close()
		
		# Read in the 
		fh = open(testFile, 'rb')
		fakeData = s60Reader.readChunk(fh, Chunk=2202)
		fh.close()

		# Test the array size and values
		self.assertEqual(len(fakeData), len(origData))
		for i in range(2202):
			self.assertAlmostEqual(fakeData[i].real, origData[i].real, 4)
			self.assertAlmostEqual(fakeData[i].imag, origData[i].imag, 4)

	def test_frame_data_errors(self):
		"""Test the different error scenarios in s60.frame2frame."""
		
		# Read in a S60 frame from the test file
		fh = open(s60File, 'rb')
		origData = s60Reader.readChunk(fh, Chunk=2202)
		fh.close()
		
		# Try to write a chunk using frame2frame
		self.assertRaises(errors.invalidDataSize, s60Writer.frame2frame, origData)
		
		# Try to write using only real data with frame2frame
		origData = origData.real
		self.assertRaises(errors.invalidDataType, s60Writer.frame2frame, origData)
	
	def test_chunk_data_errors(self):
		"""Test the different error scenarios in s60.chunk2frame."""

		# Read in a S60 frame from the test file
		fh = open(s60File, 'rb')
		origData = s60Reader.readChunk(fh, Chunk=2202)
		fh.close()
		
		# Try to write using only real data with chunk2frame
		origData = origData.real
		self.assertRaises(errors.invalidDataType, s60Writer.chunk2frame, origData)

	def tearDown(self):
		"""Remove the test path directory and its contents"""

		tempFiles = os.listdir(self.testPath)
		for tempFile in tempFiles:
			os.unlink(os.path.join(self.testPath, tempFile))
		os.rmdir(self.testPath)
		self.testPath = None


class fake_TBW_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.tbw
	module."""

	testPath = None

	def setUp(self):
		"""Turn off all numpy warnings and create the temporary file directory."""

		numpy.seterr(all='ignore')
		self.testPath = tempfile.mkdtemp(prefix='test-fakedata-', suffix='.tmp')

	def test_sim_frame(self):
		"""Test the tbw.SimFrame class."""
		
		# Read in a TBW frame from the test file
		fh = open(tbwFile, 'rb')
		origFrame = tbwReader.readFrame(fh)
		fh.close()
		
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(origFrame)
		# Test the validity of the SimFrame
		self.assertTrue(fakeFrame.isValid())

	def test_write_frame(self):
		"""Test that the TBW data writer works."""

		testFile = os.path.join(self.testPath, 'tbw-test-W.dat')

		# Read in a TBW frame from the test file
		fh = open(tbwFile, 'rb')
		origFrame = tbwReader.readFrame(fh)
		fh.close()
		
		# Write the data to a TBW test frame
		fh = open(testFile, 'wb')
		rawFrame = tbwWriter.frame2frame(origFrame)
		rawFrame.tofile(fh)
		fh.close()
		
		# Read in the 
		fh = open(testFile, 'rb')
		fakeFrame = tbwReader.readFrame(fh)
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
		for i in range(400):
			self.assertEqual(fakeFrame.data.xy[0,i], origFrame.data.xy[0,i])
			self.assertEqual(fakeFrame.data.xy[1,i], origFrame.data.xy[1,i])

	def test_frame_data_errors(self):
		"""Test the data error scenarios when validating a TBW SimFrame ."""
		
		# Read in a TBW frame from the test file
		fh = open(tbwFile, 'rb')
		origFrame = tbwReader.readFrame(fh)
		fh.close()
		
		# Try to validate frame with the wrong data type
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.xy = fakeFrame.data.xy.astype(numpy.complex64)
		self.assertRaises(errors.invalidDataType, fakeFrame.isValid, raiseErrors=True)
		
		# Try to validate frame with the wrong data size
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.xy = None
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.xy = fakeFrame.data.xy[0,:]
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.xy = fakeFrame.data.xy[:,0:50]
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.dataBits = 4
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		
	def test_frame_header_errors(self):
		"""Test the header error scenarios when validating a TBW SimFrame."""
		
		# Read in a TBW frame from the test file
		fh = open(tbwFile, 'rb')
		origFrame = tbwReader.readFrame(fh)
		fh.close()
		
		# Try to validate frame with the wrong stand number
		fakeFrame = tbwWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.stand = 300
		self.assertRaises(errors.invalidStand, fakeFrame.isValid, raiseErrors=True)

	def tearDown(self):
		"""Remove the test path directory and its contents"""

		tempFiles = os.listdir(self.testPath)
		for tempFile in tempFiles:
			os.unlink(os.path.join(self.testPath, tempFile))
		os.rmdir(self.testPath)
		self.testPath = None


class fake_TBN_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.tbn
	module."""

	testPath = None

	def setUp(self):
		"""Turn off all numpy warnings and create the temporary file directory."""

		numpy.seterr(all='ignore')
		self.testPath = tempfile.mkdtemp(prefix='test-fakedata-', suffix='.tmp')

	def test_sim_frame(self):
		"""Test the tbn.SimFrame class."""
		
		# Read in a TBN frame from the test file
		fh = open(tbnFile, 'rb')
		origFrame = tbnReader.readFrame(fh)
		fh.close()
		
		fakeFrame = tbnWriter.SimFrame()
		fakeFrame.loadFrame(origFrame)
		# Test the validity of the SimFrame
		self.assertTrue(fakeFrame.isValid())

	def test_write_frame(self):
		"""Test that the TBN data writer works."""
		
		testFile = os.path.join(self.testPath, 'tbn-test-W.dat')

		# Read in a TBN frame from the test file
		fh = open(tbnFile, 'rb')
		origFrame = tbnReader.readFrame(fh)
		fh.close()
		
		# Write the data to a TBN test frame
		fh = open(testFile, 'wb')
		rawFrame = tbnWriter.frame2frame(origFrame)
		rawFrame.tofile(fh)
		fh.close()
		
		# Read in the 
		fh = open(testFile, 'rb')
		fakeFrame = tbnReader.readFrame(fh)
		fh.close()

		# Test values returned by info functions
		self.assertEqual(fakeFrame.parseID()[0], origFrame.parseID()[0])
		self.assertEqual(fakeFrame.parseID()[1], origFrame.parseID()[1])
		
		# Test raw header values
		self.assertTrue(fakeFrame.header.isTBN())
		self.assertEqual(fakeFrame.header.frameCount, origFrame.header.frameCount)
		self.assertEqual(fakeFrame.header.tuningWord, origFrame.header.tuningWord)
		
		# Test raw data values
		self.assertEqual(fakeFrame.data.timeTag, origFrame.data.timeTag)
		for i in range(512):
			self.assertEqual(fakeFrame.data.iq[i].real, origFrame.data.iq[i].real)
			self.assertEqual(fakeFrame.data.iq[i].imag, origFrame.data.iq[i].imag)

	def test_frame_data_errors(self):
		"""Test the data error scenarios when validating a TBN SimFrame."""
		
		# Read in a TBN frame from the test file
		fh = open(tbnFile, 'rb')
		origFrame = tbnReader.readFrame(fh)
		fh.close()
		
		# Try to validate frame with the wrong data type
		fakeFrame = tbnWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.iq = fakeFrame.data.iq.real
		self.assertRaises(errors.invalidDataType, fakeFrame.isValid, raiseErrors=True)
		
		# Try to validate frame with the wrong data size
		fakeFrame = tbnWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.iq = None
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		fakeFrame = tbnWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.iq = fakeFrame.data.iq[0:50]
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		
	def test_frame_header_errors(self):
		"""Test the header error scenarios when validating a TBN SimFrame."""
		
		# Read in a TBN frame from the test file
		fh = open(tbnFile, 'rb')
		origFrame = tbnReader.readFrame(fh)
		fh.close()
		
		# Try to validate frame with the wrong stand number
		fakeFrame = tbnWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.stand = 300
		self.assertRaises(errors.invalidStand, fakeFrame.isValid, raiseErrors=True)

	def tearDown(self):
		"""Remove the test path directory and its contents"""

		tempFiles = os.listdir(self.testPath)
		for tempFile in tempFiles:
			os.unlink(os.path.join(self.testPath, tempFile))
		os.rmdir(self.testPath)
		self.testPath = None


class fake_DRX_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.sim.drx
	module."""

	testPath = None

	def setUp(self):
		"""Turn off all numpy warnings and create the temporary file directory."""

		numpy.seterr(all='ignore')
		self.testPath = tempfile.mkdtemp(prefix='test-fakedata-', suffix='.tmp')

	def test_sim_frame(self):
		"""Test the drx.SimFrame class."""
		
		# Read in a DRX frame from the test file
		fh = open(drxFile, 'rb')
		origFrame = drxReader.readFrame(fh)
		fh.close()
		
		fakeFrame = drxWriter.SimFrame()
		fakeFrame.loadFrame(origFrame)
		# Test the validity of the SimFrame
		self.assertTrue(fakeFrame.isValid())

	def test_write_frame(self):
		"""Test that the DRX data writer works."""
		
		testFile = os.path.join(self.testPath, 'drx-test-W.dat')

		# Read in a TBN frame from the test file
		fh = open(drxFile, 'rb')
		origFrame = drxReader.readFrame(fh)
		fh.close()
		
		# Write the data to a TBN test frame
		fh = open(testFile, 'wb')
		rawFrame = drxWriter.frame2frame(origFrame)
		rawFrame.tofile(fh)
		fh.close()
		
		# Read in the 
		fh = open(testFile, 'rb')
		fakeFrame = drxReader.readFrame(fh)
		fh.close()

		# Test values returned by info functions
		self.assertEqual(fakeFrame.parseID()[0], origFrame.parseID()[0])
		self.assertEqual(fakeFrame.parseID()[1], origFrame.parseID()[1])
		self.assertEqual(fakeFrame.parseID()[2], origFrame.parseID()[2])
		self.assertAlmostEqual(fakeFrame.getSampleRate(), origFrame.getSampleRate(), 4)
		
		# Test raw header values
		self.assertEqual(fakeFrame.header.secondsCount, origFrame.header.secondsCount)
		self.assertEqual(fakeFrame.header.decimation, origFrame.header.decimation)
		self.assertEqual(fakeFrame.header.timeOffset, origFrame.header.timeOffset)
		
		# Test raw data values
		self.assertEqual(fakeFrame.data.timeTag, origFrame.data.timeTag)
		self.assertEqual(fakeFrame.data.flags, origFrame.data.flags)
		for i in range(4096):
			self.assertEqual(fakeFrame.data.iq[i].real, origFrame.data.iq[i].real)
			self.assertEqual(fakeFrame.data.iq[i].imag, origFrame.data.iq[i].imag)

	def test_frame_data_errors(self):
		"""Test the data error scenarios when validating a DRX SimFrame."""
		
		# Read in a DRX frame from the test file
		fh = open(drxFile, 'rb')
		origFrame = drxReader.readFrame(fh)
		fh.close()
		
		# Try to validate frame with the wrong data type
		fakeFrame = drxWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.iq = fakeFrame.data.iq.real
		self.assertRaises(errors.invalidDataType, fakeFrame.isValid, raiseErrors=True)
		
		# Try to validate frame with the wrong data size
		fakeFrame = drxWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.iq = None
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		fakeFrame = drxWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.iq = fakeFrame.data.iq[0:50]
		self.assertRaises(errors.invalidDataSize, fakeFrame.isValid, raiseErrors=True)
		
	def test_frame_header_errors(self):
		"""Test the header error scenarios when validating a DRX SimFrame."""
		
		# Read in a DRX frame from the test file
		fh = open(drxFile, 'rb')
		origFrame = drxReader.readFrame(fh)
		fh.close()
		
		# Try to validate frame with the wrong beam number
		fakeFrame = drxWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.beam = 5
		self.assertRaises(errors.invalidBeam, fakeFrame.isValid, raiseErrors=True)

		# Try to validate frame with the wrong tuning number
		fakeFrame = drxWriter.SimFrame()
		fakeFrame.loadFrame(copy.deepcopy(origFrame))
		fakeFrame.tune = 3
		self.assertRaises(errors.invalidTune, fakeFrame.isValid, raiseErrors=True)

	def tearDown(self):
		"""Remove the test path directory and its contents"""

		tempFiles = os.listdir(self.testPath)
		for tempFile in tempFiles:
			os.unlink(os.path.join(self.testPath, tempFile))
		os.rmdir(self.testPath)
		self.testPath = None


class fakedata_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(fake_S60_tests))
		self.addTests(loader.loadTestsFromTestCase(fake_TBW_tests))
		self.addTests(loader.loadTestsFromTestCase(fake_TBN_tests))
		self.addTests(loader.loadTestsFromTestCase(fake_DRX_tests))


if __name__ == '__main__':
	unittest.main()

