# -*- coding: utf-8 -*-

"""Unit test for lsl.reader modules"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import s60
from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import drx
from lsl.reader import errors
from lsl.reader import _gofast


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


s60File = os.path.join(dataPath, 'tests', 's60-test.dat')
tbwFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')


class reader_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.reader
	modules."""

	### S60 ###

	def test_s60_read(self):
		"""Test reading in a frame from a S60 file."""

		fh = open(s60File, 'rb')
		# First frame makes it in with the correct number of elements
		frame1 = s60.readFrame(fh)
		self.assertEqual(frame1.shape[0], 734)

		# Next try a chunk that is equivalent to 3 frames
		chunk1 = s60.readChunk(fh, Chunk=2202)
		self.assertEqual(chunk1.shape[0], 2202)
		fh.close()

	def test_s60_errors(self):
		"""Test reading in all frames from a truncated S60 file."""

		fh = open(s60File, 'rb')
		# Frames 1 through 6
		for i in range(1,7):
			frame = s60.readFrame(fh)

		# Last frame should be an error (errors.numpyError)
		self.assertRaises(errors.numpyError, s60.readFrame, fh)
		fh.close()

	### TBW ###

	def test_tbw_read(self):
		"""Test reading in a frame from a TBW file."""

		fh = open(tbwFile, 'rb')
		# First frame is really TBW and stores the correct stand ID
		frame1 = tbw.readFrame(fh)
		self.assertTrue(frame1.header.isTBW())
		self.assertEqual(frame1.parseID(), 2)
		# Second frame
		frame2 = tbw.readFrame(fh)
		self.assertTrue(frame2.header.isTBW())
		self.assertEqual(frame2.parseID(), 1)
		fh.close()
		
	def test_tbw_read_gofast(self):
		"""Test the Go Fast! TBW reader."""
		
		fh = open(tbwFile, 'rb')
		# First frame is really TBW and stores the correct stand ID
		fileLoc = fh.tell()
		frame1 = tbw.readFrame(fh)
		fh.seek(fileLoc)
		frame1G = _gofast.readTBW(fh, tbw.Frame())
		
		self.assertEqual(frame1.parseID(), frame1G.parseID())
		self.assertEqual(frame1.header.tbwID, frame1G.header.tbwID)
		self.assertEqual(frame1.header.frameCount, frame1G.header.frameCount)
		self.assertEqual(frame1.header.secondsCount, frame1G.header.secondsCount)
		self.assertEqual(frame1.data.timeTag, frame1G.data.timeTag)
		for i in xrange(400):
			self.assertEqual(frame1.data.xy[0,i], frame1G.data.xy[0,i])
			self.assertEqual(frame1.data.xy[1,i], frame1G.data.xy[1,i])
		
		# Second frame
		fileLoc = fh.tell()
		frame2 = tbw.readFrame(fh)
		fh.seek(fileLoc)
		frame2G = _gofast.readTBW(fh, tbw.Frame())
		
		self.assertEqual(frame2.parseID(), frame2G.parseID())
		self.assertEqual(frame2.header.tbwID, frame2G.header.tbwID)
		self.assertEqual(frame2.header.frameCount, frame2G.header.frameCount)
		self.assertEqual(frame2.header.secondsCount, frame2G.header.secondsCount)
		self.assertEqual(frame2.data.timeTag, frame2G.data.timeTag)
		for i in xrange(400):
			self.assertEqual(frame2.data.xy[0,i], frame2G.data.xy[0,i])
			self.assertEqual(frame2.data.xy[1,i], frame2G.data.xy[1,i])
			
		# Third frame
		fileLoc = fh.tell()
		frame3 = tbw.readFrame(fh)
		fh.seek(fileLoc)
		frame3G = _gofast.readTBW(fh, tbw.Frame())
		
		self.assertEqual(frame3.parseID(), frame3G.parseID())
		self.assertEqual(frame3.header.tbwID, frame3G.header.tbwID)
		self.assertEqual(frame3.header.frameCount, frame3G.header.frameCount)
		self.assertEqual(frame3.header.secondsCount, frame3G.header.secondsCount)
		self.assertEqual(frame3.data.timeTag, frame3G.data.timeTag)
		for i in xrange(400):
			self.assertEqual(frame3.data.xy[0,i], frame3G.data.xy[0,i])
			self.assertEqual(frame3.data.xy[1,i], frame3G.data.xy[1,i])
			
		# Fourth frame
		fileLoc = fh.tell()
		frame4 = tbw.readFrame(fh)
		fh.seek(fileLoc)
		frame4G = _gofast.readTBW(fh, tbw.Frame())
		
		self.assertEqual(frame4.parseID(), frame4G.parseID())
		self.assertEqual(frame4.header.tbwID, frame4G.header.tbwID)
		self.assertEqual(frame4.header.frameCount, frame4G.header.frameCount)
		self.assertEqual(frame4.header.secondsCount, frame4G.header.secondsCount)
		self.assertEqual(frame4.data.timeTag, frame4G.data.timeTag)
		for i in xrange(400):
			self.assertEqual(frame4.data.xy[0,i], frame4G.data.xy[0,i])
			self.assertEqual(frame4.data.xy[1,i], frame4G.data.xy[1,i])
		
		fh.close()

	def test_tbw_bits(self):
		"""Test getting the data bits from a TBW file."""

		fh = open(tbwFile, 'rb')
		# File contains 12-bit data, two ways
		self.assertEqual(tbw.getDataBits(fh), 12)
		frame1 = tbw.readFrame(fh)
		self.assertEqual(frame1.getDataBits(), 12)
		fh.close()

	def test_tbw_errors(self):
		"""Test reading errors."""

		fh = open(tbwFile, 'rb')
		# Frames 1 through 8
		for i in range(1,9):
			frame = tbw.readFrame(fh)

		# Last frame should be an error (errors.numpyError)
		self.assertRaises(errors.numpyError, tbw.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(tbwFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, tbw.readFrame, fh)
		fh.close()
		
	def test_tbw_errors_gofast(self):
		"""Test the Go Fast! TBW reader errors."""
		
		fh = open(tbwFile, 'rb')
		# Frames 1 through 8
		for i in range(1,9):
			frame = _gofast.readTBW(fh, tbw.Frame())

		# Last frame should be an error (_gofast.eofError)
		self.assertRaises(_gofast.eofError, _gofast.readTBW, fh, tbw.Frame())
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (_gofast.syncError).
		fh = open(tbwFile, 'rb')
		fh.seek(1)
		self.assertRaises(_gofast.syncError, _gofast.readTBW, fh, tbw.Frame())
		fh.close()

	def test_tbw_comps(self):
		"""Test the TBW frame comparison operators (>, <, etc.) for time tags."""

		fh = open(tbwFile, 'rb')
		# Frames 1 through 3
		frames = []
		for i in range(1,4):
			frames.append(tbw.readFrame(fh))
		fh.close()

		self.assertTrue(0 < frames[0])
		self.assertTrue(frames[-1] >= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertFalse(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])

	def test_tbw_math(self):
		"""Test mathematical operations on TBW frame data via frames."""

		fh = open(tbwFile, 'rb')
		# Frames 1 through 3
		frames = []
		for i in range(1,4):
			frames.append(tbw.readFrame(fh))
		fh.close()

		# Multiplication
		frameT = frames[0] * 2.0
		for i in range(800):
			self.assertAlmostEqual(frameT.data.xy[i%2, i/2], 2*frames[0].data.xy[i%2, i/2], 6)
		frameT *= 2.0
		for i in range(800):
			self.assertAlmostEqual(frameT.data.xy[i%2, i/2], 4*frames[0].data.xy[i%2, i/2], 6)
		frameT = frames[0] * frames[1]
		for i in range(800):
			self.assertAlmostEqual(frameT.data.xy[i%2, i/2], frames[0].data.xy[i%2, i/2]*frames[1].data.xy[i%2, i/2], 6)
		
		# Addition
		frameA = frames[0] + 2.0
		for i in range(800):
			self.assertAlmostEqual(frameA.data.xy[i%2, i/2], 2+frames[0].data.xy[i%2, i/2], 6)
		frameA += 2.0
		for i in range(800):
			self.assertAlmostEqual(frameA.data.xy[i%2, i/2], 4+frames[0].data.xy[i%2, i/2], 6)
		frameA = frames[0] + frames[1]
		for i in range(800):
			self.assertAlmostEqual(frameA.data.xy[i%2, i/2], frames[0].data.xy[i%2, i/2]+frames[1].data.xy[i%2, i/2], 6)

	### TBN ###

	def test_tbn_read(self):
		"""Test reading in a frame from a TBN file."""

		fh = open(tbnFile, 'rb')
		# First frame is really TBW and stores the correct stand ID
		frame1 = tbn.readFrame(fh)
		stand, pol = frame1.parseID()
		self.assertEqual(stand, 1)
		self.assertEqual(pol, 0)
		# Second frame
		frame2 = tbn.readFrame(fh)
		stand, pol = frame2.parseID()
		self.assertEqual(stand, 1)
		self.assertEqual(pol, 1)
		fh.close()
		
	def test_tbn_read_gofast(self):
		"""Test the Go Fast! TBN reader."""
		
		fh = open(tbnFile, 'rb')
		# First frame is really TBN and stores the correct stand ID
		fileLoc = fh.tell()
		frame1 = tbn.readFrame(fh)
		fh.seek(fileLoc)
		frame1G = _gofast.readTBN(fh, tbn.Frame())
		
		self.assertEqual(frame1.header.tbnID, frame1G.header.tbnID)
		self.assertEqual(frame1.header.frameCount, frame1G.header.frameCount)
		self.assertEqual(frame1.header.secondsCount, frame1G.header.secondsCount)
		self.assertEqual(frame1.data.timeTag, frame1G.data.timeTag)
		for i in xrange(512):
			self.assertEqual(frame1.data.iq[i].real, frame1G.data.iq[i].real)
			self.assertEqual(frame1.data.iq[i].imag, frame1G.data.iq[i].imag)
		
		# Second frame
		fileLoc = fh.tell()
		frame2 = tbn.readFrame(fh)
		fh.seek(fileLoc)
		frame2G = _gofast.readTBN(fh, tbn.Frame())
		
		self.assertEqual(frame2.header.tbnID, frame2G.header.tbnID)
		self.assertEqual(frame2.header.frameCount, frame2G.header.frameCount)
		self.assertEqual(frame2.header.secondsCount, frame2G.header.secondsCount)
		self.assertEqual(frame2.data.timeTag, frame2G.data.timeTag)
		for i in xrange(512):
			self.assertEqual(frame2.data.iq[i].real, frame2G.data.iq[i].real)
			self.assertEqual(frame2.data.iq[i].imag, frame2G.data.iq[i].imag)
		
		# Third frame
		fileLoc = fh.tell()
		frame3 = tbn.readFrame(fh)
		fh.seek(fileLoc)
		frame3G = _gofast.readTBN(fh, tbn.Frame())
		
		self.assertEqual(frame3.header.tbnID, frame3G.header.tbnID)
		self.assertEqual(frame3.header.frameCount, frame3G.header.frameCount)
		self.assertEqual(frame3.header.secondsCount, frame3G.header.secondsCount)
		self.assertEqual(frame3.data.timeTag, frame3G.data.timeTag)
		for i in xrange(512):
			self.assertEqual(frame3.data.iq[i].real, frame3G.data.iq[i].real)
			self.assertEqual(frame3.data.iq[i].imag, frame3G.data.iq[i].imag)
			
		# Fourth frame
		fileLoc = fh.tell()
		frame4 = tbn.readFrame(fh)
		fh.seek(fileLoc)
		frame4G = _gofast.readTBN(fh, tbn.Frame())
		
		self.assertEqual(frame4.header.tbnID, frame4G.header.tbnID)
		self.assertEqual(frame4.header.frameCount, frame4G.header.frameCount)
		self.assertEqual(frame4.header.secondsCount, frame4G.header.secondsCount)
		self.assertEqual(frame4.data.timeTag, frame4G.data.timeTag)
		for i in xrange(512):
			self.assertEqual(frame4.data.iq[i].real, frame4G.data.iq[i].real)
			self.assertEqual(frame4.data.iq[i].imag, frame4G.data.iq[i].imag)
		
		fh.close()

	def test_tbn_errors(self):
		"""Test reading in all frames from a truncated TBN file."""

		fh = open(tbnFile, 'rb')
		# Frames 1 through 29
		for i in range(1,30):
			frame = tbn.readFrame(fh)

		# Last frame should be an error (errors.numpyError)
		self.assertRaises(errors.numpyError, tbn.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(tbnFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, tbn.readFrame, fh)
		fh.close()
		
	def test_tbn_errors_gofast(self):
		"""Test reading in all frames from a truncated TBN file."""

		fh = open(tbnFile, 'rb')
		# Frames 1 through 29
		for i in range(1,30):
			frame = _gofast.readTBN(fh, tbn.Frame())

		# Last frame should be an error (_gofast.eofError)
		self.assertRaises(_gofast.eofError, _gofast.readTBN, fh, tbn.Frame())
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (_gofast.syncError).
		fh = open(tbnFile, 'rb')
		fh.seek(1)
		self.assertRaises(_gofast.syncError, _gofast.readTBN, fh, tbn.Frame())
		fh.close()

	def test_tbn_block(self):
		"""Test finding out how many stands are in a file."""

		fh = open(tbnFile, 'rb')
		nx, ny = tbn.getFramesPerObs(fh)
		self.assertEqual(nx, 10)
		self.assertEqual(ny, 10)
		fh.close()

	def test_tbn_rate(self):
		"""Test finding out the sample rate of a TBN file."""

		fh = open(tbnFile, 'rb')
		rate = tbn.getSampleRate(fh)
		self.assertEqual(rate, 100000)
		code = tbn.getSampleRate(fh, FilterCode=True)
		self.assertEqual(code, 7)
		fh.close()

	def test_tbn_comps(self):
		"""Test the TBN frame comparison operators (>, <, etc.) for time tags."""

		fh = open(tbnFile, 'rb')
		# Frames 1 through 29
		frames = []
		for i in range(1,30):
			frames.append(tbn.readFrame(fh))
		fh.close()

		self.assertTrue(frames[0] > 0)
		self.assertTrue(frames[-1] >= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertFalse(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])

	def test_tbn_math(self):
		"""Test mathematical operations on TBN frame data via frames."""

		fh = open(tbnFile, 'rb')
		# Frames 1 through 29
		frames = []
		for i in range(1,30):
			frames.append(tbn.readFrame(fh))
		fh.close()

		# Multiplication
		frameT = frames[0] * 2.0
		for i in range(512):
			self.assertAlmostEqual(frameT.data.iq[i], 2*frames[0].data.iq[i], 6)
		frameT *= 2.0
		for i in range(512):
			self.assertAlmostEqual(frameT.data.iq[i], 4*frames[0].data.iq[i], 6)
		frameT = frames[0] * frames[1]
		for i in range(512):
			self.assertAlmostEqual(frameT.data.iq[i], frames[0].data.iq[i]*frames[1].data.iq[i], 6)
		
		# Addition
		frameA = frames[0] + 2.0
		for i in range(512):
			self.assertAlmostEqual(frameA.data.iq[i], 2+frames[0].data.iq[i], 6)
		frameA += 2.0
		for i in range(512):
			self.assertAlmostEqual(frameA.data.iq[i], 4+frames[0].data.iq[i], 6)
		frameA = frames[0] + frames[1]
		for i in range(512):
			self.assertAlmostEqual(frameA.data.iq[i], frames[0].data.iq[i]+frames[1].data.iq[i], 6)

	### TBW/TBN Mix-up ###

	def test_tbw_tbn_catch(self):
		"""Test that tbw will not read tbn files and vice versa."""

		fh = open(tbnFile, 'rb')
		frame1 = tbw.readFrame(fh)
		self.assertFalse(frame1.header.isTBW())
		fh.close()

		fh = open(tbwFile, 'rb')
		frame1 = tbn.readFrame(fh)
		self.assertFalse(frame1.header.isTBN())
		fh.close()
		
	### DRX ###

	# Maybe someday...


class reader_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(reader_tests)) 


if __name__ == '__main__':
	unittest.main()
