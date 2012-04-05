# -*- coding: utf-8 -*-

"""Unit test for lsl.reader modules"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import drx
from lsl.reader import drspec
from lsl.reader import errors


__revision__ = "$ Revision: 2 $"
__version__  = "0.4"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')
drxFile = os.path.join(dataPath, 'tests', 'drx-test.dat')
drspecFile = os.path.join(dataPath, 'tests', 'drspec-test.dat')


class reader_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.reader
	modules."""

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

		# Last frame should be an error (errors.eofError)
		self.assertRaises(errors.eofError, tbw.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(tbwFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, tbw.readFrame, fh)
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
		self.assertFalse(0 > frames[0])
		self.assertTrue(frames[-1] >= frames[0])
		self.assertFalse(frames[-1] <= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertFalse(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])
		
	def test_tbw_sort(self):
		"""Test sorting TBW frames by time tags."""
		
		fh = open(tbwFile, 'rb')
		# Frames 1 through 3
		frames = []
		for i in range(1,4):
			frames.append(tbw.readFrame(fh))
		fh.close()
		
		frames.sort()
		frames = frames[::-1]
		
		for i in xrange(1,len(frames)):
			self.assertTrue( frames[i-1] >= frames[i] )

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
		# First frame is really TBN and stores the correct IDs
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
		
	def test_tbn_errors(self):
		"""Test reading in all frames from a truncated TBN file."""

		fh = open(tbnFile, 'rb')
		# Frames 1 through 29
		for i in range(1,30):
			frame = tbn.readFrame(fh)

		# Last frame should be an error (errors.eofError)
		self.assertRaises(errors.eofError, tbn.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(tbnFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, tbn.readFrame, fh)
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

		self.assertTrue(0 < frames[0])
		self.assertFalse(0 > frames[0])
		self.assertTrue(frames[-1] >= frames[0])
		self.assertFalse(frames[-1] <= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertFalse(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])
		
	def test_tbn_sort(self):
		"""Test sorting TBN frames by time tags."""
		
		fh = open(tbnFile, 'rb')
		# Frames 1 through 29
		frames = []
		for i in range(1,30):
			frames.append(tbn.readFrame(fh))
		fh.close()
		
		frames.sort()
		frames = frames[::-1]
		
		for i in xrange(1,len(frames)):
			self.assertTrue( frames[i-1] >= frames[i] )

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

	def test_drx_read(self):
		"""Test reading in a frame from a DRX file."""

		fh = open(drxFile, 'rb')
		# First frame is really DRX and stores the IDs
		frame1 = drx.readFrame(fh)
		beam, tune, pol = frame1.parseID()
		self.assertEqual(beam, 2)
		self.assertEqual(tune, 2)
		self.assertEqual(pol,  1)
		# Second frame
		frame2 = drx.readFrame(fh)
		beam, tune, pol = frame2.parseID()
		self.assertEqual(beam, 2)
		self.assertEqual(tune, 1)
		self.assertEqual(pol,  0)
		fh.close()
		
	def test_drx_read_block(self):
		"""Test reading in a block of DRX frames."""
		
		fh = open(drxFile, 'rb')
		cFrame = drx.readFrame(fh)
		b,t,p = cFrame.parseID()
		while 2*(t-1)+p != 0:
			cFrame = drx.readFrame(fh)
			b,t,p = cFrame.parseID()
			
		# Read a block
		cFrames = drx.readBlock(fh)
		
		# Make sure everthing is in the right place
		self.assertEqual(cFrames.x1.parseID(), (2,1,0))
		self.assertEqual(cFrames.x2.parseID(), (2,2,0))
		self.assertEqual(cFrames.y1.parseID(), (2,1,1))
		self.assertEqual(cFrames.y2.parseID(), (2,2,1))
		
		fh.close()

	def test_drx_errors(self):
		"""Test reading in all frames from a truncated DRX file."""

		fh = open(drxFile, 'rb')
		# Frames 1 through 16
		for i in range(1,17):
			frame = drx.readFrame(fh)

		# Last frame should be an error (errors.eofError)
		self.assertRaises(errors.eofError, drx.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(drxFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, drx.readFrame, fh)
		fh.close()

	def test_drx_beam(self):
		"""Test finding out how many beams are present in a DRX file."""

		fh = open(drxFile, 'rb')
		nBeam = drx.getBeamCount(fh)
		self.assertEqual(nBeam, 1)
		fh.close()

	def test_drx_block(self):
		"""Test finding out how many tunings/pols. per beam are in a DRX file."""

		fh = open(drxFile, 'rb')
		b1, b2, b3, b4 = drx.getFramesPerObs(fh)
		self.assertEqual(b1, 0)
		self.assertEqual(b2, 4)
		self.assertEqual(b3, 0)
		self.assertEqual(b4, 0)
		fh.close()
		
	def test_drx_rate(self):
		"""Test finding out the DRX sample rate."""
		
		fh = open(drxFile, 'rb')
		cFrame = drx.readFrame(fh)
		fh.seek(0)
		
		# Sample rate
		self.assertEqual(cFrame.getSampleRate(), drx.getSampleRate(fh))
		
		# Filter code
		self.assertEqual(cFrame.getFilterCode(), drx.getSampleRate(fh, FilterCode=True))
		fh.close()

	def test_drx_comps(self):
		"""Test the DRX frame comparison operators (>, <, etc.) for time tags."""

		fh = open(drxFile, 'rb')
		# Frames 1 through 10
		frames = []
		for i in range(1,11):
			frames.append(drx.readFrame(fh))
		fh.close()

		self.assertTrue(0 < frames[0])
		self.assertFalse(0 > frames[0])
		self.assertTrue(frames[-1] >= frames[0])
		self.assertFalse(frames[-1] <= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertFalse(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])
		
	def test_drx_sort(self):
		"""Test sorting DRX frames by time tags."""
		
		fh = open(drxFile, 'rb')
		# Frames 1 through 10
		frames = []
		for i in range(1,11):
			frames.append(drx.readFrame(fh))
		fh.close()
		
		frames.sort()
		frames = frames[::-1]
		
		for i in xrange(1,len(frames)):
			self.assertTrue( frames[i-1] >= frames[i] )

	def test_drx_math(self):
		"""Test mathematical operations on DRX frame data via frames."""

		fh = open(drxFile, 'rb')
		# Frames 1 through 10
		frames = []
		for i in range(1,11):
			frames.append(drx.readFrame(fh))
		fh.close()

		# Multiplication
		frameT = frames[0] * 2.0
		for i in range(4096):
			self.assertAlmostEqual(frameT.data.iq[i], 2*frames[0].data.iq[i], 6)
		frameT *= 2.0
		for i in range(4096):
			self.assertAlmostEqual(frameT.data.iq[i], 4*frames[0].data.iq[i], 6)
		frameT = frames[0] * frames[1]
		for i in range(4096):
			self.assertAlmostEqual(frameT.data.iq[i], frames[0].data.iq[i]*frames[1].data.iq[i], 6)
		
		# Addition
		frameA = frames[0] + 2.0
		for i in range(4096):
			self.assertAlmostEqual(frameA.data.iq[i], 2+frames[0].data.iq[i], 6)
		frameA += 2.0
		for i in range(4096):
			self.assertAlmostEqual(frameA.data.iq[i], 4+frames[0].data.iq[i], 6)
		frameA = frames[0] + frames[1]
		for i in range(4096):
			self.assertAlmostEqual(frameA.data.iq[i], frames[0].data.iq[i]+frames[1].data.iq[i], 6)

	### DR Spectrometer ###

	def test_drspec_read(self):
		"""Test reading in a frame from a DR spectrometer file."""

		fh = open(drspecFile, 'rb')
		# First frame is really DR spectrometer and stores the IDs
		frame1 = drspec.readFrame(fh)
		beam = frame1.parseID()
		self.assertEqual(beam, 3)

		# Second frame
		frame2 = drspec.readFrame(fh)
		beam = frame2.parseID()
		self.assertEqual(beam, 3)
		fh.close()

	def test_drspec_errors(self):
		"""Test reading in all frames from a truncated DR spectrometer file."""

		fh = open(drspecFile, 'rb')
		# Frames 1 through 8
		for i in range(1,8):
			frame = drspec.readFrame(fh)

		# Last frame should be an error (errors.eofError)
		self.assertRaises(errors.eofError, drspec.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(drspecFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, drspec.readFrame, fh)
		fh.close()
		
	def test_drspec_rate(self):
		"""Test finding out the DR spectrometer sample rate."""
		
		fh = open(drspecFile, 'rb')
		cFrame = drspec.readFrame(fh)
		fh.seek(0)
		
		# Sample rate
		self.assertEqual(cFrame.getSampleRate(), drspec.getSampleRate(fh))
		
		# Filter code
		self.assertEqual(cFrame.getFilterCode(), drspec.getSampleRate(fh, FilterCode=True))
		fh.close()

	def test_drspec_comps(self):
		"""Test the DR spectrometer frame comparison operators (>, <, etc.) for time tags."""

		fh = open(drspecFile, 'rb')
		# Frames 1 through 7
		frames = []
		for i in range(1,8):
			frames.append(drspec.readFrame(fh))
		fh.close()

		self.assertTrue(0 < frames[0])
		self.assertFalse(0 > frames[0])
		self.assertTrue(frames[-1] >= frames[0])
		self.assertFalse(frames[-1] <= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertFalse(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])
		
	def test_drspec_sort(self):
		"""Test sorting DR spectrometer frames by time tags."""
		
		fh = open(drspecFile, 'rb')
		# Frames 1 through 7
		frames = []
		for i in range(1,8):
			frames.append(drspec.readFrame(fh))
		
		frames.sort()
		frames = frames[::-1]
		
		for i in xrange(1,len(frames)):
			self.assertTrue( frames[i-1] >= frames[i] )
		fh.close()

	def test_drspec_math(self):
		"""Test mathematical operations on DR spectrometer frame data via frames."""

		fh = open(drspecFile, 'rb')
		# Frames 1 through 7
		frames = []
		for i in range(1,8):
			frames.append(drspec.readFrame(fh))
		fh.close()

		nPts = frames[0].data.X0.size

		# Multiplication
		frameT = frames[0] * 2.0
		for i in xrange(nPts):
			self.assertAlmostEqual(frameT.data.X0[i], 2*frames[0].data.X0[i], 6)
		frameT *= 2.0
		for i in xrange(nPts):
			self.assertAlmostEqual(frameT.data.X1[i], 4*frames[0].data.X1[i], 6)
		frameT = frames[0] * frames[1]
		for i in xrange(nPts):
			self.assertAlmostEqual(frameT.data.Y0[i], frames[0].data.Y0[i]*frames[1].data.Y0[i], 6)
		
		# Addition
		frameA = frames[0] + 2.0
		for i in xrange(nPts):
			self.assertAlmostEqual(frameA.data.X0[i], 2+frames[0].data.X0[i], 6)
		frameA += 2.0
		for i in xrange(nPts):
			self.assertAlmostEqual(frameA.data.X1[i], 4+frames[0].data.X1[i], 6)
		frameA = frames[0] + frames[1]
		for i in xrange(nPts):
			self.assertAlmostEqual(frameA.data.Y0[i], frames[0].data.Y0[i]+frames[1].data.Y0[i], 6)

class reader_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(reader_tests)) 


if __name__ == '__main__':
	unittest.main()
