# -*- coding: utf-8 -*-

"""Unit test for the lsl.writer.tsfits modules."""

import os
import time
import unittest
import tempfile
import numpy
import pyfits

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import tbw
from lsl.reader import tbn
from lsl.reader import errors
from lsl.writer import vdif


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"

tbwFile = os.path.join(dataPath, 'tests', 'tbw-test.dat')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')


class vdif_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.writer.vdif
	module."""

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
				frame.setCentralFreq(40e6)
				frame.setGain(22)
				frame.setSampleRate(100000)

		return frames

	def test_vdif_real(self):
		"""Test writing real data to VDIF format."""

		# Setup the file names
		testPath = tempfile.mkdtemp(prefix='tsfits-', suffix='.tmp')
		testFile = os.path.join(testPath, 'ts-tbw-testS.fits')

		# Get some data
		frames = self.__getTBW()

		# Write the data
		fh = open(testFile, 'wb')
		for frame in frames:
			vFrame = vdif.Frame(frame.parseID(), frame.getTime(), bits=12, data=frame.data.xy[0,:])
			vFrame.writeRawFrame(fh)
		fh.close()

		# Read it back in
		fh = open(testFile, 'rb')
		junk = numpy.fromfile(fh, dtype=numpy.uint8, count=16)
		frameSize = (junk[10]<<16) | (junk[9]<<8) | (junk[8])
		dataBits = (junk[15]>>2) & 15
		dataSize = (frameSize - 4) * 8
		fh.seek(0)
		vFrames = []
		for i in range(8):
			vFrames.append( numpy.fromfile(fh, dtype=numpy.uint8, count=frameSize*8) )
		for vFrame, tFrame in zip(vFrames, frames):
			data = numpy.zeros(400, dtype=numpy.int16)
			j = 0
			for i in range(32,vFrame.shape[0],4):
				word = vFrame[i] | (vFrame[i+1]<<8) | (vFrame[i+2]<<16) | (vFrame[i+3]<<24)
				for k in range(32/dataBits):
					data[j] = ((word>>(dataBits*k)) & (2**dataBits-1)) - 2**(dataBits-1)
					j = j + 1
			
			for v,t in zip(data, tFrame.data.xy[0,:]):
				self.assertEqual(v, t)

		fh.close()
		os.unlink(testFile)
		os.rmdir(testPath)

	def test_vdif_complex(self):
		"""Test writing complex data to VIDF format."""

		# Setup the file names
		testPath = tempfile.mkdtemp(prefix='tsfits-', suffix='.tmp')
		testFile = os.path.join(testPath, 'ts-tbw-testS.fits')

		# Get some data
		frames = self.__getTBN()
		
		# Write the data
		fh = open(testFile, 'wb')
		for frame in frames:
			stand, pol = frame.parseID()
			if pol == 1:
				continue
			vFrame = vdif.Frame(stand, frame.getTime(), bits=8, data=frame.data.iq)
			vFrame.writeRawFrame(fh)
		fh.close()

		# Read it back in
		fh = open(testFile, 'rb')
		junk = numpy.fromfile(fh, dtype=numpy.uint8, count=16)
		frameSize = (junk[10]<<16) | (junk[9]<<8) | (junk[8])
		dataBits = (junk[15]>>2) & 15
		dataSize = (frameSize - 4) * 8
		fh.seek(0)
		vFrames = []
		for i in range(8):
			vFrames.append( numpy.fromfile(fh, dtype=numpy.uint8, count=frameSize*8) )
		for vFrame, tFrame in zip(vFrames, frames[::2]):
			data = numpy.zeros(512, dtype=numpy.singlecomplex)
			j = 0
			for i in range(32,vFrame.shape[0],4):
				word = vFrame[i] | (vFrame[i+1]<<8) | (vFrame[i+2]<<16) | (vFrame[i+3]<<24)
				for k in range(32/dataBits):
					if k % 2 == 0:
						data.real[j] = ((word>>(dataBits*k)) & (2**dataBits-1)) - 2**(dataBits-1)
					else:
						data.imag[j] = ((word>>(dataBits*k)) & (2**dataBits-1)) - 2**(dataBits-1)
						j = j + 1

			for v,t in zip(data, tFrame.data.iq):
				self.assertAlmostEqual(v, t, 6)
		
		fh.close()
		os.unlink(testFile)
		os.rmdir(testPath)


class  vdif_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.writer.vdif units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(vdif_tests)) 


if __name__ == '__main__':
	unittest.main()
