# -*- coding: utf-8 -*-

"""Unit test for lsl.reader modules that are for ADP stations"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import tbf
from lsl.reader import errors


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


tbfFile = os.path.join(dataPath, 'tests', 'tbf-test.dat')


class reader_adp_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.reader
	modules."""
	
	### TBF ###
	
	def test_tbf_read(self):
		"""Test reading in a frame from a TBF file."""
		
		fh = open(tbfFile, 'rb')
		# First frame is really TBF and stores the first channel
		frame1 = tbf.readFrame(fh)
		self.assertTrue(frame1.header.isTBF())
		self.assertEqual(frame1.header.firstChan, 2348)
		# Second frame
		frame2 = tbf.readFrame(fh)
		self.assertTrue(frame2.header.isTBF())
		self.assertEqual(frame2.header.firstChan, 2360)
		fh.close()
		
	def test_tbf_errors(self):
		"""Test reading errors."""
		
		fh = open(tbfFile, 'rb')
		# Frames 1 through 5
		for i in range(1,6):
			frame = tbf.readFrame(fh)
			
		# Last frame should be an error (errors.eofError)
		self.assertRaises(errors.eofError, tbf.readFrame, fh)
		fh.close()
		
		# If we offset in the file by 1 byte, we should be a 
		# sync error (errors.syncError).
		fh = open(tbfFile, 'rb')
		fh.seek(1)
		self.assertRaises(errors.syncError, tbf.readFrame, fh)
		fh.close()
		
	def test_tbf_comps(self):
		"""Test the TBF frame comparison operators (>, <, etc.) for time tags."""
		
		fh = open(tbfFile, 'rb')
		# Frames 1 through 4
		frames = []
		for i in range(1,5):
			frames.append(tbf.readFrame(fh))
		fh.close()
		
		self.assertTrue(0 < frames[0])
		self.assertFalse(0 > frames[0])
		self.assertTrue(frames[-1] >= frames[0])
		self.assertTrue(frames[-1] <= frames[0])
		self.assertTrue(frames[0] == frames[0])
		self.assertTrue(frames[0] == frames[-1])
		self.assertFalse(frames[0] != frames[0])
		
	def test_tbf_sort(self):
		"""Test sorting TBF frames by time tags."""
		
		fh = open(tbfFile, 'rb')
		# Frames 1 through 3
		frames = []
		for i in range(1,4):
			frames.append(tbf.readFrame(fh))
		fh.close()
		
		frames.sort()
		frames = frames[::-1]
		
		for i in xrange(1,len(frames)):
			self.assertTrue( frames[i-1] >= frames[i] )
			
	def test_tbf_math(self):
		"""Test mathematical operations on TBF frame data via frames."""
		
		fh = open(tbfFile, 'rb')
		# Frames 1 through 3
		frames = []
		for i in range(1,4):
			frames.append(tbf.readFrame(fh))
		fh.close()
		
		# Multiplication
		frameT = frames[0] * 2.0
		for i in range(12*256*2):
			c = i / 2 / 256
			s = i / 2 % 256
			p = i % 2
			self.assertAlmostEqual(frameT.data.fDomain[c,s,p], 2*frames[0].data.fDomain[c,s,p], 2)
		frameT *= 2.0
		for i in range(12*256*2):
			c = i / 2 / 256
			s = i / 2 % 256
			p = i % 2
			self.assertAlmostEqual(frameT.data.fDomain[c,s,p], 4*frames[0].data.fDomain[c,s,p], 2)
		frameT = frames[0] * frames[1]
		for i in range(12*256*2):
			c = i / 2 / 256
			s = i / 2 % 256
			p = i % 2
			self.assertAlmostEqual(frameT.data.fDomain[c,s,p], frames[0].data.fDomain[c,s,p]*frames[1].data.fDomain[c,s,p], 2)
			
		# Addition
		frameA = frames[0] + 2.0
		for i in range(800):
			c = i / 2 / 256
			s = i / 2 % 256
			p = i % 2
			self.assertAlmostEqual(frameA.data.fDomain[c,s,p], 2+frames[0].data.fDomain[c,s,p], 2)
		frameA += 2.0
		for i in range(800):
			c = i / 2 / 256
			s = i / 2 % 256
			p = i % 2
			self.assertAlmostEqual(frameA.data.fDomain[c,s,p], 4+frames[0].data.fDomain[c,s,p], 2)
		frameA = frames[0] + frames[1]
		for i in range(800):
			c = i / 2 / 256
			s = i / 2 % 256
			p = i % 2
			self.assertAlmostEqual(frameA.data.fDomain[c,s,p], frames[0].data.fDomain[c,s,p]+frames[1].data.fDomain[c,s,p], 2)


class reader_adp_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader units 
	tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(reader_adp_tests)) 


if __name__ == '__main__':
	unittest.main()
