# -*- coding: utf-8 -*-

"""
Unit test for the lsl.reader.buffer module.
"""

import os
import unittest

from lsl.common.paths import dataBuild as dataPath
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader import buffer


__revision__ = "$ Revision: 1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"


tbnFile = os.path.join(dataPath, 'tests', 'tbn-test.dat')


class buffer_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.reader.buffer
	module."""
	
	def test_tbn_default(self):
		"""Test the TBN ring buffer with the default values."""
		
		fh = open(tbnFile, 'rb')
		nFpO = tbn.getFramesPerObs(fh)
		nFpO = nFpO[0] + nFpO[1]
		
		# Create the FrameBuffer instance
		frameBuffer = buffer.TBNFrameBuffer(stands=range(1,nFpO/2+1), pols=[0, 1])
		
		# Go
		while True:
			try:
				cFrame = tbn.readFrame(fh)
			except errors.eofError:
				break
			except errors.syncError:
				continue

			frameBuffer.append(cFrame)
			cFrames = frameBuffer.get()

			if cFrames is None:
				continue
			
		fh.close()
		
		# Make sure we have the right number of frames in the buffer
		nFrames = 0
		for key in frameBuffer.buffer.keys():
			nFrames = nFrames + len(frameBuffer.buffer[key])
		self.assertEqual(nFrames, 29)
		
		# Make sure nothing has happened that shouldn't have
		self.assertEqual(frameBuffer.full,    0)
		self.assertEqual(frameBuffer.partial, 0)
		self.assertEqual(frameBuffer.missing, 0)
		self.assertEqual(frameBuffer.dropped, 0)
		
		# Make sure we have the right keys
		for key in frameBuffer.buffer.keys():
			self.assertTrue(key in (119196674956800, 119196675960320))
		
		# Make sure the buffer keys have the right sizes
		self.assertEqual(len(frameBuffer.buffer[119196674956800]), 20)
		self.assertEqual(len(frameBuffer.buffer[119196675960320]),  9)
	
	
	def test_buffer_small(self):
		"""Test a small version of the TBN ring buffer."""
		
		fh = open(tbnFile, 'rb')
		nFpO = tbn.getFramesPerObs(fh)
		nFpO = nFpO[0] + nFpO[1]
		
		# Create the FrameBuffer instance
		frameBuffer = buffer.TBNFrameBuffer(stands=range(1,nFpO/2+1), pols=[0, 1], nSegments=1)
		
		# Go
		while True:
			try:
				cFrame = tbn.readFrame(fh)
			except errors.eofError:
				break
			except errors.syncError:
				continue

			frameBuffer.append(cFrame)
			cFrames = frameBuffer.get()

			if cFrames is None:
				continue
			
			# Make sure the dump has one of the expected time tags
			self.assertTrue(cFrames[0].data.timeTag in (119196674956800, 119196675960320))
			
			# Make sure it has the right number of frames
			self.assertEqual(len(cFrames), nFpO)
		
		fh.close()
		
	def test_buffer_reorder(self):
		"""Test the reorder function of the TBN ring buffer."""
		
		fh = open(tbnFile, 'rb')
		nFpO = tbn.getFramesPerObs(fh)
		nFpO = nFpO[0] + nFpO[1]
		
		# Create the FrameBuffer instance
		frameBuffer = buffer.TBNFrameBuffer(stands=range(1,nFpO/2+1), pols=[0, 1], nSegments=1, ReorderFrames=True)
		
		# Go
		while True:
			try:
				cFrame = tbn.readFrame(fh)
			except errors.eofError:
				break
			except errors.syncError:
				continue

			frameBuffer.append(cFrame)
			cFrames = frameBuffer.get()

			if cFrames is None:
				continue
			
			# Make sure the dump has one of the expected time tags
			self.assertTrue(cFrames[0].data.timeTag in (119196674956800, 119196675960320))
			
			# Make sure it has the right number of frames
			self.assertEqual(len(cFrames), nFpO)
			
			# Check the order
			for i in xrange(1, len(cFrames)):
				pS, pP = cFrames[i-1].parseID()
				cS, cP = cFrames[i].parseID()
				
				pID = pS*2 + pP
				cID = cS*2 + cP
				self.assertTrue(cID > pID)
		
		fh.close()


class buffer_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.reader.buffer 
	unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(buffer_tests)) 


if __name__ == '__main__':
	unittest.main()
