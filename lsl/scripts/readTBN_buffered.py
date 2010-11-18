#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in TBW data and writing it to a TS-FITS file.
This version differs from the regular readTBN script in that it uses a frame
buffer to reorder out-of-order packets and dropped frames."""

import os
import sys
import time
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer as FrameBuffer
from lsl.writer import tsfits

import matplotlib.pyplot as plt

def main(args):
	nSamples = os.path.getsize(args[0]) / tbn.FrameSize
	print "Samples in file: ", nSamples
	fh = open(args[0], "rb", buffering=tbn.FrameSize)
	
	test = tbn.readFrame(fh)
	print "TBN Data:  %s" % test.header.isTBN()
	if not test.header.isTBN():
		raise notTBNError()
	fh.seek(0)

	nFpO = tbn.getFramesPerObs(fh)
	print "Samples per observations: %i in x pol., %i in y pol." % nFpO
	nFpO = nFpO[0] + nFpO[1]

	sampleRate = tbn.getSampleRate(fh, nFrames=nFpO)
	print "Filter code is: %i" % tbn.getSampleRate(fh, nFrames=nFpO, FilterCode=True)
	print "Sampling rate is: %i Hz" % sampleRate

	tStart = time.time()

	# Create the FrameBuffer instance
	buffer = FrameBuffer(stands=range(nFpO/2), pols=[0, 1])

	# Create a new FITS file with the name 'tbw.fits'
	fitsFile = tsfits.TBN('tbn-tsfits-test.fits')

	nSamples = 340000

	blocks = []
	count = {}
	syncCount = 0
	masterCount = 0
	for i in range(nSamples):
		try:
			cFrame = tbn.readFrame(fh, SampleRate=sampleRate)
		except errors.eofError:
			break
		except errors.syncError:
			syncCount = syncCount + 1
			continue
		except errors.numpyError:
			break
		buffer.append(cFrame)
		cFrames = buffer.get()
		if cFrames is None:
			continue
		blocks.append(cFrames)

		for frame in blocks[-1]:
			if frame is None:
				print "sync error"
				continue
			stand, pol = frame.parseID()
			if stand not in count.keys():
				count[stand] = 0
			count[stand] = count[stand] + 1

		for cFrame in blocks[-1]:
			fitsFile.addStandData(cFrame)

		masterCount = masterCount + 1
		#print cFrame.header.parseID(),cFrame.data.timeTag

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (nFpO*nSamples, (tEnd-tStart), nFpO*nSamples/(tEnd-tStart))

	fh.close()
	fitsFile.close()
	fitsFile.info()
	buffer.status()

	# Summary information about the file that was just read in
	print "Summary:"
	for stand in sorted(count.keys()):
		print "Stand: %2i, Frames: %5i" % (stand, count[stand])
	print "Sync Errors: %5i" % syncCount


if __name__ == "__main__":
	main(sys.argv[1:])
