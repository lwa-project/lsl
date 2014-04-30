#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in TBN data and writing it to a TS-FITS file."""

import os
import sys
import time
from lsl.reader import tbn
from lsl.reader import errors
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
	fh.seek(-tbn.FrameSize, 1)

	nFpO = tbn.getFramesPerObs(fh)
	print "Samples per observations: %i in x pol., %i in y pol." % nFpO
	nFpO = nFpO[0] + nFpO[1]

	sampleRate = tbn.getSampleRate(fh, nFrames=nFpO)
	print "Filter code is: %i" % tbn.getSampleRate(fh, nFrames=nFpO, FilterCode=True)
	print "Sampling rate is: %i Hz" % sampleRate

	tStart = time.time()

	# Create a new FITS file with the name 'tbn-tsfits.fits'
	fitsFile = tsfits.TBN('tbn-tsfits.fits')

	nSamples = 340000

	blocks = []
	count = {}
	syncCount = 0
	masterCount = 0
	for i in range(nSamples):
		try:
			cFrame = tbn.readBlock(fh, nFrames=nFpO, SampleRate=sampleRate)
		except errors.eofError:
			break
		except errors.syncError:
			syncCount = syncCount + 1
			continue
		except errors.numpyError:
			break
		blocks.append(cFrame)

		for frame in blocks[-1].x:
			if frame is None:
				print "sync error"
				continue
			stand, pol = frame.parseID()
			if stand not in count.keys():
				count[stand] = 0
			count[stand] = count[stand] + 1

		for xFrame, yFrame in zip(cFrame.x, cFrame.y):
			if xFrame is not None:
				fitsFile.addStandData(xFrame)
			if yFrame is not None:
				fitsFile.addStandData(yFrame)

		masterCount = masterCount + 1
		#print cFrame.header.parseID(),cFrame.data.timeTag

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (nFpO*nSamples, (tEnd-tStart), nFpO*nSamples/(tEnd-tStart))

	fh.close()
	fitsFile.close()
	fitsFile.info()

	# Summary information about the file that was just read in
	print "Summary:"
	for stand in sorted(count.keys()):
		print "Stand: %2i, Frames: %5i" % (stand, count[stand])
	print "Sync Errors: %5i" % syncCount


if __name__ == "__main__":
	main(sys.argv[1:])
