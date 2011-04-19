#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in DRX data and writing it to a TS-FITS file."""

import os
import sys
import time
import ephem

from lsl.reader import drx
from lsl.writer import tsfits
from lsl.astro import unix_to_utcjd, DJD_OFFSET

def main(args):
	fh = open(args[0], "rb")
	nFramesFile = os.path.getsize(args[0]) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	
	fh.seek(0)
	srate = junkFrame.getSampleRate()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol
	
	# Date
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	
	# File summary
	print "Filename: %s" % args[0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"

	tStart = time.time()
	
	# Create a new FITS file with the name 'drx-tsfits.fits'
	fitsFile = tsfits.TBN('drx-tsfits.fits')
	
	nSamples = 3400

	count = {}
	syncCount = 0
	masterCount = 0
	for i in range(0,nSamples):
		currBlock = drx.readBlock(fh)

		for attr in ['x1', 'y1', 'x2', 'y2']:
			frame = getattr(currBlock, attr)
			if frame is None:
				syncCount += 1
				print "sync error"
				continue
			else:
				beam, pol, tune = frame.parseID()
				try:
					count[beam] += 1
				except KeyError:
					count[beam] = 1

				fitsFile.addStandData(frame)

		masterCount = masterCount + 1

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (beampols*nSamples, (tEnd-tStart), beampols*nSamples/(tEnd-tStart))

	fh.close()
	fitsFile.close()
	fitsFile.info()

	# Summary information about the file that was just read in
	print "Summary:"
	for beam in sorted(count.keys()):
		print "Beam: %2i, Frames: %5i" % (beam, count[beam])
	print "Sync Errors: %5i" % syncCount


if __name__ == "__main__":
	main(sys.argv[1:])
