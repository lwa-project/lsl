#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in DRX data and writing it to a SD-FITS file."""

import os
import sys
from lsl.reader import drx
from lsl.writer import sdfits

import matplotlib.pyplot as plt

def main(args):
	nSamples = os.path.getsize(args[0]) / drx.FrameSize
	print "Samples in file: ", nSamples
	fh = open(args[0], "rb", buffering=drx.FrameSize)
	nFpO = getFramesPerObs(fh)
	nBeams = getBeamCount(fh)
	print "Beams: ", nBeams
	print "Frames per Observations: ", nFpO
	blockBuffer = []
	blocks = []

	tStart = time.time()

	nSamples = (nSamples/4/16)*16

	fig = plt.figure()

	for i in range(0,nSamples):
		currBlock = readDRXBlock(fh)
		blockBuffer.append(currBlock)

		if len(blockBuffer) == 16:
			avgBlock = averageObservations(blockBuffer)
			#avgBlock = averageObservations2(blockBuffer, timeAvg=16, chanAvg=2)
			blocks.append(avgBlock)
			blockBuffer = []
	
	nChan = blocks[0].x1.data.iq.shape[0]
	outSpec = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	outTime = numpy.zeros(nSamples/16)
	for row,block in zip(range(nSamples),blocks):
		outSpec[row,:] = block.x1.data.iq
		outTime[row] = block.x1.data.timeTag
	outSpec2 = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec2[row,:] = block.y1.data.iq
	outSpec3 = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec3[row,:] = block.x2.data.iq
	outSpec4 = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec4[row,:] = block.y2.data.iq

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (4*nSamples, (tEnd-tStart), 4*nSamples/(tEnd-tStart))
	
	writefits(outSpec, outTime)
	readfits('test-sdfits.fits')

	ax = fig.add_subplot(221)
	dB = outSpec - outSpec.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 1, Pol. 0')
	ax.set_ylabel('Channel')
	ax.axis('auto')
	
	ax = fig.add_subplot(222)
	dB = outSpec2 - outSpec2.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 1, Pol. 1')
	ax.axis('auto')

	ax = fig.add_subplot(223)
	dB = outSpec3 - outSpec3.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 2, Pol. 0')
	ax.set_xlabel('Time')
	ax.set_ylabel('Channel')
	ax.axis('auto')

	ax = fig.add_subplot(224)
	dB = outSpec4 - outSpec4.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 2, Pol. 1')
	ax.set_xlabel('Time')
	ax.axis('auto')

	plt.show()
	fig.savefig("readDRX.png")

	fh.close()

if __name__ == "__main__":
	main(sys.argv[1:])
