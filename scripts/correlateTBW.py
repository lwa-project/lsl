#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
import ephem
import numpy
import getopt
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.common import stations
from lsl.common import dp as dp_common
from lsl.statistics import robust
from lsl.reader import tbw
from lsl.reader import errors
from lsl.correlator import uvUtils
from lsl.correlator import fx as fxc
from lsl.writer import fitsidi

from matplotlib import pyplot as plt


class UTC(tzinfo):
    """tzinfo object for UTC time."""

    def utcoffset(self, dt):
        return timedelta(0)

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return timedelta(0)


def usage(exitCode=None):
	print """correlateTBW.py - cross-correlate data in a TBW file

Usage: correlateTBW.py [OPTIONS] file

Options:
-h, --help             Display this help information
-f, --fft-length       Set FFT length (default = 512)
-s, --samples          Number of average visibilities to generate
                       (default = 10)
-q, --quiet            Run correlateTBW in silent mode
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['avgTime'] = 6
	config['LFFT'] = 512
	config['cFreq'] = 38.0e6
	config['samples'] = 10
	config['offset'] = 0
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hql:s:", ["help", "quiet", "fft-length=", "samples="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-q', '--quiet'):
			config['verbose'] = False
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-s', '--samples'):
			config['samples'] = int(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config


def processChunk(fh, site, stands, filename, LFFT=64, Overlap=1, SampleRate=dp_common.fS, ChunkSize=300, dataSize=400):
	"""Given a filehandle pointing to some TBN data and various parameters for
	the cross-correlation, write cross-correlate the data and save it to a file."""

	nFrames = 30000
	
	refTime = 0.0
	setTime = 0.0
	wallTime = time.time()
	for s in list(range(ChunkSize)):
		count = {}
		masterCount = 0
		iTime = 0
		data = numpy.zeros((20,12000000), dtype=numpy.int16)
		for i in range(10*nFrames):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbw.readFrame(fh)
			except errors.eofError:
				break
			except errors.syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FrameSize-1)
				continue
			except errors.numpyError:
				break

			stand = cFrame.parseID()
			if i == 0:
				setTime = cFrame.getTime()
				if s == 0:
					refTime = setTime

			if stand not in count.keys():
				count[stand] = 0
			data[2*(stand-1),  count[stand]*dataSize:(count[stand]+1)*dataSize] = numpy.squeeze(cFrame.data.xy[0,:])
			data[2*(stand-1)+1,count[stand]*dataSize:(count[stand]+1)*dataSize] = numpy.squeeze(cFrame.data.xy[1,:])
		
			count[stand] = count[stand] + 1
			masterCount = masterCount + 1

		setDT = datetime.utcfromtimestamp(setTime)
		setDT.replace(tzinfo=UTC())
		print "Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-refTime), setDT.strftime("%Y/%m/%d %H:%M:%S.%f"))
		freqYY, outYY = fxc.FXMaster(data[1::2,:], stands, LFFT=LFFT, Overlap=1, IncludeAuto=True, SampleRate=SampleRate)
		blList = uvUtils.getBaselines(stands, IncludeAuto=True, Indicies=False)

		toUse = numpy.where( (freqYY>10.0e6) & (freqYY<88.0e6) )
		toUse = toUse[0]

		if s  == 0:
			fits = fitsidi.IDI(filename, refTime=refTime)
			fits.setStokes(['xx'])
			fits.setFrequency(freqYY[toUse])
			fits.setGeometry(site, stands)

		obsTime = astro.unix_to_taimjd(setTime)
		fits.addDataSet(obsTime, 400*nFrames/SampleRate, blList, outYY[:,toUse])
		print "->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)/(s+1))

	fits.write()
	fits.close()
	del(fits)
	del(data)
	del(outYY)
	return True		


def main(args):
	# Parse command line options
	config = parseConfig(args)
	filename = config['args'][0]

	# Length of the FFT
	LFFT = config['LFFT']

	# Setup the LWA station information
	lwa1 = stations.lwa1

	fh = open(filename, "rb", buffering=tbw.FrameSize*10000)
	test = tbw.readFrame(fh)
	if not test.header.isTBW():
		raise errors.notTBWError()
	fh.seek(0)

	jd = astro.unix_to_utcjd(test.getTime())
	date = str(ephem.Date(jd - astro.DJD_OFFSET))
	sampleRate = dp_common.fS
	nInts = os.path.getsize(filename) / tbw.FrameSize / 300000

	# Get stands
	stands = lwa1.getStands(date)
	print stands

	# Number of frames to read in at once and average
	nFrames = 30000
	nSets = os.path.getsize(filename) / tbw.FrameSize / 300000

	print "TBW Data:  %s" % test.header.isTBW()
	print "Captures in file: %i (%.3f s)" % (nInts, nInts*30000*400/sampleRate)
	print "=="
	print "Station: %s" % lwa1.name
	print "Date observed: %s" % date
	print "Julian day: %.5f" % jd
	print "Integration Time: %.3f s" % (400*nFrames/sampleRate)
	print "Number of integrations in file: %i" % nSets
	print "=="

	# Align the file with the first full capture
	fh = open(filename, "rb")
	frame = tbw.readFrame(fh)
	skip = 0
	while frame.parseID() != 1 and frame.header.frameCount != 1:
		frame = tbw.readFrame(fh)
		skip += 1
	fh.seek(fh.tell() - tbw.FrameSize)
	print "Skipped %i frames" % skip

	# Get the data frame size
	dataSize = 400
	if frame.getDataBits() == 4:
		dataSize = 1200

	
	if config['samples'] > nSets:
		config['samples'] = nSets

	s = 0
	leftToDo = config['samples']
	while leftToDo > 0:
		fitsFilename = "TEST.FITS_%i" % (s+1)
		if leftToDo > 300:
			chunk = 300
		else:
			chunk = leftToDo
		processChunk(fh, lwa1, stands, fitsFilename, LFFT=config['LFFT'], Overlap=1, SampleRate=sampleRate,
				ChunkSize=chunk, dataSize=dataSize)

		s = s + 1
		leftToDo = leftToDo - chunk

	fh.close()


if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])

