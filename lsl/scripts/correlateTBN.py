#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script that reads in TBN data and runs a cross-correlation on it.  
The results are saved in the Miriad UV format."""

import os
import re
import sys
import time
import ephem
import numpy
import getopt
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.common import stations as lwa_common
from lsl.common import dp as dp_common
from lsl.statistics import robust
from lsl.reader import tbn
from lsl.reader import errors
from lsl.correlator import uvUtils
from lsl.correlator import fx as fxc
from lsl.writer import fitsidi

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


class UTC(tzinfo):
    """tzinfo object for UTC time."""

    def utcoffset(self, dt):
        return timedelta(0)

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return timedelta(0)



def usage(exitCode=None):
	print """correlateTBN.py - cross-correlate data in a TBN file

Usage: correlateTBN.py [OPTIONS] file

Options:
-h --help             Display this help information
-c --central-freq     Central frequency of the observations in MHz
-f --fft-length       Set FFT length (default = 512)
-t --avg-time         Window to average visibilities in time (seconds; 
                      default = 6 s)
-s --samples          Number of average visibilities to generate
                      (default = 10)
-o --offset           Seconds to skip from the beginning of the file
-q --quiet            Run correlateTBN in silent mode
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
		opts, arg = getopt.getopt(args, "hqc:l:t:s:o:", ["help", "quiet", "central-freq=", "fft-length=", "avg-time=", "samples=", "offset="])
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
		elif opt in ('-c', '--central-freq'):
			config['cFreq'] = float(value)*1e6
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-t', '--avg-time'):
			config['avgTime'] = int(value)
		elif opt in ('-s', '--samples'):
			config['samples'] = int(value)
		elif opt in ('-o', '--offset'):
			config['offset'] = int(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config


def processChunk(fh, site, stands, filename, intTime=6.0, LFFT=64, Overlap=1, CentralFreq=49.0e6, SampleRate=dp_common.fS, ChunkSize=300):
	"""Given a filehandle pointing to some TBN data and various parameters for
	the cross-correlation, write cross-correlate the data and save it to a file."""

	nFrames = int(intTime*SampleRate/512)

	refTime = 0.0
	setTime = 0.0
	wallTime = time.time()
	for s in list(range(ChunkSize)):
		count = {}
		masterCount = 0
		iTime = 0
		data = numpy.zeros((20,nFrames*512), dtype=numpy.complex64)
		for i in range(20*nFrames):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbn.readFrame(fh, Verbose=False, SampleRate=SampleRate, CentralFreq=CentralFreq, Gain=22)
			except errors.eofError:
				break
			except errors.syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFrameSize-1)
				continue
			except errors.numpyError:
				break

			stand,pol = cFrame.header.parseID()
			aStand = 2*(stand-1)+pol
			if i == 0:
				setTime = cFrame.getTime()
				if s == 0:
					refTime = setTime
			
			if aStand not in count.keys():
				count[aStand] = 0

			data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
			
			count[aStand] = count[aStand] + 1
			masterCount = masterCount + 1

		setDT = datetime.utcfromtimestamp(setTime)
		setDT.replace(tzinfo=UTC())
		print "Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-refTime), setDT.strftime("%Y/%m/%d %H:%M:%S.%f"))
		freqYY, outYY = fxc.FXMaster(data, stands, LFFT=LFFT, Overlap=Overlap, IncludeAuto=True, verbose=False,  SampleRate=SampleRate, CentralFreq=CentralFreq)
		blList = uvUtils.getBaselines(stands, IncludeAuto=True, Indicies=False)

		toUse = numpy.where( (freqYY>10.0e6) & (freqYY<88.0e6) )
		toUse = toUse[0]

		if s  == 0:
			fits = fitsidi.IDI(filename, refTime=refTime)
			fits.setStokes(['xx'])
			fits.setFrequency(freqYY[toUse])
			fits.setGeometry(site, stands)

		obsTime = astro.unix_to_taimjd(setTime)
		fits.addDataSet(obsTime, 512*nFrames/SampleRate, blList, outYY[:,toUse])
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
	lwa1 = lwa_common.lwa1()

	fh = open(filename, "rb", buffering=tbn.FrameSize*10000)
	test = tbn.readFrame(fh)
	if not test.header.isTBN():
		raise errors.notTBNError()
	fh.seek(0)

	jd = astro.unix_to_utcjd(test.getTime())
	date = str(ephem.Date(jd - astro.DJD_OFFSET))
	nFpO = tbn.getFramesPerObs(fh)
	nFpO = nFpO[0] + nFpO[1]
	sampleRate = tbn.getSampleRate(fh, nFrames=nFpO)
	nInts = os.path.getsize(filename) / tbn.FrameSize / nFpO

	# Get stands
	stands = lwa1.getStands(date)
	print stands

	# Number of frames to read in at once and average
	nFrames = int(config['avgTime']*sampleRate/512)
	nSkip = int(config['offset']*sampleRate/512)
	fh.seek(nSkip*20*tbn.FrameSize)
	nSets = os.path.getsize(filename) / tbn.FrameSize / nFpO / nFrames
	nSets = nSets - nSkip / nFrames

	print "TBN Data:  %s" % test.header.isTBN()
	print "Samples per observations: %i per pol." % (nFpO/2)
	print "Filter code is: %i" % tbn.getSampleRate(fh, nFrames=nFpO, FilterCode=True)
	print "Sampling rate is: %i Hz" % sampleRate
	print "Captures in file: %i (%.1f s)" % (nInts, nInts*512 / sampleRate)
	print "=="
	print "Station: %s" % lwa1.name
	print "Date observed: %s" % date
	print "Julian day: %.5f" % jd
	print "Integration Time: %.3f s" % (512*nFrames/sampleRate)
	print "Number of integrations in file: %i" % nSets

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
		processChunk(fh, lwa1, stands, fitsFilename, intTime=config['avgTime'], LFFT=config['LFFT'], Overlap=1, CentralFreq=config['cFreq'], 
					SampleRate=sampleRate, ChunkSize=chunk)

		s = s + 1
		leftToDo = leftToDo - chunk

	fh.close()

if __name__ == "__main__":
	main(sys.argv[1:])
