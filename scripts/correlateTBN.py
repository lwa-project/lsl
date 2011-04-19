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
from lsl.common import stations
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
-h, --help             Display this help information
-c, --central-freq     Central frequency of the observations in MHz
-f, --fft-length       Set FFT length (default = 512)
-t, --avg-time         Window to average visibilities in time (seconds; 
                       default = 6 s)
-s, --samples          Number of average visibilities to generate
                       (default = 10)
-o, --offset           Seconds to skip from the beginning of the file
-q, --quiet            Run correlateTBN in silent mode
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


def processChunk(fh, site, filename, intTime=6.0, LFFT=64, Overlap=1, CentralFreq=49.0e6, SampleRate=dp_common.fS, ChunkSize=300, pol=0):
	"""Given a filehandle pointing to some TBN data and various parameters for
	the cross-correlation, write cross-correlate the data and save it to a file."""

	# Get antennas
	antennas = site.getAntennas()
	antpols = len(antennas)

	# Create the FrameBuffer instance
	buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])

	# Create the list of good digitizers
	goodDigs = []
	dig2ant = {}
	for ant in antennas:
		if ant.pol != pol:
			continue
		goodDigs.append(ant.digitizer)
		dig2ant[ant.digitizer] = ant
		
	if pol == 0:
		polCode = 'yy'
	else:
		polCode = 'xx'

	# Find out how many frames to work with at a time
	nFrames = int(intTime*SampleRate/512)

	k = 0
	mapper = []
	mapper2 = []
	refTime = 0.0
	setTime = 0.0
	wallTime = time.time()
	for s in xrange(ChunkSize):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = ChunkSize*nFrames*antpols - k
		if framesRemaining > nFrames*antpols:
			framesWork =  nFrames*antpols
			data = numpy.zeros((len(goodDigs), framesWork*512/antpols), dtype=numpy.csingle)
		else:
			framesWork = framesRemaining + antpols*buffer.nSegments
			data = numpy.zeros((len(goodDigs), framesWork/antpols*512), dtype=numpy.csingle)
			framesWork = framesRemaining
		
		count = [0 for a in xrange(len(goodDigs))]
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break
		
		j = 0
		fillsWork = framesWork / antpols
		# Inner loop that actually reads the frames into the data array
		while j < fillsWork:
			try:
				cFrame = tbn.readFrame(fh)
				k = k + 1
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
				continue
					
			buffer.append(cFrame)
			cFrames = buffer.get()

			if cFrames is None:
				continue
			
			valid = reduce(lambda x,y: x+int(y.valid), cFrames, 0)
			if valid != antpols:
				print "WARNING: frame count %i at %i missing %.2f%% of frames" % (cFrames[0].header.frameCount, cFrames[0].data.timeTag, float(antpols - valid)/antpols*100)
			
			for cFrame in cFrames:
				stand,pol = cFrame.header.parseID()
				
				# In the current configuration, stands start at 1 and go up to 260.  So, we
				# can use this little trick to populate the data array
				aStand = 2*(stand-1)+pol
				
				# Skip stands of the wrong polarization
				if aStand not in goodDigs:
					continue
				
				# Get the set time and reference time if we are in the right 
				if j == 0:
					setTime = cFrame.getTime()
					if s == 0:
						refTime = setTime
				
				data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
				
				# Update the counters so that we can average properly later on
				count[aStand] = count[aStand] + 1
				masterCount[aStand] = masterCount[aStand] + 1
			
			j += 1

		setDT = datetime.utcfromtimestamp(setTime)
		setDT.replace(tzinfo=UTC())
		print "Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-refTime), setDT.strftime("%Y/%m/%d %H:%M:%S.%f"))
		freqVis, outVis = fxc.FXMaster(data, stands, LFFT=LFFT, Overlap=Overlap, IncludeAuto=True, verbose=False,  SampleRate=SampleRate, CentralFreq=CentralFreq)
		blList = uvUtils.getBaselines(stands, IncludeAuto=True, Indicies=False)

		toUse = numpy.where( (freqVis>10.0e6) & (freqVis<88.0e6) )
		toUse = toUse[0]

		if s  == 0:
			fits = fitsidi.IDI(filename, refTime=refTime)
			fits.setStokes([polCode])
			fits.setFrequency(freqVis[toUse])
			fits.setGeometry(site, stands)

		obsTime = astro.unix_to_taimjd(setTime)
		fits.addDataSet(obsTime, 512*nFrames/SampleRate, blList, outVis[:,toUse])
		print "->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)/(s+1))
		
	# Empty the remaining portion of the buffer and integrate what's left
	for cFrames in buffer.flush():
		# Inner loop that actually reads the frames into the data array
		valid = reduce(lambda x,y: x+int(y.valid), cFrames, 0)
		if valid != antpols:
			print "WARNING: frame count %i at %i missing %.2f%% of frames" % (cFrames[0].header.frameCount, cFrames[0].data.timeTag, float(antpols - valid)/antpols*100)
		
		for cFrame in cFrames:
			stand,pol = cFrame.header.parseID()
			# In the current configuration, stands start at 1 and go up to 10.  So, we
			# can use this little trick to populate the data array
			aStand = 2*(stand-1)+pol
			if j == 0:
				setTime = cFrame.getTime()
				if s == 0:
					refTime = setTime
			
			data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
				
			# Update the counters so that we can average properly later on
			count[aStand] = count[aStand] + 1
			masterCount[aStand] = masterCount[aStand] + 1
			
		j += 1
		
	setDT = datetime.utcfromtimestamp(setTime)
	setDT.replace(tzinfo=UTC())
	print "Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-refTime), setDT.strftime("%Y/%m/%d %H:%M:%S.%f"))
	freqVis, outVis = fxc.FXMaster(data, stands, LFFT=LFFT, Overlap=Overlap, IncludeAuto=True, verbose=False,  SampleRate=SampleRate, CentralFreq=CentralFreq)
	blList = uvUtils.getBaselines(stands, IncludeAuto=True, Indicies=False)

	toUse = numpy.where( (freqVis>10.0e6) & (freqVis<88.0e6) )
	toUse = toUse[0]

	if s  == 0:
		fits = fitsidi.IDI(filename, refTime=refTime)
		fits.setStokes([polCode])
		fits.setFrequency(freqVis[toUse])
		fits.setGeometry(site, stands)

	obsTime = astro.unix_to_taimjd(setTime)
	fits.addDataSet(obsTime, 512*nFrames/SampleRate, blList, outVis[:,toUse])
	print "->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)/(s+1))

	fits.write()
	fits.close()
	del(fits)
	del(data)
	del(freqVis)
	del(outVis)
	return True


def main(args):
	# Parse command line options
	config = parseConfig(args)
	filename = config['args'][0]

	# Length of the FFT
	LFFT = config['LFFT']

	# Setup the LWA station information
	station = stations.lwa1
	antennas = station.getAntennas()

	fh = open(filename, "rb", buffering=tbn.FrameSize*10000)
	test = tbn.readFrame(fh)
	if not test.header.isTBN():
		raise errors.notTBNError()
	fh.seek(0)

	jd = astro.unix_to_utcjd(test.getTime())
	date = str(ephem.Date(jd - astro.DJD_OFFSET))
	#nFpO = tbn.getFramesPerObs(fh)
	#nFpO = nFpO[0] + nFpO[1]
	nFpO = len(antennas)
	sampleRate = tbn.getSampleRate(fh)
	nInts = os.path.getsize(filename) / tbn.FrameSize / nFpO

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
		processChunk(fh, station, fitsFilename, intTime=config['avgTime'], LFFT=config['LFFT'], Overlap=1, CentralFreq=config['cFreq'], 
					SampleRate=sampleRate, ChunkSize=chunk)

		s = s + 1
		leftToDo = leftToDo - chunk

	fh.close()

if __name__ == "__main__":
	main(sys.argv[1:])
