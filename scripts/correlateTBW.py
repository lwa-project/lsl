#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example script that reads in TBW data and runs a cross-correlation on it.  
The results are saved in the FITS IDI format.
"""

import os
import sys
import time
import ephem
import numpy
import getopt
import random
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.reader.ldp import LWA1DataFile
from lsl.common import stations, metabundle
from lsl.statistics import robust
from lsl.correlator import uvUtils
from lsl.correlator import fx as fxc
from lsl.writer import fitsidi
from lsl.common.progress import ProgressBar


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
-m, --metadata         Name of SSMIF or metadata tarball file to use for 
                       mappings
-l, --fft-length       Set FFT length (default = 2048)
-q, --quiet            Run correlateTBW in silent mode
-x, --xx               Compute only the XX polarization product (default)
-y, --yy               Compute only the YY polarization product
-2, --two-products     Compute both the XX and YY polarization products
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['metadata'] = ''
	config['LFFT'] = 2048
	config['samples'] = 1
	config['offset'] = 0
	config['verbose'] = True
	config['products'] = ['xx',]
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hm:ql:xy2", ["help", "metadata=", "quiet", "fft-length=", "xx", "yy", "two-products"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
		
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-m', '--metadata'):
			config['metadata'] = value
		elif opt in ('-q', '--quiet'):
			config['verbose'] = False
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-x', '--xx'):
			config['products'] = ['xx',]
		elif opt in ('-y', '--yy'):
			config['products'] = ['yy',]
		elif opt in ('-2', '--two-products'):
			config['products'] = ['xx', 'yy']
		else:
			assert False
			
	# Add in arguments
	config['args'] = arg
	
	# Return configuration
	return config


def processChunk(idf, site, good, filename, LFFT=64, Overlap=1, pols=['xx','yy']):
	"""
	Given an lsl.reader.ldp.TBWFile instances and various parameters for the 
	cross-correlation, write cross-correlate the data and save it to a file.
	"""
	
	# Get antennas
	antennas = site.getAntennas()
	
	# Get the metadata
	sampleRate = idf.getInfo('sampleRate')
	
	# Create the list of good digitizers and a digitizer to Antenna instance mapping.  
	# These are:
	#  toKeep  -> mapping of digitizer number to array location
	#  mapper -> mapping of Antenna instance to array location
	toKeep = [antennas[i].digitizer-1 for i in good]
	mapper = [antennas[i] for i in good]
	
	wallTime = time.time()
	readT, t, data = idf.read()
	setTime = t
	refTime = t
	
	setDT = datetime.utcfromtimestamp(setTime)
	setDT.replace(tzinfo=UTC())
	print "Working on set #1 (%.3f seconds after set #1 = %s)" % ((setTime-refTime), setDT.strftime("%Y/%m/%d %H:%M:%S.%f"))
	
	# In order for the TBW stuff to actaully run, we need to run in with sub-
	# integrations.  8 sub-integrations (61.2 ms / 8 = 7.7 ms per section) 
	# seems to work ok with a "reasonable" number of channels.
	nSec = 8
	secSize = data.shape[1]/nSec
	
	# Loop over polarizations (there should be only 1)
	for pol in pols:
		print "-> %s" % pol
		try:
			tempVis *= 0
		except:
			pass
			
		# Set up the progress bar so we can keep up with how the sub-integrations 
		# are progressing
		pb = ProgressBar(max=nSec)
		sys.stdout.write(pb.show()+'\r')
		sys.stdout.flush()
		
		# Loop over sub-integrations (set by nSec)
		for k in xrange(nSec):
			blList, freq, vis = fxc.FXMaster(data[toKeep,k*secSize:(k+1)*secSize], mapper, LFFT=LFFT, Overlap=Overlap, IncludeAuto=True, verbose=False, SampleRate=sampleRate, CentralFreq=0.0, Pol=pol, ReturnBaselines=True, GainCorrect=True)
			
			toUse = numpy.where( (freq>=5.0e6) & (freq<=93.0e6) )
			toUse = toUse[0]
			
			try:
				tempVis += vis
			except:
				tempVis = vis
				
			pb.inc(amount=1)
			sys.stdout.write(pb.show()+'\r')
			sys.stdout.flush()
			
		# Average the sub-integrations together
		vis = tempVis / float(nSec)
		
		# Set up the FITS IDI file is we need to
		if pol == pols[0]:
			pol1, pol2 = fxc.pol2pol(pol)
			
			fits = fitsidi.IDI(filename, refTime=refTime)
			fits.setStokes(pols)
			fits.setFrequency(freq[toUse])
			fits.setGeometry(site, [a for a in mapper if a.pol == pol1])
			
		# Add the visibilities
		obsTime = astro.unix_to_taimjd(setTime)
		fits.addDataSet(obsTime, readT, blList, vis[:,toUse], pol=pol)
		sys.stdout.write(pb.show()+'\r')
		sys.stdout.write('\n')
		sys.stdout.flush()
	print "->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime))
	
	fits.write()
	fits.close()
	del(fits)
	del(data)
	del(vis)
	return True


def main(args):
	# Parse command line options
	config = parseConfig(args)
	filename = config['args'][0]

	# Length of the FFT
	LFFT = config['LFFT']

	# Setup the LWA station information
	if config['metadata'] != '':
		try:
			station = stations.parseSSMIF(config['metadata'])
		except ValueError:
			station = metabundle.getStation(config['metadata'], ApplySDM=True)
	else:
		station = stations.lwa1
	antennas = station.getAntennas()
	
	idf = LWA1DataFile(filename)
	
	jd = astro.unix_to_utcjd(idf.getInfo('tStart'))
	date = str(ephem.Date(jd - astro.DJD_OFFSET))
	sampleRate = idf.getInfo('sampleRate')
	nInts = idf.getInfo('nFrames') / (30000 * len(antennas) / 2)
	
	# Get valid stands for both polarizations
	goodX = []
	goodY = []
	for i in xrange(len(antennas)):
		ant = antennas[i]
		if ant.getStatus() != 33:
			pass
		else:
			if ant.pol == 0:
				goodX.append(ant)
			else:
				goodY.append(ant)
				
	# Select which polarization to use
	good = []
	for antX in goodX:
		for antY in goodY:
			if antX.stand.id == antY.stand.id:
				good.append(antX.digitizer-1)
				good.append(antY.digitizer-1)
				break
				
	# Report on the valid stands found.  This is a little verbose,
	# but nice to see.
	print "Found %i good stands to use" % (len(good)/2,)
	for i in good:
		print "%3i, %i" % (antennas[i].stand.id, antennas[i].pol)
		
	# Number of frames to read in at once and average
	nFrames = 30000
	nSets = idf.getInfo('nFrames') / (30000*len(antennas)/2)
	
	print "Data type:  %s" % type(idf)
	print "Captures in file: %i (%.3f s)" % (nInts, nInts*30000*400/sampleRate)
	print "=="
	print "Station: %s" % station.name
	print "Date observed: %s" % date
	print "Julian day: %.5f" % jd
	print "Integration Time: %.3f s" % (400*nFrames/sampleRate)
	print "Number of integrations in file: %i" % nSets
	print "=="
	
	if config['samples'] > nSets:
		config['samples'] = nSets
		
	leftToDo = config['samples']
	basename = os.path.split(filename)[1]
	basename, ext = os.path.splitext(basename)
	
	fitsFilename = "%s.FITS_1" % basename
	processChunk(idf, station, good, fitsFilename, LFFT=config['LFFT'], Overlap=1, pols=config['products'])
	
	idf.close()


if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])
