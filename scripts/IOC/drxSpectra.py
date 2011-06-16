#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a DRX file, plot the time averaged spectra for each beam output."""

import os
import re
import sys
import math
import numpy
import ephem
import getopt

import lsl.common.metabundle as mcsMB
import lsl.reader.drx as drx
import lsl.reader.errors as errors
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


filenameRE = re.compile(r'.*(?P<mjd>\d{6})_(?P<opcode>\d{9})_(?P<dr>\d)_(?P<mode>[A-Z]{3}).dat')


def usage(exitCode=None):
	print """drxSpectra.py - Read in DRX files and create a collection of 
time-averaged spectra.

Usage: drxSpectra.py [OPTIONS] metaData data

Options:
-h, --help                  Display this help information
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to average for spectra 
                            (default = 10)
-q, --quiet                 Run drxSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
-d, --disable-chunks        Display plotting chunks in addition to the global 
                            average
-o, --output                Output file name for spectra image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['average'] = 10.0
	config['LFFT'] = 4096
	config['maxFrames'] = 19144
	config['window'] = fxc.noWindow
	config['output'] = None
	config['displayChunks'] = True
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:o:s:a:d", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "output=", "skip=", "average=", "disable-chunks"])
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
		elif opt in ('-t', '--bartlett'):
			config['window'] = numpy.bartlett
		elif opt in ('-b', '--blackman'):
			config['window'] = numpy.blackman
		elif opt in ('-n', '--hanning'):
			config['window'] = numpy.hanning
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-a', '--average'):
			config['average'] = float(value)
		elif opt in ('-d', '--disable-chunks'):
			config['displayChunks'] = False
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def bestFreqUnits(freq):
	"""Given a numpy array of frequencies in Hz, return a new array with the
	frequencies in the best units possible (kHz, MHz, etc.)."""

	# Figure out how large the data are
	scale = int(math.log10(freq.max()))
	if scale >= 9:
		divis = 1e9
		units = 'GHz'
	elif scale >= 6:
		divis = 1e6
		units = 'MHz'
	elif scale >= 3:
		divis = 1e3
		units = 'kHz'
	else:
		divis = 1
		units = 'Hz'

	# Convert the frequency
	newFreq = freq / divis

	# Return units and freq
	return (newFreq, units)


def main(args):
	# Parse command line options
	config = parseOptions(args)

	# Length of the FFT
	LFFT = config['LFFT']

	metaDataFile = config['args'][0]
	dataFile = config['args'][1]

	# First, deal with the metadata
	## Read in the meta-data bundle and split up the data filename into MJD, op 
	## code, etc.
	project = mcsMB.getSessionDefinition(metaDataFile)
	mtch = filenameRE.match('filename')
	
	## Loop over the observation defined in the meta-data bundle and set cObs to 
	## the observation with the right MJD and op code. 
	cObs = None
	for o in project.sessions[0].observations:
		if int(mtch.group('mjd')) == o.mjd and int(mtch.group('opcode')) == o.opcode:
			cObs = o
			break
	if cObs is None:
		print "WARNING:  observation entry for file '%s' not found in '%s'" % (dataFile, metaDataFile)

	fh = open(dataFile, "rb")
	nFramesFile = os.path.getsize(dataFile) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	
	fh.seek(0)
	srate = junkFrame.getSampleRate()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(config['offset'] * srate / 4096 * beampols)
	offset = int(1.0 * offset / beampols) * beampols
	config['offset'] = 1.0 * offset / beampols * 4096 / srate
	fh.seek(offset*drx.FrameSize)

	# Make sure that the file chunk size contains is an integer multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of beampols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(1.0*config['maxFrames']/beampols*4096/float(LFFT))*LFFT/4096*beampols

	# Number of frames to integrate over
	nFrames = int(config['average'] * srate / 4096 * beampols)
	nFrames = int(1.0 * nFrames / beampols*4096/float(LFFT))*LFFT/4096*beampols
	config['average'] = 1.0 * nFrames / beampols * 4096 / srate

	# Number of remaining chunks
	nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))
	
	# Date
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

	# File summary
	print "Filename: %s" % dataFile
	print "Date of First Frame: %s" % str(beginDate)
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	if cObs is not None:
		print "Tunings:  1 @ %.4f MHz, 2 @ %.4f MHz" % (cObs.frequency1/1e6, cObs.frequency2/1e6)
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames / beampols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if offset > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requested integration time+offset is greater than file length")

	# Master loop over all of the file chunks
	masterCount = {}
	standMapper = []
	masterWeight = numpy.zeros((nChunks, beampols, LFFT-1))
	masterSpectra = numpy.zeros((nChunks, beampols, LFFT-1))
	for i in range(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - i*maxFrames
		if framesRemaining > maxFrames:
			framesWork = maxFrames
		else:
			framesWork = framesRemaining
		print "Working on chunk %i, %i frames remaining" % (i, framesRemaining)
		
		count = {}
		data = numpy.zeros((beampols,framesWork*4096/beampols), dtype=numpy.csingle)
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break

		# Inner loop that actually reads the frames into the data array
		print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)

		for j in range(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FrameSize-1)
				continue

			beam,tune,pol = cFrame.parseID()
			aStand = 4*(beam-1) + 2*(tune-1) + pol
			if aStand not in standMapper:
				standMapper.append(aStand)
				oStand = 1*aStand
				aStand = standMapper.index(aStand)
				print "Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand)
			else:
				aStand = standMapper.index(aStand)

			if aStand not in count.keys():
				count[aStand] = 0
				masterCount[aStand] = 0
			#if cFrame.header.frameCount % 10000 == 0 and config['verbose']:
			#	print "%2i,%1i,%1i -> %2i  %5i  %i" % (beam, tune, pol, aStand, j/4, cFrame.data.timeTag)

			data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
			# Update the counters so that we can average properly later on
			count[aStand] = count[aStand] + 1
			masterCount[aStand] = masterCount[aStand] + 1

		## Calculate the data mean for each signal
		#for stand in range(data.shape[0]):
		#	print "Stand %i:  mean is %.4f + %.4f j" % (stand, data[stand,:].mean().real, data[stand,:].mean().imag)

		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate)
		for stand in count.keys():
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = count[stand]

		# We don't really need the data array anymore, so delete it
		del(data)

	# Now that we have read through all of the chunks, perform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
	
	# Apply what we know from the meta-data bundle to the frequencies
	if cObs is not None:
		freq1 = freq + cObs.frequency1
		freq2 = freq + cObs.frequency2
	else:
		freq1 = freq
		freq2 = freq

	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(beampols)))
	figsY = beampols / figsX
	# Put the frequencies in the best units possible
	freq1, units1 = bestFreqUnits(freq1)
	freq2, units2 = bestFreqUnits(freq2)

	sortedMapper = sorted(standMapper)
	for k, aStand in enumerate(sortedMapper):
		i = standMapper.index(aStand)
		
		if standMapper[i]%4/2+1 == 1:
			freq = freq1
			units = units1
		else:
			freq = freq2
			units = units2

		ax = fig.add_subplot(figsX,figsY,k+1)
		currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
		ax.plot(freq, currSpectra, label='%i (avg)' % (i+1))

		# If there is more than one chunk, plot the difference between the global 
		# average and each chunk
		if nChunks > 1 and config['displayChunks']:
			for j in range(nChunks):
				# Some files are padded by zeros at the end and, thus, carry no 
				# weight in the average spectra.  Skip over those.
				if masterWeight[j,i,:].sum() == 0:
					continue

				# Calculate the difference between the spectra and plot
				subspectra = numpy.squeeze( numpy.log10(masterSpectra[j,i,:])*10.0 )
				diff = subspectra - currSpectra
				ax.plot(freq, diff, label='%i' % j)

		ax.set_title('Beam %i, Tune. %i, Pol. %i' % (standMapper[i]/4+1, standMapper[i]%4/2+1, standMapper[i]%2))
		ax.set_xlabel('Frequency Offset [%s]' % units)
		ax.set_ylabel('P.S.D. [dB/RBW]')
		ax.set_xlim([freq.min(), freq.max()])
		ax.legend(loc=0)
		
		print "For beam %i, tune. %i, pol. %i maximum in PSD at %.3f %s" % (standMapper[i]/4+1, standMapper[i]%4/2+1, standMapper[i]%2, freq[numpy.where( spec[i,:] == spec[i,:].max() )][0], units)

	print "RBW 1: %.4f %s" % ((freq1[1]-freq1[0]), units1)
	print "RBW 2: %.4f %s" % ((freq2[1]-freq2[0]), units2)
	plt.subplots_adjust(hspace=0.35, wspace=0.30)
	plt.show()

	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])

if __name__ == "__main__":
	main(sys.argv[1:])
