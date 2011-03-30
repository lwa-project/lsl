#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a DRX file, plot the time averaged spectra for each beam output."""

import os
import sys
import math
import time
import numpy
import getopt

import lsl.reader.drx as drx
import lsl.reader.errors as errors

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxTimeseries.py - Read in DRX files and create a collection of 
timeseries (I/Q) plots.

Usage: drxTimeseries.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-p, --plot-range            Number of seconds of data to show in the I/Q plots
                            (default = 0.0001)
-q, --quiet                 Run drxSpectra in silent mode
-o, --output                Output file name for time series image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['average'] = 0.0001
	config['maxFrames'] = 19144
	config['output'] = None
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqo:s:p:", ["help", "quiet", "output=", "skip=", "plot-range="])
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
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-p', '--plot-range'):
			config['average'] = float(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	# Parse command line options
	config = parseOptions(args)
	
	fh = open(config['args'][0], "rb")
	nFramesFile = os.path.getsize(config['args'][0]) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	fh.seek(0)
	srate = junkFrame.getSampleRate()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(round(config['offset'] * srate / 4096 * beampols))
	offset = int(1.0 * offset / beampols) * beampols
	config['offset'] = 1.0 * offset / beampols * 4096 / srate
	fh.seek(offset*drx.FrameSize)

	# Make sure that the file chunk size contains is an intger multiple
	# of the beampols.
	maxFrames = int(config['maxFrames']/beampols)*beampols

	# Number of frames to integrate over
	toClip = False
	oldAverage = config['average']
	if config['average'] < 4096/srate:		
		toClip = True
		config['average'] = 4096/srate
	nFrames = int(config['average'] * srate / 4096 * beampols)
	nFrames = int(1.0 * nFrames / beampols) * beampols
	config['average'] = 1.0 * nFrames / beampols * 4096 / srate

	# Number of remaining chunks
	nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Plot time: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames / beampols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if offset > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requested integration time+offset is greater than file length")

	# Master loop over all of the file chuncks
	standMapper = []
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
		
		# Inner loop that actually reads the frames into the data array
		print "Working on %.1f ms of data" % ((framesWork*4096/beampols/srate)*1000.0)
		t0 = time.time()
		
		for j in xrange(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = drx.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/drx.FrameSize-1)
				continue
			except errors.numpyError:
				break
			
			beam,tune,pol = cFrame.parseID()
			aStand = 4*(beam-1) + 2*(tune-1) + pol
			#print aStand, beam, tune, pol
			if aStand not in standMapper:
				standMapper.append(aStand)
				oStand = 1*aStand
				aStand = standMapper.index(aStand)
				print "Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand)
			else:
				aStand = standMapper.index(aStand)

			if aStand not in count.keys():
				count[aStand] = 0
			#if cFrame.header.frameCount % 10000 == 0 and config['verbose']:
			#	print "%2i,%1i,%1i -> %2i  %5i  %i" % (beam, tune, pol, aStand, cFrame.header.frameCount, cFrame.data.timeTag)

			#print data.shape, count[aStand]*4096, (count[aStand]+1)*4096, cFrame.data.iq.shape
			data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
			# Update the counters so that we can average properly later on
			count[aStand] += 1
			
		# The plots:  This is setup for the current configuration of 20 beampols
		fig = plt.figure()
		figsX = int(round(math.sqrt(beampols)))
		figsY = beampols / figsX

		samples = int(oldAverage * srate)
		if toClip:
			print "Plotting only the first %i samples (%.3f ms) of data" % (samples, oldAverage*1000.0)

		sortedMapper = sorted(standMapper)
		for k, aStand in enumerate(sortedMapper):
			i = standMapper.index(aStand)

			ax = fig.add_subplot(figsX,figsY,k+1)
			if toClip:
				ax.plot(numpy.arange(0,samples)/srate, data[i,0:samples].real, label='%i (real)' % (i+1))
				ax.plot(numpy.arange(0,samples)/srate, data[i,0:samples].imag, label='%i (imag)' % (i+1))
			else:
				ax.plot(numpy.arange(0,data.shape[1])/srate, data[i,:].real, label='%i (real)' % (i+1))
				ax.plot(numpy.arange(0,data.shape[1])/srate, data[i,:].imag, label='%i (imag)' % (i+1))
			ax.set_ylim([-8, 7])
			ax.legend(loc=0)
			
			ax.set_title('Beam %i, Tune. %i, Pol. %i' % (standMapper[i]/4+1, standMapper[i]%4/2+1, standMapper[i]%2))
			ax.set_xlabel('Time [seconds]')
			ax.set_ylabel('Output Level')
		plt.show()

		# Save image if requested
		if config['output'] is not None:
			fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
