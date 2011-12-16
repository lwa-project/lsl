#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a TBN file, plot the time averaged spectra for each digitizer input."""

import os
import sys
import math
import numpy
import ephem
import getopt

from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.correlator import fx as fxc
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt

from collections import deque


def usage(exitCode=None):
	print """tbnSpectra.py - Read in TBN files and create a collection of 
time-averaged spectra.

Usage: tbnSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-m, --metadata              Name of SSMIF file to use for mappings
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to average for spectra 
                            (default = 10)
-q, --quiet                 Run tbnSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
-d, --disable-chunks        Display plotting chunks in addition to the global 
                            average
-k, --keep                  Only display the following comma-seperated list of 
                            stands (default = show all 260 dual pol)
-o, --output                Output file name for spectra image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['SSMIF'] = ''
	config['offset'] = 0.0
	config['average'] = 10.0
	config['LFFT'] = 4096
	config['maxFrames'] = 2*260*1000
	config['window'] = fxc.noWindow
	config['applyGain'] = False
	config['output'] = None
	config['displayChunks'] = True
	config['verbose'] = True
	config['keep'] = None
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hm:qtbnl:go:s:a:dk:", ["help", "metadata=", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "gain-correct", "output=", "skip=", "average=", "disable-chunks", "keep="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-m', '--metadata'):
			config['SSMIF'] = value
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
		elif opt in ('-g', '--gain-correct'):
			# This will need to be changed if we ever get frequency information
			# into the frame headers/add an option to input the frequency
			config['applyGain'] = False
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-a', '--average'):
			config['average'] = float(value)
		elif opt in ('-d', '--disable-chunks'):
			config['displayChunks'] = False
		elif opt in ('-k', '--keep'):
			config['keep'] = [int(i) for i in value.split(',')]
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
	
	# Set the station
	if config['SSMIF'] != '':
		station = stations.parseSSMIF(config['SSMIF'])
	else:
		station = stations.lwa1
	antennas = station.getAntennas()

	# Length of the FFT
	LFFT = config['LFFT']

	fh = open(config['args'][0], "rb")
	nFramesFile = os.path.getsize(config['args'][0]) / tbn.FrameSize
	srate = tbn.getSampleRate(fh)
	#antpols = tbn.getFramesPerObs(fh)
	antpols = len(antennas)

	# Offset in frames for beampols beam/tuning/pol. sets
	offset = int(config['offset'] * srate / 512 * antpols)
	offset = int(1.0 * offset / antpols) * antpols
	config['offset'] = 1.0 * offset / antpols * 512 / srate
	fh.seek(offset*tbn.FrameSize)

	# Make sure that the file chunk size contains is an integer multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of antpols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(config['maxFrames']/antpols*512/float(LFFT))*LFFT/512*antpols

	# Number of frames to integrate over
	nFrames = int(config['average'] * srate / 512 * antpols)
	nFrames = int(1.0 * nFrames / antpols*512/float(LFFT))*LFFT/512*antpols
	config['average'] = 1.0 * nFrames / antpols * 512 / srate

	# Number of remaining chunks
	nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbn.readFrame(fh)
	fh.seek(0)
	centralFreq = junkFrame.getCentralFreq()
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz" % centralFreq
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], offset)
	print "Integration: %.3f s (%i frames; %i frames per stand/pol)" % (config['average'], nFrames, nFrames / antpols)
	print "Chunks: %i" % nChunks

	# Sanity check
	if offset > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - offset):
		raise RuntimeError("Requested integration time+offset is greater than file length")

	# Create the FrameBuffer instance
	buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])

	# Master loop over all of the file chunks
	masterWeight = numpy.zeros((nChunks, antpols, LFFT-1))
	masterSpectra = numpy.zeros((nChunks, antpols, LFFT-1))

	k = 0
	for i in xrange(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFrames - k
		if framesRemaining > maxFrames:
			framesWork = maxFrames
			data = numpy.zeros((antpols, framesWork*512/antpols), dtype=numpy.csingle)
		else:
			framesWork = framesRemaining + antpols*buffer.nSegments
			data = numpy.zeros((antpols, framesWork/antpols*512), dtype=numpy.csingle)
			framesWork = framesRemaining
			print "Padding from %i to %i frames" % (framesRemaining, framesWork)
		print "Working on chunk %i, %i frames remaining" % (i, framesRemaining)
		
		count = [0 for a in xrange(len(antennas))]
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break
		
		j = 0
		fillsWork = framesWork / antpols
		# Inner loop that actually reads the frames into the data array
		while j < fillsWork:
			# Read in the next frame and anticipate any problems that could occur
			cFrames = deque()
			for l in xrange(520):
				try:
					cFrames.append( tbn.readFrame(fh) )
					k = k + 1
				except errors.eofError:
					break
				except errors.syncError:
					#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
					continue
					
			buffer.append(cFrames)
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
				
				data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
				count[aStand] += 1
			
			j += 1
		
		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate)
		for stand in xrange(len(count)):
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = int(count[stand] * 512 / LFFT)
	
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
			
			data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
			count[aStand] += 1
			
	# Calculate the spectra for this block of data and then weight the results by 
	# the total number of frames read.  This is needed to keep the averages correct.
	freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate)
	for stand in xrange(len(count)):
		masterSpectra[i,stand,:] = tempSpec[stand,:]
		masterWeight[i,stand,:] = int(count[stand] * 512 / LFFT)

	# We don't really need the data array anymore, so delete it
	del(data)

	# Apply the cable loss corrections, if requested
	if config['applyGain']:
		for s in range(masterSpectra.shape[1]):
			currGain = antennas[s].cable.gain(freq)
			for c in range(masterSpectra.shape[0]):
				masterSpectra[c,s,:] /= currGain

	# Now that we have read through all of the chunks, perform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )

	# Put the frequencies in the best units possible
	freq += centralFreq
	freq, units = bestFreqUnits(freq)
	
	# Deal with the `keep` options
	if config['keep'] is None:
		antpolsDisp = int(numpy.ceil(antpols/20))
		js = [i for i in xrange(antpols)]
	else:
		antpolsDisp = int(numpy.ceil(len(config['keep'])*2/20))
		if antpolsDisp < 1:
			antpolsDisp = 1
		
		js = []
		for k in config['keep']:
			for i,ant in enumerate(antennas):
				if ant.stand.id == k:
					js.append(i)

	for i in xrange(antpolsDisp):
		# Normal plotting
		fig = plt.figure()
		figsY = 4
		figsX = 5
		fig.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.20, hspace=0.50)
		for k in xrange(i*20, i*20+20):
			try:
				j = js[k]
				currSpectra = numpy.squeeze( numpy.log10(spec[j,:])*10.0 )
			except IndexError:
				break
			ax = fig.add_subplot(figsX, figsY, (k%20)+1)
			ax.plot(freq, currSpectra, label='Stand: %i, Pol: %i (Dig: %i)' % (antennas[j].stand.id, antennas[j].pol, antennas[j].digitizer))

			# If there is more than one chunk, plot the difference between the global 
			# average and each chunk
			if nChunks > 1 and config['displayChunks']:
				for k in xrange(nChunks):
					# Some files are padded by zeros at the end and, thus, carry no 
					# weight in the average spectra.  Skip over those.
					if masterWeight[k,j,:].sum() == 0:
						continue

					# Calculate the difference between the spectra and plot
					subspectra = numpy.squeeze( numpy.log10(masterSpectra[k,j,:])*10.0 )
					diff = subspectra - currSpectra
					ax.plot(freq, diff)

			ax.set_title('Stand: %i (%i); Dig: %i [%i]' % (antennas[j].stand.id, antennas[j].pol, antennas[j].digitizer, antennas[j].getStatus()))
			ax.set_xlabel('Frequency [%s]' % units)
			ax.set_ylabel('P.S.D. [dB/RBW]')
			ax.set_ylim([-10, 30])
			
		# Save spectra image if requested
		if config['output'] is not None:
			base, ext = os.path.splitext(config['output'])
			outFigure = "%s-%02i%s" % (base, i+1, ext)
			fig.savefig(outFigure)
			
		plt.draw()
			
	print "RBW: %.4f %s" % ((freq[1]-freq[0]), units)
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
