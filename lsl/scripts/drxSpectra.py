#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a TBN file, plot the time averaged spectra for each beam output."""

import os
import sys
import math
import numpy
import getopt

import lsl.reader.drx as drx
import lsl.reader.errors as errors
import lsl.correlator.fx as fxc

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxSpectra.py - Read in DRX files and create a collection of 
time- averaged spectra.

Usage: drxSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-q, --quiet                 Run tbwSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
-o, --output                Output file name for spectra image
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['LFFT'] = 4096
	config['maxFrames'] = 50000
	config['window'] = fxc.noWindow
	config['output'] = None
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:o:", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "output="])
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

	fh = open(config['args'][0], "rb")
	nFrames = os.path.getsize(config['args'][0]) / drx.FrameSize
	junkFrame = drx.readFrame(fh)
	print junkFrame
	fh.seek(0)
	srate = junkFrame.getSampleRate()
	beams = drx.getBeamCount(fh)
	tunepols = drx.getFramesPerObs(fh)
	tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
	beampols = tunepol

	# Make sure that the file chunk size contains is an intger multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of beampols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(config['maxFrames']/beampols*4096/float(LFFT))*LFFT/4096*beampols

	nChunks = int(math.ceil(1.0*nFrames/maxFrames))

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Beams: %i" % beams
	print "Tune/Pols: %i %i %i %i" % tunepols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i" % nFrames
	print "Chunks: %i" % nChunks

	# Master loop over all of the file chuncks
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
		for j in range(framesWork):
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
			if aStand not in standMapper:
				standMapper.append(aStand)
				oStand = 1*aStand
				aStand = standMapper.index(aStand)
				print "Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand)

			if aStand not in count.keys():
				count[aStand] = 0
				masterCount[aStand] = 0
			if cFrame.header.frameCount % 10000 == 0 and config['verbose']:
				print "%2i,%1i,%1i -> %2i  %5i  %i" % (beam, tune, pol, aStand, cFrame.header.frameCount, cFrame.data.timeTag)

			# Additional check on the data array bounds so that we don't overflow it.  
			# This check may be redundant...
			if (count[aStand]+1)*4096 >= data.shape[1]:
				continue
			data[aStand, count[aStand]*4096:(count[aStand]+1)*4096] = cFrame.data.iq
			# Update the counters so that we can average properly later on
			count[aStand] = count[aStand] + 1
			masterCount[aStand] = masterCount[aStand] + 1

		## Calculate the data mean for each signal
		#for stand in range(data.shape[0]):
		#	print "Stand %i:  mean is %.4f + %.4f j" % (stand, data[stand,:].mean().real, data[stand,:].mean().imag)

		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec = fxc.calcSpectra(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate)
		for stand in count.keys():
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = count[stand]

		# We don't really need the data array anymore, so delete it
		del(data)

	# Now that we have read through all of the chunks, peform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )

	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(beampols)))
	figsY = beampols / figsX
	# Put the freqencies in the best units possible
	freq, units = bestFreqUnits(freq)

	sortedMapper = sorted(standMapper)
	for k, aStand in enumerate(sortedMapper):
		i = standMapper.index(aStand)

		ax = fig.add_subplot(figsX,figsY,k+1)
		currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
		ax.plot(freq, currSpectra, label='%i' % (i+1))

		# If there is more than one chunk, plot the difference between the global 
		# average and each chunk
		if nChunks > 1:
			for j in range(nChunks):
				# Some files are padded by zeros at the end and, thus, carry no 
				# weight in the average spectra.  Skip over those.
				if masterWeight[j,i,:].sum() == 0:
					continue

				# Calculate the difference between the spectra and plot
				subspectra = numpy.squeeze( numpy.log10(masterSpectra[j,i,:])*10.0 )
				diff = subspectra - currSpectra
				ax.plot(freq, diff)

		ax.set_title('Beam %i, Tune. %i, Pol. %i' % (standMapper[i]/4+1, standMapper[i]%4/2+1, standMapper[i]%2))
		ax.set_xlabel('Frequency Offset [%s]' % units)
		ax.set_ylabel('P.S.D. [dB/RBW]')
		ax.set_xlim([freq.min(), freq.max()])

	print "RBW: %.4f %s" % ((freq[1]-freq[0]), units)
	plt.show()

	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])

if __name__ == "__main__":
	main(sys.argv[1:])
