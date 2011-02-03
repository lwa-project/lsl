#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a TBN file, plot the time averaged spectra for each digitizer input."""

import os
import sys
import math
import numpy
import ephem
import getopt

import lsl.reader.tbn as tbn
import lsl.reader.errors as errors
import lsl.correlator.fx as fxc
from lsl.correlator.uvUtils import CableCache
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """tbnSpectra.py - Read in TBN files and create a collection of 
time- averaged spectra.

Usage: tbnSpectra.py [OPTIONS] file

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
	config['maxFrames'] = 400000
	config['window'] = fxc.noWindow
	config['applyGain'] = False
	config['output'] = None
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:go:", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "gain-correct", "output="])
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
		elif opt in ('-g', '--gain-correct'):
			# This will need to be changed if we ever get frequency information
			# into the frame headers/add an option to input the frequency
			config['applyGain'] = False
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
	nFrames = os.path.getsize(config['args'][0]) / tbn.FrameSize
	srate = tbn.getSampleRate(fh)
	antpols = tbn.getFramesPerObs(fh)
	antpols = antpols[0]+antpols[1]

	# Make sure that the file chunk size contains is an intger multiple
	# of the FFT length so that no data gets dropped.  This needs to
	# take into account the number of antpols in the data, the FFT length,
	# and the number of samples per frame.
	maxFrames = int(config['maxFrames']/antpols*512/float(LFFT))*LFFT/512*antpols

	nChunks = int(math.ceil(1.0*nFrames/maxFrames))

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbn.readFrame(fh)
	fh.seek(0)
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	stands = stations.lwa1().getStands(beginDate)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Rate: %i Hz" % srate
	print "Frames: %i" % nFrames
	print "Chunks: %i" % nChunks

	# Master loop over all of the file chuncks
	masterCount = {}
	standMapper = []
	masterWeight = numpy.zeros((nChunks, antpols, LFFT-1))
	masterSpectra = numpy.zeros((nChunks, antpols, LFFT-1))
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
		data = numpy.zeros((antpols,framesWork*512/antpols), dtype=numpy.csingle)
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < LFFT:
			break
		# Inner loop that actually reads the frames into the data array
		for j in range(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbn.readFrame(fh, SampleRate=srate, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
				continue
			except errors.numpyError:
				break

			stand,pol = cFrame.header.parseID()
			# In the current configuration, stands start at 1 and go up to 10.  So, we
			# can use this little trick to populate the data array
			aStand = 2*(stand-1)+pol
			if aStand not in standMapper:
				standMapper.append(aStand)
				oStand = 1*aStand
				aStand = standMapper.index(aStand)
				print "Mapping stand %3i, pol. %1i (%3i) to array index %3i" % (stand, pol, oStand, aStand)

			if aStand not in count.keys():
				count[aStand] = 0
				masterCount[aStand] = 0
			if cFrame.header.frameCount % 10000 == 0 and config['verbose']:
				print "%2i,%1i -> %2i  %5i  %i" % (stand, pol, aStand, cFrame.header.frameCount, cFrame.data.timeTag)

			# Additional check on the data array bounds so that we don't overflow it.  
			# This check may be redundant...
			if (count[aStand]+1)*512 >= data.shape[1]:
				continue
			data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
			# Update the counters so that we can average properly later on
			count[aStand] = count[aStand] + 1
			masterCount[aStand] = masterCount[aStand] + 1

		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec = fxc.calcSpectra(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate)
		for stand in count.keys():
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = count[stand]

		# We don't really need the data array anymore, so delete it
		del(data)

	# Apply the cable loss corrections, if requested
	if config['applyGain']:
		cbl = CableCache(freq=freq)
		for s in range(masterSpectra.shape[1]):
			for c in range(masterSpectra.shape[0]):
				masterSpectra[c,s,:] *= cbl.cableGain(stands[s])

	# Now that we have read through all of the chunks, peform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )

	# The plots:  This is setup for the current configuration of 20 antpols
	fig = plt.figure()
	figsX = int(round(math.sqrt(antpols)))
	figsY = antpols / figsX
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

		ax.set_title('Stand %i, Pol. %i' % (standMapper[i]/2+1, standMapper[i]%2))
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
