#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a TBW file, plot the time averaged spectra for each digitizer input."""

import os
import sys
import math
import numpy
import ephem
import getopt

import lsl.reader.tbw as tbw
import lsl.reader.errors as errors
import lsl.correlator.fx as fxc
from lsl.correlator.uvUtils import CableCache
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """tbwSpectra.py - Read in TBW files and create a collection of 
time-averaged spectra.

Usage: tbwSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-q, --quiet                 Run tbwSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
-g, --gain-correct          Correct signals for the cable losses
-s, --stack                 Stack spectra in groups of 6 (if '-g' is enabled only)
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
	config['LFFT'] = 4096
	config['maxFrames'] = 300000
	config['window'] = fxc.noWindow
	config['applyGain'] = False
	config['stack'] = False
	config['output'] = None
	config['displayChunks'] = True
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqtbnl:gso:d", ["help", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "gain-correct", "stack", "output=", "disable-chunks"])
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
			config['applyGain'] = True
		elif opt in ('-s', '--stack'):
			config['stack'] = True
		elif opt in ('-o', '--output'):
			config['output'] = value
		elif opt in ('-d', '--disable-chunks'):
			config['displayChunks'] = False
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config

def main(args):
	# Parse command line options
	config = parseOptions(args)

	# Length of the FFT
	LFFT = config['LFFT']

	# Make sure that the file chunk size contains is an intger multiple
	# of the FFT length so that no data gets dropped
	maxFrames = int(config['maxFrames']/float(LFFT))*LFFT
	# It seems like that would be a good idea, however...  TBW data comes one
	# capture at a time so doing something like this actually truncates data 
	# from the last set of stands for the first integration.  So, we really 
	# should stick with
	maxFrames = config['maxFrames']

	fh = open(config['args'][0], "rb")
	nFrames = os.path.getsize(config['args'][0]) / tbw.FrameSize
	dataBits = tbw.getDataBits(fh)
	# The number of ant/pols in the file is hard coded because I cannot figure out 
	# a way to get this number in a systematic fashion
	antpols = 20 #2*getFramesPerObs(fh
	nChunks = int(math.ceil(1.0*nFrames/maxFrames))
	if dataBits == 12:
		nSamples = 400
	else:
		nSamples = 1200

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbw.readFrame(fh)
	fh.seek(0)
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	stands = stations.lwa1().getStands(beginDate)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Length: %i-bit" % dataBits
	print "Frames: %i" % nFrames
	print "Chunks: %i" % nChunks

	# Master loop over all of the file chuncks
	masterCount = {}
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
		data = numpy.zeros((antpols,framesWork*nSamples/antpols), dtype=numpy.int16)
		# If there are fewer frames than we need to fill an FFT, skip this chunk
		if data.shape[1] < 2*LFFT:
			break
		# Inner loop that actually reads the frames into the data array
		for j in range(framesWork):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbw.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FrameSize-1)
				continue
			except errors.numpyError:
				break

			stand = cFrame.header.parseID()
			# In the current configuration, stands start at 1 and go up to 10.  So, we
			# can use this little trick to populate the data array
			aStand = 2*(stand-1)
			if aStand not in count.keys():
				count[aStand] = 0
				masterCount[aStand] = 0
			if cFrame.header.frameCount % 10000 == 0 and config['verbose']:
				print "%2i -> %2i  %6.3f  %5i  %i" % (stand, aStand, cFrame.getTime(), cFrame.header.frameCount, cFrame.data.timeTag)

			# Additional check on the data array bounds so that we don't overflow it.  
			# This check may be redundant...
			if (count[aStand]+1)*nSamples >= data.shape[1]:
				continue

			# Actaully load the data.  x pol goes into the even numbers, y pol into the 
			# odd numbers
			data[aStand, count[aStand]*nSamples:(count[aStand]+1)*nSamples] = numpy.squeeze(cFrame.data.xy[0,:])
			data[aStand+1, count[aStand]*nSamples:(count[aStand]+1)*nSamples] = numpy.squeeze(cFrame.data.xy[1,:])
			# Update the counters so that we can average properly later on
			count[aStand] = count[aStand] + 1
			masterCount[aStand] = masterCount[aStand] + 1

		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		# NB:  The weighting is the same for the x and y polarizations because of how 
		# the data are packed in TBW
		freq, tempSpec = fxc.calcSpectra(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'])
		for stand in count.keys():
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterSpectra[i,stand+1,:] = tempSpec[stand+1,:]
			masterWeight[i,stand,:] = count[stand]
			masterWeight[i,stand+1,:] = count[stand]

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
	if config['applyGain'] & config['stack']:
		# Stacked spectra - only if cable loss corrections are to be applied
		colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 
			'purple', 'salmon', 'olive', 'maroon', 'saddlebrown', 'yellowgreen', 
			'teal', 'steelblue', 'seagreen', 'slategray', 'mediumorchid', 'lime', 
			'dodgerblue', 'darkorange']
		fig = plt.figure()
		ax1 = fig.add_subplot(1, 1, 1)
		for i in range(antpols):
			currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
			ax1.plot(freq/1e6, currSpectra, label='%i' % stands[i], color=colors[i])

		ax1.set_xlabel('Frequency [MHz]')
		ax1.set_ylabel('P.S.D. [dB/RBW]')
		ax1.set_xlim([20,88])
		#ax1.set_ylim([10,90])
		leg = ax1.legend(loc=0, ncol=3)
		for l in leg.get_lines():
			l.set_linewidth(1.7)  # the legend line width
	else:
		# Normal plotting
		fig = plt.figure()
		figsY = int(round(math.sqrt(antpols)))
		figsX = antpols / figsY
		for i in range(antpols):
			ax = fig.add_subplot(figsX,figsY,i+1)
			currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
			ax.plot(freq/1e6, currSpectra, label='%i' % (i+1))

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
					ax.plot(freq/1e6, diff)

			ax.set_title('Input: %i' % (i+1))
			ax.set_xlabel('Frequency [MHz]')
			ax.set_ylabel('P.S.D. [dB/RBW]')
			ax.set_xlim([10,90])

	print "RBW: %.1f Hz" % (freq[1]-freq[0])
	plt.show()

	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])

