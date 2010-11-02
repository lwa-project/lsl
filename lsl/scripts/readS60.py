#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python script to read in a S60 file and average it in time.  The output is a 
npz file of the time-averaged spectra and a PNG of the bandpass/waterfall diagram."""

import os
import sys
import numpy
import getopt

from lsl.common.paths import data as dataPath
from lsl.reader import s60

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """readS60.py - Read in raw S60 files and create a collection of 
time- averaged spectra that are saved to NPZ files.

Usage: readS60.py [OPTIONS] file

Options:
-h --help                  Display this help information
-e --enable-model          Use the CFTOOL bandpass model if 
                            it is present in the current 
                            directory
-q --quiet                 Run readS60 in silent mode
-l --fft-length            Set FFT length (default = 4096)
-t --avg-time              Window to average spectra in time
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['avgTime'] = 10
	config['LFFT'] = 4096
	config['verbose'] = True
	config['enableModel'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "heql:t:", ["help", "enable-model", "quiet", "fft-length=", "avg-time="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-e', '--enable-model'):
			config['enableModel'] = False
		elif opt in ('-q', '--quiet'):
			config['verbose'] = False
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-t', '--avg-time'):
			config['avgTime'] = int(value)
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

	# Primary integration time in seconds.  This sets how large of chunks to read in 
	# from the file at one time.  This is limited by the ammount of memory avaliable 
	# on the host machine and the upper limit on fornax is around 20-30 seconds.  The
	# integration time should be a multiple of two to ensure that the number of 
	# samples read in divisable by 4096.
	integrationTime = config['avgTime']

	for filename in config['args']:
		# The total number of samples in a file is its size divided by 2 bytes per sample
		nSamples = os.path.getsize(filename) / 2
		# The number of seconds in a file is its size divided by the number of bytes per
		# packet, divided by the number of packets per second
		nSeconds = os.path.getsize(filename) / 2 / s60.SampleRate
		print "Working on file: '%s'" % filename
		print "Samples in file: %i" % nSamples
		print "Seconds in file: %i" % nSeconds

		# Open file for reading in binary mode
		fh = open(filename, "rb", buffering=7340000)

		# Frequency array for LFFT channels and a 3.75 MHz sample rate
		w = numpy.arange(-1.875, 1.875, 3.50/LFFT)

		# Work on the data in ``integrationTime" second segments
		timeBlocks = numpy.zeros((nSeconds/integrationTime, LFFT))
		for i in range(timeBlocks.shape[0]):
			if config['verbose']:
				print "Working on time: %i to %i... %6.2f%%" % (i*integrationTime, (i+1)*integrationTime, (100.0*i/timeBlocks.shape[0]))
			data = s60.readChunk(fh, Chunk=integrationTime*5109*734)
			# Work on computing an integrated spectrum.  This is done by reading in a chunk 
			# of length LFFT, computing is power spectrum, and adding it to all of the pre-
			# vious samples.  The min, mean, and max powers are printed out at every 1,000-
			# th step
			avgCount = 0
			integratedPS = numpy.zeros(LFFT)
			for j in range(data.shape[0]/LFFT):
				avgCount = avgCount + 1

				oneChunk = data[j*LFFT:(j+1)*LFFT]
				oneChunk = oneChunk - oneChunk.mean()
				chunkFFT = numpy.fft.fftshift(numpy.fft.fft(oneChunk))
				chunkPS = numpy.abs(chunkFFT)**2.0
				
				# This middle point of the power spectrum is always zero and it
				# always causes problems with the fitting later on.  So, replace
				# it with the average of its neighbors.
				chunkPS[LFFT/2] = 0.5*chunkPS[LFFT/2-1] + 0.5*chunkPS[LFFT/2+1]

				integratedPS = integratedPS + chunkPS
				
				if j % 1000 == 0 and config['verbose']:
					print " %6.2f%%" % (100.0*j*LFFT/data.shape[0])

			# Print out the last integrated power spectrum to the master array
			timeBlocks[i,:] = 10.0*numpy.log10(integratedPS/avgCount)
		fh.close()

		# Save the data so we don't have to go through all that again
		basename = filename[filename.rfind('/')+1:]
		output = basename.replace('.txt', '.npz')
		output = output.replace('.dat', '.npz')
		numpy.savez(output, w=w, timeBlocks=timeBlocks, integrationTime=integrationTime)

		# Bandpass fitting can either be done with Joe's model from the wiki (if in the
		# current directory) or with a 17-th order polynomial.  In general, I found that
		# the higher order works better on some of the bandpasses.
		if config['enableModel']:
			bandpass = s60.getBandpassModel()
			bandpass = bandpass[300:3700]
			bandpass = bandpass - bandpass.mean() + timeBlocks[:,300:3700].mean()
		else:
			bandpassCoeff = numpy.polyfit(w[300:3700], numpy.squeeze(timeBlocks[:,300:3700].mean(axis=0)), 
									17)
			print "Bandpass Polynomial Coefficients: ", bandpassCoeff
			bandpass = numpy.polyval(bandpassCoeff, w[300:3700])

		# Apply the bandpass correction
		for i in range(timeBlocks.shape[0]):
			timeBlocks[i,300:3700] = timeBlocks[i,300:3700] - bandpass

		# Display and write the results to a PNG
		fig = plt.figure()
		ax1 = fig.add_subplot(211)
		ax1.plot(w[300:3700], numpy.squeeze(timeBlocks[:,300:3700]).mean(axis=0)+bandpass, color='green', 
				linestyle=' ', marker='x', markersize=1.0)
		ax1.plot(w[300:3700], bandpass, color='blue', linestyle='-')
		ax1.set_xlabel('Frequncy Bin')
		ax1.set_ylabel('Bandpass')
		ax1.set_xlim([w[300], w[3700]])
		ax1.axis('auto')
		ax2 = fig.add_subplot(212)
		ax2.imshow(timeBlocks[:,300:3700], origin='lower', vmin=0, vmax=0.3, 
					extent=(w[300], w[3700], 0, (timeBlocks.shape[0]+1)*integrationTime))
		ax2.set_xlabel('Frequncy Bin')
		ax2.set_ylabel('Time [seconds]')
		ax2.axis('auto')

		basename = filename[filename.rfind('/')+1:]
		output = basename.replace('.txt', '.png')
		output = output.replace('.dat', '.png')
		fig.savefig(output)
	

if __name__ == "__main__":
	main(sys.argv[1:])
