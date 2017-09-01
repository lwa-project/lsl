#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a DRX file, plot the time averaged spectra for each beam output."""

import sys
import math
import numpy
import ephem
import getopt

import lsl.correlator.fx as fxc
from lsl.reader.ldp import LWA1DataFile
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drxSpectra.py - Read in DRX files and create a collection of 
time-averaged spectra.

Usage: drxSpectra.py [OPTIONS] file

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
	
	idf = LWA1DataFile(config['args'][0])
	
	nFramesFile = idf.getInfo('nFrames')
	srate = idf.getInfo('sampleRate')
	beampols = idf.getInfo('beampols')
	
	# Offset in frames for beampols beam/tuning/pol. sets
	config['offset'] = idf.offset(config['offset'])
	
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
	
	# Date & Central Frequnecy
	beginDate = ephem.Date(unix_to_utcjd(idf.getInfo('tStart')) - DJD_OFFSET)
	centralFreq1 = idf.getInfo('freq1')
	centralFreq2 = idf.getInfo('freq2')
	beam = idf.getInfo('beam')
	
	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Beam: %i" % beam
	print "Tune/Pols: %i" % beampols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], config['offset']*srate*beampols/4096)
	print "Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames / beampols)
	print "Chunks: %i" % nChunks
	
	# Sanity check
	if config['offset']*srate*beampols/4096 > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - config['offset']*srate*beampols/4096):
		raise RuntimeError("Requested integration time+offset is greater than file length")
		
	# Master loop over all of the file chunks
	standMapper = [4*(beam-1) + i for i in xrange(4)]
	masterWeight = numpy.zeros((nChunks, 4, LFFT))
	masterSpectra = numpy.zeros((nChunks, 4, LFFT))
	for i in range(nChunks):
		print "Working on chunk #%i of %i" % (i+1, nChunks)
		
		try:
			readT, t, data = idf.read(config['average']/nChunks)
		except Exception, e:
			print "Error: %s" % str(e)
			continue
			
		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate, ClipLevel=0)
		for stand in xrange(tempSpec.shape[0]):
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = int(readT*srate/LFFT)
			
		# We don't really need the data array anymore, so delete it
		del(data)
		
	# Now that we have read through all of the chunks, perform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
	
	# Frequencies
	freq1 = freq + centralFreq1
	freq2 = freq + centralFreq2
	
	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(4)))
	figsY = 4 / figsX
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
		ax.set_xlabel('Frequency [%s]' % units)
		ax.set_ylabel('P.S.D. [dB/RBW]')
		ax.set_xlim([freq.min(), freq.max()])
		ax.legend(loc=0)
		
		print "For beam %i, tune. %i, pol. %i maximum in PSD at %.3f %s" % (standMapper[i]/4+1, standMapper[i]%4/2+1, standMapper[i]%2, freq[numpy.where( spec[i,:] == spec[i,:].max() )][0], units)
		
	print "RBW: %.4f %s" % ((freq[1]-freq[0]), units)
	plt.subplots_adjust(hspace=0.35, wspace=0.30)
	plt.show()
	
	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
