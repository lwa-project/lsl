#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a DR spectrometer file, plot the time averaged spectra for each 
polarization product."""

import sys
import math
import numpy
import ephem
import getopt

from lsl.reader.ldp import LWA1DataFile
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """drSpecSpectra.py - Read in DR spectrometer files and create a collection of 
time-averaged spectra.

Usage: drSpecSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-a, --average               Number of seconds of data to average for spectra 
                            (default = 10)
-q, --quiet                 Run drSpecSpectra in silent mode
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
	config['maxFrames'] = 10000
	config['average'] = 10.0
	config['output'] = None
	config['displayChunks'] = True
	config['verbose'] = True
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hqo:s:a:d", ["help", "quiet", "output=", "skip=", "average=", "disable-chunks"])
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
	
	idf = LWA1DataFile(config['args'][0])
	
	# Basic file informaiton
	nFramesFile = idf.getInfo('nFrames')
	srate = idf.getInfo('sampleRate')
	beam = idf.getInfo('beam')
	beampols = idf.getInfo('beampols')
	tInt = idf.getInfo('tInt')
	LFFT = idf.getInfo('LFFT')
	products = idf.getInfo('dataProducts')
	
	# Offset in frames for beampols beam/tuning/pol. sets
	config['offset'] = idf.offset(config['offset'])
	
	# Number of frames to integrate over
	maxFrames = config['maxFrames']
	nFrames = int(config['average'] / tInt)
	config['average'] = nFrames * tInt
	
	# Number of remaining chunks
	maxFramesTime = maxFrames*tInt
	nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))
	
	# Date & Central Frequnecy
	beginDate = ephem.Date(unix_to_utcjd(idf.getInfo('tStart')) - DJD_OFFSET)
	centralFreq1 = idf.getInfo('freq1')
	centralFreq2 = idf.getInfo('freq2')
	freq = numpy.fft.fftfreq(LFFT, d=1.0/srate)
	freq = numpy.fft.fftshift(freq)
	
	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Beam: %i" % beam
	print "Tune/Pols: %i" % beampols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile*tInt)
	print "---"
	print "Transform Length: %i channels" % LFFT
	print "Integration Time: %.3f s" % tInt
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], config['offset']*srate*beampols/4096)
	print "Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (config['average'], nFrames, nFrames)
	print "Chunks: %i" % nChunks
	
	# Sanity check
	if config['offset']/tInt > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - config['offset']/tInt):
		raise RuntimeError("Requested integration time+offset is greater than file length")
		
	# Master loop over all of the file chunks
	masterWeight = numpy.zeros((nChunks, 2*len(products), LFFT))
	masterSpectra = numpy.zeros((nChunks, 2*len(products), LFFT))
	for i in xrange(nChunks):
		print "Working on chunk #%i of %i" % (i+1, nChunks)
		
		try:
			readT, t, data = idf.read(config['average']/nChunks)
		except Exception, e:
			print "Error: %s" % str(e)
			continue
			
		## Integrate up the chunck
		data = data.mean(axis=1)
		
		## Save
		for stand in xrange(data.shape[0]):
			masterSpectra[i,stand,:] = data[stand,:]
			masterWeight[i,stand,:] = int(readT*srate/LFFT)
			
		## We don't really need the data array anymore, so delete it
		del(data)
		
	# Now that we have read through all of the chunks, perform the final averaging by
	# dividing by all of the chunks
	spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
	
	# Frequencies
	freq1 = freq + centralFreq1
	freq2 = freq + centralFreq2
	
	# The plots:  This is setup for the current configuration of 20 beampols
	fig = plt.figure()
	figsX = int(round(math.sqrt(2*len(products))))
	figsY = 2*len(products) / figsX
	# Put the frequencies in the best units possible
	freq1, units1 = bestFreqUnits(freq1)
	freq2, units2 = bestFreqUnits(freq2)
	
	for i in xrange(masterSpectra.shape[1]):
		if i/len(products) == 0:
			freq = freq1
			units = units1
		else:
			freq = freq2
			units = units2
			
		ax = fig.add_subplot(figsX,figsY,i+1)
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
				
		ax.set_title('Beam %i, Tune. %i, %s' % (beam, i/len(products), products[i % len(products)]))
		ax.set_xlabel('Frequency [%s]' % units)
		ax.set_ylabel('P.S.D. [dB/RBW]')
		ax.set_xlim([freq.min(), freq.max()])
		ax.legend(loc=0)
		
		print "For beam %i, tune. %i, %s maximum in PSD at %.3f %s" % (beam, i/len(products), products[i % len(products)], freq[numpy.where( spec[i,:] == spec[i,:].max() )][0], units)
		
	print "RBW: %.4f %s" % ((freq[1]-freq[0]), units)
	plt.subplots_adjust(hspace=0.35, wspace=0.30)
	plt.show()
	
	# Save spectra image if requested
	if config['output'] is not None:
		fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
