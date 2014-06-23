#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a TBN file, plot the time averaged spectra for each digitizer input."""

import os
import sys
import math
import numpy
import ephem
import getopt

from lsl.common import stations
from lsl.reader.ldp import LWA1DataFile
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """tbnSpectra.py - Read in TBN files and create a collection of 
time-averaged spectra.

Usage: tbnSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-m, --metadata              Name of SSMIF or metadata tarball file to use for 
                            mappings
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
	config['metadata'] = ''
	config['offset'] = 0.0
	config['average'] = 10.0
	config['LFFT'] = 4096
	config['maxFrames'] = 2*260*750
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
			config['metadata'] = value
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
			config['keep'] = []
			for entry in value.split(','):
				if entry.find('-') != -1:
					start, stop = entry.split('-', 1)
					config['keep'].extend( range(int(start), int(stop)+1) )
				else:
					config['keep'].append( int(entry) )
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
	
	# Setup the LWA station information
	if config['metadata'] != '':
		try:
			station = stations.parseSSMIF(config['metadata'])
		except ValueError:
			station = metabundle.getStation(config['metadata'], ApplySDM=True)
	else:
		station = stations.lwa1
	antennas = station.getAntennas()
	
	# Length of the FFT
	LFFT = config['LFFT']
	
	idf = LWA1DataFile(config['args'][0])
	
	nFramesFile = idf.getInfo('nFrames')
	srate = idf.getInfo('sampleRate')
	antpols = len(antennas)
	
	# Offset in frames for beampols beam/tuning/pol. sets
	config['offset'] = idf.offset(config['offset'])
	
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
	beginDate = ephem.Date(unix_to_utcjd(idf.getInfo('tStart')) - DJD_OFFSET)
	centralFreq = idf.getInfo('freq1')
	
	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz" % centralFreq
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
	print "---"
	print "Offset: %.3f s (%i frames)" % (config['offset'], config['offset']*srate*antpols/512)
	print "Integration: %.3f s (%i frames; %i frames per stand/pol)" % (config['average'], nFrames, nFrames / antpols)
	print "Chunks: %i" % nChunks
	
	# Sanity check
	if config['offset']*srate*antpols/512 > nFramesFile:
		raise RuntimeError("Requested offset is greater than file length")
	if nFrames > (nFramesFile - config['offset']*srate*antpols/512):
		raise RuntimeError("Requested integration time+offset is greater than file length")
		
	# Master loop over all of the file chunks
	masterWeight = numpy.zeros((nChunks, antpols, LFFT))
	masterSpectra = numpy.zeros((nChunks, antpols, LFFT))
	
	for i in xrange(nChunks):
		print "Working on chunk #%i of %i" % (i+1, nChunks)
		
		try:
			readT, t, data = idf.read(config['average']/nChunks)
		except Exception, e:
			print "Error: %s" % str(e)
			continue
			
		# Calculate the spectra for this block of data and then weight the results by 
		# the total number of frames read.  This is needed to keep the averages correct.
		freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'], SampleRate=srate)
		for stand in xrange(tempSpec.shape[0]):
			masterSpectra[i,stand,:] = tempSpec[stand,:]
			masterWeight[i,stand,:] = int(readT*srate/LFFT)
			
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
					
	nPlot = len(js)
	if nPlot < 20:
		if nPlot % 4 == 0 and nPlot != 4:
			figsY = 4
		else:
			figsY = 2
		figsX = int(numpy.ceil(1.0*nPlot/figsY))
	else:
		figsY = 4
		figsX = 5
	figsN = figsX*figsY
	for i in xrange(antpolsDisp):
		# Normal plotting
		fig = plt.figure()
		for k in xrange(i*figsN, i*figsN+figsN):
			try:
				j = js[k]
				currSpectra = numpy.squeeze( numpy.log10(spec[j,:])*10.0 )
			except IndexError:
				break
			ax = fig.add_subplot(figsX, figsY, (k%figsN)+1)
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
