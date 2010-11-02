#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script that reads in TBN data and runs a cross-correlation on it.  
The results are saved in the Miriad UV format."""

import os
import re
import sys
import aipy
import numpy
import getopt

import lsl.common.stations as lwa_common
import lsl.statistics.robust as robust
import lsl.reader.tbn as tbn
import lsl.reader.errors as errors
import lsl.correlator.uvUtils as uvUtils
import lsl.correlator.fx as fxc

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def isValidDate(date):
	dateRE = re.compile(r'(?P<valid>\d{4}/\d{1,2}/\d{1,2}[T ]\d{1,2}:\d{1,2}:\d{1,2})')

	try:
		mtch = re.match(dateRE, date)
		mtch.group('valid')
		return True
	except:
		return False


def usage(exitCode=None):
	print """correlateTBN.py - cross-correlate data in a TBN file

Usage: readS60.py [OPTIONS] file

Options:
-h --help             Display this help information
-c --central-freq     Central frequency of the observations in MHz
-d --date             Specify the date when the data were obtained
                      (YYYY/MM/DD HH:MM:SS)
-f --fft-length       Set FFT length (default = 512)
-t --avg-time         Window to average visibilities in time (seconds; 
                      default = 6 s)
-s --samples          Number of average visibilities to generate
                      (default = 10)
-q --quiet            Run correlateTBN in silent mode
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['avgTime'] = 6
	config['LFFT'] = 512
	config['cFreq'] = 38.0e6
	config['date'] = '2010/01/01 00:00:00'
	config['samples'] = 10
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hqc:d:l:t:s:", ["help", "quiet", "central-freq=", "date=", "fft-length=", "avg-time=", "samples="])
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
		elif opt in ('-c', '--central-freq'):
			config['cFreq'] = float(value)*1e6
		elif opt in ('-d', '--date'):
			if not isValidDate(value):
				print "Invalid Date: %s" % value
				usage(exitCode=1)
			else:
				config['date'] = value
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		elif opt in ('-t', '--avg-time'):
			config['avgTime'] = int(value)
		elif opt in ('-s', '--samples'):
			config['samples'] = int(samples)
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config


def main(args):
	# Parse command line options
	config = parseConfig(args)
	filename = config['args'][0]

	# Length of the FFT
	LFFT = config['LFFT']

	# Number of frames to read in at once and average
	nFrames = int(391*config['avgTime']/2)

	# Setup the LWA station information and get the JD of observations
	lwa1 = lwa_common.lwa1()
	stands = lwa1.getStands(config['date'])
	lwa1_OO = lwa1.getObserver(date=config['date'], JD=False)
	jd = float(lwa1_OO.date) + 2415020.0

	fh = open(filename, "rb", buffering=tbn.FrameSize)
	test = tbn.readFrame(fh)
	if not test.header.isTBN():
		raise errors.notTBNError()
	fh.seek(0)
	nFpO = tbn.getFramesPerObs(fh)
	nFpO = nFpO[0] + nFpO[1]
	sampleRate = tbn.getSampleRate(fh, nFrames=nFpO)
	nInts = os.path.getsize(filename) / tbn.FrameSize / nFpO
	nSets = os.path.getsize(filename) / tbn.FrameSize / nFpO / nFrames

	print "TBN Data:  %s" % test.header.isTBN()
	print "Samples per observations: %i per pol." % (nFpO/2)
	print "Filter code is: %i" % tbn.getSampleRate(fh, nFrames=nFpO, FilterCode=True)
	print "Sampling rate is: %i Hz" % sampleRate
	print "Captures in file: %i (%.1f s)" % (nInts, nInts*512 / sampleRate)
	print "=="
	print "Station: %s" % lwa1.name
	print "Date observed: %s" % config['date']
	print "Julian day: %.5f" % jd
	print "Integration Time: %.2f s" % (512*nFrames/sampleRate)
	print "Number of averaged samples: %i" % nSets

	if config['samples'] > nSets:
		config['samples'] = nSets

	refTime = 0.0
	setTime = 0.0
	for s in range(config['samples']):
		count = {}
		masterCount = 0
		iTime = 0
		data = numpy.zeros((20,nFrames*512), dtype=numpy.complex64)
		for i in range(20*nFrames):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbn.readFrame(fh, Verbose=False, SampleRate=sampleRate, CentralFreq=config['cFreq'], Gain=22)
			except errors.eofError:
				break
			except errors.syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFrameSize-1)
				continue
			except errors.numpyError:
				break

			stand,pol = cFrame.header.parseID()
			aStand = 2*(stand-1)+pol
			if i == 0:
				if s == 0:
					refTime = 0
				setTime = s*512*nFrames/sampleRate
			if cFrame.header.frameCount % 500 == 0:
				print "%2i,%1i  %6.3f  %5i  %5i" % (stand, pol, cFrame.getTime(), cFrame.header.frameCount, cFrame.header.secondsCount)
			
			if aStand not in count.keys():
				count[aStand] = 0

			data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
			
			count[aStand] = count[aStand] + 1
			masterCount = masterCount + 1

		print "Working on set #%i (%i seconds after set #1)" % ((s+1), (setTime-refTime))
		freqYY, outYY = fxc.FXCorrelator(data, stands, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=True,  SampleRate=sampleRate, CentralFreq=config['cFreq'])
		uvw = uvUtils.computeUVW(stands, HA=((setTime-refTime)/3600.0), freq=freqYY, IncludeAuto=True)
		ccList = uvUtils.getBaselines(stands, IncludeAuto=True, Indicies=True)

		import pylab
		for i in [0, 1]:
			print freqYY.shape, outYY.shape, outYY[i,:].shape
			#pylab.plot(freqYY/1e6, numpy.log10(numpy.abs(outYY[i,:]))*10, label='%i' % i)
			pylab.plot(freqYY/1e6, numpy.angle(outYY[i,:]), label='%i' % i)
		pylab.legend(loc=0)
		pylab.show()

		toUse = numpy.where( (freqYY>10.0e6) & (freqYY<88.0e6) )
		toUse = toUse[0]

		if s == 0:
			uvoutname = os.path.basename(filename)
			uvoutname = uvoutname.replace('.dat', '.uv')
			if uvoutname.find('.uv') == -1:
				uvoutname = "%s.uv" % uvoutname

			uv = aipy.miriad.UV(uvoutname, status='new')
			tags = ['nants', 'nchan', 'npol', 'pol', 'sfreq', 'sdf', 'longitu', 'latitud', 'dec', 'ra', 'source']
			codes = ['i', 'i', 'i', 'i', 'r', 'r', 'r', 'r', 'r', 'r', 'a']
			for tag,code in zip(tags, codes):
				print "%s -> %s" % (tag, code)
				uv.add_var(tag, code)
			uv['source'] = 'zenith'
			uv['longitu'] = -107.628
			uv['latitud'] = 34.070
			uv['dec'] = 34.070
			uv['nants'] = stands.shape[0]
			uv['nchan'] = toUse.shape[0]
			uv['npol'] = 1
			uv['sfreq'] = freqYY[toUse[0]]
			uv['sdf'] = freqYY[toUse[1]]-freqYY[toUse[0]]

		uv['pol'] = -6
		for i in range(uvw.shape[0]):
			preamble = (numpy.squeeze(uvw[i,:,toUse[0]]), jd+(setTime-refTime)/3600.0/24.0, ccList[i])

			current = numpy.squeeze( outYY[i,toUse] )
			mask = numpy.zeros_like(current)

			uv.write(preamble, current, mask)
		del(data)

	fh.close()

if __name__ == "__main__":
	main(sys.argv[1:])

	import aipy
	import robust



	# Pre 8/25/10
	#stands = numpy.array([214, 212, 228, 204, 225, 210, 187, 174, 123, 125])
	# Post 8/25/10
	stands = numpy.array([214, 212, 228, 206, 127, 157, 187, 208, 123, 125])

	#fh = open("/home/jdowell/055424_740_threeTBW.dat", "rb")
	#fh = open("/home/jdowell/triple12M_TBW_Midnight_Sept_3.dat", "rb")
	#fh = open("/home/jdowell/12M_TBW_430pm_Sept_8.dat", "rb")
	fh = open("/home/jdowell/12M_TBW_900pm_Sept_8.dat", "rb")
	jd = 2455448.62500
	#fh = open("/home/jdowell/tbw_test1_091310.dat", "rb")
	#jd = 2455453.31531
	nSamples = 300000

	refTime = 0.0
	setTime = 0.0

	for s in range(1):
		count = {}
		masterCount = 0
		iTime = 0
		data = numpy.zeros((20,12000000), dtype=numpy.int16)
		for i in range(nSamples):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = readTBWFrame(fh, Verbose=False)
			except eofError:
				break
			except syncError:
				print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFrameSize-1)
				continue
			except numpyError:
				break

			stand = cFrame.header.parseID()
			if i == 0:
				if s == 0:
					refTime = cFrame.header.secondsCount
				setTime = cFrame.header.secondsCount
			print "%2i  %6.3f  %5i  %5i" % (stand, cFrame.getTime(), cFrame.header.frameCount, cFrame.header.secondsCount)
			if stand not in count.keys():
				count[stand] = 0

			if stand < 11:
				if cFrame.getDataBits() == 12:
					data[2*(stand-1),  count[stand]*400:(count[stand]+1)*400] = numpy.squeeze(cFrame.data.xy[0,:])
					data[2*(stand-1)+1,count[stand]*400:(count[stand]+1)*400] = numpy.squeeze(cFrame.data.xy[1,:])
				else:
					data[2*(stand-1),  count[stand]*1200:(count[stand]+1)*1200] = numpy.squeeze(cFrame.data.xy[0,:])
					data[2*(stand-1)+1,count[stand]*1200:(count[stand]+1)*1200] = numpy.squeeze(cFrame.data.xy[1,:])

			count[stand] = count[stand] + 1
			masterCount = masterCount + 1

		LFFT = 64
		freqY, specY = calcSpectra(data[0::2,:], LFFT=LFFT, PolyphaseFilter=True, verbose=True)
		freqX, specX = calcSpectra(data[1::2,:], LFFT=LFFT, PolyphaseFilter=True, verbose=True)

		rfiClip = 5

		rfiX = {}
		bandpassX = {}
		rfiY = {}
		bandpassY = {}

		valid = numpy.where( (freqX>20.0e6) & (freqX<88.0e6) )
		valid = valid[0]
		for i in range(specX.shape[0]):
			sm = robust.robustMean(specX[i,valid])
			ss = robust.robustSigma(specX[i,valid])

			rfiX[i] = freqX[valid[numpy.where( specX[i,valid] > sm+rfiClip*ss )]]
			bandpassX[i] = numpy.zeros_like(valid)*sm

		valid = numpy.where( (freqY>20.0e6) & (freqY<88.0e6) )
		valid = valid[0]
		for i in range(specY.shape[0]):
			sm = robust.robustMean(specY[i,valid])
			ss = robust.robustSigma(specY[i,valid])
			
			rfiY[i] = freqY[valid[numpy.where( specY[i,valid] > sm+rfiClip*ss )]]
			bandpassY[i] = numpy.zeros_like(valid)*sm

		print "Working on set #%i (%i seconds after set #1)" % ((s+1), (setTime-refTime))
		freqYY, outYY = FXCorrelator(data[0::2,:], stands, LFFT=LFFT, PolyphaseFilter=True, verbose=True)
		freqXX, outXX = FXCorrelator(data[1::2,:], stands, LFFT=LFFT, PolyphaseFilter=True, verbose=True)
		freqYX, outYX = FXCorrelator(data[0::2,:], stands, LFFT=LFFT, PolyphaseFilter=True, verbose=True, CrossPol=data[1::2,:])
		freqXY, outXY = FXCorrelator(data[1::2,:], stands, LFFT=LFFT, PolyphaseFilter=True, verbose=True, CrossPol=data[0::2,:])
		uvw = uvUtils.computeUVW(stands, HA=((setTime-refTime)/3600.0), freq=freqXX)
		ccList = uvUtils.getBaselines(stands, Indicies=True)

		toUse = numpy.where( (freqYY>20.0e6) & (freqYY<88.0e6) )
		toUse = toUse[0]

		if s == 0:
			uv = aipy.miriad.UV('test-20100908-2100.uv', status='new')
			tags = ['nants', 'nchan', 'npol', 'pol', 'sfreq', 'sdf', 'longitu', 'latitud', 'dec', 'ra', 'source']
			codes = ['i', 'i', 'i', 'i', 'r', 'r', 'r', 'r', 'r', 'r', 'a']
			for tag,code in zip(tags, codes):
				print "%s -> %s" % (tag, code)
				uv.add_var(tag, code)
			uv['source'] = 'zenith'
			uv['longitu'] = -107.628
			uv['latitud'] = 34.070
			uv['dec'] = 34.070
			uv['nants'] = stands.shape[0]
			uv['nchan'] = toUse.shape[0]
			uv['npol'] = 4
			uv['sfreq'] = freqYY[toUse[0]]
			uv['sdf'] = freqYY[toUse[1]]-freqYY[toUse[0]]

		uv['pol'] = -6
		for i in range(uvw.shape[0]):
			preamble = (numpy.squeeze(uvw[i,:,toUse[0]]), jd+(setTime-refTime)/3600.0/24.0, ccList[i])

			current = numpy.squeeze( outYY[i,toUse] )
			#current = current / bandpassY[ccList[i][0]] / bandpassY[ccList[i][1]]

			mask = numpy.zeros_like(current)
			for cf in rfiY[ccList[i][0]]:
				bad = numpy.where( cf == freqYY[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1
			for cf in rfiY[ccList[i][1]]:
				bad = numpy.where( cf == freqYY[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1

			uv.write(preamble, current, mask)

		uv['pol'] = -5
		for i in range(uvw.shape[0]):
			preamble = (numpy.squeeze(uvw[i,:,toUse[0]]), jd+(setTime-refTime)/3600.0/24.0, ccList[i])

			current = numpy.squeeze( outXX[i,toUse] )
			#current = current / bandpassX[ccList[i][0]] / bandpassX[ccList[i][1]]

			mask = numpy.zeros_like(current)
			for cf in rfiX[ccList[i][0]]:
				bad = numpy.where( cf == freqXX[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1
			for cf in rfiX[ccList[i][1]]:
				bad = numpy.where( cf == freqXX[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1

			uv.write(preamble, current, mask)

		uv['pol'] = -8
		for i in range(uvw.shape[0]):
			preamble = (numpy.squeeze(uvw[i,:,toUse[0]]), jd+(setTime-refTime)/3600.0/24.0, ccList[i])

			current = numpy.squeeze( outYX[i,toUse] )
			#current = current / bandpassY[ccList[i][0]] / bandpassX[ccList[i][1]]

			mask = numpy.zeros_like(current)
			for cf in rfiY[ccList[i][0]]:
				bad = numpy.where( cf == freqYY[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1
			for cf in rfiX[ccList[i][1]]:
				bad = numpy.where( cf == freqXX[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1

			uv.write(preamble, current, mask)

		uv['pol'] = -7
		for i in range(uvw.shape[0]):
			preamble = (numpy.squeeze(uvw[i,:,toUse[0]]), jd+(setTime-refTime)/3600.0/24.0, ccList[i])

			current = numpy.squeeze( outXY[i,toUse] )
			#current = current / bandpassX[ccList[i][0]] / bandpassY[ccList[i][1]]

			mask = numpy.zeros_like(current)
			for cf in rfiX[ccList[i][0]]:
				bad = numpy.where( cf == freqXX[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1
			for cf in rfiY[ccList[i][1]]:
				bad = numpy.where( cf == freqYY[toUse] )
				if len(bad[0]) != 0:
					mask[bad] = 1

			uv.write(preamble, current, mask)

		del(data)

	fh.close()
	
