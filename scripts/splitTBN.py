#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import getopt
from datetime import datetime

from lsl.reader import tbn
from lsl.common.progress import ProgressBar


def usage(exitCode=None):
	print """splitTBN.py - split a TBN file containing multiple seconds into
several files

Usage: splitTBN.py [OPTIONS] file

Options:
-h, --help             Display this help information
-c, --count            Number of seconds to keep
-o, --offset           Number of seconds to skip before splitting
-d, --date             Label the split files with a date rather than a 
                       sequence number
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True

def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0
	config['count'] = 0
	config['date'] = False

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hc:o:d", ["help", "count=", "offset=", "date"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-c', '--count'):
			config['count'] = float(value)
		elif opt in ('-o', '--offset'):
			config['offset'] = float(value)
		elif opt in ('-d', '--date'):
			config['date'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config


def main(args):
	config = parseConfig(args)
	filename = config['args'][0]

	sizeB = os.path.getsize(filename)

	# Open the file and get some basic info about the data contained
	fh = open(filename, 'rb')
	sampleRate = tbn.getSampleRate(fh)
	nFramesX, nFramesY = tbn.getFramesPerObs(fh)

	nCaptures = sizeB / tbn.FrameSize / (nFramesX + nFramesY)

	print "Filename:    %s" % filename
	print "Size:        %.1f MB" % (float(sizeB)/1024/1024)
	print "Captures:    %i" % nCaptures
	print "Stands:      %i (%i x pol., %i y pol.)" % ((nFramesX+nFramesY), nFramesX, nFramesY)
	print "Sample Rate: %.2f kHz" % (sampleRate/1000.0)
	print "==="

	if config['count'] > 0:
		nCaptures = config['count'] * sampleRate / 512
	else:
		config['count'] = nCaptures * 512 / sampleRate
	nSkip = config['offset'] * sampleRate / 512

	print "Seconds to Skip:  %.2f (%i captures)" % (config['offset'], nSkip)
	print "Seconds to Split: %.2f (%i captures)" % (config['count'], nCaptures)

	# Make sure that the first frame in the file is the first frame if a capture 
	# (stand 1, pol 0).  If not, read in as many frames as necessary to get to 
	# the beginning of a complete capture.
	frame = tbn.readFrame(fh)
	stand, pol = frame.parseID()

	skip = 0
	while (2*(stand-1)+pol) != 0:
		frame = tbn.readFrame(fh)
		stand, pol = frame.parseID()
		skip += 1
	fh.seek(fh.tell() - tbn.FrameSize)

	if skip != 0:
		print "Skipped %i frames at the beginning of the file" % skip
	
	for c in list(range(nSkip)):
		if c < nSkip:
			fh.seek(fh.tell() + tbn.FrameSize*(nFramesX+nFramesY))
			continue

	if config['date']:
		filePos = fh.tell()
		junkFrame = tbn.readFrame(fh)
		fh.seek(filePos)

		dt = datetime.utcfromtimestamp(junkFrame.getTime())
		captFilename = "%s_%s.dat" % (os.path.splitext(os.path.basename(filename))[0], dt.isoformat())
	else:
		captFilename = "%s_s%04i.dat" % (os.path.splitext(os.path.basename(filename))[0], config['count'])

	print "Writing %.2f s to file '%s'" % (nCaptures*512/sampleRate, captFilename)
	fhOut = open(captFilename, 'wb')
	pb = ProgressBar(max=nCaptures)
	for c in list(range(int(nCaptures))):
		for i in list(range(nFramesX+nFramesY)):
			cFrame = fh.read(tbn.FrameSize)
			fhOut.write(cFrame)

		pb.inc(amount=1)
		if c != 0 and c % 100 == 0:
			sys.stdout.write(pb.show()+'\r')
			sys.stdout.flush()

	sys.stdout.write(pb.show()+'\r')
	sys.stdout.write('\n')
	sys.stdout.flush()
	fhOut.close()

	fh.close()
	
	

if __name__ == "__main__":
	main(sys.argv[1:])
