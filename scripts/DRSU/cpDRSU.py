#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import time
import getopt

from datetime import datetime

from lsl.reader import drsu
from lsl.common.progress import ProgressBar


def usage(exitCode=None):
	print """cpDRSU.py - Simple `cp` type script for DRSUs.

Usage: cpDRSU.py md_device tag dest

Options:
-h, --help              Display this help information
-v, --verbose           Be verbose
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	# Build up the configuration
	config = {}
	config['verbose'] = False
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hv", ["help", "verbose"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-v', '--verbose'):
			config['verbose'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	device = config['args'][0]
	pattern = config['args'][1]
	dest = config['args'][2]

	# Get files to copy off the DRSU
	toCopy = drsu.globFiles(device, pattern)

	# Make sure we can get the file
	if not len(toCopy):
		print "ERROR:  cannot find '%s' on device '%s'" % (pattern, device)
		sys.exit(1)
	
	# Make sure we have a place to go to
	if len(toCopy) > 1 and not os.path.isdir(dest):
		print "ERROR:  target '%s' is not a directory" % dest
	
	tStart = time.time()
	copiedBytes = 0
	for fileObject in toCopy:
		# Deal with a destination directory or filename
		if os.path.isdir(dest):
			outname = os.path.join(dest, '%s_%s.dat' % (fileObject.name, fileObject.mode))
		else:
			outname = dest
	
		# Verbosity
		if config['verbose']:
			print "%s on %s -> %s" % (fileObject.name, device, outname)
	
		# Actually perform the copy
		nSections = fileObject.size / fileObject.chunkSize

		pb = ProgressBar(max=nSections)
		sys.stdout.write(pb.show()+'\r')
		sys.stdout.flush()

		fileObject.open()
		ofh = open(outname, 'wb')
		for c in xrange(nSections):
			part = fileObject.fh.read(fileObject.chunkSize)
			ofh.write(part)
			ofh.flush()

			pb.inc(amount=1)
			if c % 10 == 0:
				sys.stdout.write(pb.show()+'\r')
				sys.stdout.flush()
				
		leftover = fileObject.size - nSections*fileObject.chunkSize
		if leftover > 0:
			part = fileObject.fh.read(fileObject.chunkSize)
			ofh.write(part[:leftover])
			ofh.flush()
		
		fileObject.close()
		ofh.close()

		copiedBytes += fileObject.size

		sys.stdout.write(pb.show()+'\r')
		sys.stdout.write('\n')
		sys.stdout.flush()

	tStop = time.time()

	print "Finished in %.1f seconds (%.1f MB/s)" % ((tStop - tStart), float(copiedBytes) / 1024**2 / (tStop - tStart))
	
	# Done


if __name__ == "__main__":
	main(sys.argv[1:])
	
