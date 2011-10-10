#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import getopt
from hashlib import md5

from lsl.reader import drsu


def md5sum(f, size, buf_size=8192):
	"""
	Compute the MD5 sum, binary style.
	"""
	
	m = md5()
	
	count = 0
	while count < size:
		data = f.read(buf_size)
		count += buf_size
		
		m.update(data)
		
	# We return the md5 hash in hexadecimal format
	return m.hexdigest()


def usage(exitCode=None):
	print """md5sumDRSU.py - Simple `md5sum` type script for DRSUs.

Usage: md5sumDRSU.py md_device tag [tag [tag [...]]]

Options:
-h, --help              Display this help information
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
		opts, args = getopt.getopt(args, "h", ["help",])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	
	device = config['args'][0]
	patterns = config['args'][1:]

	files = []
	for pattern in patterns:
		files.extend( drsu.globFiles(device, pattern) )
	
	for fileObject in files:
		fileObject.open()
		checksum = md5sum(fileObject.fh, fileObject.size, buf_size=fileObject.chunkSize)
		print "%s  %s on %s" % (str(checksum), fileObject.name, device)
	
	
if __name__ == "__main__":
	main(sys.argv[1:])
