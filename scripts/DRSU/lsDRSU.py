#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import getopt

from datetime import datetime

from lsl.reader import drsu


def usage(exitCode=None):
	print """lsDRSU.py - Simple `ls` type script for DRSUs.

Usage: lsDRSU.py [OPTIONS] md_device

Options:
-h, --help              Display this help information
-l                      Long listing format of size, mtime, and name
-t                      Sort by modification time
-1                      List one file per line
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	# Build up the configuration
	config = {}
	config['long'] = False
	config['mtimes'] = False
	config['onefile'] = False
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hlt1", ["help",])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-l'):
			config['long'] = True
		elif opt in ('-t'):
			config['mtimes'] = True
		elif opt in ('-1'):
			config['onefile'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def sortFilesMTime(x, y):
	"""
	Sort drsu.File objects by modicication time.
	"""
	
	if x.mtime < y.mtime:
		return 1
	elif x.mtime > y.mtime:
		return -1
	else:
		return 0


def main(args):
	config = parseOptions(args)
	device = config['args'][0]
	patterns = config['args'][1:]
	if len(patterns) == 0:
		patterns = ['*',]

	files = []
	for pattern in patterns:
		files.extend( drsu.globFiles(device, pattern) )
	if config['mtimes']:
		files.sort(sortFilesMTime)
	
	dates = []
	for f in files:
		dates.append( datetime.utcfromtimestamp(f.mtime) )
		
	if config['long']:
		for f,d in zip(files, dates):
			print f.mode, f.size, d.strftime("%b %d %H:%M"), f.name
	elif config['onefile']:
		for f in files:
			print f.name
	else:
		nCol = int(math.ceil(len(files) / 4.0))
		for i in xrange(nCol):
			out = files[i].name
			for j in xrange(1, 4):
				try:
					out = '%s\t%s' % (out, files[nCol*j+i].name)
				except IndexError:
					pass
			print out


if __name__ == "__main__":
	main(sys.argv[1:])
