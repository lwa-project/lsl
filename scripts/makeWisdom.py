#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import getopt

from lsl.misc import wisdom


def usage(exitCode=None):
	print """makeWisdom.py - Build LSL-specific wisdom files for FFTW and PyFFTW

Usage: makeWisdom.py [OPTIONS]

Options:
-h, --help             Display this help information
-f, --fftw-only        Build/show only FFTW wisdom only (default = no)
-p, --pyfftw-only      Build/show only PyFFTW wisdom only (default = no)
-s, --show             Show information about the avaliable wisdom
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['show'] = False
	config['fftw'] = True
	config['pyfftw'] = True

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hfps", ["help", "fftw-only", "pyfftw-only", "show"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-f', '--fftw-only'):
			config['fftw'] = True
			config['pyfftw'] = False
		elif opt in ('-p', '--pyfftw-only'):
			config['fftw'] = False
			config['pyfftw'] = True
		elif opt in ('-s', '--show'):
			config['show'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config

def main(args):
	config = parseConfig(args)
	
	if config['show']:
		wisdom.show(FFTW=config['fftw'], PyFFTW=config['pyfftw'])
	else:
		wisdom.make(FFTW=config['fftw'], PyFFTW=config['pyfftw'])


if __name__ == "__main__":
	main(sys.argv[1:])
