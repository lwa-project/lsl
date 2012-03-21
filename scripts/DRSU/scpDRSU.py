#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import time
import getopt
import getpass
import paramiko

from datetime import datetime

from lsl.reader import drsu
from lsl.common.progress import ProgressBar


def usage(exitCode=None):
	print """scpDRSU.py - Simple `scp` type script for DRSUs.

Usage: scpDRSU.py md_device tag user@host:dest

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
	tag = config['args'][1]
	dest = config['args'][2]
	user, dest = dest.split('@', 1)
	hostname, dest = dest.split(':', 1)

	# Make sure we can get the file
	fileObject = drsu.getFileByName(device, tag)
	if fileObject is None:
		print "ERROR:  cannot find '%s' on device '%s'" % (tag, device)
		sys.exit(1)
	
	ssh = paramiko.SSHClient()
	ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
	
	password = getpass.getpass("%s@%s's password: " % (user, hostname))
	ssh.connect(hostname, username=user, password=password)
	
	# Deal with a destination directory or filename
	pyin, pyout, pyerr = ssh.exec_command("python -c 'import os,sys; print os.path.isdir(sys.argv[2]);' - %s" % dest)
	result = pyout.readline()
	if result.find('True') != -1:
		outname = os.path.join(dest, '%s_%s.dat' % (fileObject.name, fileObject.mode))
	else:
		outname = dest
	
	# Verbosity
	if config['verbose']:
		print "%s on %s -> %s on %s" % (tag, device, outname, hostname)
	
	# Actually perform the copy
	tStart = time.time()
	nSections = fileObject.size / fileObject.chunkSize

	pb = ProgressBar(max=nSections)
	sys.stdout.write(pb.show()+'\r')
	sys.stdout.flush()

	fileObject.open()
	scpin, scpout, scperr = ssh.exec_command("dd bs=%i of=%s" % (fileObject.chunkSize, outname), bufsize=10*fileObject.chunkSize)
	for c in xrange(nSections):
		part = fileObject.fh.read(fileObject.chunkSize)
		scpin.write(part)
		scpin.flush()

		pb.inc(amount=1)
		if c % 10 == 0:
			sys.stdout.write(pb.show()+'\r')
			sys.stdout.flush()
			
	leftover = fileObject.size - nSections*fileObject.chunkSize
	if leftover > 0:
		part = fileObject.fh.read(fileObject.chunkSize)
		scpin.write(part[:leftover])
		scpin.flush()
		
	fileObject.close()
	ssh.close()

	sys.stdout.write(pb.show()+'\r')
	sys.stdout.write('\n')
	sys.stdout.flush()

	tStop = time.time()

	print "Finished in %.1f seconds (%.1f MB/s)" % ((tStop - tStart), float(fileObject.size) / 1024**2 / (tStop - tStart))
	
	# Done


if __name__ == "__main__":
	main(sys.argv[1:])
	
