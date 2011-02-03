#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import getopt
from datetime import datetime

from lsl.reader import tbw

try:
	from lsl.common.progress import ProgressBar
except ImportError:
	class ProgressBar(object):
		def __init__(self,max=100):
			self.amount = 0
			self.max = max
			self.span = 70
		
		def inc(self, amount=1):
			self.amount += amount

		def show(self):
			barSpan = self.span - 9
			nMarks = int(round(float(self.amount)/self.max * barSpan))
			bar = '=' * nMarks
			bar = bar+(' ' * (barSpan-nMarks))
			nte = "%5.1f%%" % (float(self.amount)/self.max*100)
			return "|%s| %s" % (bar, nte)


def usage(exitCode=None):
	print """splitTBW.py - split a TBW file containing multiple captures into
several single capture files

Usage: splitTBW.py [OPTIONS] file

Options:
-h, --help             Display this help information
-c, --count            Number of capturs to split off
-o, --offset           Number of captures to skip before splitting
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
			config['count'] = int(value)
		elif opt in ('-o', '--offset'):
			config['offset'] = int(value)
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
	nCaptures = sizeB / tbw.FrameSize / 300000

	print "Filename: %s" % filename
	print "Size:     %.1f MB" % (float(sizeB)/1024/1024)
	print "Captures: %i" % nCaptures
	print "==="

	if config['count'] > 0:
		nCaptures = config['count']
	nSkip = config['offset']

	print "Captures to Skip:  %i" % nSkip
	print "Captures to Split: %i" % nCaptures

	if config['date']:
		baseFilename = "%s_%s.dat" % (os.path.splitext(os.path.basename(filename))[0], '%s')
	else:
		baseFilename = "%s_s%s.dat" % (os.path.splitext(os.path.basename(filename))[0], '%04i')
	
	fh = open(filename, 'rb')
	frame = tbw.readFrame(fh)

	skip = 0
	while frame.parseID() != 1:
		frame = tbw.readFrame(fh)
		skip += 1
	fh.seek(fh.tell() - tbw.FrameSize)

	if skip != 0:
		print "Skipped %i frames at the beginning of the file" % skip

	for c in list(range(nCaptures+nSkip)):
		if c < nSkip:
			fh.seek(fh.tell() + tbw.FrameSize*300000)
			continue

		if config['date']:
			filePos = fh.tell()
			junkFrame = tbw.readFrame(fh)
			fh.seek(filePos)

			dt = datetime.utcfromtimestamp(junkFrame.getTime())
			captFilename = baseFilename % dt.isoformat()
		else:
			captFilename = baseFilename % (c+1)

		print "Writing capture #%i to file '%s'" % ((c+1-nSkip), captFilename)
		fhOut = open(captFilename, 'wb')
		pb = ProgressBar(max=300000)
		for i in list(range(300000)):
			cFrame = fh.read(tbw.FrameSize)
			fhOut.write(cFrame)

			pb.inc(amount=1)
			if i != 0 and i % 5000 == 0:
				sys.stdout.write(pb.show()+'\r')
				sys.stdout.flush()

		sys.stdout.write(pb.show()+'\r')
		sys.stdout.write('\n')
		sys.stdout.flush()
		fhOut.close()

	fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])