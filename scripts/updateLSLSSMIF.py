#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import getopt
import urllib
import hashlib
from datetime import datetime

from lsl.common.paths import data as dataPath


# Regular expression for finding the SSMIF version from the file contents
versionRE = re.compile(r'(?P<year>\d{4}) (?P<month>[a-zA-Z]{3,4}) (?P<day>\d{1,2}) by (?P<author>.*)')


def usage(exitCode=None):
	print """updateSSMIF.py - Update the internal LWA1/LWA-SV SSMIFs used by LSL.

Usage: updateSSMIF.py [OPTIONS]

Options:
-h, --help          Display this help information
-s, --lwasv         Update LWA-SV instead of LWA1
-u, --update        Update the default LWA1 SSMIF
-r, --revert        Revert the default LWA1 SSMIF to an older version
-f, --file          Update the default LWA1 SSMIF using the specified file

NOTE:  The -s/--lwasv option changes the behavior of all other options, e.g., 
       -u/--update updates the LWA-SV SSMIF.
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['site'] = 'lwa1'
	config['show'] = True
	config['update'] = False
	config['revert'] = False
	config['filename'] = None
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hsurf:", ["help", "lwasv", "update", "revert", "file="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
		
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-s', '--lwasv'):
			config['site'] = 'lwasv'
		elif opt in ('-u', '--update'):
			config['update'] = True
			
			config['revert'] = False
			config['filename'] = None
		elif opt in ('-r', '--revert'):
			config['revert'] = True
			
			config['update'] = False
			config['filename'] = None
		elif opt in ('-f', '--file'):
			config['filename'] = value
			
			config['update'] = False
			config['revert'] = False
		else:
			assert False
			
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def parseIndex(index):
	"""
	Parse the archive listing of SSMIF version and return a list of 
	filename/date tuples.
	"""
	
	from xml.etree import ElementTree as ET
	from BeautifulSoup import BeautifulSoup
	
	# Find the table
	start = index.find('<table>')
	stop  = index.find('</table>')
	index = index[start:stop+8]
	
	# Clean it up in such a way that ElementTree can parse it
	myMassage = [(re.compile('<!([^--].*)>'), lambda match: '<!--' + match.group(1) + '-->'), 
			   (re.compile('<hr>'), lambda match: ''), 
			   (re.compile('&nbsp;'), lambda match: ' '), 
			   (re.compile('<a.*?>(.*)</a>'), lambda mtch: mtch.group(1))]
	soup = BeautifulSoup(index, markupMassage=myMassage)
	index = soup.prettify()
	
	# Parse it
	table = ET.XML(index)
	rows = iter(table)
	
	# Extract the SSMIF entries
	versions = []
	for row in rows:
		values = [col.text for col in row]
		if len(values) != 5:
			 continue
			
		filename = values[1].lstrip().rstrip()
		if filename[:5] != 'SSMIF':
			 continue
		if filename.find('CURRENT') != -1:
			 continue
			
		date = filename.split('_', 1)[1]
		date = date.split('.', 1)[0]
		versions.append( (filename, date) )
		
	# Done
	return versions


def getMD5(filename, blockSize=262144):
	"""
	Compute the MD5 checksum of a file.
	"""
	
	fh = open(filename, 'r')
	
	m = hashlib.md5()
	while True:
		block = fh.read(blockSize)
		if len(block) == 0:
			break
		m.update(block)
		
	fh.close()
	
	return m.hexdigest()


def main(args):
	# Parse the options
	config = parseOptions(args)
	
	# Current LSL SSMIF
	if config['site'] == 'lwa1':
		_ssmif = os.path.join(dataPath, 'lwa1-ssmif.txt')
		_url = "http://lda10g.alliance.unm.edu/metadata/lwa1/ssmif/"
	elif config['site'] == 'lwasv':
		_ssmif = os.path.join(dataPath, 'lwa1-ssmif.txt')
		_url = "http://lda10g.alliance.unm.edu/metadata/lwasv/ssmif/"
	else:
		raise RuntimeError("Unknown site name: %s" % config['site'])
		
	urlToDownload = None
	if config['revert']:
		# Revert and Upgrade URL
		
		try:
			## Retrieve the list
			ah = urllib.urlopen(_url)
			index = ah.read()
			ah.close()
			
			## Parse
			versions = parseIndex(index)
			
			## Prompt the user for the version to revert to
			for i,(filename,date) in enumerate(versions):
				print "%i: %s" % (i, filename)
			i = -1
			while i not in range(0, len(versions)):
				i = raw_input("Enter SSMIF to revert to: ")
				try:
					i = int(i)
				except ValueError:
					print "-> Invalid value"
					i = -1
			print " "
			
			## Build the URL
			urlToDownload = "%s/%s" % (_url, versions[i][0])
		except Exception, e:
			print "Error:  Cannot process reversion, %s" % str(e)
			
	elif config['update']:
		# Update to the latest version
		
		urlToDownload = "%s/SSMIF_CURRENT.txt" % _url
		
	elif config['filename'] is not None:
		# Use the specified file
		
		urlToDownload = os.path.abspath(config['filename'])
		
	# Put the new SSMIF in place
	if urlToDownload is not None:
		## Retrieve
		try:
			ah = urllib.urlopen(urlToDownload)
			newSSMIF = ah.read()
			ah.close()
		except Exception, e:
			print "Error:  Cannot download SSMIF, %s" % str(e)
			
		## Save
		try:
			fh = open(_ssmif, 'wb')
			fh.write(newSSMIF)
			fh.close()
		except Exception, e:
			print "Error:  Cannot %s SSMIF, %s" % ('update' if config['update'] else 'revert', str(e))
			
	# Summarize the SSMIF
	if config['show']:
		## Filesystem information
		size = os.path.getsize(_ssmif)
		mtime = datetime.utcfromtimestamp(os.stat(_ssmif)[8])
		age = datetime.utcnow() - mtime
		
		## MD5 checksum
		md5 = getMD5(_ssmif)
		
		## SSMIF version (date)
		fh = open(_ssmif, 'r')
		lines = [fh.readline() for i in xrange(10)]
		fh.close()
		
		version = None
		for line in lines:
			mtch = versionRE.search(line)
			if mtch is not None:
				try:
					version = datetime.strptime("%s %s %s" % (mtch.group('year'), mtch.group('month'), mtch.group('day')), "%Y %b %d")
				except:
					version = datetime.strptime("%s %s %s" % (mtch.group('year'), mtch.group('month'), mtch.group('day')), "%Y %B %d")
				break
				
		print "LSL %s SSMIF%s:" % (config['site'].upper(), ' (updated)' if config['update'] else '',)
		print "  Size: %i bytes" % size
		print "  SSMIF Version: %s" % version.strftime("%Y %b %d")
		print "  File Last Modified: %s (%i day%s ago)" % (mtime.strftime("%Y-%m-%d %H:%M:%S UTC"), age.days, 's' if age.days != 1 else '')
		print "  MD5 Sum: %s" % md5


if __name__ == "__main__":
	main(sys.argv[1:])