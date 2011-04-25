# -*- coding: utf-8 -*-
#
# Supplemental configuration file for the LWA Software Library documentation 
# build configuration file (conf.py in this directory).  This file should set 
# the version (long and short), year, and release notes information/files.

import sys, os, re


def version2version(version):
	try:
		version, junk = version.split('-', 1)
	except ValueError:
		pass
	parts = version.split('.')
	
	value = 0.0
	for i,p in enumerate(parts):
		value += float(p) / 10**i
		
	return value


def isDevelopment():
	fh = open(os.path.join('..', '..', 'setup.cfg'))
	
	devel = False
	for line in fh:
		if line.find('tag_build =') == 0:
			keyword, value = line.split(' = ')
			if len(value) > 3:
				devel = True
			break
	fh.close()
	
	return devel


def getLongVersion():
	fh = open(os.path.join('..', '..', 'VERSION'))
	longVersion = fh.read().strip()
	fh.close()
	
	if isDevelopment():
		longVersion = "%s-dev" % longVersion
		
	return longVersion


def getShortVersion():
	longVersion = getLongVersion()
	shortVersion = '.'.join(longVersion.split('.')[0:2])
	
	return shortVersion


def updateReleaseNotes():
	refVersion = version2version(getShortVersion())
	inOlder = False
	
	try:
		fh = open(os.path.join('..', '..', 'CHANGELOG'))
		
		versionRE = re.compile(r'^\d+\.\d+(\.\d+)?')
		
		output = "Release Notes\n=============\n\n"
		for line in fh:
			if len(line) == 0:
				continue
				
			mtch = versionRE.match(line)
			if mtch is not None:
				lineChar = '-'
				if version2version(line) < refVersion:
					lineChar = '+'
					if not inOlder:
						output += "Older Versions\n--------------\n"
						inOlder = True
				output += "Version %s%s\n" % (line, lineChar*(len(line)+8))
			else:
				output += "%s" % line
				
		fh.close()
		
		fh = open('rnotes.rst', 'w')
		fh.write(output)
		fh.close()
	except:
		return False
	
	return True
		