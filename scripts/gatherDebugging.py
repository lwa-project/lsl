#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to gather information about the Python interpreter, modules, 
C libraries, numpy installation, and LSL installation to help with 
debugging and install issues.
"""

import os
import sys
import platform
import subprocess


def main(args):
	#
	# Clean the path
	#
	sys.path = sys.path[1:]

	#
	# Python interpreter
	#
	print "Executable path: %s" % sys.executable
	print "Platform: %s" % sys.platform
	print "Version: %s" % sys.version
	print "API: %s" % sys.api_version
	print "Bits: %s\nLinkage: %s" % platform.architecture()

	print " "
	
	#
	# Python Module Check
	#
	for mod in ('numpy', 'scipy', 'pyfits', 'ephem', 'aipy', 'jinja2'):
		try:
			exec "import %s" % mod
		except ImportError, e:
			notFound += 1
			if (str(e)).find('not found') != -1:
				print "%s: not found" % mod
			else:
				print "%s: import error '%s'" % (mod, str(e))
		else:
			try:
				version = eval("%s.__version__" % mod)
			except AttributeError:
				version = "unknown"
			print "%s:  version %s" % (mod.capitalize(), version)
			
	print " "
	
	#
	# C library check (linux only)
	#
	if platform.system() == 'Linux':
		p = subprocess.Popen(['ldconfig', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		o, e = p.communicate()
		o = o.split('\n')
		
		for lib in ('libatlas', 'libcblas', 'libfftw3', 'libnova', 'libgdbm'):
			found = False
			currPath = None
			
			for line in o:
				if len(line) == 0:
					continue				
				elif line[0] != '\t':
					currPath, junk = line.split(':', 1)
					continue
				elif line.find(lib) != -1:
					found = True
					break
					
			print "%s: %s" % (lib, "found in %s" % currPath if found else "not found")
			
		print " "
	
	#
	# Numpy
	#
	try:
		import numpy
		nfp,junk = os.path.split(numpy.__file__)
		nfp = os.path.join(nfp, 'core', 'umath.so')
		nfp = os.path.realpath(nfp)

		p = subprocess.Popen(['file', nfp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		o, e = p.communicate()
		junk, numpyLinkage = o.split(None, 1)
		numpyLinkage = numpyLinkage.replace('\n', '')

		nfp, junk = os.path.split(numpy.__file__)
		print "Numpy Path: %s" % nfp
		print "Numpy Version: %s" % numpy.version.version
		print "Numpy Linkage: %s" % numpyLinkage
	except ImportError, e:
		print "Numpy Import Error: %s" % str(e)

	print " "

	#
	# LSL
	#
	try:
		import lsl
		lfp,junk = os.path.split(lsl.__file__)
		lfp = os.path.join(lfp, 'correlator', '_core.so')
		lfp = os.path.realpath(lfp)

		p = subprocess.Popen(['file', lfp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		o, e = p.communicate()
		junk, lslLinkage = o.split(None, 1)
		lslLinkage = lslLinkage.replace('\n', '')

		lfp, junk = os.path.split(lsl.__file__)
		print "LSL Path: %s" % lfp
		print "LSL Base Version: %s" % lsl.__version__
		print "LSL Linkage: %s" % lslLinkage
	except ImportError, e:
		print "LSL Import Error: %s" % str(e)


if __name__ == "__main__":
	main(sys.argv[1:])

