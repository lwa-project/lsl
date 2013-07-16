#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to gather information about the Python interpreter, modules, 
C libraries, numpy installation, and LSL installation to help with 
debugging and install issues.
"""

import os
import re
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
	notFound = 0
	for mod in ('numpy', 'scipy', 'pyfits', 'ephem', 'aipy'):
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
				try:	
					versionRE = re.compile(r'%s-(?P<version>[\d\.]+)-py.*' % mod)
					mtch = versionRE.search(eval("%s.__file__" % mod))
					version = mtch.group('version')
				except:
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
		
		for lib in ('libatlas', 'libcblas', 'libfftw3', 'libgdbm'):
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
	# Compiler check
	#
	import shutil
	import tempfile
	from distutils import sysconfig
	from distutils import ccompiler
	compiler = ccompiler.new_compiler()
	sysconfig.customize_compiler(compiler)
	cc = compiler.compiler
	
	print "Compiler: %s" % cc[0]
	
	tmpdir = tempfile.mkdtemp()
	curdir = os.getcwd()
	os.chdir(tmpdir)
	
	fh = open('test.c', 'w')
	fh.write(r"""#include <omp.h>
#include <stdio.h>
int main() {
#pragma omp parallel
printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}
""")
	fh.close()
	
	cmd = cc
	cmd.extend(['-fopenmp', 'test.c', '-o', 'test', '-lgomp'])
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	o, e = p.communicate()
	print "Compiler OpenMP Support: %s" % ("Yes" if p.returncode == 0 else "No",)
	if p.returncode != 0:
		o = o.split('\n')[:-1]
		for i in xrange(len(o)):
			o[i] = '  %s' % o[i]
		o = '\n'.join(o)
		e = e.split('\n')[:-1]
		for i in xrange(len(e)):
			e[i] = '  %s' % e[i]
		e = '\n'.join(e)
		
		print "Compiler OpenMP Test Command:"
		print "  %s" % ' '.join(cmd)
		print "Compiler OpenMP Test Output:"
		print o
		print "Compiler OpenMP Test Errors:"
		print e
		
	os.chdir(curdir)
	shutil.rmtree(tmpdir)
	
	p = subprocess.Popen([cc[0], '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	o, e = p.communicate()
	e = e.split('\n')[:-1]
	for i in xrange(len(e)):
		e[i] = '  %s' % e[i]
	e = '\n'.join(e)
	print "Compiler Version:"
	print e
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
		import lsl, lsl.version
		lfp,junk = os.path.split(lsl.__file__)
		lfp = os.path.join(lfp, 'correlator', '_core.so')
		lfp = os.path.realpath(lfp)

		p = subprocess.Popen(['file', lfp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		o, e = p.communicate()
		junk, lslLinkage = o.split(None, 1)
		lslLinkage = lslLinkage.replace('\n', '')

		lfp, junk = os.path.split(lsl.__file__)
		print "LSL Path: %s" % lfp
		print "LSL Version: %s" % lsl.version.version
		print "LSL Linkage: %s" % lslLinkage
	except ImportError, e:
		print "LSL Import Error: %s" % str(e)


if __name__ == "__main__":
	main(sys.argv[1:])

