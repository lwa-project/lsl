# -*- coding: utf-8 -*

import ez_setup
ez_setup.use_setuptools()

import os
import imp
import glob
import unittest
import commands

from setuptools import setup, Extension, Distribution, find_packages
try:
	import numpy
except ImportError:
	pass


def get_version():
	"""Read the VERSION file and return the version number as a string."""

	return open('VERSION').read().strip()

def get_description(filename):
	"""Read in a README-type file and return the contents of the DESCRIPTION
	section."""

	desc = ''
	fh = open(filename, 'r')
	lines = fh.readlines()
	fh.close()

	inDescription = False
	for line in lines:
		line = line.replace('\n', '')
		line = line.replace('\t', '')
		if line.find('DESCRIPTION') == 0:
			inDescription = True
			continue
		if line.find('REQUIREMENTS') == 0:
			inDescription = False
			break
		if inDescription:
			desc = ' '.join([desc, line])

	return desc


def get_fftw():
	"""Use pkg-config (if installed) to figure out the C flags and linker flags
	needed to compile a C program with FFTW3.  If FFTW3 cannot be found via
	pkg-config, some 'sane' values are returned."""

	status, output = commands.getstatusoutput('pkg-config fftw3 --exists')
	if status == 0:
		configCommand = 'pkg-config fftw3'
		outCFLAGS = os.popen('%s --cflags' % configCommand, 'r').readline().rstrip().split()
		outLIBS = os.popen('%s --libs' % configCommand, 'r').readline().rstrip().split()

	else:
		print "WARNING:  FFTW3 cannot be found, using defaults"
		outCFLAGS = []
		outLIBS = ['-lfftw3', '-lm']

	return outCFLAGS, outLIBS
		

class LSLDist(Distribution):
	"""Sub-class of setupuptools (distutils.core) Distribution class that fixes
	problems building the libnova extension."""

	def parse_config_files(self, filenames=None):
		# parse cfg file, but remove any pacakge specific options
		# otherwise distutils complains
		
		Distribution.parse_config_files(self, filenames)
		self.get_option_dict('build_ext').pop('libnova_prefix')
		try:
			self.get_option_dict('egg_info').pop('libnova_prefix')
		except:
			pass

# Get the FFTW flags/libs and manipulate the flags and libraries for 
# correlator._core appropriately.  This will, hopefully, fix the build
# problems on Mac
cflags, libs = get_fftw()
coreExtraFlags = ['-fopenmp',]
coreExtraFlags.extend(cflags)
coreExtraLibs = ['-fopenmp', '-lcblas', '-latlas']
coreExtraLibs.extend(libs)

setup(
	distclass = LSLDist, 
	name = "lsl", 
	version = get_version(), 
	description = "LWA Software Library", 
	author = "Jayce Dowell", 
	author_email = "jdowell@unm.edu", 
	url = "http://fornax.phys.unm.edu/lwa/trac/", 
	long_description = get_description('README'), 
	classifiers = ['Development Status :: 4 - Beta',
			'Intended Audience :: Science/Research',
			'Topic :: Scientific/Engineering :: Astronomy'],
	packages = find_packages(), 
	scripts = glob.glob('lsl/scripts/*.py'), 
	setup_requires = ['numpy>=1.2'], 
	install_requires = ['pyfits>=2.1', 'numpy>=1.2', 'scipy>=0.7', 'pyephem>=3.7.3', 'aipy>=0.9.1'], 
	dependency_links = ['http://www.stsci.edu/resources/software_hardware/pyfits/Download'], 
	include_package_data = True,  
	ext_package = 'lsl', 
	ext_modules = [Extension('_libnova', ['lsl/libnova.i']), 
				Extension('astro_array', ['lsl/astro_array.c'], include_dirs=[numpy.get_include()]),
				Extension('reader._gofast', ['lsl/reader/gofast.c'], include_dirs=[numpy.get_include()], extra_compile_args=['-funroll-loops']),
				Extension('correlator._core', ['lsl/correlator/core.c'], include_dirs=[numpy.get_include()], libraries=['m'], extra_compile_args=coreExtraFlags, extra_link_args=coreExtraLibs)], 
	zip_safe = False,  
	test_suite = "lsl.tests.test_lsl.lsl_tests"
) 
