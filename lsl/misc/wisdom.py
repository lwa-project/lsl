# -*- coding: utf-8 -*-
"""
Module for building and saving LSL-specific FFTW wisdom.

.. versionadded:: 1.0.1
"""

import os
import time

from lsl.common.paths import data as dataPath
from lsl.common.busy import BusyIndicator

import _wisdom


__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ["make", "show", "__version__", "__revision__", "__all__"]


# Path to the LSL-specific FFTW wisdom file
_wisdomFilename = os.path.join(dataPath, 'fftw_wisdom.txt')


def make():
	"""
	Build a new set of LSL-specific FFTW wisdom.
	"""
	
	
	bi = BusyIndicator(message="Building FFTW wisdom")
	bi.start()
	
	t0 = time.time()
	_wisdom.buildWisdom(_wisdomFilename)
	t1 = time.time()
	
	bi.stop()
	
	print "Finished in %.1f s" % (t1-t0)


def show():
	"""
	List information about the current LSL-specific FFTW wisdom.
	"""
	
	if not os.path.exists(_wisdomFilename):
		print "No LSL-specific wisdom file found, consider running make()"
		return False
		
	print "Listing '%s':" % _wisdomFilename
	
	fh = open(_wisdomFilename, 'r')
	for line in fh:
		print "  %s" % line.rstrip()
	fh.close()
	
	return True
