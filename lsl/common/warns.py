# -*- coding: utf-8 -*-

"""
Module that contains the warning classes for the DRX, TBN, and TBW 
readers.  The two warnings defined here are:
  * warnDeprecated and
  * warnExperimental.

These are for warning users about deprecated functions and experimental 
features, respectively.
"""

import sys
import warnings as sysWarnings
import linecache

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['warnDeprecated', 'warnExperimental', '__version__', '__revision__', '__all__']


def __modifiedFormatWarning(message, category, filename, lineno, line=None):
	"""
	For some reason, the default warning.formatwarning function does not 
	support the 'line' argument and all warning calls throw their own warning
	calls.  This function is meant to fix that.
	"""

	s =  "%s:%s: %s: %s\n" % (filename, lineno, category.__name__, message)
	if line is None:
		line = linecache.getline(filename, lineno)
	if line:
		line = line.strip()
		s += "  %s\n" % line
	return s


def __modifiedShowWarning(message, category, filename, lineno, file=None, line=None):
	"""
	This is a new warning.showwarning function to address the same problem as 
	seen with the warning.formatwarning default function.
	"""

	s =  __modifiedFormatWarning(message, category, filename, lineno, line=line)
	
	if file is None:
		sys.stderr.write(s)
	else:
		file.write(s)


def warnDeprecated(name, memo=None):
	"""
	Wrapper around warnins.warn to fix the problems with 
	[(format)(show)]warnings functions described above.
	"""

	sysWarnings.showwarning = __modifiedShowWarning
	sysWarnings.formatwarning = __modifiedFormatWarning
	
	if memo is None:
		sysWarnings.warn("'%s' is deprecated" % name, DeprecationWarning, 3)
	else:
		sysWarnings.warn("'%s' is deprecated: %s" % (name, memo), DeprecationWarning, 3)


def warnExperimental(name, memo=None):
	"""
	Wrapper around warnins.warn to fix the problems with 
	[(format)(show)]warnings functions described above.
	"""

	sysWarnings.showwarning = __modifiedShowWarning
	sysWarnings.formatwarning = __modifiedFormatWarning
	
	if memo is None:
		sysWarnings.warn("'%s' is experimental" % name, ImportWarning, 3)
	else:
		sysWarnings.warn("'%s' is experimental: %s" % (name, memo), ImportWarning, 3)
