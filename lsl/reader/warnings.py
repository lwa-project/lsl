# -*- coding: utf-8 -*-

"""Module that contains the warning classes for the DRX, TBN, and TBW readers."""

import sys
import warnings
import linecache


__version__ = '0.1'
__revision__ = '$ Revision: 1 $'
__all__ = ['warnDeprecated', '__version__', '__revision__', '__all__']


def __modifiedFormatWarning(message, category, filename, lineno, line=None):
	"""For some reason, the default warning.formatwarning function does not 
	support the 'line' argument and all warning calls throw their own warning
	calls.  This function is meant to fix that."""

	s =  "%s:%s: %s: %s\n" % (filename, lineno, category.__name__, message)
	if line is None:
		line = linecache.getline(filename, lineno)
	if line:
		line = line.strip()
		s += "  %s\n" % line
	return s


def __modifiedShowWarning(message, category, filename, lineno, file=None, line=None):
	"""This is a new warning.showwarning function to address the same problem as 
	seen with the warning.formatwarning default function."""

	s =  __modifiedFormatWarning(message, category, filename, lineno, line=line)
	
	if file is None:
		sys.stderr.write(s)
	else:
		file.write(s)


def warnDeprecated(name, memo=None):
	"""Wrapper around warnins.warn to fix the problems with 
	[(format)(show)]warnings functions described above."""

	warnings.showwarning = __modifiedShowWarning
	warnings.formatwarning = __modifiedFormatWarning
	
	if memo is None:
		warnings.warn("'%s' is deprecated" % name, DeprecationWarning, 3)
	else:
		warnings.warn("'%s' is deprecated: %s" % (name, memo), DeprecationWarning, 3)
	
	