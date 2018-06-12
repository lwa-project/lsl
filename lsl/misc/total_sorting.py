# -*- coding: utf-8 -*-

"""
Module to add all six sorting comparison methods to a class that only defines
the  __cmp__() method.

This module is based on the total_ordering decorator from functools at:  
	https://github.com/python/cpython/blob/master/Lib/functools.py
	
.. versionadded:: 1.2.1
"""

__version__   = '0.1'
__revision__ = '$Rev$'
__all__ = ['cmp_to_total', '__version__', '__revision__', '__all__']


def _cmp_to_lt(self, other):
	"""
	Return a < b.  Compute by @cmp_to_total from __cmp__
	"""
	
	return True if self.__cmp__(other) < 0 else False


def _cmp_to_le(self, other):
	"""
	Return a <= b.  Compute by @cmp_to_total from __cmp__
	"""
	
	return True if self.__cmp__(other) <= 0 else False


def _cmp_to_gt(self, other):
	"""
	Return a > b.  Compute by @cmp_to_total from __cmp__
	"""
	
	return True if self.__cmp__(other) > 0 else False


def _cmp_to_ge(self, other):
	"""
	Return a >= b.  Compute by @cmp_to_total from __cmp__
	"""
	
	return True if self.__cmp__(other) >= 0 else False


def _cmp_to_eq(self, other):
	"""
	Return a == b.  Compute by @cmp_to_total from __cmp__
	"""
	
	return True if self.__cmp__(other) == 0 else False


def _cmp_to_ne(self, other):
	"""
	Return a != b.  Compute by @cmp_to_total from __cmp__
	"""
	
	return True if self.__cmp__(other) != 0 else False


def cmp_to_total(cls):
	"""
	Decorator to define the six comparison operators that are needed
	for total ordering that can be used in all versions of Python 
	from the __cmp__() method.
	"""
	
	names = ['__lt__', '__le__', '__gt__', '__ge__', '__eq__', '__ne__']
	funcs = [_cmp_to_lt, _cmp_to_le, _cmp_to_gt, _cmp_to_ge, _cmp_to_eq, _cmp_to_ne]
	
	for name,func in zip(names, funcs):
		# Is it defined?
		if name not in dir(cls):
			func.__name__ = name
			setattr(cls, name, func)
			
	return cls
