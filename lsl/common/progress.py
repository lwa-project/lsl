# -*- coding: utf-8 -*-

"""
Module to make an ASCII progress bar.
"""

import copy

__version__ = '0.1'
__revision__ = '$ Revision: 3 $'
__all__ = ['ProgressBar', '__version__', '__revision__', '__all__']


class ProgressBar(object):
	"""
	Object to make a ASCII progress bar for use with various long-
	run tasks.

	Example Usage:
		>>> import sys
		>>> from progess import ProgressBar
		>>> pb = ProgressBar()
		>>> pb.inc()
		>>> sys.stdout.write(pb.show())
		>>> sys.stdout.flush()
	"""

	def __init__(self, max=100, span=70, sym='=', printP=True):
		"""
		Initialize the ProgressBar class with various parameters:
		  * max: maximum count for the progress bar (default: 100)
		  * span: width in characters to make the bar (default: 70)
		  * sym: character to use in the progress bar (default: '=')
		  * printP: whether or not to print the percentage in addition to
		    the bar or not (default: True)
		"""

		self.amount = 0
		self.max = max
		self.span = span
		self.sym = sym
		self.printP = printP
	
	def inc(self, amount=1):
		"""
		Increment the progress bar's internal counter by some amount.  The
		default is one.
		"""

		self.__iadd__(amount)

	def dec(self, amount=1):
		"""
		Decrement the progress bar's internal counter by some amount.  The
		default is one.
		"""
			
		self.__isub__(amount)

	def show(self):
		"""
		Build a string representation of the progress bar and return it.
		"""

		if self.printP:
			# If we want the percentage also displayed, trim a little 
			# more from the progress bar's wdith
			barSpan = self.span - 9
			nMarks = int(round(float(self.amount)/self.max * barSpan))
			bar = self.sym * nMarks
			bar = bar+(' ' * (barSpan-nMarks))
			nte = "%5.1f%%" % (float(self.amount)/self.max*100)

			out = "|%s| %s" % (bar, nte)
		else:
			# Progress bar only
			barSpan = self.span - 2
			nMarks = int(round(float(self.amount)/self.max * barSpan))
			bar = self.sym * nMarks
			bar = bar+(' ' * (barSpan-nMarks))

			out = "|%s|" % bar

		return out
		
	def __add__(self, amount):
		"""
		Increment the internal counter by a certain amount, return a new
		ProgressBar object.
		"""
		
		newBar = copy.deepcopy(self)
		newBar += amount
		return newBar

	def __iadd__(self, amount):
		"""
		Increment the internal counter by a certain amount.
		"""

		self.amount += amount
		return self
		
	def __sub__(self, amount):
		"""
		Decrement the internal counter by a certain amount, return a new
		ProgressBar object.
		"""
		
		newBar = copy.deepcopy(self)
		if newBar.amount >= amount:
			newBar += amount
		return newBar

	def __isub__(self, amount):
		"""
		Decrement the internal counter by a certain amount.
		"""

		if self.amount >= amount:
			self.amount -= amount
		return self

	def __str__(self):
		"""
		Alternative to self.show().
		"""

		return self.show()
