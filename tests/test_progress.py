# -*- coding: utf-8 -*-

"""Unit test for regressions in the lsl.common.progress module."""

import os
import numpy
import unittest

from lsl.common import progress

__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


class progress_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the regressions in LSL."""
	
	def test_bar_default(self):
		"""Test the progress bar, default options."""
		
		pbar = progress.ProgressBar()
		for i in range(101):
			pbar.inc(2)
			pbar.dec(1)
		
	def test_bar_attributes(self):
		"""Test the progress bar's attributes."""
		
		pbar2 = progress.ProgressBar()
		for i in range(101):
			pbar2 += 2
			pbar2 -= 1
			
			pbar2 = pbar2 + 1
			pbar2 = pbar2 - 1
			
			self.assertEqual(pbar2.amount, i+1)
			
	def test_bar_show(self):
		"""Test the progress bar's rendering."""
		
		# With percentage
		pbar = progress.ProgressBar()
		for i in range(101):
			pbar.inc(1)
			pbar.show()
			
		# Without percentage
		pbar = progress.ProgressBar(printP=False)
		for i in range(101):
			pbar.inc(1)
			pbar.show()


class progress_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.progress unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(progress_tests)) 


if __name__ == '__main__':
	unittest.main()
