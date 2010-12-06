# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.fx module."""

import os
import time
import unittest
import numpy

from lsl.correlator import fx


__revision__ = "$Revision:1 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"

class fx_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx
	module."""

	def test_blackman(self):
		"""Test the Blackman window function."""

		window = {1: [1], 
				2: [-1.3878e-17, -1.3878e-17], 
				4: [-1.3878e-17, 6.3000e-1, 6.3000e-1, -1.3878e-17], 
				8: [-1.3878e-17, 9.0453e-2, 4.5918e-1, 9.2036e-1, 9.2036e-1, 4.5918e-1, 9.0453e-2, -1.3878e-17]}

		for length in window.keys():
			ref = numpy.array(window[length])
			out = fx.blackmanWindow(length)
			
			i = 0
			for r,o in zip(ref, out):
				self.assertAlmostEqual(r, o, 4, "Error on window length %i, value %i" % (length,i))
				i = i + 1

	def test_sinc(self):
		"""Test the sinc window function."""

		window = {1: [1], 
				2: [-3.8980e-17, 1.0], 
				4: [-3.8980e-17, 3.8980e-17, 1.0, 3.8980e-17]}

		for length in window.keys():
			ref = numpy.array(window[length])
			out = fx.sincWindow(length)

			i = 0
			for r,o in zip(ref, out):
				self.assertAlmostEqual(r, o, 5, "Error on window length %i, value %i" % (length,i))
				i = i + 1

	def test_spectra_single(self):
		"""Test the calcSpectra/calcSpectrum functions in single thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, DisablePool=True)

	def test_spectra_multi(self):
		"""Test the calcSpectra/calcSpectrum functions in multi-thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, DisablePool=False)

	def test_correlator_single(self):
		"""Test the FXCorrelaotr/correlate functions in single thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, cps = fx.FXCorrelator(fakeData, numpy.array([1,2,3,4]), DisablePool=True)

	def test_correlator_multi(self):
		"""Test the FXCorrelaotr/correlate functions in multi-thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, cps = fx.FXCorrelator(fakeData, numpy.array([1,2,3,4]), DisablePool=False)



class fx_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.correlator.fx
	units tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(fx_tests)) 


if __name__ == '__main__':
	unittest.main()
