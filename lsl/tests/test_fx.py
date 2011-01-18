# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.fx module."""

import os
import time
import unittest
import numpy

from lsl.correlator import fx


__revision__ = "$ Revision: 3 $"
__version__  = "0.2"
__author__    = "Jayce Dowell"

class fx_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx
	module."""

	def setUp(self):
		"""Turn off all numpy warnings."""

		numpy.seterr(all='ignore')

	def test_window(self):
		"""Test that window functions can be passed to calcSpectra/calcSpectrum."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, window=numpy.blackman, DisablePool=True)

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, window=numpy.hamming, DisablePool=True)

	def test_window_custom(self):
		"""Test that custom window functions can be passed to calcSpectra/calcSpectrum."""

		def wndw(L):
			return numpy.kaiser(L, 5)

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, window=wndw, DisablePool=True)

	def test_spectra_single(self):
		"""Test the calcSpectra/calcSpectrum functions in single thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, DisablePool=True)

	def test_spectra_multi(self):
		"""Test the calcSpectra/calcSpectrum functions in multi-thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, spectra = fx.calcSpectra(fakeData, DisablePool=False)

	def test_correlator_old_single(self):
		"""Test the FXCorrelator/correlate functions in single thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, cps = fx.FXCorrelator(fakeData, numpy.array([1,2,3,4]), DisablePool=True)

	def test_correlator_old_multi(self):
		"""Test the FXCorrelator/correlate functions in multi-thread mode."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, cps = fx.FXCorrelator(fakeData, numpy.array([1,2,3,4]), DisablePool=False)

	def test_correlator_real(self):
		"""Test the C-based correlator on real-valued data."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		freq, cps = fx.FXMaster(fakeData, numpy.array([1,2,3,4]))

	def test_correlator_complex(self):
		"""Test the C-based correlator on complex-valued data."""

		fakeData = numpy.random.rand(4,1024) + 1j*numpy.random.rand(4,1024)
		freq, cps = fx.FXMaster(fakeData, numpy.array([1,2,3,4]))



class fx_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.correlator.fx
	units tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(fx_tests)) 


if __name__ == '__main__':
	unittest.main()
