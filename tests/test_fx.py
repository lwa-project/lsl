# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.fx module."""

import os
import time
import warnings
import unittest
import numpy

from lsl.common import stations
from lsl.correlator import fx


__revision__ = "$ Revision: 4 $"
__version__  = "0.3"
__author__    = "Jayce Dowell"

class fx_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx
	module."""

	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')
	
	### SpecMaster Function ###
	
	def test_window(self):
		"""Test that window functions can be passed to SpecMaster."""
		
		fakeData = numpy.random.rand(4,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMaster(fakeData, window=numpy.blackman)

		fakeData = numpy.random.rand(4,1024) + 1j*numpy.random.rand(4,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMaster(fakeData, window=numpy.hamming)
		
	def test_window_custom(self):
		"""Test that custom window functions can be passed to SpecMaster."""

		def wndw(L):
			return numpy.kaiser(L, 5)

		fakeData = numpy.random.rand(4,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMaster(fakeData, window=wndw)
		
		def wndw2(L):
			return numpy.kaiser(L, 1)
		
		fakeData = numpy.random.rand(4,1024) + 1j*numpy.random.rand(4,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMaster(fakeData, window=wndw2)

	def test_spectra_real(self):
		"""Test the SpecMaster function on real-valued data."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMaster(fakeData)
		
	def test_spectra_complex(self):
		"""Test the SpecMaster function on complex-valued data."""
		
		fakeData = numpy.random.rand(4,1024) + 1j*numpy.random.rand(4,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMaster(fakeData, SampleRate=1e5, CentralFreq=38e6)

	### FXMaster Function ###

	def test_correlator_real(self):
		"""Test the C-based correlator on real-valued data."""

		fakeData = numpy.random.rand(4,1024) + 3.0
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:4])

	def test_correlator_complex(self):
		"""Test the C-based correlator on complex-valued data."""

		fakeData = numpy.random.rand(4,1024) + 1j*numpy.random.rand(4,1024)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:4], SampleRate=1e5, CentralFreq=38e6)
		
	def test_correlator_real_window(self):
		"""Test the C-based correlator on real-valued data window."""
		
		fakeData = numpy.random.rand(4,1024) + 3.0
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:4], 
							window=numpy.blackman)
		
	def test_correlator_complex_window(self):
		"""Test the C-based correlator on complex-valued data window."""

		fakeData = numpy.random.rand(4,1024) + 1j*numpy.random.rand(4,1024)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:4], SampleRate=1e5, CentralFreq=38e6, 
							window=numpy.blackman)
	
	#### PXMaster Function ###
	
	#def test_fb_correlator_correlator_real(self):
		#"""Test the C-based correlator on real-valued data."""

		#fakeData = numpy.random.rand(4,4096) + 3.0
		#freq, cps = fx.PXMaster(fakeData, numpy.array([1,2,3,4]))

	#def test_fb_correlator_complex(self):
		#"""Test the C-based correlator on complex-valued data."""

		#fakeData = numpy.random.rand(4,4096) + 1j*numpy.random.rand(4,4096)
		#freq, cps = fx.PXMaster(fakeData, numpy.array([1,2,3,4]), SampleRate=1e5, CentralFreq=38e6)
		
	#def test_fb_correlator_real_window(self):
		#"""Test the C-based correlator on real-valued data window."""
		
		#fakeData = numpy.random.rand(4,4096) + 3.0
		#freq, cps = fx.PXMaster(fakeData, numpy.array([1,2,3,4]), 
							#window=numpy.blackman)
		
	#def test_fb_correlator_complex_window(self):
		#"""Test the C-based correlator on complex-valued data."""

		#fakeData = numpy.random.rand(4,4096) + 1j*numpy.random.rand(4,4096)
		#freq, cps = fx.PXMaster(fakeData, numpy.array([1,2,3,4]), SampleRate=1e5, CentralFreq=38e6, 
							#window=numpy.blackman)


class fx_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.correlator.fx
	units tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(fx_tests)) 


if __name__ == '__main__':
	unittest.main()
