# -*- coding: utf-8 -*-

"""Unit test for lsl.correlator.filterbank module."""

import os
import numpy
import unittest

from lsl.correlator import filterbank


__revision__ = "$ Revision: 2 $"
__version__  = "0.2"
__author__    = "Jayce Dowell"


class filterbank_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.filterbank
	module."""
	
	def test_filterbank2(self):
		"""Test that the 2-tap filter band works"""
		
		data = numpy.random.rand(4096)
		filterbank.fft2(data, 256)
		
	def test_filterbank4(self):
		"""Test that the 4-tap filter band works"""
		
		data = numpy.random.rand(4096)
		filterbank.fft4(data, 256)
		
	def test_filterbank8(self):
		"""Test that the 8-tap filter band works"""
		
		data = numpy.random.rand(4096)
		filterbank.fft8(data, 256)
		
	def test_filterbank16(self):
		"""Test that the 16-tap filter band works"""
		
		data = numpy.random.rand(4096)
		filterbank.fft16(data, 256)
		
	def test_filterbankN(self):
		"""Test that a N-tap filter bank works"""
		
		data = numpy.random.rand(4096)
		filterbank.fft(data, 256, P=2)
		
	def test_window(self):
		"""Test that a window function can be applied to the data"""
		
		data = numpy.random.rand(4096)
		
		# N-taps
		filterbank.fft(data, 256, P=2, window=numpy.blackman)
		
		# 4-taps
		filterbank.fft4(data, 256,window=numpy.bartlett)
		
		# 16-taps
		filterbank.fft16(data, 256,window=numpy.hanning)
		
	def test_filterbank2_complex64(self):
		"""Test the 2-tap filter band with complex64 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex64)
		filterbank.fft2(data, 256)
		
	def test_filterbank2_complex128(self):
		"""Test the 2-tap filter band with complex128 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex128)
		filterbank.fft2(data, 256)
		
	def test_filterbank4_complex64(self):
		"""Test the 4-tap filter band with complex64 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex64)
		filterbank.fft4(data, 256)
		
	def test_filterbank4_complex128(self):
		"""Test the 4-tap filter band with complex128 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex128)
		filterbank.fft4(data, 256)
		
	def test_filterbank8_complex64(self):
		"""Test the 8-tap filter band with complex64 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex64)
		filterbank.fft8(data, 256)
		
	def test_filterbank8_complex128(self):
		"""Test the 8-tap filter band with complex128 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex128)
		filterbank.fft8(data, 256)
		
	def test_filterbank16_complex64(self):
		"""Test the 16-tap filter band with complex64 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex64)
		filterbank.fft16(data, 256)
		
	def test_filterbank16_complex128(self):
		"""Test the 16-tap filter band with complex128 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex128)
		filterbank.fft16(data, 256)
		
	def test_window_complex64(self):
		"""Test that a window function can be applied to complex64 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex64)
		
		# N-taps
		filterbank.fft(data, 256, P=2, window=numpy.blackman)
		
		# 4-taps
		filterbank.fft4(data, 256,window=numpy.bartlett)
		
		# 16-taps
		filterbank.fft16(data, 256,window=numpy.hanning)
		
	def test_window_complex128(self):
		"""Test that a window function can be applied to complex128 data"""
		
		data = numpy.random.rand(4096)
		data = data.astype(numpy.complex128)
		
		# N-taps
		filterbank.fft(data, 256, P=2, window=numpy.blackman)
		
		# 4-taps
		filterbank.fft4(data, 256,window=numpy.bartlett)
		
		# 16-taps
		filterbank.fft16(data, 256,window=numpy.hanning)


class filterbank_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.misc.geodesy
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(filterbank_tests)) 


if __name__ == '__main__':
	unittest.main()

