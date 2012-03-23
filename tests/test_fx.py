# -*- coding: utf-8 -*-

"""Unit test for the lsl.correlator.fx module."""

import os
import time
import warnings
import unittest
import numpy

from lsl.common import stations
from lsl.correlator import fx


__version__  = "0.4"
__revision__ = "$ Revision: 6 $"
__author__    = "Jayce Dowell"

class SpecMaster_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx.SpecMaster
	function."""

	nAnt = 8

	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')
		
	def test_window(self):
		"""Test that window functions can be passed to SpecMaster."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMaster(fakeData, window=numpy.blackman)

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMaster(fakeData, window=numpy.hamming)
		
	def test_window_custom(self):
		"""Test that custom window functions can be passed to SpecMaster."""

		def wndw(L):
			return numpy.kaiser(L, 5)

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMaster(fakeData, window=wndw)
		
		def wndw2(L):
			return numpy.kaiser(L, 1)
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMaster(fakeData, window=wndw2)

	def test_spectra_real(self):
		"""Test the SpecMaster function on real-valued data."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMaster(fakeData)
		
	def test_spectra_complex(self):
		"""Test the SpecMaster function on complex-valued data."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMaster(fakeData, SampleRate=1e5, CentralFreq=38e6)


class SpecMasterP_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx.SpecMasterP
	function."""
	
	nAnt = 8

	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')
	
	def test_window_p(self):
		"""Test that window functions can be passed to SpecMasterP."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMasterP(fakeData, LFFT=128, window=numpy.blackman)

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMasterP(fakeData, LFFT=128, window=numpy.hamming)
		
	def test_window_custom_p(self):
		"""Test that custom window functions can be passed to SpecMasterP."""

		def wndw(L):
			return numpy.kaiser(L, 5)

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMasterP(fakeData, LFFT=128, window=wndw)
		
		def wndw2(L):
			return numpy.kaiser(L, 1)
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMasterP(fakeData, LFFT=128, window=wndw2)

	def test_spectra_real_p(self):
		"""Test the SpecMasterP function on real-valued data."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.SpecMasterP(fakeData, LFFT=128)
		
	def test_spectra_complex_p(self):
		"""Test the SpecMasterP function on complex-valued data."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.SpecMasterP(fakeData, LFFT=128, SampleRate=1e5, CentralFreq=38e6)


class StokesMaster_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx.StokesMaster
	function."""

	nAnt = 8

	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')
		
	def test_window(self):
		"""Test that window functions can be passed to StokesMaster."""
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=numpy.blackman)

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=numpy.hamming)
		
	def test_window_custom(self):
		"""Test that custom window functions can be passed to StokesMaster."""

		def wndw(L):
			return numpy.kaiser(L, 5)
			
		station = stations.lwa1
		antennas = station.getAntennas()

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=wndw)
		
		def wndw2(L):
			return numpy.kaiser(L, 1)
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=wndw2)

	def test_spectra_real(self):
		"""Test the StokesMaster function on real-valued data."""
		
		station = stations.lwa1
		antennas = station.getAntennas()

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt])
		
	def test_spectra_complex(self):
		"""Test the StokesMaster function on complex-valued data."""
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
		fakeData = fakeData.astype(numpy.csingle)
		freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6)


class FXMaster_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx.FXMaster
	function."""
	
	nAnt = 8

	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')

	def test_correlator_real(self):
		"""Test the C-based correlator on real-valued data."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt])

	def test_correlator_complex(self):
		"""Test the C-based correlator on complex-valued data."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6)
		
	def test_correlator_real_window(self):
		"""Test the C-based correlator on real-valued data window."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], 
							window=numpy.blackman)
		
	def test_correlator_complex_window(self):
		"""Test the C-based correlator on complex-valued data window."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
							window=numpy.blackman)

	def test_correlator_gaincorrect(self):
		"""Test appling gain correction to the correlator output."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
							GainCorrect=True)
							
	def test_correlator_baselines(self):
		"""Test that the ReturnBaselines keyword works."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
									ReturnBaselines=True)
									
	def test_correlator_pol(self):
		"""Test various correlator polarization settings."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		## XX
		blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
									ReturnBaselines=True, Pol='XX')
		for (ant1,ant2) in blList:
			self.assertEqual(ant1.pol, 0)
			self.assertEqual(ant2.pol, 0)
		
		## YY
		blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
									ReturnBaselines=True, Pol='YY')
		for (ant1,ant2) in blList:
			self.assertEqual(ant1.pol, 1)
			self.assertEqual(ant2.pol, 1)
		
		## XY
		blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
									ReturnBaselines=True, Pol='XY')
		for (ant1,ant2) in blList:
			self.assertEqual(ant1.pol, 0)
			self.assertEqual(ant2.pol, 1)
			
		## YX
		blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
									ReturnBaselines=True, Pol='YX')
		for (ant1,ant2) in blList:
			self.assertEqual(ant1.pol, 1)
			self.assertEqual(ant2.pol, 0)


class FXStokes_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.correlator.fx.FXStokes
	function."""
	
	nAnt = 8

	def setUp(self):
		"""Turn off all numpy and python warnings."""

		numpy.seterr(all='ignore')
		warnings.simplefilter('ignore')

	def test_correlator_real(self):
		"""Test the C-based correlator on real-valued data."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt])

	def test_correlator_complex(self):
		"""Test the C-based correlator on complex-valued data."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6)
		
	def test_correlator_real_window(self):
		"""Test the C-based correlator on real-valued data window."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 3.0
		fakeData = fakeData.astype(numpy.int16)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], 
							window=numpy.blackman)
		
	def test_correlator_complex_window(self):
		"""Test the C-based correlator on complex-valued data window."""

		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
							window=numpy.blackman)

	def test_correlator_gaincorrect(self):
		"""Test appling gain correction to the correlator output."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
							GainCorrect=True)
							
	def test_correlator_baselines(self):
		"""Test that the ReturnBaselines keyword works."""
		
		fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
		fakeData = fakeData.astype(numpy.csingle)
		
		station = stations.lwa1
		antennas = station.getAntennas()
		
		blList, freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], SampleRate=1e5, CentralFreq=38e6, 
									ReturnBaselines=True)


class fx_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.correlator.fx
	units tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(SpecMaster_tests))
		self.addTests(loader.loadTestsFromTestCase(SpecMasterP_tests))
		self.addTests(loader.loadTestsFromTestCase(StokesMaster_tests))
		self.addTests(loader.loadTestsFromTestCase(FXMaster_tests))
		self.addTests(loader.loadTestsFromTestCase(FXStokes_tests))


if __name__ == '__main__':
	unittest.main()
