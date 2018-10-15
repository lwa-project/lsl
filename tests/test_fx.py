# -*- coding: utf-8 -*-

# Python3 compatibility
from __future__ import division
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for the lsl.correlator.fx module."""

import os
import time
import warnings
import unittest
import numpy

from lsl.common.paths import DATA_BUILD
from lsl.common import stations
from lsl.correlator import fx

_SSMIF = os.path.join(DATA_BUILD, 'lwa1-ssmif.txt')

__version__  = "0.5"
__revision__ = "$Rev$"
__author__    = "Jayce Dowell"

class SpecMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.SpecMaster
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def test_window(self):
        """Test that window functions can be passed to SpecMaster."""
        
        #
        # Real
        #
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, window=numpy.blackman)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            wndw = numpy.blackman(2*LFFT)
            for i in xrange(self.nAnt):
                for j in xrange(nFFT):
                    spectra2[i,:] += (numpy.abs( numpy.fft.fft(fakeData[i,j*2*LFFT:(j+1)*2*LFFT]*wndw) )**2)[:LFFT]
            spectra2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
        #
        # Complex
        #
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, window=numpy.hamming)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            ewndw = numpy.hamming(LFFT)
            for i in xrange(self.nAnt):
                for j in xrange(nFFT):
                    spectra2[i,:] += numpy.fft.fftshift( numpy.abs( numpy.fft.fft(fakeData[i,j*LFFT:(j+1)*LFFT]*ewndw) )**2 )
            spectra2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_window_custom(self):
        """Test that custom window functions can be passed to SpecMaster."""
        
        #
        # Real
        #
        def wndw(L):
            return numpy.kaiser(L, 5)
            
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, window=wndw)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            ewndw = wndw(2*LFFT)
            for i in xrange(self.nAnt):
                for j in xrange(nFFT):
                    spectra2[i,:] += (numpy.abs( numpy.fft.fft(fakeData[i,j*2*LFFT:(j+1)*2*LFFT]*ewndw) )**2)[:LFFT]
            spectra2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
        #
        # Complex
        #
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, window=wndw2)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            ewndw = wndw2(LFFT)
            for i in xrange(self.nAnt):
                for j in xrange(nFFT):
                    spectra2[i,:] += numpy.fft.fftshift( numpy.abs( numpy.fft.fft(fakeData[i,j*LFFT:(j+1)*LFFT]*ewndw) )**2 )
            spectra2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_spectra_real(self):
        """Test the SpecMaster function on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            for i in xrange(self.nAnt):
                for j in xrange(nFFT):
                    spectra2[i,:] += (numpy.abs( numpy.fft.fft(fakeData[i,j*2*LFFT:(j+1)*2*LFFT]) )**2)[:LFFT]
            spectra2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_spectra_complex(self):
        """Test the SpecMaster function on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, sample_rate=1e5, central_freq=38e6)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            for i in xrange(self.nAnt):
                for j in xrange(nFFT):
                    spectra2[i,:] += numpy.fft.fftshift( numpy.abs( numpy.fft.fft(fakeData[i,j*LFFT:(j+1)*LFFT]) )**2 )
            spectra2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_spectra_odd_complex(self):
        """Test the SpecMaster function on odd-sized complex transforms."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,10) + 3.0
            fakeData = fakeData + 0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, LFFT=9, sample_rate=1e5, central_freq=38e6)
            
            for i in xrange(spectra.shape[0]):
                self.assertTrue(numpy.abs(spectra[i,0]-spectra[i,-1]) < 1e-6*spectra[i,:].max())
                
            def wndw2(L):
                return numpy.kaiser(L, 1)
                
            freq, spectra = fx.SpecMaster(fakeData, LFFT=9, sample_rate=1e5, central_freq=38e6, window=wndw2)
            
            for i in xrange(spectra.shape[0]):
                self.assertTrue(numpy.abs(spectra[i,0]-spectra[i,-1]) < 1e-6*spectra[i,:].max())


class StokesMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.StokesMaster
    function."""

    nAnt = 8

    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def test_window(self):
        """Test that window functions can be passed to StokesMaster."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        #
        # Real
        #
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=numpy.blackman)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            wndw = numpy.blackman(2*LFFT)
            for i in xrange(self.nAnt//2):
                for j in xrange(nFFT):
                    xF = numpy.fft.fft(fakeData[2*i+0,j*2*LFFT:(j+1)*2*LFFT]*wndw)[:LFFT]
                    yF = numpy.fft.fft(fakeData[2*i+1,j*2*LFFT:(j+1)*2*LFFT]*wndw)[:LFFT]
                    
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
        #
        # Complex
        #
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=numpy.hamming)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            wndw = numpy.hamming(LFFT)
            for i in xrange(self.nAnt//2):
                for j in xrange(nFFT):
                    xF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+0,j*LFFT:(j+1)*LFFT]*wndw) )
                    yF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+1,j*LFFT:(j+1)*LFFT]*wndw) )
                    
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_window_custom(self):
        """Test that custom window functions can be passed to StokesMaster."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        #
        # Real
        #
        
        def wndw(L):
            return numpy.kaiser(L, 5)
            
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=wndw)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            ewndw = wndw(2*LFFT)
            for i in xrange(self.nAnt//2):
                for j in xrange(nFFT):
                    xF = numpy.fft.fft(fakeData[2*i+0,j*2*LFFT:(j+1)*2*LFFT]*ewndw)[:LFFT]
                    yF = numpy.fft.fft(fakeData[2*i+1,j*2*LFFT:(j+1)*2*LFFT]*ewndw)[:LFFT]
                    
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
        #
        # Complex
        #
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], window=wndw2)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            ewndw = wndw2(LFFT)
            for i in xrange(self.nAnt//2):
                for j in xrange(nFFT):
                    xF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+0,j*LFFT:(j+1)*LFFT]*ewndw) )
                    yF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+1,j*LFFT:(j+1)*LFFT]*ewndw) )
                    
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_spectra_real(self):
        """Test the StokesMaster function on real-valued data."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt])
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            for i in xrange(self.nAnt//2):
                for j in xrange(nFFT):
                    xF = numpy.fft.fft(fakeData[2*i+0,j*2*LFFT:(j+1)*2*LFFT])[:LFFT]
                    yF = numpy.fft.fft(fakeData[2*i+1,j*2*LFFT:(j+1)*2*LFFT])[:LFFT]
                    
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_spectra_complex(self):
        """Test the StokesMaster function on complex-valued data."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6)
            
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            for i in xrange(self.nAnt//2):
                for j in xrange(nFFT):
                    xF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+0,j*LFFT:(j+1)*LFFT]) )
                    yF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+1,j*LFFT:(j+1)*LFFT]) )
                    
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(spectra-spectra2).max() < 1e-6*spectra2.max())
            
    def test_spectra_odd_complex(self):
        """Test the SpecMaster function on odd-sized complex transforms."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,10) + 3.0
            fakeData = fakeData + 0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], LFFT=9, sample_rate=1e5, central_freq=38e6)
            
            for i in xrange(spectra.shape[0]):
                self.assertTrue(numpy.abs(spectra[0,i,0]-spectra[0,i,-1]) < 1e-6*spectra[0,i,:].max())
                
            def wndw2(L):
                return numpy.kaiser(L, 1)
                
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], LFFT=9, sample_rate=1e5, central_freq=38e6, window=wndw2)
            
            for i in xrange(spectra.shape[0]):
                self.assertTrue(numpy.abs(spectra[0,i,0]-spectra[0,i,-1]) < 1e-6*spectra[0,i,:].max())


class FXMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.FXMaster
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def test_correlator_real(self):
        """Test the C-based correlator on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt])
            
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt])
            
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            blc = 0
            for i in xrange(0, self.nAnt):
                if antennas[i].pol != 0:
                    continue
                for j in xrange(i+1, self.nAnt):
                    if antennas[j].pol != 0:
                        continue
                        
                    for k in xrange(nFFT):
                        f1 = numpy.fft.fft(fakeData[i,k*2*LFFT:(k+1)*2*LFFT])[:LFFT]
                        f2 = numpy.fft.fft(fakeData[j,k*2*LFFT:(k+1)*2*LFFT])[:LFFT]
                        
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps[:,1:]-cps2[:,1:])).max() < 1e-6*numpy.abs(numpy.abs(cps2[:,1:])).max())
            
    def test_correlator_complex(self):
        """Test the C-based correlator on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6)
            
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6)
            
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            blc = 0
            for i in xrange(0, self.nAnt):
                if antennas[i].pol != 0:
                    continue
                for j in xrange(i+1, self.nAnt):
                    if antennas[j].pol != 0:
                        continue
                        
                    for k in xrange(nFFT):
                        f1 = numpy.fft.fftshift( numpy.fft.fft(fakeData[i,k*LFFT:(k+1)*LFFT]) )
                        f2 = numpy.fft.fftshift( numpy.fft.fft(fakeData[j,k*LFFT:(k+1)*LFFT]) )
                        
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps-cps2)).max() < 1e-6*numpy.abs(cps2).max())
            
    def test_correlator_real_window(self):
        """Test the C-based correlator on real-valued data window."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], 
                                window=numpy.blackman)
                                
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], 
                                window=numpy.blackman)
                                
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            wndw = numpy.blackman(2*LFFT)
            blc = 0
            for i in xrange(0, self.nAnt):
                if antennas[i].pol != 0:
                    continue
                for j in xrange(i+1, self.nAnt):
                    if antennas[j].pol != 0:
                        continue
                        
                    for k in xrange(nFFT):
                        f1 = numpy.fft.fft(fakeData[i,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                        f2 = numpy.fft.fft(fakeData[j,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                        
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps[:,1:]-cps2[:,1:])).max() < 1e-6*numpy.abs(numpy.abs(cps2[:,1:])).max())
            
    def test_correlator_complex_window(self):
        """Test the C-based correlator on complex-valued data window."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                window=numpy.blackman)
                                
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                window=numpy.blackman)
                                
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            ewndw = numpy.blackman(LFFT)
            blc = 0
            for i in xrange(0, self.nAnt):
                if antennas[i].pol != 0:
                    continue
                for j in xrange(i+1, self.nAnt):
                    if antennas[j].pol != 0:
                        continue
                        
                    for k in xrange(nFFT):
                        f1 = numpy.fft.fftshift( numpy.fft.fft(fakeData[i,k*LFFT:(k+1)*LFFT]*ewndw) )
                        f2 = numpy.fft.fftshift( numpy.fft.fft(fakeData[j,k*LFFT:(k+1)*LFFT]*ewndw) )
                        
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps-cps2)).max() < 1e-6*numpy.abs(cps2).max())
            
    def test_correlator_gaincorrect(self):
        """Test appling gain correction to the correlator output."""
        
        fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(numpy.csingle)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                            gain_correct=True)
                            
    def test_correlator_baselines(self):
        """Test that the return_baselines keyword works."""
        
        fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(numpy.csingle)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True)
                                    
    def test_correlator_pol(self):
        """Test various correlator polarization settings."""
        
        fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(numpy.csingle)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        ## XX
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, Pol='XX')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 0)
            self.assertEqual(ant2.pol, 0)
        
        ## YY
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, Pol='YY')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 1)
            self.assertEqual(ant2.pol, 1)
        
        ## XY
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, Pol='XY')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 0)
            self.assertEqual(ant2.pol, 1)
            
        ## YX
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, Pol='YX')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 1)
            self.assertEqual(ant2.pol, 0)
            
    def test_correlator_odd_complex(self):
        """Test the FXMaster function on odd-sized complex transforms."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,10) + 3.0
            fakeData = fakeData + 0j
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            for ant in antennas:
                ant.stand.x = 0.0
                ant.stand.y = 0.0
                ant.stand.z = 0.0
                ant.cable.length = 0.0
                
            blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=9, sample_rate=1e5, central_freq=38e6, 
                                        return_baselines=True, Pol='XX')
            
            for i in xrange(cps.shape[0]):
                self.assertTrue(numpy.abs(cps[i,0]-cps[i,-1].conj()) < 1e-6*numpy.abs(cps[i,:]).max())
                
            def wndw2(L):
                return numpy.kaiser(L, 1)
                
            blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=9, sample_rate=1e5, central_freq=38e6, 
                                        return_baselines=True, Pol='XX', window=wndw2)
            
            for i in xrange(cps.shape[0]):
                self.assertTrue(numpy.abs(cps[i,0]-cps[i,-1].conj()) < 1e-6*numpy.abs(cps[i,:]).max())


class FXStokes_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.FXStokes
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def test_correlator_real(self):
        """Test the C-based correlator on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt])
            
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt])
            
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            blc = 0
            for i in xrange(0, self.nAnt//2):
                for j in xrange(i+1, self.nAnt//2):
                    for k in xrange(nFFT):
                        f1X = numpy.fft.fft(fakeData[2*i+0,k*2*LFFT:(k+1)*2*LFFT])[:LFFT]
                        f1Y = numpy.fft.fft(fakeData[2*i+1,k*2*LFFT:(k+1)*2*LFFT])[:LFFT]
                        f2X = numpy.fft.fft(fakeData[2*j+0,k*2*LFFT:(k+1)*2*LFFT])[:LFFT]
                        f2Y = numpy.fft.fft(fakeData[2*j+1,k*2*LFFT:(k+1)*2*LFFT])[:LFFT]
                        
                        cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                        cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                        cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                        cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                    blc += 1
            cps2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps[:,:,1:]-cps2[:,:,1:])).max() < 1e-6*numpy.abs(cps2[:,:,1:]).max())
            
    def test_correlator_complex(self):
        """Test the C-based correlator on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6)
            
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6)
            
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            blc = 0
            for i in xrange(0, self.nAnt//2):
                for j in xrange(i+1, self.nAnt//2):
                    for k in xrange(nFFT):
                        f1X = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+0,k*LFFT:(k+1)*LFFT]) )
                        f1Y = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+1,k*LFFT:(k+1)*LFFT]) )
                        f2X = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*j+0,k*LFFT:(k+1)*LFFT]) )
                        f2Y = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*j+1,k*LFFT:(k+1)*LFFT]) )
                        
                        cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                        cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                        cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                        cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                    blc += 1
            cps2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps-cps2)).max() < 1e-6*numpy.abs(cps2).max())
            
    def test_correlator_real_window(self):
        """Test the C-based correlator on real-valued data window."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024) + 3.0
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], 
                                window=numpy.blackman)
                                
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt],
                                window=numpy.blackman)
                                
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            wndw = numpy.blackman(2*LFFT)
            blc = 0
            for i in xrange(0, self.nAnt//2):
                for j in xrange(i+1, self.nAnt//2):
                    for k in xrange(nFFT):
                        f1X = numpy.fft.fft(fakeData[2*i+0,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                        f1Y = numpy.fft.fft(fakeData[2*i+1,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                        f2X = numpy.fft.fft(fakeData[2*j+0,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                        f2Y = numpy.fft.fft(fakeData[2*j+1,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                        
                        cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                        cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                        cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                        cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                    blc += 1
            cps2 /= (2*LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps[:,:,1:]-cps2[:,:,1:])).max() < 1e-6*numpy.abs(cps2[:,:,1:]).max())
            
    def test_correlator_complex_window(self):
        """Test the C-based correlator on complex-valued data window."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                window=numpy.blackman)
                                
            # Numpy comparison
            for i in xrange(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
                
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                window=numpy.blackman)
                                
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            wndw = numpy.blackman(LFFT)
            blc = 0
            for i in xrange(0, self.nAnt//2):
                for j in xrange(i+1, self.nAnt//2):
                    for k in xrange(nFFT):
                        f1X = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+0,k*LFFT:(k+1)*LFFT]*wndw) )
                        f1Y = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+1,k*LFFT:(k+1)*LFFT]*wndw) )
                        f2X = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*j+0,k*LFFT:(k+1)*LFFT]*wndw) )
                        f2Y = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*j+1,k*LFFT:(k+1)*LFFT]*wndw) )
                        
                        cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                        cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                        cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                        cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                    blc += 1
            cps2 /= (LFFT * nFFT)
            self.assertTrue(numpy.abs(numpy.abs(cps-cps2)).max() < 1e-6*numpy.abs(cps2).max())
            
    def test_correlator_gaincorrect(self):
        """Test appling gain correction to the correlator output."""
        
        fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(numpy.csingle)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                            gain_correct=True)
                            
    def test_correlator_baselines(self):
        """Test that the return_baselines keyword works."""
        
        fakeData = numpy.random.rand(self.nAnt,1024) + 1j*numpy.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(numpy.csingle)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        blList, freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True)
        
    def test_correlator_odd_complex(self):
        """Test the FXStokes function on odd-sized complex transforms."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,10) + 3.0
            fakeData = fakeData + 0j
            fakeData = fakeData.astype(dtype)
            
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
            for ant in antennas:
                ant.stand.x = 0.0
                ant.stand.y = 0.0
                ant.stand.z = 0.0
                ant.cable.length = 0.0
                    
            blList, freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=9, sample_rate=1e5, central_freq=38e6, 
                                        return_baselines=True)
            
            for i in xrange(cps.shape[0]):
                self.assertTrue(numpy.abs(cps[0,i,0]-cps[0,i,-1].conj()) < 1e-6*cps[0,i,:].max())
                
            def wndw2(L):
                return numpy.kaiser(L, 1)
                
            blList, freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=9, sample_rate=1e5, central_freq=38e6, 
                                        return_baselines=True, window=wndw2)
            
            for i in xrange(cps.shape[0]):
                self.assertTrue(numpy.abs(cps[0,i,0]-cps[0,i,-1].conj()) < 1e-6*cps[0,i,:].max())


class fx_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.correlator.fx
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(SpecMaster_tests))
        self.addTests(loader.loadTestsFromTestCase(StokesMaster_tests))
        self.addTests(loader.loadTestsFromTestCase(FXMaster_tests))
        self.addTests(loader.loadTestsFromTestCase(FXStokes_tests))


if __name__ == '__main__':
    unittest.main()
