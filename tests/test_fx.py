"""
Unit test for the lsl.correlator.fx module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import time
import warnings
import unittest
import numpy
from scipy.signal import get_window as scipy_get_window

from lsl.common.paths import DATA_BUILD
from lsl.common import stations
from lsl.correlator import fx
import lsl.testing

_SSMIF = os.path.join(DATA_BUILD, 'lwa1-ssmif.txt')

__version__  = "0.7"
__author__    = "Jayce Dowell"


def _pfb_filter_coeff(LFFT, ntaps):
    """
    Private function to generate the filter bank coefficients for LFFT
    channels using ntaps taps.
    """

    t = numpy.arange(LFFT*ntaps)
    return numpy.sinc((t - LFFT*ntaps/2.0 + 0.5)/LFFT)


def _pfb(data, start, LFFT, ntaps=4):
    """
    Private function to compute a PFB at the specified location in the data.
    """
    
    sub = numpy.zeros(LFFT*ntaps, dtype=data.dtype)
    for i in range(ntaps):
        j = start//LFFT - (ntaps-1) + i
        if j < 0:
            continue
        sub[i*LFFT:(i+1)*LFFT] = data[j*LFFT:(j+1)*LFFT]
    sub = sub*_pfb_filter_coeff(LFFT, ntaps)*scipy_get_window("hamming", LFFT*ntaps)
    pfb  = numpy.fft.fft(sub[0*LFFT:1*LFFT])
    for i in range(1, ntaps):
        pfb += numpy.fft.fft(sub[i*LFFT:(i+1)*LFFT])
    return pfb


class SpecMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.SpecMaster
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def run_specmaster_test_real(self, dtype, nchan=256, window=fx.null_window):
        fakeData = 10.0*numpy.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        freq, spectra = fx.SpecMaster(fakeData, LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = numpy.zeros_like(spectra)
        LFFT = spectra2.shape[1]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        for i in range(self.nAnt):
            for j in range(nFFT):
                spectra2[i,:] += (numpy.abs( numpy.fft.fft(fakeData[i,j*2*LFFT:(j+1)*2*LFFT]*wndw) )**2)[:LFFT]
        spectra2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def run_specmaster_test_complex(self, dtype, nchan=256, window=fx.null_window):
        fakeData = numpy.random.rand(self.nAnt,nchan*4) + 1j*numpy.random.rand(self.nAnt,nchan*4) + 3.0 + 3.0j
        fakeData = fakeData.astype(dtype)
        freq, spectra = fx.SpecMaster(fakeData, LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = numpy.zeros_like(spectra)
        LFFT = spectra2.shape[1]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        for i in range(self.nAnt):
            for j in range(nFFT):
                spectra2[i,:] += numpy.fft.fftshift( numpy.abs( numpy.fft.fft(fakeData[i,j*LFFT:(j+1)*LFFT]*wndw) )**2 )
        spectra2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def test_spectra_real(self):
        """Test the SpecMaster function on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_real(dtype)
            
    def test_spectra_complex(self):
        """Test the SpecMaster function on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_complex(dtype)
                
    def test_window(self):
        """Test that window functions can be passed to SpecMaster."""
        
        #
        # Real
        #
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_real(dtype, window=numpy.blackman)
                
        #
        # Complex
        #
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_complex(dtype, window=numpy.hamming)
            
    def test_window_custom(self):
        """Test that custom window functions can be passed to SpecMaster."""
        
        #
        # Real
        #
        def wndw(L):
            return numpy.kaiser(L, 5)
            
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_real(dtype, window=wndw)
                
        #
        # Complex
        #
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_complex(dtype, window=wndw2)
                
    def test_spectra_odd_complex(self):
        """Test the SpecMaster function on odd-sized complex transforms."""
        
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype, window='none'):
                self.run_specmaster_test_complex(dtype, nchan=259)
            
            with self.subTest(dtype=dtype, window='custom'):
                self.run_specmaster_test_complex(dtype, nchan=259, window=wndw2)
                
    def test_spectra_real_pfb(self):
        """Test the PFB version of the SpecMaster function on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, pfb=True)
        
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            for i in range(self.nAnt):
                for j in range(nFFT):
                    spectra2[i,:] += (numpy.abs( _pfb(fakeData[i,:], 2*j*LFFT, 2*LFFT) )**2)[:LFFT]
            spectra2 /= (2*LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)
        
    def test_spectra_complex_pfb(self):
        """Test the PFB version of the SpecMaster function on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024*4) + 1j*numpy.random.rand(self.nAnt,1024*4) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, pfb=True, sample_rate=1e5, central_freq=38e6)
        
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            for i in range(self.nAnt):
                for j in range(nFFT): 
                    spectra2[i,:] += numpy.fft.fftshift( numpy.abs( _pfb(fakeData[i,:], j*LFFT, LFFT) )**2 )
            spectra2 /= (LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)


class StokesMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.StokesMaster
    function."""

    nAnt = 8

    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def run_stokesmaster_test_real(self, dtype, nchan=256, window=fx.null_window):
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        fakeData = 10.0*numpy.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = numpy.zeros_like(spectra)
        LFFT = spectra2.shape[2]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        for i in range(self.nAnt//2):
            for j in range(nFFT):
                xF = numpy.fft.fft(fakeData[2*i+0,j*2*LFFT:(j+1)*2*LFFT]*wndw)[:LFFT]
                yF = numpy.fft.fft(fakeData[2*i+1,j*2*LFFT:(j+1)*2*LFFT]*wndw)[:LFFT]
                
                spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                spectra2[2,i,:] += 2*(xF*yF.conj()).real
                spectra2[3,i,:] += 2*(xF*yF.conj()).imag
        spectra2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def run_stokesmaster_test_complex(self, dtype, nchan=256, window=fx.null_window):
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        fakeData = numpy.random.rand(self.nAnt,nchan*4) + 1j*numpy.random.rand(self.nAnt,nchan*4) + 3.0 + 3.0j
        fakeData = fakeData.astype(dtype)
        freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = numpy.zeros_like(spectra)
        LFFT = spectra2.shape[2]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        for i in range(self.nAnt//2):
            for j in range(nFFT):
                xF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+0,j*LFFT:(j+1)*LFFT]*wndw) )
                yF = numpy.fft.fftshift( numpy.fft.fft(fakeData[2*i+1,j*LFFT:(j+1)*LFFT]*wndw) )
                
                spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                spectra2[2,i,:] += 2*(xF*yF.conj()).real
                spectra2[3,i,:] += 2*(xF*yF.conj()).imag
        spectra2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def test_spectra_real(self):
        """Test the StokesMaster function on real-valued data."""

        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_real(dtype)
                
    def test_spectra_complex(self):
        """Test the StokesMaster function on complex-valued data."""

        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_complex(dtype)
                
    def test_window(self):
        """Test that window functions can be passed to StokesMaster."""
        
        #
        # Real
        #
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_real(dtype, window=numpy.blackman)
                
        #
        # Complex
        #
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_complex(dtype, window=numpy.hamming)
                
    def test_window_custom(self):
        """Test that custom window functions can be passed to StokesMaster."""
        
        #
        # Real
        #
        
        def wndw(L):
            return numpy.kaiser(L, 5)
            
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_real(dtype, window=wndw)
                
        #
        # Complex
        #
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_complex(dtype, window=wndw2)
                
    def test_spectra_odd_complex(self):
        """Test the SpecMaster function on odd-sized complex transforms."""
        
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype, window='none'):
                self.run_stokesmaster_test_complex(dtype, nchan=259)
            
            with self.subTest(dtype=dtype, window='custom'):
                self.run_stokesmaster_test_complex(dtype, nchan=259, window=wndw2)
                
    def test_spectra_real_pfb(self):
        """Test the PFB version of the StokesMaster function on real-valued data."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], pfb=True)
        
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            for i in range(self.nAnt//2):
                for j in range(nFFT):
                    xF = _pfb(fakeData[2*i+0,:], 2*j*LFFT, 2*LFFT)[:LFFT]
                    yF = _pfb(fakeData[2*i+1,:], 2*j*LFFT, 2*LFFT)[:LFFT]
                
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (2*LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)
            
    def test_spectra_complex_pfb(self):
        """Test the PFB version of the StokesMaster function on complex-valued data."""
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024*4) + 1j*numpy.random.rand(self.nAnt,1024*4) + 3.0 + 3.0j
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
            # Numpy comparison
            spectra2 = numpy.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            for i in range(self.nAnt//2):
                for j in range(nFFT):
                    xF = numpy.fft.fftshift( _pfb(fakeData[2*i+0,:], j*LFFT, LFFT) )
                    yF = numpy.fft.fftshift( _pfb(fakeData[2*i+1,:], j*LFFT, LFFT) )
                
                    spectra2[0,i,:] += numpy.abs(xF)**2 + numpy.abs(yF)**2
                    spectra2[1,i,:] += numpy.abs(xF)**2 - numpy.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)


class FXMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.FXMaster
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def run_correlator_test_real(self, dtype, nchan=256, window=fx.null_window):
        fakeData = 10.0*numpy.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        cps2 = numpy.zeros_like(cps)
        LFFT = cps.shape[1]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        blc = 0
        for i in range(0, self.nAnt):
            if antennas[i].pol != 0:
                continue
            for j in range(i+1, self.nAnt):
                if antennas[j].pol != 0:
                    continue
                    
                for k in range(nFFT):
                    f1 = numpy.fft.fft(fakeData[i,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    f2 = numpy.fft.fft(fakeData[j,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    
                    cps2[blc,:] += f1*f2.conj()
                blc += 1
        cps2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def run_correlator_test_complex(self, dtype, nchan=256, window=fx.null_window):
        fakeData = numpy.random.rand(self.nAnt,nchan*4) + 1j*numpy.random.rand(self.nAnt,nchan*4)
        fakeData = fakeData.astype(dtype)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan,
                                sample_rate=1e5, central_freq=38e6,
                                window=window)
        
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, 
                                sample_rate=1e5, central_freq=38e6,
                                window=window)
        
        cps2 = numpy.zeros_like(cps)
        LFFT = cps.shape[1]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        blc = 0
        for i in range(0, self.nAnt):
            if antennas[i].pol != 0:
                continue
            for j in range(i+1, self.nAnt):
                if antennas[j].pol != 0:
                    continue
                    
                for k in range(nFFT):
                    f1 = numpy.fft.fftshift( numpy.fft.fft(fakeData[i,k*LFFT:(k+1)*LFFT]*wndw) )
                    f2 = numpy.fft.fftshift( numpy.fft.fft(fakeData[j,k*LFFT:(k+1)*LFFT]*wndw) )
                    
                    cps2[blc,:] += f1*f2.conj()
                blc += 1
        cps2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def test_correlator_real(self):
        """Test the C-based correlator on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype)
                
    def test_correlator_complex(self):
        """Test the C-based correlator on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype)
                
    def test_correlator_real_window(self):
        """Test the C-based correlator on real-valued data window."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype, window=numpy.blackman)
                
    def test_correlator_complex_window(self):
        """Test the C-based correlator on complex-valued data window."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype, window=numpy.hamming)
                
    def test_correlator_real_pfb(self):
        """Test the C-based PFB version of the correlator on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
        
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
        
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True)
        
            # Numpy comparison
            for i in range(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True)
        
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            blc = 0
            for i in range(0, self.nAnt):
                if antennas[i].pol != 0:
                    continue
                for j in range(i+1, self.nAnt):
                    if antennas[j].pol != 0:
                        continue
                    
                    for k in range(nFFT):
                        f1 = _pfb(fakeData[i,:], 2*LFFT*k, 2*LFFT)[:LFFT]
                        f2 = _pfb(fakeData[j,:], 2*LFFT*k, 2*LFFT)[:LFFT]
                    
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (2*LFFT * nFFT)
            lsl.testing.assert_allclose(cps, cps2)
            
    def test_correlator_complex_pfb(self):
        """Test the C-based PFB version of the correlator on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            fakeData = numpy.random.rand(self.nAnt,1024*4) + 1j*numpy.random.rand(self.nAnt,1024*4)
            fakeData = fakeData.astype(dtype)
        
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
        
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
            # Numpy comparison
            for i in range(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            blc = 0
            for i in range(0, self.nAnt):
                if antennas[i].pol != 0:
                    continue
                for j in range(i+1, self.nAnt):
                    if antennas[j].pol != 0:
                        continue
                    
                    for k in range(nFFT):
                        f1 = numpy.fft.fftshift( _pfb(fakeData[i,:], LFFT*k, LFFT) )
                        f2 = numpy.fft.fftshift( _pfb(fakeData[j,:], LFFT*k, LFFT) )
                    
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (LFFT * nFFT)
            lsl.testing.assert_allclose(cps, cps2)
        
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
                                    return_baselines=True, pol='XX')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 0)
            self.assertEqual(ant2.pol, 0)
        
        ## YY
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, pol='YY')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 1)
            self.assertEqual(ant2.pol, 1)
        
        ## XY
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, pol='XY')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 0)
            self.assertEqual(ant2.pol, 1)
            
        ## YX
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True, pol='YX')
        for (ant1,ant2) in blList:
            self.assertEqual(ant1.pol, 1)
            self.assertEqual(ant2.pol, 0)
            
    def test_correlator_odd_complex(self):
        """Test the FXMaster function on odd-sized complex transforms."""
        
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype, window='none'):
                self.run_correlator_test_complex(dtype, nchan=259)
            
            with self.subTest(dtype=dtype, window='custom'):
                self.run_correlator_test_complex(dtype, nchan=259, window=wndw2)


class FXStokes_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.FXStokes
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        numpy.random.seed(1234)
        
    def run_correlator_test_real(self, dtype, nchan=256, window=fx.null_window):
        fakeData = 10.0*numpy.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
                            
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
                            
        cps2 = numpy.zeros_like(cps)
        LFFT = cps.shape[2]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        blc = 0
        for i in range(0, self.nAnt//2):
            for j in range(i+1, self.nAnt//2):
                for k in range(nFFT):
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
        lsl.testing.assert_allclose(cps, cps2)
        
    def run_correlator_test_complex(self, dtype, nchan=256, window=fx.null_window):
        fakeData = numpy.random.rand(self.nAnt,nchan*4) + 1j*numpy.random.rand(self.nAnt,nchan*4)
        fakeData = fakeData.astype(dtype)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan,
                                sample_rate=1e5, central_freq=38e6, 
                                window=window)
                            
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan,
                                sample_rate=1e5, central_freq=38e6, 
                                window=window)
                            
        cps2 = numpy.zeros_like(cps)
        LFFT = cps.shape[2]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        blc = 0
        for i in range(0, self.nAnt//2):
            for j in range(i+1, self.nAnt//2):
                for k in range(nFFT):
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
        lsl.testing.assert_allclose(cps, cps2)
        
    def test_correlator_real(self):
        """Test the C-based correlator on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype)
                
    def test_correlator_complex(self):
        """Test the C-based correlator on complex-valued data."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype)
                
    def test_correlator_real_window(self):
        """Test the C-based correlator on real-valued data window."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype, window=numpy.blackman)
            
    def test_correlator_complex_window(self):
        """Test the C-based correlator on complex-valued data window."""
        
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype, window=numpy.hamming)
            
    def test_correlator_real_pfb(self):
        """Test the C-based PFB version of the correlator on real-valued data."""
        
        for dtype in (numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.float32, numpy.float64):
            fakeData = 10.0*numpy.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
        
            station = stations.parse_ssmif(_SSMIF)
            antennas = station.antennas
        
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True)
        
            # Numpy comparison
            for i in range(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
            
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True)
        
            cps2 = numpy.zeros_like(cps)
            LFFT = cps.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            blc = 0
            for i in range(0, self.nAnt//2):
                for j in range(i+1, self.nAnt//2):
                    for k in range(nFFT):
                        f1X = _pfb(fakeData[2*i+0,:], k*2*LFFT, 2*LFFT)[:LFFT]
                        f1Y = _pfb(fakeData[2*i+1,:], k*2*LFFT, 2*LFFT)[:LFFT]
                        f2X = _pfb(fakeData[2*j+0,:], k*2*LFFT, 2*LFFT)[:LFFT]
                        f2Y = _pfb(fakeData[2*j+1,:], k*2*LFFT, 2*LFFT)[:LFFT]
                    
                        cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                        cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                        cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                        cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                    blc += 1
            cps2 /= (2*LFFT * nFFT)
            lsl.testing.assert_allclose(cps, cps2)
            
    def test_correlator_complex_pfb(self):
        """Test the C-based PFB version of the correlator on complex-valued data."""
        
        fakeData = numpy.random.rand(self.nAnt,1024*4) + 1j*numpy.random.rand(self.nAnt,1024*4)
        fakeData = fakeData.astype(numpy.csingle)
        
        station = stations.parse_ssmif(_SSMIF)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
        cps2 = numpy.zeros_like(cps)
        LFFT = cps.shape[2]
        nFFT = fakeData.shape[1]//LFFT
        blc = 0
        for i in range(0, self.nAnt//2):
            for j in range(i+1, self.nAnt//2):
                for k in range(nFFT):
                    f1X = numpy.fft.fftshift( _pfb(fakeData[2*i+0,:], k*LFFT, LFFT) )
                    f1Y = numpy.fft.fftshift( _pfb(fakeData[2*i+1,:], k*LFFT, LFFT) )
                    f2X = numpy.fft.fftshift( _pfb(fakeData[2*j+0,:], k*LFFT, LFFT) )
                    f2Y = numpy.fft.fftshift( _pfb(fakeData[2*j+1,:], k*LFFT, LFFT) )
                    
                    cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                    cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                    cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                    cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                blc += 1
        cps2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
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
        
        def wndw2(L):
            return numpy.kaiser(L, 1)
            
        for dtype in (numpy.complex64, numpy.complex128):
            with self.subTest(dtype=dtype, window='none'):
                self.run_correlator_test_complex(dtype, nchan=259)
            
            with self.subTest(dtype=dtype, window='custom'):
                self.run_correlator_test_complex(dtype, nchan=259, window=wndw2)


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
