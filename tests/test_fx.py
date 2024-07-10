"""
Unit test for the lsl.correlator.fx module.
"""

import os
import time
import warnings
import unittest
import ephem
import numpy as np
from scipy.signal import get_window as scipy_get_window

from astropy.coordinates import SkyCoord, AltAz

from lsl.common.data_access import DataAccess
from lsl.common import stations
from lsl.correlator import fx
from lsl.reader.base import CI8
import lsl.testing

_SSMIF = 'lwa1-ssmif.txt'

__version__  = "0.7"
__author__    = "Jayce Dowell"


def _make_complex_data(shape, scale=1.0, offset=0.0, dtype=np.complex64):
    """
    Private function to make a complex data set and return it along with a
    version that is suitable for numpy.  This is useful for testing LSL's
    ci8 return data.
    """
    
    i = np.random.rand(*shape)*scale + offset
    q = np.random.rand(*shape)*scale + offset
    data = i + 1j*q
    if dtype == CI8:
        data = data.astype(np.complex64)
        data = data.view(np.float32)
        new_shape = list(data.shape)
        new_shape[-1] = new_shape[-1] // 2
        new_shape.append(2)
        data.shape = tuple(new_shape)
        data = np.round(data)
        data = data.astype(np.int8)
        data_comp = data[...,0] + 1j*data[...,1]
        data_comp = data_comp.astype(np.complex64)
        
        data = data.view(CI8)
        data.shape = data.shape[:-1]
    else:
        data = data.astype(dtype)
        data_comp = data.copy()
    return data, data_comp


def _pfb_filter_coeff(LFFT, ntaps):
    """
    Private function to generate the filter bank coefficients for LFFT
    channels using ntaps taps.
    """

    t = np.arange(LFFT*ntaps)
    return np.sinc((t - LFFT*ntaps/2.0 + 0.5)/LFFT)


def _pfb(data, start, LFFT, ntaps=4):
    """
    Private function to compute a PFB at the specified location in the data.
    """
    
    sub = np.zeros(LFFT*ntaps, dtype=data.dtype)
    for i in range(ntaps):
        j = start//LFFT - (ntaps-1) + i
        if j < 0:
            continue
        sub[i*LFFT:(i+1)*LFFT] = data[j*LFFT:(j+1)*LFFT]
    sub = sub*_pfb_filter_coeff(LFFT, ntaps)*scipy_get_window("hamming", LFFT*ntaps)
    pfb  = np.fft.fft(sub[0*LFFT:1*LFFT])
    for i in range(1, ntaps):
        pfb += np.fft.fft(sub[i*LFFT:(i+1)*LFFT])
    return pfb


class SpecMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.SpecMaster
    function."""
    
    nAnt = 8
    
    def setUp(self):
        """Turn off all numpy and python warnings."""
        
        np.seterr(all='ignore')
        warnings.simplefilter('ignore')
        np.random.seed(1234)
        
    def run_specmaster_test_real(self, dtype, nchan=256, window=fx.null_window):
        fakeData = 10.0*np.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        freq, spectra = fx.SpecMaster(fakeData, LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = np.zeros_like(spectra)
        LFFT = spectra2.shape[1]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        for i in range(self.nAnt):
            for j in range(nFFT):
                spectra2[i,:] += (np.abs( np.fft.fft(fakeData[i,j*2*LFFT:(j+1)*2*LFFT]*wndw) )**2)[:LFFT]
        spectra2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def run_specmaster_test_complex(self, dtype, nchan=256, window=fx.null_window):
        fakeData, fakeData_np = _make_complex_data((self.nAnt,nchan*4), scale=16, offset=3+3j, dtype=dtype)
        freq, spectra = fx.SpecMaster(fakeData, LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = np.zeros_like(spectra)
        LFFT = spectra2.shape[1]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        for i in range(self.nAnt):
            for j in range(nFFT):
                spectra2[i,:] += np.fft.fftshift( np.abs( np.fft.fft(fakeData_np[i,j*LFFT:(j+1)*LFFT]*wndw) )**2 )
        spectra2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def test_spectra_real(self):
        """Test the SpecMaster function on real-valued data."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_real(dtype)
            
    def test_spectra_complex(self):
        """Test the SpecMaster function on complex-valued data."""
        
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_complex(dtype)
                
    def test_window(self):
        """Test that window functions can be passed to SpecMaster."""
        
        #
        # Real
        #
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_real(dtype, window=np.blackman)
                
        #
        # Complex
        #
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_complex(dtype, window=np.hamming)
            
    def test_window_custom(self):
        """Test that custom window functions can be passed to SpecMaster."""
        
        #
        # Real
        #
        def wndw(L):
            return np.kaiser(L, 5)
            
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_real(dtype, window=wndw)
                
        #
        # Complex
        #
        def wndw2(L):
            return np.kaiser(L, 1)
            
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_specmaster_test_complex(dtype, window=wndw2)
                
    def test_spectra_odd_complex(self):
        """Test the SpecMaster function on odd-sized complex transforms."""
        
        def wndw2(L):
            return np.kaiser(L, 1)
            
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype, window='none'):
                self.run_specmaster_test_complex(dtype, nchan=259)
            
            with self.subTest(dtype=dtype, window='custom'):
                self.run_specmaster_test_complex(dtype, nchan=259, window=wndw2)
                
    def test_spectra_real_pfb(self):
        """Test the PFB version of the SpecMaster function on real-valued data."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            fakeData = 10.0*np.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.SpecMaster(fakeData, pfb=True)
        
            # Numpy comparison
            spectra2 = np.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//2//LFFT
            for i in range(self.nAnt):
                for j in range(nFFT):
                    spectra2[i,:] += (np.abs( _pfb(fakeData[i,:], 2*j*LFFT, 2*LFFT) )**2)[:LFFT]
            spectra2 /= (2*LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)
        
    def test_spectra_complex_pfb(self):
        """Test the PFB version of the SpecMaster function on complex-valued data."""
        
        for dtype in (np.complex64, np.complex128):
            fakeData, fakeData_np = _make_complex_data((self.nAnt,1024*4), scale=16, offset=3+3j, dtype=dtype)
            freq, spectra = fx.SpecMaster(fakeData, pfb=True, sample_rate=1e5, central_freq=38e6)
        
            # Numpy comparison
            spectra2 = np.zeros_like(spectra)
            LFFT = spectra2.shape[1]
            nFFT = fakeData.shape[1]//LFFT
            for i in range(self.nAnt):
                for j in range(nFFT): 
                    spectra2[i,:] += np.fft.fftshift( np.abs( _pfb(fakeData_np[i,:], j*LFFT, LFFT) )**2 )
            spectra2 /= (LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)


class StokesMaster_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.fx.StokesMaster
    function."""

    nAnt = 8

    def setUp(self):
        """Turn off all numpy and python warnings."""

        np.seterr(all='ignore')
        warnings.simplefilter('ignore')
        np.random.seed(1234)
        
    def run_stokesmaster_test_real(self, dtype, nchan=256, window=fx.null_window):
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        fakeData = 10.0*np.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = np.zeros_like(spectra)
        LFFT = spectra2.shape[2]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        for i in range(self.nAnt//2):
            for j in range(nFFT):
                xF = np.fft.fft(fakeData[2*i+0,j*2*LFFT:(j+1)*2*LFFT]*wndw)[:LFFT]
                yF = np.fft.fft(fakeData[2*i+1,j*2*LFFT:(j+1)*2*LFFT]*wndw)[:LFFT]
                
                spectra2[0,i,:] += np.abs(xF)**2 + np.abs(yF)**2
                spectra2[1,i,:] += np.abs(xF)**2 - np.abs(yF)**2
                spectra2[2,i,:] += 2*(xF*yF.conj()).real
                spectra2[3,i,:] += 2*(xF*yF.conj()).imag
        spectra2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def run_stokesmaster_test_complex(self, dtype, nchan=256, window=fx.null_window):
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        fakeData, fakeData_np = _make_complex_data((self.nAnt,nchan*4), scale=16, offset=3+3j, dtype=dtype)
        freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        # Numpy comparison
        spectra2 = np.zeros_like(spectra)
        LFFT = spectra2.shape[2]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        for i in range(self.nAnt//2):
            for j in range(nFFT):
                xF = np.fft.fftshift( np.fft.fft(fakeData_np[2*i+0,j*LFFT:(j+1)*LFFT]*wndw) )
                yF = np.fft.fftshift( np.fft.fft(fakeData_np[2*i+1,j*LFFT:(j+1)*LFFT]*wndw) )
                
                spectra2[0,i,:] += np.abs(xF)**2 + np.abs(yF)**2
                spectra2[1,i,:] += np.abs(xF)**2 - np.abs(yF)**2
                spectra2[2,i,:] += 2*(xF*yF.conj()).real
                spectra2[3,i,:] += 2*(xF*yF.conj()).imag
        spectra2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(spectra, spectra2)
        
    def test_spectra_real(self):
        """Test the StokesMaster function on real-valued data."""

        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_real(dtype)
                
    def test_spectra_complex(self):
        """Test the StokesMaster function on complex-valued data."""

        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_complex(dtype)
                
    def test_window(self):
        """Test that window functions can be passed to StokesMaster."""
        
        #
        # Real
        #
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_real(dtype, window=np.blackman)
                
        #
        # Complex
        #
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_complex(dtype, window=np.hamming)
                
    def test_window_custom(self):
        """Test that custom window functions can be passed to StokesMaster."""
        
        #
        # Real
        #
        
        def wndw(L):
            return np.kaiser(L, 5)
            
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_real(dtype, window=wndw)
                
        #
        # Complex
        #
        def wndw2(L):
            return np.kaiser(L, 1)
            
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_stokesmaster_test_complex(dtype, window=wndw2)
                
    def test_spectra_odd_complex(self):
        """Test the SpecMaster function on odd-sized complex transforms."""
        
        def wndw2(L):
            return np.kaiser(L, 1)
            
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype, window='none'):
                self.run_stokesmaster_test_complex(dtype, nchan=259)
                
            with self.subTest(dtype=dtype, window='custom'):
                self.run_stokesmaster_test_complex(dtype, nchan=259, window=wndw2)
                
    def test_spectra_real_pfb(self):
        """Test the PFB version of the StokesMaster function on real-valued data."""
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            fakeData = 10.0*np.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], pfb=True)
        
            # Numpy comparison
            spectra2 = np.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//2//LFFT
            for i in range(self.nAnt//2):
                for j in range(nFFT):
                    xF = _pfb(fakeData[2*i+0,:], 2*j*LFFT, 2*LFFT)[:LFFT]
                    yF = _pfb(fakeData[2*i+1,:], 2*j*LFFT, 2*LFFT)[:LFFT]
                
                    spectra2[0,i,:] += np.abs(xF)**2 + np.abs(yF)**2
                    spectra2[1,i,:] += np.abs(xF)**2 - np.abs(yF)**2
                    spectra2[2,i,:] += 2*(xF*yF.conj()).real
                    spectra2[3,i,:] += 2*(xF*yF.conj()).imag
            spectra2 /= (2*LFFT * nFFT)
            lsl.testing.assert_allclose(spectra, spectra2)
            
    def test_spectra_complex_pfb(self):
        """Test the PFB version of the StokesMaster function on complex-valued data."""
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        for dtype in (CI8, np.complex64, np.complex128):
            fakeData, fakeData_np = _make_complex_data((self.nAnt,1024*4), scale=16, offset=3+3j, dtype=dtype)
            freq, spectra = fx.StokesMaster(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
            # Numpy comparison
            spectra2 = np.zeros_like(spectra)
            LFFT = spectra2.shape[2]
            nFFT = fakeData.shape[1]//LFFT
            for i in range(self.nAnt//2):
                for j in range(nFFT):
                    xF = np.fft.fftshift( _pfb(fakeData_np[2*i+0,:], j*LFFT, LFFT) )
                    yF = np.fft.fftshift( _pfb(fakeData_np[2*i+1,:], j*LFFT, LFFT) )
                
                    spectra2[0,i,:] += np.abs(xF)**2 + np.abs(yF)**2
                    spectra2[1,i,:] += np.abs(xF)**2 - np.abs(yF)**2
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
        
        np.seterr(all='ignore')
        warnings.simplefilter('ignore')
        np.random.seed(1234)
        
    def run_correlator_test_real(self, dtype, nchan=256, window=fx.null_window, phase_center='z'):
        fakeData = 10.0*np.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window, phase_center=phase_center)
        
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
        
        cps2 = np.zeros_like(cps)
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
                    f1 = np.fft.fft(fakeData[i,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    f2 = np.fft.fft(fakeData[j,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    
                    cps2[blc,:] += f1*f2.conj()
                blc += 1
        cps2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def run_correlator_test_complex(self, dtype, nchan=256, window=fx.null_window, phase_center='z'):
        fakeData, fakeData_np = _make_complex_data((self.nAnt,1024), scale=16, offset=3+3j, dtype=dtype)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan,
                                sample_rate=1e5, central_freq=38e6,
                                window=window, phase_center=phase_center)
        
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], LFFT=nchan, 
                                sample_rate=1e5, central_freq=38e6,
                                window=window)
        
        cps2 = np.zeros_like(cps)
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
                    f1 = np.fft.fftshift( np.fft.fft(fakeData_np[i,k*LFFT:(k+1)*LFFT]*wndw) )
                    f2 = np.fft.fftshift( np.fft.fft(fakeData_np[j,k*LFFT:(k+1)*LFFT]*wndw) )
                    
                    cps2[blc,:] += f1*f2.conj()
                blc += 1
        cps2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def test_phase_centers(self):
        """Test the C-base correlator with multiple phase center types."""
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
            
        with self.subTest(type='str'):
            self.run_correlator_test_real(np.int8, phase_center='z')
            
            self.assertRaises(ValueError, self.run_correlator_test_real, np.int8, phase_center='notgoingtowork')
            
        with self.subTest(type='ephem.Body'):
            bdy = ephem.FixedBody()
            bdy._ra = '1:02:03'
            bdy._dec = '+89:00:00'
            bdy._epoch = ephem.J2000
            bdy.compute(station)
            
            self.run_correlator_test_real(np.int8, phase_center=bdy)
            
        with self.subTest(type='astropy.coordinates.SkyCoord'):
            sc = SkyCoord('1h02m03s', '+89d00m00s', frame='fk5', equinox='J2000')
            aa = AltAz(location=station.earth_location, obstime='J2000')
            sc = sc.transform_to(aa)
            
            self.run_correlator_test_real(np.int8, phase_center=sc)
            
        with self.subTest(type='tuple'):
            self.run_correlator_test_real(np.int8, phase_center=(0,90))
            
    def test_correlator_real(self):
        """Test the C-based correlator on real-valued data."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype)
                
    def test_correlator_complex(self):
        """Test the C-based correlator on complex-valued data."""
        
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype)
                
    def test_correlator_real_window(self):
        """Test the C-based correlator on real-valued data window."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype, window=np.blackman)
                
    def test_correlator_complex_window(self):
        """Test the C-based correlator on complex-valued data window."""
        
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype, window=np.hamming)
                
    def test_correlator_real_pfb(self):
        """Test the C-based PFB version of the correlator on real-valued data."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            fakeData = 10.0*np.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
        
            with DataAccess.open(_SSMIF, 'r') as fh:
                station = stations.parse_ssmif(fh)
            antennas = station.antennas
        
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True)
        
            # Numpy comparison
            for i in range(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True)
        
            cps2 = np.zeros_like(cps)
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
        
        for dtype in (CI8, np.complex64, np.complex128):
            fakeData, fakeData_np = _make_complex_data((self.nAnt,1024*4), scale=16, offset=3+3j, dtype=dtype)
            
            with DataAccess.open(_SSMIF, 'r') as fh:
                station = stations.parse_ssmif(fh)
            antennas = station.antennas
        
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
            # Numpy comparison
            for i in range(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
            
            freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
            cps2 = np.zeros_like(cps)
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
                        f1 = np.fft.fftshift( _pfb(fakeData_np[i,:], LFFT*k, LFFT) )
                        f2 = np.fft.fftshift( _pfb(fakeData_np[j,:], LFFT*k, LFFT) )
                    
                        cps2[blc,:] += f1*f2.conj()
                    blc += 1
            cps2 /= (LFFT * nFFT)
            lsl.testing.assert_allclose(cps, cps2)
        
    def test_correlator_gaincorrect(self):
        """Test appling gain correction to the correlator output."""
        
        fakeData = np.random.rand(self.nAnt,1024) + 1j*np.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(np.csingle)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                            gain_correct=True)
                            
    def test_correlator_baselines(self):
        """Test that the return_baselines keyword works."""
        
        fakeData = np.random.rand(self.nAnt,1024) + 1j*np.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(np.csingle)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        blList, freq, cps = fx.FXMaster(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True)
                                    
    def test_correlator_pol(self):
        """Test various correlator polarization settings."""
        
        fakeData = np.random.rand(self.nAnt,1024) + 1j*np.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(np.csingle)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
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
            return np.kaiser(L, 1)
            
        for dtype in (CI8, np.complex64, np.complex128):
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
        
        np.seterr(all='ignore')
        warnings.simplefilter('ignore')
        np.random.seed(1234)
        
    def run_correlator_test_real(self, dtype, nchan=256, window=fx.null_window, phase_center='z'):
        fakeData = 10.0*np.random.rand(self.nAnt,nchan*8) + 3.0
        fakeData = fakeData.astype(dtype)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window, phase_center=phase_center)
                            
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan, window=window)
                            
        cps2 = np.zeros_like(cps)
        LFFT = cps.shape[2]
        nFFT = fakeData.shape[1]//2//LFFT
        wndw = window(2*LFFT)
        blc = 0
        for i in range(0, self.nAnt//2):
            for j in range(i+1, self.nAnt//2):
                for k in range(nFFT):
                    f1X = np.fft.fft(fakeData[2*i+0,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    f1Y = np.fft.fft(fakeData[2*i+1,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    f2X = np.fft.fft(fakeData[2*j+0,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    f2Y = np.fft.fft(fakeData[2*j+1,k*2*LFFT:(k+1)*2*LFFT]*wndw)[:LFFT]
                    
                    cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                    cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                    cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                    cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                blc += 1
        cps2 /= (2*LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def run_correlator_test_complex(self, dtype, nchan=256, window=fx.null_window, phase_center='z'):
        fakeData, fakeData_np = _make_complex_data((self.nAnt,1024), scale=16, offset=3+3j, dtype=dtype)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan,
                                sample_rate=1e5, central_freq=38e6, 
                                window=window, phase_center=phase_center)
                            
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], LFFT=nchan,
                                sample_rate=1e5, central_freq=38e6, 
                                window=window)
                            
        cps2 = np.zeros_like(cps)
        LFFT = cps.shape[2]
        nFFT = fakeData.shape[1]//LFFT
        wndw = window(LFFT)
        blc = 0
        for i in range(0, self.nAnt//2):
            for j in range(i+1, self.nAnt//2):
                for k in range(nFFT):
                    f1X = np.fft.fftshift( np.fft.fft(fakeData_np[2*i+0,k*LFFT:(k+1)*LFFT]*wndw) )
                    f1Y = np.fft.fftshift( np.fft.fft(fakeData_np[2*i+1,k*LFFT:(k+1)*LFFT]*wndw) )
                    f2X = np.fft.fftshift( np.fft.fft(fakeData_np[2*j+0,k*LFFT:(k+1)*LFFT]*wndw) )
                    f2Y = np.fft.fftshift( np.fft.fft(fakeData_np[2*j+1,k*LFFT:(k+1)*LFFT]*wndw) )
                    
                    cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                    cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                    cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                    cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                blc += 1
        cps2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def test_phase_centers(self):
        """Test the C-base correlator with multiple phase center types."""
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
            
        with self.subTest(type='str'):
            self.run_correlator_test_real(np.int8, phase_center='z')
            
            self.assertRaises(ValueError, self.run_correlator_test_real, np.int8, phase_center='notgoingtowork')
            
        with self.subTest(type='ephem.Body'):
            bdy = ephem.FixedBody()
            bdy._ra = '1:02:03'
            bdy._dec = '+89:00:00'
            bdy._epoch = ephem.J2000
            bdy.compute(station)
            
            self.run_correlator_test_real(np.int8, phase_center=bdy)
            
        with self.subTest(type='astropy.coordinates.SkyCoord'):
            sc = SkyCoord('1h02m03s', '+89d00m00s', frame='fk5', equinox='J2000')
            aa = AltAz(location=station.earth_location, obstime='J2000')
            sc = sc.transform_to(aa)
            
            self.run_correlator_test_real(np.int8, phase_center=sc)
            
        with self.subTest(type='tuple'):
            self.run_correlator_test_real(np.int8, phase_center=(0,90))
            
    def test_correlator_real(self):
        """Test the C-based correlator on real-valued data."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype)
                
    def test_correlator_complex(self):
        """Test the C-based correlator on complex-valued data."""
        
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype)
                
    def test_correlator_real_window(self):
        """Test the C-based correlator on real-valued data window."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_real(dtype, window=np.blackman)
            
    def test_correlator_complex_window(self):
        """Test the C-based correlator on complex-valued data window."""
        
        for dtype in (CI8, np.complex64, np.complex128):
            with self.subTest(dtype=dtype):
                self.run_correlator_test_complex(dtype, window=np.hamming)
            
    def test_correlator_real_pfb(self):
        """Test the C-based PFB version of the correlator on real-valued data."""
        
        for dtype in (np.int8, np.int16, np.int32, np.int64, np.float32, np.float64):
            fakeData = 10.0*np.random.rand(self.nAnt,1024*4) + 3.0
            fakeData = fakeData.astype(dtype)
        
            with DataAccess.open(_SSMIF, 'r') as fh:
                station = stations.parse_ssmif(fh)
            antennas = station.antennas
        
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True)
        
            # Numpy comparison
            for i in range(self.nAnt):
                antennas[i].stand.x = 0.0
                antennas[i].stand.y = 0.0
                antennas[i].stand.z = 0.0
                antennas[i].cable.length = 0.0
            
            freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True)
        
            cps2 = np.zeros_like(cps)
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
        
        fakeData = np.random.rand(self.nAnt,1024*4) + 1j*np.random.rand(self.nAnt,1024*4)
        fakeData = fakeData.astype(np.csingle)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
        # Numpy comparison
        for i in range(self.nAnt):
            antennas[i].stand.x = 0.0
            antennas[i].stand.y = 0.0
            antennas[i].stand.z = 0.0
            antennas[i].cable.length = 0.0
            
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], pfb=True, sample_rate=1e5, central_freq=38e6)
        
        cps2 = np.zeros_like(cps)
        LFFT = cps.shape[2]
        nFFT = fakeData.shape[1]//LFFT
        blc = 0
        for i in range(0, self.nAnt//2):
            for j in range(i+1, self.nAnt//2):
                for k in range(nFFT):
                    f1X = np.fft.fftshift( _pfb(fakeData[2*i+0,:], k*LFFT, LFFT) )
                    f1Y = np.fft.fftshift( _pfb(fakeData[2*i+1,:], k*LFFT, LFFT) )
                    f2X = np.fft.fftshift( _pfb(fakeData[2*j+0,:], k*LFFT, LFFT) )
                    f2Y = np.fft.fftshift( _pfb(fakeData[2*j+1,:], k*LFFT, LFFT) )
                    
                    cps2[0,blc,:] += f1X*f2X.conj() + f1Y*f2Y.conj()
                    cps2[1,blc,:] += f1X*f2X.conj() - f1Y*f2Y.conj()
                    cps2[2,blc,:] += f1X*f2Y.conj() + f1X.conj()*f2Y
                    cps2[3,blc,:] += (f1X*f2Y.conj() - f1X.conj()*f2Y)/1j
                blc += 1
        cps2 /= (LFFT * nFFT)
        lsl.testing.assert_allclose(cps, cps2)
        
    def test_correlator_gaincorrect(self):
        """Test appling gain correction to the correlator output."""
        
        fakeData = np.random.rand(self.nAnt,1024) + 1j*np.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(np.csingle)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                            gain_correct=True)
                            
    def test_correlator_baselines(self):
        """Test that the return_baselines keyword works."""
        
        fakeData = np.random.rand(self.nAnt,1024) + 1j*np.random.rand(self.nAnt,1024)
        fakeData = fakeData.astype(np.csingle)
        
        with DataAccess.open(_SSMIF, 'r') as fh:
            station = stations.parse_ssmif(fh)
        antennas = station.antennas
        
        blList, freq, cps = fx.FXStokes(fakeData, antennas[:self.nAnt], sample_rate=1e5, central_freq=38e6, 
                                    return_baselines=True)
        
    def test_correlator_odd_complex(self):
        """Test the FXStokes function on odd-sized complex transforms."""
        
        def wndw2(L):
            return np.kaiser(L, 1)
            
        for dtype in (CI8, np.complex64, np.complex128):
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
