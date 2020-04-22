"""
Unit test for the lsl.misc.dedispersion module.
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

from lsl.misc import dedispersion


__version__  = "0.1"
__author__    = "Jayce Dowell"

class dedispersion_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.dedispersion
    module."""
    
    def test_delay(self):
        """Test calculating the dispersion delay"""
        
        d0 = dedispersion.delay(74e6, 12.455)
        d1 = dedispersion.delay([74e6,], 12.455)
        d2 = dedispersion.delay([74e6, 80e6], 12.455)
        d3 = dedispersion.delay(numpy.array([74e6, 80e6]), 12.455)
        self.assertAlmostEqual(d0, d1, 6)
        self.assertAlmostEqual(d2[0], d3[0], 6)
        
    def _create_pulse(self, freq, bw, dm):
        data = numpy.random.randn(4096*1024) + numpy.random.randn(4096*1024)*1j
        data = data.astype(numpy.complex64)
        data[1000*1024] += 100
        
        freq = numpy.fft.fftfreq(data.size, d=1.0/bw) + freq
        inv_chirp = dedispersion._chirp(freq, -dm)
        data = numpy.fft.fft(data)
        data *= inv_chirp
        data = numpy.fft.ifft(data)
        return data
        
    def test_incoherent(self):
        """Test incoherent dedispersion"""
        
        data = self._create_pulse(74e6, 1e6, 12.455)
        data = data.reshape(-1, 1024)
        data = numpy.abs(numpy.fft.fft(data, axis=1))**2
        data = numpy.fft.fftshift(data, axes=1)
        
        
        tint = data.shape[1]/1e6
        freq = numpy.fft.fftshift(numpy.fft.fftfreq(data.shape[1], d=1.0/1e6) + 74e6)
        data1 = dedispersion.incoherent(freq, data, tint, 12.455)
        data2 = dedispersion.incoherent(freq, data, tint,  2.455)
        self.assertEqual(numpy.argmax(data1.mean(axis=1)), 1000-124)   # 124 is half the delay
        self.assertTrue(data2[1000-124,:].mean() < data1[1000-124,:].mean())
        
    def test_coherent(self):
        """Test coherent dedispersion"""
        
        data = self._create_pulse(74e6, 1e6, 12.455)
        t = numpy.arange(data.size)/1e6
        
        t, data = dedispersion.coherent(t, data, 74e6, 1e6, 12.455)
        data = data.reshape(-1, 1024)
        data = numpy.abs(numpy.fft.fft(data, axis=1))**2
        data = numpy.fft.fftshift(data, axes=1)
        
        data = data.mean(axis=1)
        self.assertEqual(numpy.argmax(data), 1000)


class dedispersion_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.dedispersion 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(dedispersion_tests)) 


if __name__ == '__main__':
    unittest.main()
