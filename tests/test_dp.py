# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for lsl.common.dp module."""

import warnings
import unittest
import numpy

from lsl.common import dp
from lsl.common import stations


__revision__ = "$Rev$"
__version__  = "0.2"
__author__    = "Jayce Dowell"


class dp_bandpass_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.dp
    module function for the bandpass."""
    
    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
    
    def test_tbn_bandpass(self):
        """Test that the TBN bandpass generator actually runs."""
        
        fnc = dp.tbn_filter(sample_rate=1e5, npts=256)
        junk = fnc(1e3)

    def test_drx_bandpass(self):
        """Test that the DRX bandpass generator actually runs."""
        
        fnc = dp.drx_filter(sample_rate=19.6e6, npts=256)
        junk = fnc(1e3)


class dp_software_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.dp
    module function for SoftwareDP instance."""
    
    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        
    def test_modes(self):
        """Test that various SoftwareDP modes are recognized."""
        
        sdp = dp.SoftwareDP()
        
        sdp.set_mode("DRX")
        sdp.set_mode("TBN")
        self.assertRaises(ValueError, sdp.set_mode, "TBW")
        
    def test_filters(self):
        """Test that various SoftwareDP filters work."""
        
        sdp = dp.SoftwareDP()
        
        sdp.set_mode("DRX")
        for i in xrange(7,0,-1):
            if i >= 3:
                sdp.set_filter(i)
            else:
                self.assertRaises(ValueError, sdp.set_filter, i)
                
        sdp.set_mode("TBN")
        for i in xrange(7,0,-1):
            if i >= 5:
                sdp.set_filter(i)
            else:
                self.assertRaises(ValueError, sdp.set_filter, i)
                
    def test_frequency(self):
        """Test that SoftwareDP the tuning works."""
        
        sdp = dp.SoftwareDP()
        
        for f in xrange(5, 95, 5):
            if f < 10 or f > 88:
                self.assertRaises(ValueError, sdp.set_tuning_freq, f*1e6)
            else:
                sdp.set_tuning_freq(f*1e6)
                
    def test_input(self):
        """Test the SoftwareDP filtering on some data."""
        
        sdp = dp.SoftwareDP()
        
        npts = 10000
        time = numpy.arange(npts)
        data = numpy.random.rand(npts)
        
        sdp.set_mode("DRX")
        sdp.set_filter(7)
        sdp.set_tuning_freq(40e6)
        
        output = sdp.apply_filter(time, data)
        self.assertEqual(output.size, npts/10)
        
    def test_beam(self):
        """Test the SoftwareDP beamformer on some data."""
        
        sdp = dp.SoftwareDP()
        
        antennas = stations.lwa1.get_antennas()
        
        npts = 10000
        time = numpy.arange(npts)
        data = numpy.random.rand(len(antennas), npts)*2048 - 1024
        data = data.astype(numpy.int16)
        
        course = numpy.zeros(data.shape[0])
        course[0] = 2
        fine = numpy.zeros(data.shape[0])
        gains = numpy.zeros((data.shape[0]/2, 4))
        gains[0,0] = 1
        gains[1,1] = 1
        
        sdp.set_mode("DRX")
        sdp.set_filter(7)
        sdp.set_tuning_freq(40e6)
        
        beamX, beamY = sdp.form_beam(antennas, time, data, course, fine, gains)
        beamX = numpy.round(beamX).astype(numpy.int16)
        beamY = numpy.round(beamY).astype(numpy.int16)
        
        self.assertEqual(beamX[0], data[0,0])
        self.assertEqual(beamY[0], data[2,2])


class dp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.dp
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(dp_bandpass_tests)) 
        self.addTests(loader.loadTestsFromTestCase(dp_software_tests))


if __name__ == '__main__':
    unittest.main()
