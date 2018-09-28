# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for the lsl.correlator.uvUtils module."""

import warnings
import unittest
import numpy

from lsl.correlator import uvUtils
from lsl.common import stations


__version__  = "0.5"
__revision__ = "$Rev$"
__author__    = "Jayce Dowell"


class uvUtils_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.correlator.uvUtils
    module."""
    
    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')

    def test_baseline_gen(self):
        """Test that the generated baselines contain the correct numbers of elements."""

        standList = numpy.array([100, 101, 102, 103])

        bl = uvUtils.getBaselines(standList, include_auto=False, indicies=False)
        self.assertEqual(len(bl), 6)
        bl = uvUtils.getBaselines(standList, include_auto=True, indicies=False)
        self.assertEqual(len(bl), 10)

    def test_baseline_ind(self):
        """Test that the baselines generated with indicies do return indicies and vice
        versa."""

        standList = numpy.array([100, 101, 102, 103])

        bl = uvUtils.getBaselines(standList, include_auto=False, indicies=False)
        bl = numpy.array(bl)
        self.assertTrue(bl.min() == 100)
        bl = uvUtils.getBaselines(standList, include_auto=True, indicies=False)
        bl = numpy.array(bl)
        self.assertTrue(bl.min() == 100)

        bl = uvUtils.getBaselines(standList, include_auto=False, indicies=True)
        bl = numpy.array(bl)
        self.assertTrue(bl.max() < 100)
        bl = uvUtils.getBaselines(standList, include_auto=True, indicies=True)
        bl = numpy.array(bl)
        self.assertTrue(bl.max() < 100)
        
    def test_antenna_lookup(self):
        """Test baseline number to antenna lookup function."""
        
        standList = numpy.array([100, 101, 102, 103])

        bl = uvUtils.getBaselines(standList, include_auto=False, indicies=False)
        ind = uvUtils.baseline2antenna(0, standList)
        self.assertEqual(ind[0], 100)
        self.assertEqual(ind[1], 101)
        
        ind = uvUtils.baseline2antenna(1, standList, baseline_list=bl)
        self.assertEqual(ind[0], 100)
        self.assertEqual(ind[1], 102)
        
    def test_baseline_lookup(self):
        """Test antennas to baseline lookup function."""
        
        standList = numpy.array([100, 101, 102, 103])
        bl = uvUtils.getBaselines(standList, include_auto=False, indicies=False)
        
        ind = uvUtils.antenna2baseline(100, 101, standList, include_auto=False, indicies=False)
        self.assertEqual(ind, 0)
        
        ind = uvUtils.antenna2baseline(100, 102, standList, baseline_list=bl)
        self.assertEqual(ind, 1)
        
        ind = uvUtils.antenna2baseline(0, 3, standList, include_auto=False, indicies=True)
        self.assertEqual(ind, 2)
        
    def test_computeUVW(self):
        """Test the the computeUVW function runs."""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        
        # Frequency is a scalar
        freq = 45e6
        out = uvUtils.computeUVW(antennas[0:60:2], freq=freq)
        self.assertEqual(len(out.shape), 3)
        self.assertEqual(out.shape[-1], 1)
        
        # Frequency is a list
        freq = [45e6, 60e6]
        out = uvUtils.computeUVW(antennas[0:60:2], freq=freq)
        self.assertEqual(len(out.shape), 3)
        self.assertEqual(out.shape[-1], 2)

        # Frequency is an array
        ## 1-D
        freq = numpy.linspace(45e6, 60e6, 1024)
        out0 = uvUtils.computeUVW(antennas[0:60:2], freq=freq)
        
        ## 2-D
        freq.shape = (512, 2)
        out1 = uvUtils.computeUVW(antennas[0:60:2], freq=freq)
        
        ## 3-D
        freq.shape = (128, 4, 2)
        out2 = uvUtils.computeUVW(antennas[0:60:2], freq=freq)
        
        shape0 = (out0.shape[0], 3, 1024)
        shape1 = (out0.shape[0], 3, 512, 2)
        shape2 = (out0.shape[0], 3, 128, 4, 2)
        
        # Make sure we have the right dimensions
        for i in xrange(len(shape0)):
            self.assertEqual(out0.shape[i], shape0[i])
        for i in xrange(len(shape1)):
            self.assertEqual(out1.shape[i], shape1[i])
        for i in xrange(len(shape2)):
            self.assertEqual(out2.shape[i], shape2[i])
            
        # Make sure the values are the same
        out1.shape = shape0
        out2.shape = shape0
        diff01 = ((out0 - out1)**2).sum()
        diff02 = ((out0 - out2)**2).sum()
        self.assertAlmostEqual(diff01, 0.0, 6)
        self.assertAlmostEqual(diff02, 0.0, 6)
        
    def test_computeUVTrack(self):
        """Test that the computeUVTrack function runs."""
        
        station = stations.lwa1
        antennas = station.getAntennas()
        
        out = uvUtils.computeUVTrack(antennas[0:60:2])
        
        # Make sure we have the right dimensions
        self.assertEqual(out.shape, (435,2,512))
        
        
class uvUtils_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(uvUtils_tests)) 


if __name__ == '__main__':
    unittest.main()
