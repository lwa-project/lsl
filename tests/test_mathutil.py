# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for lsl.misc.mathutil module."""

import unittest
import logging
import math
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


import numpy
from scipy.special import sph_harm

from lsl.misc import mathutil


__revision__  = "$Rev$"
__version__   = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class mathutil_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.mathutil
    module."""
    
    def setUp(self):
        """Setup for unit tests."""
        
        # disable logger since we intentionally generate errors 
        logging.basicConfig(stream = StringIO())
    
    def test_downsample(self):
        """Test mathutil.downsample() function."""
        
        x = numpy.array((0, 1, 2, 3))
        
        y = mathutil.downsample(x, 2, False)
        self.assertEqual(len(y), 2)
        self.assertAlmostEqual(y[0], 1.0)
        self.assertAlmostEqual(y[1], 5.0)
        
        y = mathutil.downsample(x, 2, True)
        self.assertEqual(len(y), 2)
        self.assertAlmostEqual(y[0], 0.5)
        self.assertAlmostEqual(y[1], 2.5)
        
        y = mathutil.downsample(x, 3)
        self.assertEqual(len(y), 1)
        self.assertAlmostEqual(y[0], 1.0)
        
    def test_cmagnitde(self):
        """Test mathutil.cmagnitude() function."""
        
        x = complex(1, 0)
        self.assertAlmostEqual(mathutil.cmagnitude(x), 1.0)
        
        x = complex(0, 1)
        self.assertAlmostEqual(mathutil.cmagnitude(x), 1.0)
        
        x = complex(-1, 0)
        self.assertAlmostEqual(mathutil.cmagnitude(x), 1.0)
        
        x = complex(0, -1)
        self.assertAlmostEqual(mathutil.cmagnitude(x), 1.0)
        
        x = complex(3, 4)
        self.assertAlmostEqual(mathutil.cmagnitude(x), 5.0)
        
    def test_cphase(self):
        """Test mathutil.cphase() function."""
        
        x = complex(1, 0)
        self.assertAlmostEqual(mathutil.cphase(x), math.radians(0.0))
        
        x = complex(0, 1)
        self.assertAlmostEqual(mathutil.cphase(x), math.radians(90.0))
        
        x = complex(-1, 0)
        self.assertAlmostEqual(mathutil.cphase(x), math.radians(180.0))
        
        x = complex(0, -1)
        self.assertAlmostEqual(mathutil.cphase(x), math.radians(-90.0))
        
    def test_cpolar(self):
        """Test mathutil.cpolar() function."""
        
        x = numpy.array((0+1j, 1+0j, 1+1j, -1+0j))
        mag = (1.0, 1.0, math.sqrt(2.0), 1.0)
        phase = (math.radians(90.0), 0.0, math.radians(45.0), math.radians(180.0))
        im = iter(mag)
        ip = iter(phase)
        
        pol = mathutil.cpolar(x)
        for p in pol:
            self.assertAlmostEqual(p[0], next(im))
            self.assertAlmostEqual(p[1], next(ip))
            
    def test_creal(self):
        """Test mathtutil.creal() function."""
        
        x = (1, math.radians(0.0))
        self.assertAlmostEqual(mathutil.creal(x), 1.0)
        
        x = (1, math.radians(90.0))
        self.assertAlmostEqual(mathutil.creal(x), 0.0)
        
        x = (1, math.radians(180.0))
        self.assertAlmostEqual(mathutil.creal(x), -1.0)
        
        x = (1, math.radians(-90.0))
        self.assertAlmostEqual(mathutil.creal(x), 0.0)
            
    def test_cimag(self):
        """Test mathtutil.cimag() function."""
        
        x = (1, math.radians(0.0))
        self.assertAlmostEqual(mathutil.cimag(x), 0.0)
        
        x = (1, math.radians(90.0))
        self.assertAlmostEqual(mathutil.cimag(x), 1.0)
        
        x = (1, math.radians(180.0))
        self.assertAlmostEqual(mathutil.cimag(x), 0.0)
        
        x = (1, math.radians(-90.0))
        self.assertAlmostEqual(mathutil.cimag(x), -1.0)
        
    def test_crect(self):
        """Test mathutil.crect() function."""
        
        x = numpy.array(((1.0, math.radians(90.0)),
                    (1.0, 0.0),
                    (math.sqrt(2.0), math.radians(45.0)),
                    (1.0, math.radians(180.0))))
        
        c = (0+1j, 1+0j, 1+1j, -1+0j)
        ic = iter(c)
        
        rect = mathutil.crect(x)
        for r in rect:
            co = next(ic)
            self.assertAlmostEqual(r.real, co.real)
            self.assertAlmostEqual(r.imag, co.imag)
            
    def test_regrid(self):
        """Test mathutil.regrid() function."""
        
        yout = (0.0, 2.0, 4.0, 6.0, 8.0)
        
        x = numpy.arange(0, 10, dtype = numpy.float32)
        y = numpy.arange(0, 10, dtype = numpy.float32)
        xnew = numpy.arange(0, 10, 2, dtype = numpy.float32)
        
        ynew = mathutil.regrid(x, y, xnew, method = 'spline')
        iy = iter(yout)
        for yn in ynew:
            self.assertAlmostEqual(yn, next(iy))
            
        ynew = mathutil.regrid(x, y, xnew, method = 'linear')
        iy = iter(yout)
        for yn in ynew:
            self.assertAlmostEqual(yn, next(iy))
        
        xnew = numpy.arange(-2, 10, 2, dtype = numpy.float32)    
        self.assertRaises(ValueError, mathutil.regrid, x, y, xnew)
        
        xnew = numpy.arange(0, 12, 2, dtype = numpy.float32)    
        self.assertRaises(ValueError, mathutil.regrid, x, y, xnew)
        
    def test_smooth(self):
        """Test mathutil.smooth() function."""
        
        x = numpy.arange(0, 100, dtype = numpy.float32)
        mathutil.smooth(x)
        mathutil.smooth(x, window='flat')
        
    def test_to_dB(self):
        """Test mathutil.to_dB() function."""
        
        x = numpy.random.randn(100) + 1000.0
        mathutil.to_dB(x)

    def test_from_dB(self):
        """Test mathutil.from_dB function."""
        
        x = numpy.arange(1, 100, dtype = numpy.float32)
        mathutil.from_dB(x)
        
    def test_savitzky_golay(self):
        """Test Savitzky-Golay smoothing."""
        
        # Actual smoothing
        data = numpy.random.rand(4096)
        out = mathutil.savitzky_golay(data, 15, 3)
        
        # Derivative
        out = mathutil.savitzky_golay(data, 15, 3, deriv=1)
        
    def test_gaussian_gen(self):
        """Test 1-D and 2-D Gaussisan generating functions."""
        
        # 1-D
        height = 1
        center = 5.0
        width = 2.1
        gauFnc = mathutil.gaussian1d(height, center, width)
        value = gauFnc(numpy.arange(0, 100))
        
        # 2-D
        centerX = center
        centerY = -centerX
        widthX = width
        widthY = widthX/2
        gauFnc = mathutil.gaussian2d(height, centerX, centerY, widthX, widthY)
        value = gauFnc(numpy.arange(0, 100), numpy.arange(0,100))
        
    def test_gaussian_par(self):
        """Test 1-D and 2-D Gaussisan parameter estimation."""
        
        # 1-D
        height = 1.5
        center = 50.0
        width = 2.1
        gauFnc = mathutil.gaussian1d(height, center, width)
        
        x = numpy.arange(0, 100)
        value = gauFnc(x)
        
        params = mathutil.gaussparams(value)
        params = mathutil.gaussparams(value, x=x)
        self.assertAlmostEqual(height, params[0], 1)
        self.assertAlmostEqual(center, params[1], 1)
        self.assertAlmostEqual(width,  params[2], 1)
        
        # 2-D
        centerX = center
        centerY = center - 20.0
        widthX = width
        widthY = widthX/2.0
        gauFnc = mathutil.gaussian2d(height, centerX, centerY, widthX, widthY)
        
        x = numpy.zeros((100,100))
        y = numpy.zeros_like(x)
        for i in range(100):
            x[i,:] = i
            y[:,i] = i
        value = gauFnc(x, y)
        
        params = mathutil.gaussparams(value)
        params = mathutil.gaussparams(value, x=x, y=y)
        self.assertAlmostEqual(height,  params[0], 1)
        self.assertAlmostEqual(centerX, params[1], 1)
        self.assertAlmostEqual(centerY, params[2], 1)
        self.assertAlmostEqual(widthX,  params[3], 1)
        self.assertAlmostEqual(widthY,  params[4], 1)
    
    def test_sphval(self):
        """Test that the sphval() function runs."""
        
        terms = numpy.array([1, 0.5, 0.4, 0.01, -0.02, -0.005])
        az  = numpy.zeros((180,45))
        alt = numpy.zeros((180,45))
        for i in xrange(180):
            az[i,:] = 2*i
        for i in xrange(45):
            alt[:,i] = 2*i
        
        out = mathutil.sphval(terms, az, alt, degrees=True, real_only=True)
        
        # Make sure we have the right output shape
        self.assertEqual(out.shape, az.shape)
    
    def test_sphfit(self):
        """Test the sphfit() function."""
        
        az  = numpy.zeros((180,45))
        alt = numpy.zeros((180,45))
        for i in xrange(180):
            az[i,:] = 2*i
        for i in xrange(45):
            alt[:,i] = 2*i
        
        # Setup a nice, easy problem
        out = 10*sph_harm(-1, 2, az*numpy.pi/180.0, alt*numpy.pi/180.0 + numpy.pi/2)

        # Evaluate the fit
        terms = mathutil.sphfit(az, alt, out, lmax=2, degrees=True)
        
        # Make sure the term with the most power correspond to the 2,-1 mode
        terms = (numpy.abs(terms) / numpy.abs(terms).max())**2
        self.assertEqual(numpy.where(terms == 1)[0][0], 5)

    
class mathutil_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lwa_user.mathutil
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(mathutil_tests))        
        
        
if __name__ == '__main__':
    unittest.main()
