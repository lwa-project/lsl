"""
Unit test for the lsl.misc.mathutils module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import unittest
import math
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


import numpy
from scipy.special import sph_harm

from lsl.misc import mathutils


__version__   = "0.2"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class mathutils_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.mathutils
    module."""
    
    def test_downsample(self):
        """Test mathutils.downsample() function."""
        
        x = numpy.array((0, 1, 2, 3))
        
        y = mathutils.downsample(x, 2, False)
        self.assertEqual(len(y), 2)
        self.assertAlmostEqual(y[0], 1.0)
        self.assertAlmostEqual(y[1], 5.0)
        
        y = mathutils.downsample(x, 2, True)
        self.assertEqual(len(y), 2)
        self.assertAlmostEqual(y[0], 0.5)
        self.assertAlmostEqual(y[1], 2.5)
        
        y = mathutils.downsample(x, 3)
        self.assertEqual(len(y), 1)
        self.assertAlmostEqual(y[0], 1.0)
        
    def test_cmagnitde(self):
        """Test mathutils.cmagnitude() function."""
        
        x = complex(1, 0)
        self.assertAlmostEqual(mathutils.cmagnitude(x), 1.0)
        
        x = complex(0, 1)
        self.assertAlmostEqual(mathutils.cmagnitude(x), 1.0)
        
        x = complex(-1, 0)
        self.assertAlmostEqual(mathutils.cmagnitude(x), 1.0)
        
        x = complex(0, -1)
        self.assertAlmostEqual(mathutils.cmagnitude(x), 1.0)
        
        x = complex(3, 4)
        self.assertAlmostEqual(mathutils.cmagnitude(x), 5.0)
        
    def test_cphase(self):
        """Test mathutils.cphase() function."""
        
        x = complex(1, 0)
        self.assertAlmostEqual(mathutils.cphase(x), math.radians(0.0))
        
        x = complex(0, 1)
        self.assertAlmostEqual(mathutils.cphase(x), math.radians(90.0))
        
        x = complex(-1, 0)
        self.assertAlmostEqual(mathutils.cphase(x), math.radians(180.0))
        
        x = complex(0, -1)
        self.assertAlmostEqual(mathutils.cphase(x), math.radians(-90.0))
        
    def test_cpolar(self):
        """Test mathutils.cpolar() function."""
        
        x = numpy.array((0+1j, 1+0j, 1+1j, -1+0j))
        mag = (1.0, 1.0, math.sqrt(2.0), 1.0)
        phase = (math.radians(90.0), 0.0, math.radians(45.0), math.radians(180.0))
        im = iter(mag)
        ip = iter(phase)
        
        pol = mathutils.cpolar(x)
        for p in pol:
            self.assertAlmostEqual(p[0], next(im))
            self.assertAlmostEqual(p[1], next(ip))
            
    def test_creal(self):
        """Test mathtutil.creal() function."""
        
        x = (1, math.radians(0.0))
        self.assertAlmostEqual(mathutils.creal(x), 1.0)
        
        x = (1, math.radians(90.0))
        self.assertAlmostEqual(mathutils.creal(x), 0.0)
        
        x = (1, math.radians(180.0))
        self.assertAlmostEqual(mathutils.creal(x), -1.0)
        
        x = (1, math.radians(-90.0))
        self.assertAlmostEqual(mathutils.creal(x), 0.0)
            
    def test_cimag(self):
        """Test mathtutil.cimag() function."""
        
        x = (1, math.radians(0.0))
        self.assertAlmostEqual(mathutils.cimag(x), 0.0)
        
        x = (1, math.radians(90.0))
        self.assertAlmostEqual(mathutils.cimag(x), 1.0)
        
        x = (1, math.radians(180.0))
        self.assertAlmostEqual(mathutils.cimag(x), 0.0)
        
        x = (1, math.radians(-90.0))
        self.assertAlmostEqual(mathutils.cimag(x), -1.0)
        
    def test_crect(self):
        """Test mathutils.crect() function."""
        
        x = numpy.array(((1.0, math.radians(90.0)),
                    (1.0, 0.0),
                    (math.sqrt(2.0), math.radians(45.0)),
                    (1.0, math.radians(180.0))))
        
        c = (0+1j, 1+0j, 1+1j, -1+0j)
        ic = iter(c)
        
        rect = mathutils.crect(x)
        for r in rect:
            co = next(ic)
            self.assertAlmostEqual(r.real, co.real)
            self.assertAlmostEqual(r.imag, co.imag)
            
    def test_regrid(self):
        """Test mathutils.regrid() function."""
        
        yout = (0.0, 2.0, 4.0, 6.0, 8.0)
        
        x = numpy.arange(0, 10, dtype = numpy.float32)
        y = numpy.arange(0, 10, dtype = numpy.float32)
        xnew = numpy.arange(0, 10, 2, dtype = numpy.float32)
        
        ynew = mathutils.regrid(x, y, xnew, method = 'spline')
        iy = iter(yout)
        for yn in ynew:
            self.assertAlmostEqual(yn, next(iy))
            
        ynew = mathutils.regrid(x, y, xnew, method = 'linear')
        iy = iter(yout)
        for yn in ynew:
            self.assertAlmostEqual(yn, next(iy))
        
        xnew = numpy.arange(-2, 10, 2, dtype = numpy.float32)    
        self.assertRaises(ValueError, mathutils.regrid, x, y, xnew)
        
        xnew = numpy.arange(0, 12, 2, dtype = numpy.float32)    
        self.assertRaises(ValueError, mathutils.regrid, x, y, xnew)
        
    def test_smooth(self):
        """Test mathutils.smooth() function."""
        
        x = numpy.arange(0, 100, dtype = numpy.float32)
        mathutils.smooth(x)
        mathutils.smooth(x, window='flat')
        
    def test_to_dB(self):
        """Test mathutils.to_dB() function."""
        
        x = numpy.random.randn(100) + 1000.0
        mathutils.to_dB(x)

    def test_from_dB(self):
        """Test mathutils.from_dB function."""
        
        x = numpy.arange(1, 100, dtype = numpy.float32)
        mathutils.from_dB(x)
        
    def test_gaussian_gen(self):
        """Test 1-D and 2-D Gaussisan generating functions."""
        
        # 1-D
        height = 1
        center = 5.0
        width = 2.1
        gauFnc = mathutils.gaussian1d(height, center, width)
        value = gauFnc(numpy.arange(0, 100))
        
        # 2-D
        centerX = center
        centerY = -centerX
        widthX = width
        widthY = widthX/2
        gauFnc = mathutils.gaussian2d(height, centerX, centerY, widthX, widthY)
        value = gauFnc(numpy.arange(0, 100), numpy.arange(0,100))
        
    def test_gaussian_par(self):
        """Test 1-D and 2-D Gaussisan parameter estimation."""
        
        # 1-D
        height = 1.5
        center = 50.0
        width = 2.1
        gauFnc = mathutils.gaussian1d(height, center, width)
        
        x = numpy.arange(0, 100)
        value = gauFnc(x)
        
        params = mathutils.gaussparams(value)
        params = mathutils.gaussparams(value, x=x)
        self.assertAlmostEqual(height, params[0], 1)
        self.assertAlmostEqual(center, params[1], 1)
        self.assertAlmostEqual(width,  params[2], 1)
        
        # 2-D
        centerX = center
        centerY = center - 20.0
        widthX = width
        widthY = widthX/2.0
        gauFnc = mathutils.gaussian2d(height, centerX, centerY, widthX, widthY)
        
        x = numpy.zeros((100,100))
        y = numpy.zeros_like(x)
        for i in range(100):
            x[i,:] = i
            y[:,i] = i
        value = gauFnc(x, y)
        
        params = mathutils.gaussparams(value)
        params = mathutils.gaussparams(value, x=x, y=y)
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
        for i in range(180):
            az[i,:] = 2*i
        for i in range(45):
            alt[:,i] = 2*i
        
        out = mathutils.sphval(terms, az, alt, degrees=True, real_only=True)
        
        # Make sure we have the right output shape
        self.assertEqual(out.shape, az.shape)
    
    def test_sphfit(self):
        """Test the sphfit() function."""
        
        az  = numpy.zeros((180,45))
        alt = numpy.zeros((180,45))
        for i in range(180):
            az[i,:] = 2*i
        for i in range(45):
            alt[:,i] = 2*i
        
        # Setup a nice, easy problem
        out = 10*sph_harm(-1, 2, az*numpy.pi/180.0, alt*numpy.pi/180.0 + numpy.pi/2)

        # Evaluate the fit
        terms = mathutils.sphfit(az, alt, out, lmax=2, degrees=True)
        
        # Make sure the term with the most power correspond to the 2,-1 mode
        terms = (numpy.abs(terms) / numpy.abs(terms).max())**2
        self.assertEqual(numpy.where(terms == 1)[0][0], 5)

    
class mathutils_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lwa_user.mathutils
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(mathutils_tests))        
        
        
if __name__ == '__main__':
    unittest.main()
