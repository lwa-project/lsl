"""
Unit test for the lsl.statistics.robust module.
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

from lsl.statistics import robust


__version__  = "0.2"
__author__    = "Jayce Dowell"

class robust_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.statistics.robust
    module."""
    
    def setUp(self):
        """Turn off all numpy and python warnings."""

        numpy.seterr(all='ignore')
        warnings.simplefilter('ignore')
        
        self.a = numpy.array([-0.76150888, 10.41767353,  0.02535887,  0.06452704,  0.39077118,\
                              -0.46012172, -1.29228544, -0.98805374, -0.80730252,  0.40390247,\
                              -0.47337081,  0.26920613, 19.80545004,  2.63821775, -0.05748825,\
                               0.25774227, -1.00710488,  0.19595035, -0.48490421, -0.84487291,\
                               2.98662808, -1.28819736,  0.95018522,  0.12127676,  0.72756194,\
                              -0.19603773, -0.08496962, -0.06449421, -1.33233888, -5.84305279])
        self.x = numpy.array([ 1.0,          2.0,          3.0,          4.0,\
                               5.0,          6.0,          7.0,          8.0,\
                               9.0,         10.0,         11.0,         12.0,\
                              13.0,         14.0,         15.0])
        self.y = numpy.array([ 0.61325352,  -0.49450263,   4.77442597,   5.12496658,\
                               7.56993015,   7.87082987,   9.74932012,  10.22709808,\
                              13.04686648,  43.84141192,  15.9691922 ,  18.58200951,\
                              20.01579363,  20.2176033 ,  17.67552437])
        self.y1 = numpy.array([3.30994402e-01,   8.53591441e+00,   1.98878210e+01,\
                               3.77759480e+01,   5.85981538e+01,   9.70162465e+01,\
                               1.00401918e+02,   1.56681195e+02,   1.97736474e+02,\
                               2.45063947e+02,   2.96371189e+02,   3.54253174e+02,\
                               4.14467552e+02,   4.83910450e+02,   5.20380581e+02])
        self.x2 = numpy.arange(12000)*0.01
        self.y2 = self.x2**2 - 0.5*self.x2 + 1
        
    def test_biweight(self):
        """Test the biweight mean function."""
        
        """
        IDL:
        .comp biweight_mean
        .comp resistant_mean
        .comp robust_sigma
        .comp robust_linefit
        .comp robust_poly_fit

        a = [-0.76150888, 10.41767353,  0.02535887,  0.06452704,  0.39077118,$
             -0.46012172, -1.29228544, -0.98805374, -0.80730252,  0.40390247,$
             -0.47337081,  0.26920613, 19.80545004,  2.63821775, -0.05748825,$
              0.25774227, -1.00710488,  0.19595035, -0.48490421, -0.84487291,$
              2.98662808, -1.28819736,  0.95018522,  0.12127676,  0.72756194,$
             -0.19603773, -0.08496962, -0.06449421, -1.33233888, -5.84305279]

        print, '+++ Full +++'
        m = biweight_mean(a)
        print, 'biweight_mean:', m

        resistant_mean, a, 2, m, s, n
        print, 'resistant_mean @ 2:', m
        resistant_mean, a, 3, m, s, n
        print, 'resistant_mean @ 3:', m
        resistant_mean, a, 4, m, s, n
        print, 'resistant_mean @ 4:', m
        resistant_mean, a, 5, m, s, n
        print, 'resistant_mean @ 5:', m

        s = robust_sigma(a)
        print, 'robust_sigma:', s
        s = robust_sigma(a, /zero)
        print, 'robust_sigma @ zero:', s

        print, '+++ First Half +++'
        b = a[0:14]
        m = biweight_mean(b)
        print, 'biweight_mean:', m

        resistant_mean, b, 2, m, s, n
        print, 'resistant_mean @ 2:', m
        resistant_mean, b, 3, m, s, n
        print, 'resistant_mean @ 3:', m
        resistant_mean, b, 4, m, s, n
        print, 'resistant_mean @ 4:', m
        resistant_mean, b, 5, m, s, n
        print, 'resistant_mean @ 5:', m

        s = robust_sigma(b)
        print, 'robust_sigma:', s
        s = robust_sigma(b, /zero)
        print, 'robust_sigma @ zero:', s

        print, '+++ Second Half +++'
        b = a[15:29]
        m = biweight_mean(b)
        print, 'biweight_mean:', m

        resistant_mean, b, 2, m, s, n
        print, 'resistant_mean @ 2:', m
        resistant_mean, b, 3, m, s, n
        print, 'resistant_mean @ 3:', m
        resistant_mean, b, 4, m, s, n
        print, 'resistant_mean @ 4:', m
        resistant_mean, b, 5, m, s, n
        print, 'resistant_mean @ 5:', m

        s = robust_sigma(b)
        print, 'robust_sigma:', s
        s = robust_sigma(b, /zero)
        print, 'robust_sigma @ zero:', s

        x = [ 1.0,          2.0,          3.0,          4.0,$
              5.0,          6.0,          7.0,          8.0,$
              9.0,         10.0,         11.0,         12.0,$
             13.0,         14.0,         15.0]
        y = [ 0.61325352,  -0.49450263,   4.77442597,   5.12496658,$
              7.56993015,   7.87082987,   9.74932012,  10.22709808,$
             13.04686648,  43.84141192,  15.9691922 ,  18.58200951,$
             20.01579363,  20.2176033 ,  17.67552437]

        print, '+++ Full +++'
        c = robust_linefit(x, y, yfit, sig, coef_sig)
        print, 'robust_linefit:', c
        c = robust_linefit(x, y, yfit, sig, coef_sig, /bisect)
        print, 'robust_linefit @ bisect:', c
        c = robust_linefit(x, y, yfits, sig, coeff_sig, bisquare_limit=2)
        print, 'robust_linefit @ bisquare_limit=2:', c
        c = robust_linefit(x, y, yfits, sig, coeff_sig, close_factor=0.1)
        print, 'robust_linefit @ close_factor=0.1:', c

        print, '+++ First 10 +++'
        bx = x[0:9]
        by = y[0:9]
        c = robust_linefit(bx, by, yfit, sig, coef_sig)
        print, 'robust_linefit:', c
        c = robust_linefit(bx, by, yfit, sig, coef_sig, /bisect)
        print, 'robust_linefit @ bisect:', c
        c = robust_linefit(bx, by, yfits, sig, coeff_sig, bisquare_limit=2)
        print, 'robust_linefit @ bisquare_limit=2:', c
        c = robust_linefit(bx, by, yfits, sig, coeff_sig, close_factor=0.1)
        print, 'robust_linefit @ close_factor=0.1:', c
        
        print, '+++ First 4 +++'
        bx = x[0:4]
        by = y[0:4]
        c = robust_linefit(bx, by, yfit, sig, coef_sig)
        print, 'robust_linefit:', c
        c = robust_linefit(bx, by, yfit, sig, coef_sig, /bisect)
        print, 'robust_linefit @ bisect:', c
        c = robust_linefit(bx, by, yfits, sig, coeff_sig, bisquare_limit=2)
        print, 'robust_linefit @ bisquare_limit=2:', c
        c = robust_linefit(bx, by, yfits, sig, coeff_sig, close_factor=0.1)
        print, 'robust_linefit @ close_factor=0.1:', c
        
        y = [3.30994402e-01,   8.53591441e+00,   1.98878210e+01,$
             3.77759480e+01,   5.85981538e+01,   9.70162465e+01,$
             1.00401918e+02,   1.56681195e+02,   1.97736474e+02,$
             2.45063947e+02,   2.96371189e+02,   3.54253174e+02,$
             4.14467552e+02,   4.83910450e+02,   5.20380581e+02]

        print, '+++ Full +++'
        c = robust_poly_fit(x, y, 2)
        print, 'robust_poly_fit:', c
        c = robust_poly_fit(x, y, 2, numit=10)
        print, 'robust_poly_fit @ numit=10:', c

        print, '+++ First 10 +++'
        bx = x[0:9]
        by = y[0:9]
        c = robust_poly_fit(bx, by, 2)
        print, 'robust_poly_fit:', c
        c = robust_poly_fit(bx, by, 2, numit=10)
        print, 'robust_poly_fit @ numit=10:', c
        
        x = findgen(12000)*0.01
        y = x^2 - 0.5*x + 1

        print, '+++ Large +++'
        c = robust_poly_fit(x, y, 2)
        print, 'robust_poly_fit:', c
        
        Output:
        +++ Full +++
        biweight_mean:    -0.138530
        % Compiled module: POLY.
        resistant_mean @ 2:    -0.269463
        resistant_mean @ 3:    -0.269463
        resistant_mean @ 4:   -0.0411749
        resistant_mean @ 5:   -0.0411749
        robust_sigma:     0.929675
        robust_sigma @ zero:     0.960737
        +++ First Half +++
        biweight_mean:    -0.177742
        resistant_mean @ 2:    -0.217644
        resistant_mean @ 3:    -0.307197
        resistant_mean @ 4:   -0.0806268
        resistant_mean @ 5:   -0.0806268
        robust_sigma:     0.923335
        robust_sigma @ zero:     0.896659
        +++ Second Half +++
        biweight_mean:    -0.106250
        resistant_mean @ 2:    -0.234631
        resistant_mean @ 3:  -0.00454110
        resistant_mean @ 4:  -0.00454110
        resistant_mean @ 5:  -0.00454110
        robust_sigma:      1.11836
        robust_sigma @ zero:      1.10655
        +++ Full +++
        % Compiled module: ROB_CHECKFIT.
        robust_linefit:    -0.894907      1.49426
        robust_linefit @ bisect:     -1.32388      1.56321
        robust_linefit @ bisquare_limit=2:    -0.799888      1.54240
        robust_linefit @ close_factor=0.1:    -0.900498      1.49534
        +++ First 10 +++
        robust_linefit:     -1.31115      1.56634
        robust_linefit @ bisect:     -1.45311      1.59971
        robust_linefit @ bisquare_limit=2:    -0.349420      1.43579
        robust_linefit @ close_factor=0.1:     -1.30168      1.56505
        +++ First 5 +++
        robust_linefit:     -2.29467      1.94324      
        robust_linefit @ bisect:     -2.84038      2.11933
        robust_linefit @ bisquare_limit=2:     -2.29467      1.94324
        robust_linefit @ close_factor=0.1:     -2.29467      1.94324
        +++ Full +++
        % Compiled module: POLY_FIT.
        robust_poly_fit:     -1.49629    -0.269859      2.49115
        robust_poly_fit @ numit=10:     -1.49629    -0.269859      2.49115
        +++ First 10 +++
        robust_poly_fit:     -1.67901    -0.188028      2.48656
        robust_poly_fit @ numit=10:     -1.67901    -0.188028      2.48656
        +++ Large +++
        robust_poly_fit:       1.0005343     -0.50001344       1.0000001
        """
        
        self.assertAlmostEqual(robust.biweight_mean(self.a), -0.138530, 6)
        
        b = self.a.reshape(2,-1)
        self.assertAlmostEqual(robust.biweight_mean(b, axis=1)[0], -0.177742, 6)
        self.assertAlmostEqual(robust.biweight_mean(b, axis=1)[1], -0.106250, 6)
        
        b = numpy.ma.array(self.a, mask=numpy.zeros(self.a.size, dtype=bool))
        b.mask[:b.size//2] = False
        b.mask[b.size//2:] = True
        self.assertAlmostEqual(robust.biweight_mean(b), -0.177742, 6)
        b.mask[:b.size//2] = True
        b.mask[b.size//2:] = False
        self.assertAlmostEqual(robust.biweight_mean(b), -0.106250, 6)
        
    def test_mean(self):
        """Test the outlier-resistant mean function."""
        
        self.assertAlmostEqual(robust.mean(self.a, cut=2.0), -0.269463, 6)
        self.assertAlmostEqual(robust.mean(self.a, cut=3.0), -0.269463, 6)
        self.assertAlmostEqual(robust.mean(self.a, cut=4.0), -0.0411749, 6)
        self.assertAlmostEqual(robust.mean(self.a, cut=5.0), -0.0411749, 6)
        
        b = self.a.reshape(2,-1)
        self.assertAlmostEqual(robust.mean(b, cut=2.0, axis=1)[0], -0.217644, 6)
        self.assertAlmostEqual(robust.mean(b, cut=2.0, axis=1)[1], -0.234631, 6)
        
        b = numpy.ma.array(self.a, mask=numpy.zeros(self.a.size, dtype=bool))
        b.mask[:b.size//2] = False
        b.mask[b.size//2:] = True
        self.assertAlmostEqual(robust.mean(b, cut=2.0), -0.217644, 6)
        b.mask[:b.size//2] = True
        b.mask[b.size//2:] = False
        self.assertAlmostEqual(robust.mean(b, cut=2.0), -0.234631, 6)
        
    def test_mode(self):
        """Test the half-sample mode function."""
        
        b = numpy.random.randn(512)**2 + 5
        self.assertAlmostEqual(robust.mode(b), 5.0, 2)
        
    def test_std(self):
        """Test the outlier-resistant standard deviation function."""
        
        self.assertAlmostEqual(robust.std(self.a), 0.929675, 6)
        self.assertAlmostEqual(robust.std(self.a, zero=True), 0.960737, 6)
        
        b = self.a.reshape(2,-1)
        self.assertAlmostEqual(robust.std(b, axis=1)[0], 0.923335, 6)
        self.assertAlmostEqual(robust.std(b, axis=1)[1], 1.11836, 5)
        
        b = numpy.ma.array(self.a, mask=numpy.zeros(self.a.size, dtype=bool))
        b.mask[:b.size//2] = False
        b.mask[b.size//2:] = True
        self.assertAlmostEqual(robust.std(b), 0.923335, 6)
        b.mask[:b.size//2] = True
        b.mask[b.size//2:] = False
        self.assertAlmostEqual(robust.std(b, zero=True), 1.10655, 5)
        
    def test_linefit(self):
        """Test the outlier-resistant line fitter."""

        cc = robust.linefit(self.x, self.y)
        self.assertAlmostEqual(cc[0],  1.49426, 4)
        self.assertAlmostEqual(cc[1], -0.894907, 5)

        cc = robust.linefit(self.x, self.y, bisector=True)
        self.assertAlmostEqual(cc[0],  1.56321, 5)
        self.assertAlmostEqual(cc[1], -1.32388, 5)
        
        cc = robust.linefit(self.x, self.y, bisquare_limit=2.0)
        self.assertAlmostEqual(cc[0],  1.54240, 5)
        self.assertAlmostEqual(cc[1], -0.799888, 5)
        
        cc = robust.linefit(self.x, self.y, close_factor=0.1)
        self.assertAlmostEqual(cc[0],  1.49534, 5)
        self.assertAlmostEqual(cc[1], -0.900498, 5)
        
        b = numpy.ma.array(self.y, mask=numpy.zeros(self.y.size, dtype=bool))
        b.mask[:10] = False
        b.mask[10:] = True
        cc = robust.linefit(self.x, b)
        self.assertAlmostEqual(cc[0],  1.56634, 5)
        self.assertAlmostEqual(cc[1], -1.31115, 5)
        
        b = numpy.ma.array(self.y, mask=numpy.zeros(self.y.size, dtype=bool))
        b.mask[:5] = False
        b.mask[5:] = True
        cc = robust.linefit(self.x, b)
        self.assertAlmostEqual(cc[0],  1.94324, 5)
        self.assertAlmostEqual(cc[1], -2.29467, 5)
        
        b = numpy.ma.array(self.y, mask=numpy.zeros(self.y.size, dtype=bool))
        b.mask[:5] = False
        b.mask[5:] = True
        cc = robust.linefit(self.x, b, bisector=True)
        
        self.assertAlmostEqual(cc[0],  2.11933, 5)
        self.assertAlmostEqual(cc[1], -2.84038, 5)
        
    def test_polyfit(self):
        """Test the outlier-resistant polynomial fitter."""

        cc = robust.polyfit(self.x, self.y1, 2)
        self.assertAlmostEqual(cc[0],  2.49115, 5)
        self.assertAlmostEqual(cc[1], -0.269859, 5)
        self.assertAlmostEqual(cc[2], -1.49629, 5)
        
        cc = robust.polyfit(self.x, self.y1, 2, max_iter=10)
        self.assertAlmostEqual(cc[0],  2.49115, 5)
        self.assertAlmostEqual(cc[1], -0.269859, 5)
        self.assertAlmostEqual(cc[2], -1.49629, 5)
        
        b = numpy.ma.array(self.y1, mask=numpy.zeros(self.y1.size, dtype=bool))
        b.mask[:10] = False
        b.mask[10:] = True
        cc = robust.polyfit(self.x, b, 2)
        self.assertAlmostEqual(cc[0],  2.48656, 5)
        self.assertAlmostEqual(cc[1], -0.188028, 5)
        self.assertAlmostEqual(cc[2], -1.67901, 4)
        
        cc = robust.polyfit(self.x2, self.y2, 2)
        self.assertAlmostEqual(cc[0],  1.0000001, 5)
        self.assertAlmostEqual(cc[1], -0.50001344, 4)
        self.assertAlmostEqual(cc[2],  1.0005343, 2)


class robust_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.statistics.robust 
    units tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(robust_tests)) 


if __name__ == '__main__':
    unittest.main()
