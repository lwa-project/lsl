"""
Unit tests for the lsl.complex module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import numpy
import unittest

from lsl.complex import *


__version__  = "0.1"
__author__    = "Jayce Dowell"


def _c2c(value):
    r, i = numpy.int32(value.real), numpy.int32(value.imag)
    return numpy.complex64(r+i*1j)


class complex_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.complex
    module."""
    
    def test_casting_ci8(self):
        """Test building 4+4-bit complex integer objects."""
        
        for i in range(-8, 8):
            for j in range(-8, 8):
                self.assertEqual(complex_int8(i), complex_int8(i, 0))
                self.assertEqual(complex_int8(j*1j), complex_int8(0, j))
                
                a = complex_int8(i, j)
                self.assertEqual(i+j*1j, numpy.complex64(a))
                self.assertEqual(a.real_part(), i)
                self.assertEqual(a.imag_part(), j)
                
    def test_casting_ci16(self):
        """Test building 8+8-bit complex integer objects."""
        
        for i in range(-8, 8):
            for j in range(-8, 8):
                self.assertEqual(complex_int16(i), complex_int16(i, 0))
                self.assertEqual(complex_int16(j*1j), complex_int16(0, j))
                self.assertEqual(complex_int16(complex_int8(i, j)), complex_int16(i, j))
                
                a = complex_int16(i, j)
                self.assertEqual(i+j*1j, numpy.complex64(a))
                self.assertEqual(a.real_part(), i)
                self.assertEqual(a.imag_part(), j)
                
    def test_casting_ci32(self):
        """Test building 16+16-bit complex integer objects."""
        
        for i in range(-8, 8):
            for j in range(-8, 8):
                self.assertEqual(complex_int32(i), complex_int32(i, 0))
                self.assertEqual(complex_int32(j*1j), complex_int32(0, j))
                self.assertEqual(complex_int32(complex_int8(i, j)), complex_int32(i, j))
                self.assertEqual(complex_int32(complex_int16(i*10, j*10)), complex_int32(i*10, j*10))
                
                a = complex_int32(i, j)
                self.assertEqual(i+j*1j, numpy.complex64(a))
                self.assertEqual(a.real_part(), i)
                self.assertEqual(a.imag_part(), j)
                
    def test_math_ci8(self):
        """Test 4+4-bit complex integer math."""
        
        num1 = (1, 3)
        num2 = (0, 2)
        
        a = complex_int8(*num1)
        b = complex_int8(*num2)
        
        c = num1[0] + 1j*num1[1]
        d = num2[0] + 1j*num2[1]
        
        
        self.assertAlmostEqual(numpy.abs(a), numpy.abs(c), 6)
        self.assertAlmostEqual(numpy.angle(a), numpy.angle(c), 6)
        
        self.assertEqual(numpy.complex64(a+b), c+d)
        self.assertEqual(numpy.complex64(a+2), c+2)
        self.assertEqual(numpy.complex64(2+a), 2+c)
        
        self.assertEqual(numpy.complex64(a-b), c-d)
        self.assertEqual(numpy.complex64(a-2), c-2)
        self.assertEqual(numpy.complex64(2-a), 2-c)
        
        self.assertEqual(numpy.complex64(a*b), c*d)
        self.assertEqual(numpy.complex64(a*2), c*2)
        self.assertEqual(numpy.complex64(2*a), 2*c)
        
        self.assertEqual(numpy.complex64(a/b), _c2c(c/d))
        self.assertEqual(numpy.complex64(a/2), _c2c(c/2))
        self.assertEqual(numpy.complex64(2/a), _c2c(2/c))
        
    def test_math_ci16(self):
        """Test 8+8-bit complex integer math."""
        
        num1 = (1, 3)
        num2 = (0, 2)
        
        a = complex_int16(*num1)
        b = complex_int16(*num2)
        
        c = num1[0] + 1j*num1[1]
        d = num2[0] + 1j*num2[1]
        
        self.assertAlmostEqual(numpy.abs(a), numpy.abs(c), 6)
        self.assertAlmostEqual(numpy.angle(a), numpy.angle(c), 6)
        
        self.assertEqual(numpy.complex64(a+b), c+d)
        self.assertEqual(numpy.complex64(a+2), c+2)
        self.assertEqual(numpy.complex64(2+a), 2+c)
        
        self.assertEqual(numpy.complex64(a-b), c-d)
        self.assertEqual(numpy.complex64(a-2), c-2)
        self.assertEqual(numpy.complex64(2-a), 2-c)
        
        self.assertEqual(numpy.complex64(a*b), c*d)
        self.assertEqual(numpy.complex64(a*2), c*2)
        self.assertEqual(numpy.complex64(2*a), 2*c)
        
        self.assertEqual(numpy.complex64(a/b), _c2c(c/d))
        self.assertEqual(numpy.complex64(a/2), _c2c(c/2))
        self.assertEqual(numpy.complex64(2/a), _c2c(2/c))
        
    def test_math_ci32(self):
        """Test 16+16-bit complex integer math."""
        
        num1 = (1, 3)
        num2 = (0, 2)
        
        a = complex_int32(*num1)
        b = complex_int32(*num2)
        
        c = num1[0] + 1j*num1[1]
        d = num2[0] + 1j*num2[1]
        
        self.assertAlmostEqual(numpy.abs(a), numpy.abs(c), 6)
        self.assertAlmostEqual(numpy.angle(a), numpy.angle(c), 6)
        
        self.assertEqual(numpy.complex64(a+b), c+d)
        self.assertEqual(numpy.complex64(a+2), c+2)
        self.assertEqual(numpy.complex64(2+a), 2+c)
        
        self.assertEqual(numpy.complex64(a-b), c-d)
        self.assertEqual(numpy.complex64(a-2), c-2)
        self.assertEqual(numpy.complex64(2-a), 2-c)
        
        self.assertEqual(numpy.complex64(a*b), c*d)
        self.assertEqual(numpy.complex64(a*2), c*2)
        self.assertEqual(numpy.complex64(2*a), 2*c)
        
        self.assertEqual(numpy.complex64(a/b), _c2c(c/d))
        self.assertEqual(numpy.complex64(a/2), _c2c(c/2))
        self.assertEqual(numpy.complex64(2/a), _c2c(2/c))
        
    def test_casting_ci8_array(self):
        """Test building 4+4-bit complex integer arrays."""
        
        a = []
        c = []
        for i in range(-8, 8):
            for j in range(-8, 8):
                a.append(complex_int8(i,j))
                c.append(i+j*1j)
        a = numpy.array(a)
        c = numpy.array(c)
        
        self.assertEqual(a.dtype, complex_int8)
        for i in range(a.size):
            self.assertEqual(numpy.complex64(a[i]), c[i])
            
    def test_casting_ci16_array(self):
        """Test building 8+8-bit complex integer arrays."""
        
        a = []
        c = []
        for i in range(-8, 8):
            for j in range(-8, 8):
                a.append(complex_int16(i,j))
                c.append(i+j*1j)
        a = numpy.array(a)
        c = numpy.array(c)
        
        self.assertEqual(a.dtype, complex_int16)
        for i in range(a.size):
            self.assertEqual(numpy.complex64(a[i]), c[i])
            
    def test_casting_ci32_array(self):
        """Test building 16+16-bit complex integer arrays."""
        
        a = []
        c = []
        for i in range(-8, 8):
            for j in range(-8, 8):
                a.append(complex_int32(i,j))
                c.append(i+j*1j)
        a = numpy.array(a)
        c = numpy.array(c)
        
        self.assertEqual(a.dtype, complex_int32)
        for i in range(a.size):
            self.assertEqual(numpy.complex64(a[i]), c[i])
            
    def test_math_ci8_array(self):
        """Test 4+4-bit complex integer math on arrays."""
        
        num1 = (1, 1)
        num2 = (0, 2)
        
        a = complex_int8(*num1)
        b = complex_int8(*num2)
        b = numpy.array([b, b, b, b, a, a, a, a])
        a = numpy.array([a, a+b[0], b[0]-a, b[0], a, a+b[0], a-b[0], b[0]])
        
        c = num1[0] + 1j*num1[1]
        d = num2[0] + 1j*num2[1]
        d = numpy.array([d, d, d, d, c, c, c, c])
        c = numpy.array([c, c+d[0], d[0]-c, d[0], c, c+d[0], c-d[0], d[0]])
        
        numpy.testing.assert_almost_equal(numpy.abs(a), numpy.abs(c), 6)
        numpy.testing.assert_almost_equal(numpy.angle(a), numpy.angle(c), 6)
        
        numpy.testing.assert_equal(numpy.complex64(a+b), c+d)
        numpy.testing.assert_equal(numpy.complex64(a+2), c+2)
        numpy.testing.assert_equal(numpy.complex64(2+a), 2+c)
        
        numpy.testing.assert_equal(numpy.complex64(a-b), c-d)
        numpy.testing.assert_equal(numpy.complex64(a-2), c-2)
        numpy.testing.assert_equal(numpy.complex64(2-a), 2-c)
        
        numpy.testing.assert_equal(numpy.complex64(a*b), c*d)
        numpy.testing.assert_equal(numpy.complex64(a*2), c*2)
        numpy.testing.assert_equal(numpy.complex64(2*a), 2*c)
        
        numpy.testing.assert_equal(numpy.complex64(a/b), _c2c(c/d))
        numpy.testing.assert_equal(numpy.complex64(a/2), _c2c(c/2))
        numpy.testing.assert_equal(numpy.complex64(2/a), _c2c(2/c))
        
    def test_math_ci16_array(self):
        """Test 8+8-bit complex integer math on arrays."""
        
        num1 = (1, 1)
        num2 = (0, 2)
        
        a = complex_int16(*num1)
        b = complex_int16(*num2)
        b = numpy.array([b, b, b, b, a, a, a, a])
        a = numpy.array([a, a+b[0], b[0]-a, b[0], a, a+b[0], a-b[0], b[0]])
        
        c = num1[0] + 1j*num1[1]
        d = num2[0] + 1j*num2[1]
        d = numpy.array([d, d, d, d, c, c, c, c])
        c = numpy.array([c, c+d[0], d[0]-c, d[0], c, c+d[0], c-d[0], d[0]])
        
        numpy.testing.assert_almost_equal(numpy.abs(a), numpy.abs(c), 6)
        numpy.testing.assert_almost_equal(numpy.angle(a), numpy.angle(c), 6)
        
        numpy.testing.assert_equal(numpy.complex64(a+b), c+d)
        numpy.testing.assert_equal(numpy.complex64(a+2), c+2)
        numpy.testing.assert_equal(numpy.complex64(2+a), 2+c)
        
        numpy.testing.assert_equal(numpy.complex64(a-b), c-d)
        numpy.testing.assert_equal(numpy.complex64(a-2), c-2)
        numpy.testing.assert_equal(numpy.complex64(2-a), 2-c)
        
        numpy.testing.assert_equal(numpy.complex64(a*b), c*d)
        numpy.testing.assert_equal(numpy.complex64(a*2), c*2)
        numpy.testing.assert_equal(numpy.complex64(2*a), 2*c)
        
        numpy.testing.assert_equal(numpy.complex64(a/b), _c2c(c/d))
        numpy.testing.assert_equal(numpy.complex64(a/2), _c2c(c/2))
        numpy.testing.assert_equal(numpy.complex64(2/a), _c2c(2/c))
        
    def test_math_ci32_array(self):
        """Test 16+16-bit complex integer math on arrays."""
        
        num1 = (1, 1)
        num2 = (0, 2)
        
        a = complex_int32(*num1)
        b = complex_int32(*num2)
        b = numpy.array([b, b, b, b, a, a, a, a])
        a = numpy.array([a, a+b[0], b[0]-a, b[0], a, a+b[0], a-b[0], b[0]])
        
        c = num1[0] + 1j*num1[1]
        d = num2[0] + 1j*num2[1]
        d = numpy.array([d, d, d, d, c, c, c, c])
        c = numpy.array([c, c+d[0], d[0]-c, d[0], c, c+d[0], c-d[0], d[0]])
        
        numpy.testing.assert_almost_equal(numpy.abs(a), numpy.abs(c), 6)
        numpy.testing.assert_almost_equal(numpy.angle(a), numpy.angle(c), 6)
        
        numpy.testing.assert_equal(numpy.complex64(a+b), c+d)
        numpy.testing.assert_equal(numpy.complex64(a+2), c+2)
        numpy.testing.assert_equal(numpy.complex64(2+a), 2+c)
        
        numpy.testing.assert_equal(numpy.complex64(a-b), c-d)
        numpy.testing.assert_equal(numpy.complex64(a-2), c-2)
        numpy.testing.assert_equal(numpy.complex64(2-a), 2-c)
        
        numpy.testing.assert_equal(numpy.complex64(a*b), c*d)
        numpy.testing.assert_equal(numpy.complex64(a*2), c*2)
        numpy.testing.assert_equal(numpy.complex64(2*a), 2*c)
        
        numpy.testing.assert_equal(numpy.complex64(a/b), _c2c(c/d))
        numpy.testing.assert_equal(numpy.complex64(a/2), _c2c(c/2))
        numpy.testing.assert_equal(numpy.complex64(2/a), _c2c(2/c))


class complex_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.complex unit
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(complex_tests)) 


if __name__ == '__main__':
    unittest.main()
