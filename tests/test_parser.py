# -*- coding: utf-8 -*-

"""Unit test for lsl.misc.parser module."""


import unittest
import ephem
import numpy
import argparse
from astropy.constants import c

from lsl.misc import parser


__revision__  = "$Rev$"
__version__   = "0.1"
__author__    = "Jayce Dowell"


c = c.to('m/s').value


class parser_tests(unittest.TestCase):
    def test_positive_or_zero_int(self):
        """Test the parser.positive_or_zero_int conversion function."""
        
        self.assertEqual(parser.positive_or_zero_int('5'), 5)
        self.assertEqual(parser.positive_or_zero_int('0'), 0)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_or_zero_int, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_or_zero_int, '1.5')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_or_zero_int, '-1')
        
    def test_positive_int(self):
        """Test the parser.positive_int conversion function."""
        
        self.assertEqual(parser.positive_int('5'), 5)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_int, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_int, '1.5')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_int, '0')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_int, '-1')
        
    def test_positive_or_zero_float(self):
        """Test the parser.positive_or_zero_float conversion function."""
        
        self.assertEqual(parser.positive_or_zero_float('5'), 5.0)
        self.assertEqual(parser.positive_or_zero_float('5.1'), 5.1)
        self.assertEqual(parser.positive_or_zero_float('0'), 0.0)
        self.assertEqual(parser.positive_or_zero_float('0.0'), 0.0)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_or_zero_float, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_or_zero_float, '-1.1')
        
    def test_positive_float(self):
        """Test the parser.positive_float conversion function."""
        
        self.assertEqual(parser.positive_float('5'), 5.0)
        self.assertEqual(parser.positive_float('5.1'), 5.1)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_float, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_float, '0')
        self.assertRaises(argparse.ArgumentTypeError, parser.positive_float, '-1.1')
        
    def test_frequency(self):
        """Test the parser.frequency conversion function."""
        
        self.assertEqual(parser.frequency('45'), 45e6)
        self.assertEqual(parser.frequency('45kHz'), 45e3)
        self.assertEqual(parser.frequency('45.6 GHz'), 45.6e9)
        self.assertEqual(parser.frequency('4.4m'), c/4.4)
        self.assertEqual(parser.frequency('21cm'), c/0.21)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.frequency, '45G')
        self.assertRaises(argparse.ArgumentTypeError, parser.frequency, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.frequency, '4~5')
        
    def test_frequency_range(self):
        """Test the parser.frequency_range conversion function."""
        
        self.assertEqual(parser.frequency_range('4~5')[0], 4e6)
        self.assertEqual(parser.frequency_range('4~5')[1], 5e6)
        self.assertEqual(parser.frequency_range('4~5GHz')[0], 4e9)
        self.assertEqual(parser.frequency_range('4~5GHz')[1], 5e9)
        self.assertEqual(parser.frequency_range('4MHz~5GHz')[0], 4e6)
        self.assertEqual(parser.frequency_range('5GHz~4MHz')[1], 5e9)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.frequency_range, '45MHz~50')
        self.assertRaises(argparse.ArgumentTypeError, parser.frequency_range, '4~c')
        self.assertRaises(argparse.ArgumentTypeError, parser.frequency_range, '4')
        
    def test_wavelength(self):
        """Test the parser.wavelength conversion function."""
        
        self.assertEqual(parser.wavelength('45'), 45.0)
        self.assertEqual(parser.wavelength('45km'), 45e3)
        self.assertEqual(parser.wavelength('45.6 Gm'), 45.6e9)
        self.assertEqual(parser.wavelength('47MHz'), c/47e6)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.wavelength, '45G')
        self.assertRaises(argparse.ArgumentTypeError, parser.wavelength, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.wavelength, '4~5')
        
    def test_wavelength_range(self):
        """Test the parser.wavelength_range conversion function."""
        
        self.assertEqual(parser.wavelength_range('4~5')[0], 4.0)
        self.assertEqual(parser.wavelength_range('4~5')[1], 5.0)
        self.assertEqual(parser.wavelength_range('4~5cm')[0], 4e-2)
        self.assertEqual(parser.wavelength_range('4~5cm')[1], 5e-2)
        self.assertEqual(parser.wavelength_range('4cm~5m')[0], 4e-2)
        self.assertEqual(parser.wavelength_range('5m~4cm')[1], 5.0)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.wavelength_range, '45Mm~50')
        self.assertRaises(argparse.ArgumentTypeError, parser.wavelength_range, '4~c')
        self.assertRaises(argparse.ArgumentTypeError, parser.wavelength_range, '4')
        
    def test_date(self):
        """Test the parser.date conversion function."""
        
        self.assertEqual(parser.date('2001/1/1'), '2001/01/01')
        self.assertEqual(parser.date('2010/02/23'), '2010/02/23')
        self.assertEqual(parser.date('2008/02/29'), '2008/02/29')
        self.assertEqual(parser.date('2001/1/1'), '2001/01/01')
        self.assertEqual(parser.date('58427'), '2018/11/05')
        
        self.assertRaises(argparse.ArgumentTypeError, parser.date, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.date, '2018/1/32')
        self.assertRaises(argparse.ArgumentTypeError, parser.date, '2018/13/1')
        self.assertRaises(argparse.ArgumentTypeError, parser.date, '2018/2/29')
        
    def test_mjd(self):
        """Test the parser.mjd conversion function."""
        
        self.assertEqual(parser.mjd('56274'), 56274)
        self.assertEqual(parser.mjd('2018/11/5'), 58427)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.mjd, 'c')
        self.assertRaises(argparse.ArgumentTypeError, parser.mjd, '2018/1/32')
        self.assertRaises(argparse.ArgumentTypeError, parser.mjd, '2018/13/1')
        
    def test_time(self):
        """Test the parser.time conversion function."""
        
        self.assertEqual(parser.time('1:23:45.6'), '1:23:45.600000')
        self.assertEqual(parser.time('1:23:45'),   '1:23:45.000000')
        self.assertEqual(parser.time('23:59:59'),  '23:59:59.000000')
        
        self.assertRaises(argparse.ArgumentTypeError, parser.time, '24:58:23.5')
        self.assertRaises(argparse.ArgumentTypeError, parser.time, '-1:58:23.5')
        
    def test_hours(self):
        """Test the parser.hours conversion function."""
        
        self.assertAlmostEqual(parser.hours('5'), 5*numpy.pi/12, 15)
        self.assertAlmostEqual(parser.hours('5:01:02.3'), (5+1/60.+2.3/3600.)*numpy.pi/12, 15)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.hours, 'c')
        
    def test_csv_hours_list(self):
        """Test the parser.hours conversion function."""
        
        eh = ephem.hours
        
        cin, cout = '5,', [eh('5:00:00'),]
        cin = parser.csv_hours_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertAlmostEqual(cin[i], cout[i], 15)
            
        cin, cout = '5,6.5,7:10,8:10:20.1', [eh('5:00:00'),eh('6:30:00'),eh('7:10:00'),eh('8:10:20.1')]
        cin = parser.csv_hours_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertAlmostEqual(cin[i], cout[i], 15)
            
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_hours_list, '5,6,c')
        
    def test_degrees(self):
        """Test the parser.degrees conversion function."""
        
        self.assertAlmostEqual(parser.degrees('5'), 5*numpy.pi/180, 15)
        self.assertAlmostEqual(parser.degrees('5:01:02.3'), (5+1/60.+2.3/3600.)*numpy.pi/180, 15)
        
        self.assertRaises(argparse.ArgumentTypeError, parser.degrees, 'c')
        
    def test_csv_degrees_list(self):
        """Test the parser.degrees conversion function."""
        
        ed = ephem.degrees
        
        cin, cout = '5,', [ed('5:00:00'),]
        cin = parser.csv_degrees_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertAlmostEqual(cin[i], cout[i], 15)
            
        cin, cout = '5,-6.5,7:10,8:10:20.1', [ed('5:00:00'),ed('-6:30:00'),ed('7:10:00'),ed('8:10:20.1')]
        cin = parser.csv_degrees_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertAlmostEqual(cin[i], cout[i], 15)
            
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_degrees_list, '5,6,c')
        
    def test_csv_int_list(self):
        """Test the parser.csv_int_list conversion function."""
        
        cin, cout = '*', 'all'
        self.assertEqual(parser.csv_int_list(cin), cout)
        
        cin, cout = '', 'none'
        self.assertEqual(parser.csv_int_list(cin), cout)
        
        cin, cout = '1', [1,]
        cin = parser.csv_int_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = '1,2,', [1,2]
        cin = parser.csv_int_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = '10~12,1,2,5~8', [10,11,12,1,2,5,6,7,8]
        cin = parser.csv_int_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_int_list, '1,2,!5')
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_int_list, '1,2,c')
        
    def test_csv_baseline_list(self):
        """Test the parser.csv_baseline_list conversion function."""
        
        cin, cout = '*', 'all'
        self.assertEqual(parser.csv_baseline_list(cin), cout)
        
        cin, cout = '', 'none'
        self.assertEqual(parser.csv_baseline_list(cin), cout)
        
        cin, cout = '1-2,', [(1,2),]
        cin = parser.csv_baseline_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            for j in range(2):
                self.assertEqual(cin[i][j], cout[i][j])
                
        cin, cout = '10-12,1-2,1-3~4,2~3-3~4', [(10,12),(1,2),(1,3),(1,4),(2,3),(2,4),(3,3),(3,4)]
        cin = parser.csv_baseline_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            for j in range(2):
                self.assertEqual(cin[i][j], cout[i][j])
            
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_baseline_list, '1')
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_baseline_list, '1-2,2')
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_baseline_list, '1-2,1-!2')
        
    def test_csv_hostname_list(self):
        """Test the parser.csv_hostname_list conversion function."""
        
        cin, cout = 'host1', ['host1',]
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = 'host-1,host-2', ['host-1','host-2']
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = 'host1,host2,host5~6', ['host1','host2','host5','host6']
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_hostname_list, 'c&')
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_hostname_list, '(9')
        
        cin, cout = '192.168.1.1', ['192.168.1.1',]
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = '192.168.1.1,192.168.1.2,192.168.1.5~6', ['192.168.1.1','192.168.1.2','192.168.1.5','192.168.1.6']
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = '192.168.5~6.1', ['192.168.5.1','192.168.6.1']
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = '192.168~169.1.1~2,', ['192.168.1.1','192.168.1.2','192.169.1.1','192.169.1.2']
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        cin, cout = '192~193.168.5~6.1', ['192.168.5.1','192.168.6.1','193.168.5.1','193.168.6.1']
        cin = parser.csv_hostname_list(cin)
        self.assertEqual(len(cin), len(cout))
        for i in range(len(cout)):
            self.assertEqual(cin[i], cout[i])
            
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_baseline_list, '192.168.1')
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_baseline_list, '192.168.1.c')
        self.assertRaises(argparse.ArgumentTypeError, parser.csv_baseline_list, '192.168.1.1.1')


class parser_test_suite(unittest.TestSuite):
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(parser_tests))
                

if __name__ == '__main__':
    unittest.main()
