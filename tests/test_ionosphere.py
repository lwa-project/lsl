"""
Unit test for the lsl.misc.ionosphere module.
"""

import os
import unittest
import numpy as np

from lsl.common.stations import lwa1
from lsl.misc import ionosphere
import lsl.testing


__version__  = "0.3"
__author__    = "Jayce Dowell"


class ionosphere_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.misc.ionosphere
    module."""
    
    def test_magnetic_field(self):
        """Test the IGRF model"""
        
        """
        ******************************************************
        *              IGRF SYNTHESIS PROGRAM                *
        *                                                    *
        * A program for the computation of geomagnetic       *
        * field elements from the International Geomagnetic  *
        * Reference Field (14th generation) as revised in    *
        * December 2024 by the IAGA Working Group V-MOD.     *
        *                                                    *
        * It is valid for dates from 1900.0 to 2030.0,       *
        * values up to 2035.0 will be computed but with      *
        * reduced accuracy. Values for dates before 1945.0   *
        * and after 2020.0 are non-definitive, otherwise the *
        * values are definitive.                             *
        *                                                    *
        * Susan Macmillan, William Brown, Ciaran Beggan      *
        *                          British Geological Survey *
        *                           IAGA Working Group V-MOD *
        ******************************************************

        Enter name of output file (30 characters maximum)
        or press "Return" for output to screen

        Enter value for coordinate system:
        1 - geodetic (shape of Earth is approximated by a spheroid)
        2 - geocentric (shape of Earth is approximated by a sphere)
        1
        Choose an option:
        1 - values at one or more locations & dates
        2 - values at yearly intervals at one location
        3 - values on a latitude/longitude grid at one date
        1
        Enter value for format of latitudes and longitudes:
        1 - in degrees & minutes
        2 - in decimal degrees
        2
        Enter date in years A.D.
        2025.0
        Enter altitude in km
        50
        Enter latitude & longitude in decimal degrees
        34.0 -107.0
        Enter place name (20 characters maximum)

        2025.000 Lat  34.000 geodetic   Long -107.000    50.000 km                     
                    D =    7 deg  53 min    SV =      -5 min/yr
                    I =   60 deg  50 min    SV =      -2 min/yr
                    H =   22706 nT          SV =     -32 nT/yr
                    X =   22492 nT          SV =     -27 nT/yr
                    Y =    3112 nT          SV =     -36 nT/yr
                    Z =   40688 nT          SV =    -115 nT/yr
                    F =   46595 nT          SV =    -116 nT/yr

        Do you want values for another date & position? (y/n)
        n
        D is declination (+ve east)
        I is inclination (+ve down)
        H is horizontal intensity
        X is north component
        Y is east component
        Z is vertical component (+ve down)
        F is total intensity

        SV is secular variation (annual rate of change)
        """
        
        bx, by, bz = ionosphere.get_magnetic_field(34.0, -107.0, 50e3, mjd=60676)
        self.assertAlmostEqual(bx, 22492, 0)
        self.assertAlmostEqual(by, 3112, 0)
        self.assertAlmostEqual(bz, -40688, 0)
        
    def test_magnetic_declination(self):
        """Test computing the magnetic declination"""
        
        bx, by, bz = ionosphere.get_magnetic_field(34.0, -107.0, 50e3, mjd=60676)
        dec = ionosphere.compute_magnetic_declination(bx, by, bz)
        d = int(dec)
        m = int(round((dec-d)*60))
        self.assertEqual(d, 7)
        self.assertEqual(m, 53)
        
    def test_magnetic_inclination(self):
        """Test computing the magnetic inclination"""
        
        bx, by, bz = ionosphere.get_magnetic_field(34.0, -107.0, 50e3, mjd=60676)
        inc = ionosphere.compute_magnetic_inclination(bx, by, bz)
        d = int(inc)
        m = int(round((inc-d)*60))
        self.assertEqual(d, 60)
        self.assertEqual(m, 50)
        
    def test_pierce_point(self):
        """Test the ionospheric pierce point"""
        
        az, el, h = 45.0, 45.0, 350e3   # m
        pos = ionosphere.get_ionospheric_pierce_point(lwa1, az, el, height=h)
        
        # https://gssc.esa.int/navipedia/index.php/Klobuchar_Ionospheric_Model
        psi = 0.0137/(el/180.0+0.11) - 0.022
        lat = lwa1.lat/np.pi + psi*np.cos(az*np.pi/180)
        lat *= 180.0
        lng = lwa1.lon/np.pi + psi*np.sin(az*np.pi/180)/np.cos(lat*np.pi/180)
        lng *= 180.0
        
        self.assertAlmostEqual(pos[0], lat, 1)
        self.assertAlmostEqual(pos[1], lng, 1)
        
    def test_tec_value(self):
        """Test retrieving the TEC value at a particular location"""
        
        with self.subTest(service='IGS'):
            """
            LSL 1.2.5
            
            Python 2.7.17 (default, Nov  7 2019, 10:07:09) 
            [GCC 7.4.0] on linux2
            Type "help", "copyright", "credits" or "license" for more information.
            >>> from lsl.misc import ionosphere
            >>> ionosphere.getTECValue(58215, lat=34.0, lng=-107.0, includeRMS=True, type='IGS')
            (array([[ 14.87999992]]), array([[ 0.73999999]]))
            """
            
            tec, rms = ionosphere.get_tec_value(58215, lat=34.0, lng=-107.0, include_rms=True, type='IGS')
            self.assertAlmostEqual(tec[0][0], 14.87999992, 6)
            self.assertAlmostEqual(rms[0][0],  0.73999999, 6)
            
        with self.subTest(service='JPL'):
            """
            LSL 1.2.5
            
            Python 2.7.17 (default, Nov  7 2019, 10:07:09) 
            [GCC 7.4.0] on linux2
            Type "help", "copyright", "credits" or "license" for more information.
            >>> from lsl.misc import ionosphere
            >>> ionosphere.getTECValue(58215, lat=34.0, lng=-107.0, includeRMS=True, type='JPL')
            (array([[ 15.66000019]]), array([[ 2.5]]))
            """
            
            tec, rms = ionosphere.get_tec_value(58215, lat=34.0, lng=-107.0, include_rms=True, type='JPL')
            self.assertAlmostEqual(tec[0][0], 15.66000019, 6)
            self.assertAlmostEqual(rms[0][0],  2.50000000, 6)
            
        with self.subTest(service='EMR'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(58215, lat=34.0, lng=-107.0, include_rms=True, type='EMR')
            # self.assertAlmostEqual(tec[0][0], 15.66000019, 6)
            # self.assertAlmostEqual(rms[0][0],  2.50000000, 6)
            
        with self.subTest(service='CODE'):
            """
            LSL 1.2.5
            
            Python 2.7.17 (default, Nov  7 2019, 10:07:09) 
            [GCC 7.4.0] on linux2
            Type "help", "copyright", "credits" or "license" for more information.
            >>> from lsl.misc import ionosphere
            ionosphere.getTECValue(58215, lat=34.0, lng=-107.0, includeRMS=True, type='CODE')
            (array([[ 14.14000015]]), array([[ 0.76]]))
            """
            
            tec, rms = ionosphere.get_tec_value(58215, lat=34.0, lng=-107.0, include_rms=True, type='CODE')
            self.assertAlmostEqual(tec[0][0], 14.14000015, 6)
            self.assertAlmostEqual(rms[0][0],  0.76000000, 6)
            
        with self.subTest(service='UQR'):
            """
            LSL 1.2.5
            
            Python 2.7.17 (default, Nov  7 2019, 10:07:09) 
            [GCC 7.4.0] on linux2
            Type "help", "copyright", "credits" or "license" for more information.
            >>> from lsl.misc import ionosphere
            >>> ionosphere.getTECValue(58215, lat=34.0, lng=-107.0, includeRMS=True, type='UQR')
            (array([[ 13.25999996]]), array([[ 6.91600008]]))
            """
            
            tec, rms = ionosphere.get_tec_value(58215, lat=34.0, lng=-107.0, include_rms=True, type='UQR')
            self.assertAlmostEqual(tec[0][0], 13.25999996, 6)
            self.assertAlmostEqual(rms[0][0],  6.91600008, 6)
            
        with self.subTest(service='USTEC'):
            """
            LSL 1.2.5
            
            Python 2.7.17 (default, Nov  7 2019, 10:07:09) 
            [GCC 7.4.0] on linux2
            Type "help", "copyright", "credits" or "license" for more information.
            >>> from lsl.misc import ionosphere
            >>> ionosphere.getTECValue(58415, lat=34.0, lng=-107.0, includeRMS=True, type='USTEC')
            (array([[ 17.29999924]]), array([[ 2.5999999]]))
            """
            
            tec, rms = ionosphere.get_tec_value(58415, lat=34.0, lng=-107.0, include_rms=True, type='USTEC')
            self.assertAlmostEqual(tec[0][0], 17.29999924, 6)
            self.assertAlmostEqual(rms[0][0],  2.5999999, 6)
            
        # with self.subTest(service='GloTEC'):
        #     """
        #     $ less glotec_icao_20240927T092500Z.geojson
        #     ...
        #     {
        #         "type": "Feature",
        #         "geometry": {
        #             "type": "Point",
        #             "coordinates": [
        #                 -107.5,
        #                 33.75
        #             ]
        #         },
        #         "properties": {
        #             "tec": 15.556710199197358,
        #             "anomaly": -2.8757901895049898,
        #             "quality_flag": 5
        #         }
        #     },
        #     ...
        #     """
        # 
        #     tec, rms = ionosphere.get_tec_value(60580.39236111111, lat=33.75, lng=-107.5, include_rms=True, type='GloTEC')
        #     self.assertAlmostEqual(tec[0][0], 15.5567102, 6)
        #     self.assertAlmostEqual(rms[0][0],  2.8757901, 6)
            
    def test_tec_value_lpn(self):
        """Test retrieving the TEC value at a particular location in the era of long product names"""
        
        with self.subTest(service='IGS'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(60025, lat=34.0, lng=-107.0, include_rms=True, type='IGS')
            # self.assertAlmostEqual(tec[0][0], 14.87999992, 6)
            # self.assertAlmostEqual(rms[0][0],  0.73999999, 6)
            
        with self.subTest(service='JPL'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(60025, lat=34.0, lng=-107.0, include_rms=True, type='JPL')
            # self.assertAlmostEqual(tec[0][0], 15.66000019, 6)
            # self.assertAlmostEqual(rms[0][0],  2.50000000, 6)
            
        with self.subTest(service='EMR'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(60025, lat=34.0, lng=-107.0, include_rms=True, type='EMR')
            # self.assertAlmostEqual(tec[0][0], 15.66000019, 6)
            # self.assertAlmostEqual(rms[0][0],  2.50000000, 6)
            
        with self.subTest(service='CODE'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(60025, lat=34.0, lng=-107.0, include_rms=True, type='CODE')
            # self.assertAlmostEqual(tec[0][0], 14.14000015, 6)
            # self.assertAlmostEqual(rms[0][0],  0.76000000, 6)
            
        with self.subTest(service='UQR'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(60025, lat=34.0, lng=-107.0, include_rms=True, type='UQR')
            # self.assertAlmostEqual(tec[0][0], 13.25999996, 6)
            # self.assertAlmostEqual(rms[0][0],  6.91600008, 6)
            
        with self.subTest(service='USTEC'):
            """
            ...
            """
            
            tec, rms = ionosphere.get_tec_value(60025, lat=34.0, lng=-107.0, include_rms=True, type='USTEC')
            # self.assertAlmostEqual(tec[0][0], 17.29999924, 6)
            # self.assertAlmostEqual(rms[0][0],  2.5999999, 6)


class ionosphere_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.misc.ionosphere
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(ionosphere_tests))        


if __name__ == '__main__':
    unittest.main()
    
