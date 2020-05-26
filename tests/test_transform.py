"""
Unit test for the lsl.transform module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import math
import ephem
import unittest

from lsl import transform
from lsl.astro import DJD_OFFSET
from lsl.common.stations import lwa1


__version__  = "0.1"
__author__    = "Jayce Dowell"


class transform_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.transform
    module."""
    
    def test_time_init(self):
        """Test the transform.Time constructor."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        str(t0)
        repr(t0)
        
        self.assertEqual(t0.utc_str, '2013-01-08 01:23:45.000')
        self.assertEqual(t0.utc_mcs, (56300, 5025000))
        self.assertAlmostEqual(t0.utc_dp/196e6, 1357608225.0, 3)
        
        t1 = transform.Time((56300, 5026000), format='MCS')
        self.assertEqual(t1.utc_str, '2013-01-08 01:23:46.000')
        self.assertAlmostEqual(t1.utc_dp/196e6, 1357608226.0, 3)
        
        t2 = transform.Time(1357608224.0, format='TIMET')
        self.assertEqual(t2.utc_str, '2013-01-08 01:23:44.000')
        self.assertEqual(t2.utc_mcs, (56300, 5024000))
        
        t3 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        t3.utc_mcs = (56300, 5024000)
        self.assertEqual(t3.utc_str, '2013-01-08 01:23:44.000')
        
    def test_time_comparisons(self):
        """Test ordering transform.Time instances."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        t1 = transform.Time((56300, 5026000), format='MCS')
        t2 = transform.Time(1357608224.0, format='TIMET')
        
        self.assertTrue(t0 < t1)
        self.assertTrue(t0 > t2)
        self.assertTrue(t2 != t1)
        
    def test_planetaryposition_init(self):
        """Test the transform.PlanetaryPosition constructor."""
        
        p0 = transform.PlanetaryPosition('Jupiter')
        p1 = transform.PlanetaryPosition('Sun')
        str(p0)
        repr(p0)
        
    def test_planetaryposition_saturn(self):
        """Test the location of Saturn."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        sat = ephem.Saturn()
        sat.compute(obs)
        
        p0 = transform.PlanetaryPosition('Saturn')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sat.g_ra *180.0/math.pi, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sat.g_dec*180.0/math.pi, 4)
        
    def test_planetaryposition_jupiter(self):
        """Test the location of Jupiter."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        jove = ephem.Jupiter()
        jove.compute(obs)
        
        p0 = transform.PlanetaryPosition('Jupiter')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], jove.g_ra *180.0/math.pi, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], jove.g_dec*180.0/math.pi, 4)
        
    def test_planetaryposition_mars(self):
        """Test the location of Mars."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        mars = ephem.Mars()
        mars.compute(obs)
        
        p0 = transform.PlanetaryPosition('Mars')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], mars.g_ra *180.0/math.pi, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], mars.g_dec*180.0/math.pi, 4)
        
    def test_planetaryposition_venus(self):
        """Test the location of Venus."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        venu = ephem.Venus()
        venu.compute(obs)
        
        p0 = transform.PlanetaryPosition('Venus')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], venu.g_ra *180.0/math.pi, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], venu.g_dec*180.0/math.pi, 4)
        
    def test_planetaryposition_sun(self):
        """Test the location of the Sun."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        sol = ephem.Sun()
        sol.compute(obs)
        
        p0 = transform.PlanetaryPosition('Sun')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sol.g_ra *180.0/math.pi, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sol.g_dec*180.0/math.pi, 4)
        
    def test_planetaryposition_moon(self):
        """Test the location of the Moon."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        lun = ephem.Moon()
        lun.compute(obs)
        
        p0 = transform.PlanetaryPosition('Moon')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], lun.g_ra *180.0/math.pi, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], lun.g_dec*180.0/math.pi, 4)
        
    def test_geographicalposition_init(self):
        """Test the transform.GeographicalPosition constructor."""
        
        lon = lwa1.long * 180.0/math.pi
        lat = lwa1.lat  * 180.0/math.pi
        elv = lwa1.elev
        
        g0 = transform.GeographicalPosition([lon,lat,elv])
        str(g0)
        repr(g0)
        
    def test_geographicalposition_ecef(self):
        """Test the tranform.GeographicalPosition EC-EF transform."""
        
        lon = lwa1.long * 180.0/math.pi
        lat = lwa1.lat  * 180.0/math.pi
        elv = lwa1.elev
        
        g0 = transform.GeographicalPosition([lon,lat,elv])
        self.assertAlmostEqual(g0.ecef[0], lwa1.geocentric_location[0], 6)
        self.assertAlmostEqual(g0.ecef[1], lwa1.geocentric_location[1], 6)
        self.assertAlmostEqual(g0.ecef[2], lwa1.geocentric_location[2], 6)
        
    def test_geographicalposition_lst(self):
        """Test the tranform.GeographicalPosition sidereal time."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        lon = lwa1.long * 180.0/math.pi
        lat = lwa1.lat  * 180.0/math.pi
        elv = lwa1.elev
        obs = lwa1.get_observer()
        
        g0 = transform.GeographicalPosition([lon,lat,elv])
        
        # The astro.get_apparent_sidereal_time() function doesn't care about
        # elevation
        obs.elev = 0.0
        
        obs.date = t0.utc_str
        self.assertAlmostEqual(g0.sidereal(t0), obs.sidereal_time()*12.0/math.pi, 4)
        
    def test_pointingdirection_init(self):
        """Test the transform.PointingDirection constructor."""
        
        lon = lwa1.long * 180.0/math.pi
        lat = lwa1.lat  * 180.0/math.pi
        elv = lwa1.elev
        
        g0 = transform.GeographicalPosition([lon,lat,elv])
        p0 = transform.PlanetaryPosition('Jupiter')
        str(p0)
        repr(p0)
        
        d0 = transform.PointingDirection(p0, g0)
        
    def test_pointingdirection_azalt(self):
        """Test the transform.PointingDirection az/alt transform."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        lon = lwa1.long * 180.0/math.pi
        lat = lwa1.lat  * 180.0/math.pi
        elv = lwa1.elev
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        jove = ephem.Jupiter()
        jove.compute(obs)
        
        g0 = transform.GeographicalPosition([lon,lat,elv])
        p0 = transform.PlanetaryPosition('Jupiter')
        d0 = transform.PointingDirection(p0, g0)
        str(d0)
        repr(d0)
        
        self.assertAlmostEqual(d0.hrz(t0)[0], jove.az *180.0/math.pi, 0)
        self.assertAlmostEqual(d0.hrz(t0)[1], jove.alt*180.0/math.pi, 0)
        
    def test_pointingdirection_rst(self):
        """Test the transform.PointingDirection az/alt transform."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        lon = lwa1.long * 180.0/math.pi
        lat = lwa1.lat  * 180.0/math.pi
        elv = lwa1.elev
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        jove = ephem.Jupiter()
        jove.compute(obs)
        
        g0 = transform.GeographicalPosition([lon,lat,elv])
        p0 = transform.PlanetaryPosition('Jupiter')
        d0 = transform.PointingDirection(p0, g0)
        
        rst = d0.rst(t0)
        self.assertAlmostEqual(rst.rise, obs.next_rising(jove)*1.0+DJD_OFFSET, 2)
        self.assertAlmostEqual(rst.transit, obs.next_transit(jove)*1.0+DJD_OFFSET, 2)
        self.assertAlmostEqual(rst.set, obs.next_setting(jove)*1.0+DJD_OFFSET, 2)


class transform_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.transform units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(transform_tests)) 


if __name__ == '__main__':
    unittest.main()
