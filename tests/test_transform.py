"""
Unit test for the lsl.transform module.
"""

import os
import math
import ephem
import unittest

from astropy.time import Time as AstroTime
from astropy.coordinates import SkyCoord, PrecessedGeocentric, GeocentricTrueEcliptic, Galactic, FK4

from lsl import transform
from lsl.astro import DJD_OFFSET, equ_posn
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
        
        t4 = transform.Time(AstroTime('2013-01-08 01:23:45.000', format='iso', scale='utc'),
                            format='ASTROPY')
        t4.utc_mcs = (56300, 5024000)
        self.assertEqual(t4.utc_str, '2013-01-08 01:23:44.000')
        
    def test_time_comparisons(self):
        """Test ordering transform.Time instances."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        t1 = transform.Time((56300, 5026000), format='MCS')
        t2 = transform.Time(1357608224.0, format='TIMET')
        
        self.assertTrue(t0 < t1)
        self.assertTrue(t0 > t2)
        self.assertTrue(t2 != t1)
        
    def test_celestialposition_init(self):
        """Test the transform.CelestialPosition constructor."""
        
        sc = SkyCoord('12h34m56.7s', '+1d23m45.6s', frame='fk5', equinox='J2000')
        equ = equ_posn.from_astropy(sc)
        
        p0 = transform.CelestialPosition(equ, format='EQU', epoch='J2000', name='test')
        
    def test_celestialposition_convert(self):
        """Test conversions within transform.CelestialPosition."""
        
        sc = SkyCoord('12h34m56.7s', '+1d23m45.6s', frame='fk5', equinox='J2000')
        equ = equ_posn.from_astropy(sc)
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        p0 = transform.CelestialPosition(equ, format='EQU', epoch='J2000', name='test')
        
        sc1 = sc.transform_to(FK4(equinox='B1950'))
        
        self.assertAlmostEqual(p0.b1950_equ[0], sc1.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.b1950_equ[1], sc1.dec.to('deg').value, 4)
        
        sc1 = sc.transform_to(Galactic())
        
        self.assertAlmostEqual(p0.j2000_gal[0], sc1.l.to('deg').value, 4)
        self.assertAlmostEqual(p0.j2000_gal[1], sc1.b.to('deg').value, 4)
        
        sc1 = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc1.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc1.lat.to('deg').value, 4)
        
    def test_planetaryposition_init(self):
        """Test the transform.PlanetaryPosition constructor."""
        
        p0 = transform.PlanetaryPosition('Jupiter')
        p1 = transform.PlanetaryPosition('Sun')
        str(p0)
        repr(p0)
        
    def test_planetaryposition_saturn(self):
        """Test the location of Saturn."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        # Test set generated at https://ssd.jpl.nasa.gov/horizons/app.html#/
        sc = SkyCoord('14h32m55.72s', '-12d31m50.8s', frame='icrs')
        sc = sc.transform_to(PrecessedGeocentric(equinox=t0.astropy, obstime=t0.astropy))
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        sat = ephem.Saturn()
        sat.compute(obs)
        
        p0 = transform.PlanetaryPosition('Saturn')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sc.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sc.dec.to('deg').value, 4)
        
        sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc.lat.to('deg').value, 4)
        
    def test_planetaryposition_jupiter(self):
        """Test the location of Jupiter."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        # Test set generated at https://ssd.jpl.nasa.gov/horizons/app.html#/
        sc = SkyCoord('4h21m06.91s', '20d48m13.2s', frame='icrs')
        sc = sc.transform_to(PrecessedGeocentric(equinox=t0.astropy, obstime=t0.astropy))
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        jove = ephem.Jupiter()
        jove.compute(obs)
        
        p0 = transform.PlanetaryPosition('Jupiter')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sc.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sc.dec.to('deg').value, 4)
        
        sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc.lat.to('deg').value, 4)
        
    def test_planetaryposition_mars(self):
        """Test the location of Mars."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        # Test set generated at https://ssd.jpl.nasa.gov/horizons/app.html#/
        sc = SkyCoord('20h51m16.93s',  '-18d49m05.7s', frame='icrs')
        sc = sc.transform_to(PrecessedGeocentric(equinox=t0.astropy, obstime=t0.astropy))
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        mars = ephem.Mars()
        mars.compute(obs)
        
        p0 = transform.PlanetaryPosition('Mars')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sc.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sc.dec.to('deg').value, 4)
        
        sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc.lat.to('deg').value, 4)
        
    def test_planetaryposition_venus(self):
        """Test the location of Venus."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        # Test set generated at https://ssd.jpl.nasa.gov/horizons/app.html#/
        sc = SkyCoord('17h53m08.27s', '-23d01m37.2s', frame='icrs')
        sc = sc.transform_to(PrecessedGeocentric(equinox=t0.astropy, obstime=t0.astropy))
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        venu = ephem.Venus()
        venu.compute(obs)
        
        p0 = transform.PlanetaryPosition('Venus')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sc.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sc.dec.to('deg').value, 4)
        
        sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc.lat.to('deg').value, 4)
        
    def test_planetaryposition_sun(self):
        """Test the location of the Sun."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        # Test set generated at https://ssd.jpl.nasa.gov/horizons/app.html#/
        sc = SkyCoord('19h16m54.31s', '-22d15m40.1s', frame='icrs')
        sc = sc.transform_to(PrecessedGeocentric(equinox=t0.astropy, obstime=t0.astropy))
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        sol = ephem.Sun()
        sol.compute(obs)
        
        p0 = transform.PlanetaryPosition('Sun')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sc.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sc.dec.to('deg').value, 4)
        
        sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc.lat.to('deg').value, 4)
        
    def test_planetaryposition_moon(self):
        """Test the location of the Moon."""
        
        t0 = transform.Time('2013-01-08 01:23:45.000', format='STR')
        
        # Test set generated at https://ssd.jpl.nasa.gov/horizons/app.html#/
        sc = SkyCoord('15h32m48.21s', '-19d03m23.0s', frame='icrs')
        sc = sc.transform_to(PrecessedGeocentric(equinox=t0.astropy, obstime=t0.astropy))
        
        obs = lwa1.get_observer()
        obs.date = t0.utc_str
        lun = ephem.Moon()
        lun.compute(obs)
        
        p0 = transform.PlanetaryPosition('Moon')
        
        self.assertAlmostEqual(p0.apparent_equ(t0)[0], sc.ra.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_equ(t0)[1], sc.dec.to('deg').value, 4)
        
        sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t0.astropy))
        
        self.assertAlmostEqual(p0.apparent_ecl(t0)[0], sc.lon.to('deg').value, 4)
        self.assertAlmostEqual(p0.apparent_ecl(t0)[1], sc.lat.to('deg').value, 4)
        
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
        self.assertAlmostEqual(g0.sidereal(t0), obs.sidereal_time()*12.0/math.pi, 3)
        
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
