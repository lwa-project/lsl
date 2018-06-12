# -*- coding: utf-8 -*-

"""Unit test for lsl.astro module."""


import unittest
import math
import pickle
import operator

from lsl import astro


__revision__  = "$Rev$"
__version__   = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class astro_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.astro
	module."""

	def setUp(self):
		"""
		Setup for running unit tests.
		"""
		
		lng = astro.dms(True,  90, 0, 0)
		lat = astro.dms(False, 40, 0, 0)
		self.geo = astro.geo_posn(lng, lat)
		
		self.times = (\
			astro.date(2001,  1, 22, 0, 0, 0),
			astro.date(2001,  4, 16, 0, 0, 0),
			astro.date(2001,  8,  2, 0, 0, 0),
			astro.date(2001, 10, 10, 0, 0, 0),
			astro.date(1979,  2, 28, 0, 0, 0))
			
	def test_hms_init(self):
		"""Test astro.hms constructor."""
		
		h = astro.hms(3, 42, 20.3242)
		self.assertEqual(h.hours, 3)
		self.assertEqual(h.minutes, 42)
		self.assertAlmostEqual(h.seconds, 20.3242)
			
		self.assertRaises(ValueError, astro.hms, 30, 0, 0)
		self.assertRaises(ValueError, astro.hms, 0, 90, 0)
		self.assertRaises(ValueError, astro.hms, 0, 0, 100) 
		
	def test_hms_reduce(self):
		"""Test astro.hms.__reduce__() method."""
		
		h = astro.hms(3, 42, 20.3242)
		s = pickle.dumps(h)
		h = pickle.loads(s)
		self.assertEqual(h.hours, 3)   
		self.assertEqual(h.minutes, 42)
		self.assertAlmostEqual(h.seconds, 20.3242)
		
	def test_hms_cmp(self):
		"""Test astro.hms.__cmp__() method."""
		
		h1 = astro.hms(2, 43, 32.232)
		h2 = astro.hms(2, 43, 32.232)
		h3 = astro.hms(4, 23, 42.221)
		
		self.assertTrue(operator.eq(h1, h2))
		self.assertFalse(operator.eq(h1, h3))
		self.assertFalse(operator.ne(h1, h2))
		self.assertTrue(operator.ne(h1, h3))
		self.assertTrue(operator.lt(h1, h3))
		self.assertTrue(operator.gt(h3, h1))
		self.assertTrue(operator.eq(h1, h1.to_deg()))
		self.assertTrue(operator.eq(h2, h2.to_dms()))
		
		self.assertRaises(TypeError, operator.eq, h3, "string")
		
	def test_hms_to_deg(self):
		"""Test astro.hms_to_deg() function."""
		
		h = astro.hms(0, 0, 0)
		self.assertAlmostEqual(astro.hms_to_deg(h), 0.0)
		
		h = astro.hms(6, 0, 0)
		self.assertAlmostEqual(astro.hms_to_deg(h), 90.0)
		
		h = astro.hms(12, 0, 0)
		self.assertAlmostEqual(astro.hms_to_deg(h), 180.0)
		
		h = astro.hms(18, 0, 0)
		self.assertAlmostEqual(astro.hms_to_deg(h), 270.0)
		
	def test_deg_to_hms(self):
		"""Test astro.deg_to_hms() function."""
		
		h = astro.deg_to_hms(0.0)
		self.assertEqual(h.hours, 0)
		self.assertEqual(h.minutes, 0)
		self.assertAlmostEqual(h.seconds, 0.0)
		
		h = astro.deg_to_hms(90.0)
		self.assertEqual(h.hours, 6)
		self.assertEqual(h.minutes, 0)
		self.assertAlmostEqual(h.seconds, 0.0)
		
		h= astro.deg_to_hms(180.0)
		self.assertEqual(h.hours, 12)
		self.assertEqual(h.minutes, 0)
		self.assertAlmostEqual(h.seconds, 0.0)
		
		h = astro.deg_to_hms(270.0)
		self.assertEqual(h.hours, 18)
		self.assertEqual(h.minutes, 0)
		self.assertAlmostEqual(h.seconds, 0.0)
		
	def test_dms_init(self):
		"""Test astro.dms constructor."""
		
		d = astro.dms(True, 1, 30, 29.3245)
		self.assert_(d.neg)
		self.assertEqual(d.degrees, 1)
		self.assertEqual(d.minutes, 30)
		self.assertAlmostEqual(d.seconds, 29.3245)
		
		self.assertRaises(ValueError, astro.dms, False, 400, 0, 0)
		self.assertRaises(ValueError, astro.dms, False, 0, 80, 0)
		self.assertRaises(ValueError, astro.dms, False, 0, 0, 500)
		
	def test_dms_reduce(self):
		"""Test astro.dms.__reduce__() method."""
		
		d = astro.dms(True, 1, 30, 29.3245)
		s = pickle.dumps(d)
		d = pickle.loads(s)
		self.assert_(d.neg)
		self.assertEqual(d.degrees, 1)
		self.assertEqual(d.minutes, 30)
		self.assertAlmostEqual(d.seconds, 29.3245)
		
	def test_dms_cmp(self):
		"""Test astro.dms.__cmp__() method."""
		
		d1 = astro.dms(False, 2, 43, 32.232)
		d2 = astro.dms(False, 2, 43, 32.232)
		d3 = astro.dms(False, 4, 23, 42.221)
		
		self.assertTrue(operator.eq(d1, d2))
		self.assertFalse(operator.eq(d1, d3))
		self.assertFalse(operator.ne(d1, d2))
		self.assertTrue(operator.ne(d1, d3)) 
		self.assertTrue(operator.lt(d1, d3))
		self.assertTrue(operator.gt(d3, d1))
		self.assertTrue(operator.eq(d1, d1.to_deg())) 
		self.assertTrue(operator.eq(d2, d2.to_hms()))
		
		self.assertRaises(TypeError, operator.eq, d3, "string")
		
	def test_dms_to_deg(self):
		"""Test astro.dms_to_deg() function."""
		
		d = astro.dms(False, 0, 0, 0)
		self.assertAlmostEqual(astro.dms_to_deg(d), 0.0)
		
		d = astro.dms(False, 100, 0, 0)
		self.assertAlmostEqual(astro.dms_to_deg(d), 100.0)
		
		d = astro.dms(True, 100, 0, 0)
		self.assertAlmostEqual(astro.dms_to_deg(d), -100.0)
		
	def test_deg_to_dms(self):
		"""Test astro.deg_to_dms() function."""
		
		d = astro.deg_to_dms(0.0)
		self.assertEqual(d.neg, 0)
		self.assertEqual(d.degrees, 0)
		self.assertEqual(d.minutes, 0)
		self.assertAlmostEqual(d.seconds, 0.0)
		
		d = astro.deg_to_dms(100.0)
		self.assertEqual(d.neg, 0)
		self.assertEqual(d.degrees, 100)
		self.assertEqual(d.minutes, 0)
		self.assertAlmostEqual(d.seconds, 0.0)
		
		d = astro.deg_to_dms(-100.0)
		self.assertEqual(d.neg, 1)
		self.assertEqual(d.degrees, 100)
		self.assertEqual(d.minutes, 0)
		self.assertAlmostEqual(d.seconds, 0.0)
		
	def test_dir_cos(self):
		"""Test astro.dir_cos() function."""
		
		h = astro.hrz_posn(0.0, 0.0)
		(l,m,n) = astro.dir_cos(h)
		self.assertAlmostEqual(l, 0.0)
		self.assertAlmostEqual(m, 1.0)
		self.assertAlmostEqual(n, 0.0)
		
		h = astro.hrz_posn(90.0, 0.0)
		(l,m,n) = astro.dir_cos(h)
		self.assertAlmostEqual(l, 1.0)
		self.assertAlmostEqual(m, 0.0)
		self.assertAlmostEqual(n, 0.0)
		
		h = astro.hrz_posn(180.0, 0.0)
		(l,m,n) = astro.dir_cos(h)
		self.assertAlmostEqual(l, 0.0)
		self.assertAlmostEqual(m, -1.0)
		self.assertAlmostEqual(n, 0.0)
		
		h = astro.hrz_posn(270.0, 0.0)
		(l,m,n) = astro.dir_cos(h)
		self.assertAlmostEqual(l, -1.0)
		self.assertAlmostEqual(m, 0.0)
		self.assertAlmostEqual(n, 0.0)
		
		h = astro.hrz_posn(0.0, 90.0)
		(l,m,n) = astro.dir_cos(h)
		self.assertAlmostEqual(l, 0.0)
		self.assertAlmostEqual(m, 0.0)
		self.assertAlmostEqual(n, 1.0)
		
		h = astro.hrz_posn(0.0, -90.0)
		(l,m,n) = astro.dir_cos(h)
		self.assertAlmostEqual(l, 0.0)
		self.assertAlmostEqual(m, 0.0)
		self.assertAlmostEqual(n, -1.0)
		
	def test_range_degrees(self):
		"""Test astro.range_degrees() function."""
		
		self.assertAlmostEqual(astro.range_degrees(20.0), 20.0)
		self.assertAlmostEqual(astro.range_degrees(370.0), 10.0)
		self.assertAlmostEqual(astro.range_degrees(-10.0), 350.0)
		
	def test_range_hours(self):
		"""Test astro.range_hours() function."""
		
		self.assertAlmostEqual(astro.range_hours(10.0), 10.0)
		self.assertAlmostEqual(astro.range_hours(30.0), 6.0)
		self.assertAlmostEqual(astro.range_hours(-6.0), 18.0)
		
	def test_date_init(self):
		"""Test astro.date constructor."""
		
		d = astro.date(2000, 4, 28, 21, 49, 13.0238)
		self.assertEqual(d.years, 2000)
		self.assertEqual(d.months, 4)
		self.assertEqual(d.days, 28)
		self.assertEqual(d.hours, 21)
		self.assertEqual(d.minutes, 49)
		self.assertAlmostEqual(d.seconds, 13.0238)
		
		self.assertRaises(ValueError, astro.date, 2000, 30, 1, 0, 0, 0)
		self.assertRaises(ValueError, astro.date, 2000, 1, 48, 0, 0, 0)
		self.assertRaises(ValueError, astro.date, 2000, 1, 1, 39, 0, 0)
		self.assertRaises(ValueError, astro.date, 2000, 1, 1, 0, 69, 0)
		self.assertRaises(ValueError, astro.date, 2000, 1, 1, 0, 0, 73)
		
	def test_date_reduce(self):
		"""Test astro.date.__reduce__() method."""
		
		d = astro.date(2000, 4, 28, 21, 49, 13.0238)
		s = pickle.dumps(d)
		d = pickle.loads(s)
		self.assertEqual(d.years, 2000)
		self.assertEqual(d.months, 4)
		self.assertEqual(d.days, 28)
		self.assertEqual(d.hours, 21)
		self.assertEqual(d.minutes, 49)
		self.assertAlmostEqual(d.seconds, 13.0238)
		
	def test_date_cmp(self):
		"""Test astro.date.__cmp__() method."""
		
		d1 = astro.date(2004, 5, 9, 12, 34, 22)
		d2 = astro.date(2009, 2, 6,  9,  2,  6)
		d3 = astro.date(2009, 2, 6,  9,  2,  6)
		self.assertTrue(operator.lt(d1, d2))
		self.assertTrue(operator.gt(d2, d1))
		self.assertTrue(operator.eq(d2, d3))
		self.assertTrue(operator.ne(d1, d3))
		
	def test_date_load(self):
		"""Test astro.date.load() method."""
		
		d = astro.date()
		d.load('2000-4-28', '21:49:13.0238')
		self.assertEqual(d.years, 2000)
		self.assertEqual(d.months, 4)
		self.assertEqual(d.days, 28)
		self.assertEqual(d.hours, 21)
		self.assertEqual(d.minutes, 49)
		self.assertAlmostEqual(d.seconds, 13.0238)
		
	def test_zonedate_init(self):
		"""Test astro.zonedate constructor."""
		
		d = astro.zonedate(2000, 4, 28, 21, 49, 13.0238, -3600 * 5)
		self.assertEqual(d.years, 2000)
		self.assertEqual(d.months, 4)
		self.assertEqual(d.days, 28)
		self.assertEqual(d.hours, 21)
		self.assertEqual(d.minutes, 49)
		self.assertAlmostEqual(d.seconds, 13.0238)
		self.assertEqual(d.gmtoff, -3600 * 5)
		
		self.assertRaises(ValueError, astro.zonedate, 2000, 30, 1, 0, 0, 0, 0)
		self.assertRaises(ValueError, astro.zonedate, 2000, 1, 48, 0, 0, 0, 0)
		self.assertRaises(ValueError, astro.zonedate, 2000, 1, 1, 39, 0, 0, 0)
		self.assertRaises(ValueError, astro.zonedate, 2000, 1, 1, 0, 69, 0, 0)
		self.assertRaises(ValueError, astro.zonedate, 2000, 1, 1, 0, 0, 73, 0)
		
	def test_zonedate_reduce(self):
		"""Test astro.zonedate reduce method."""
		
		d = astro.zonedate(2000, 4, 28, 21, 49, 13.0238, -3600 * 5)
		s = pickle.dumps(d)
		d = pickle.loads(s)
		self.assertEqual(d.years, 2000)
		self.assertEqual(d.months, 4)
		self.assertEqual(d.days, 28)
		self.assertEqual(d.hours, 21)
		self.assertEqual(d.minutes, 49)
		self.assertAlmostEqual(d.seconds, 13.0238)
		self.assertEqual(d.gmtoff, -3600 * 5)
		
	def test_date_to_zonedate(self):
		"""Test astro.date_to_zonedate() function."""
		
		ACCURACY = 6
		
		d = astro.date(2000, 4, 28, 21, 49, 13.0238)
		z = astro.date_to_zonedate(d, -18000)
		self.assertEqual(z.years, 2000)
		self.assertEqual(z.months, 4)
		self.assertEqual(z.days, 28)
		self.assertEqual(z.hours, 16)
		self.assertEqual(z.minutes, 49)
		self.assertAlmostEqual(z.seconds, 13.0238, ACCURACY)
		self.assertEqual(z.gmtoff, -18000)
		
	def test_zonedate_to_date(self):
		"""Test astro.zonedate_to_date() function."""
		
		ACCURACY = 4
		
		z = astro.zonedate(2000, 4, 28, 16, 49, 13.0238, -18000)
		d = astro.zonedate_to_date(z)
		self.assertEqual(d.years, 2000)
		self.assertEqual(d.months, 4)
		self.assertEqual(d.days, 28)
		self.assertEqual(d.hours, 21)
		self.assertEqual(d.minutes, 49)
		self.assertAlmostEqual(d.seconds, 13.0238, ACCURACY)
		
	def test_get_date(self):
		"""Test astro.get_date() function."""
		
		year_out = (2001, 2001, 2001, 2001, 1979)
		month_out = (1, 4, 8, 10, 2)
		day_out = (22, 16, 2, 10, 28)
		
		iyear = iter(year_out)
		imonth = iter(month_out)
		iday = iter(day_out)
		for t in self.times:
			j = t.to_jd()
			d = astro.get_date(j)
			self.assertEqual(d.years, next(iyear))
			self.assertEqual(d.months, next(imonth))
			self.assertEqual(d.days, next(iday))
			self.assertEqual(d.hours, 0)
			self.assertEqual(d.minutes, 0)
			self.assertAlmostEqual(d.seconds, 0.0)
			
	def test_utc_to_tai(self):
		"""Test astro.utc_to_tai() function."""
		
		ACCURACY = 4
		
		jd_aa = (\
			astro.date(1972, 3, 1, 0, 0, 0).to_jd(),
			astro.date(1978, 6, 1, 0, 0, 0).to_jd(),
			astro.date(1983, 10, 1, 0, 0, 0).to_jd(),
			astro.date(1992, 2, 1, 0, 0, 0).to_jd())
			
		at_aa = (\
			astro.sec_to_jd(10.0), 
			astro.sec_to_jd(17.0), 
			astro.sec_to_jd(22.0),
			astro.sec_to_jd(27.0))
			
		iat = iter(at_aa)
		for t in jd_aa:
			tai = astro.utc_to_tai(t)
			self.assertAlmostEqual(tai - t, next(iat), ACCURACY) 
			
	def test_tai_to_utc(self):
		"""Test astro.tai_to_utc() function."""
		
		ACCURACY = 4
		
		jd_aa = (\
			astro.date(1972, 3, 1, 0, 0, 0).to_jd(),
			astro.date(1978, 6, 1, 0, 0, 0).to_jd(),
			astro.date(1983, 10, 1, 0, 0, 0).to_jd(),
			astro.date(1992, 2, 1, 0, 0, 0).to_jd())
			
		at_aa = (\
			-astro.sec_to_jd(10.0), 
			-astro.sec_to_jd(17.0), 
			-astro.sec_to_jd(22.0),
			-astro.sec_to_jd(27.0))
			
		iat = iter(at_aa)
		for t in jd_aa:
			utc = astro.tai_to_utc(t)
			self.assertAlmostEqual(utc - t, next(iat), ACCURACY)            
			
	def test_get_geo_from_rect(self):
		"""Test astro.get_geo_from_rect() function."""
		
		rect = astro.rect_posn(-1602196.60, -5042313.47, 3553971.51)
		geo = astro.get_geo_from_rect(rect)
		self.assertAlmostEqual(geo.lng, 252.3723299788652)
		self.assertAlmostEqual(geo.lat, 34.068899974626)
		self.assertAlmostEqual(geo.elv, 2126.9966418296099)
		
		rect = astro.rect_posn(1117051.64, -4848777.05, 3976930.95) 
		geo = astro.get_geo_from_rect(rect)
		self.assertAlmostEqual(geo.lng, 282.97333338968951)
		self.assertAlmostEqual(geo.lat, 38.821666658918112)
		self.assertAlmostEqual(geo.elv, 29.999061864800751)
		
	def test_get_rect_from_geo(self):
		"""Test astro.get_rect_from_geo() function."""
		
		geo = astro.geo_posn(252.3723299788652, 34.068899974626, 2126.9966418296099)
		rect = astro.get_rect_from_geo(geo)
		self.assertAlmostEqual(rect.X, -1602196.60)
		self.assertAlmostEqual(rect.Y, -5042313.47)
		self.assertAlmostEqual(rect.Z, 3553971.51)
		
		geo = astro.geo_posn(282.97333338968951, 38.821666658918112, 29.999061864800751)
		rect = astro.get_rect_from_geo(geo)
		self.assertAlmostEqual(rect.X, 1117051.64)
		self.assertAlmostEqual(rect.Y, -4848777.05)
		self.assertAlmostEqual(rect.Z, 3976930.95)          
		
	def test_get_julian_day(self):
		"""Test astro.get_julian_day() function and astro.date.to_jd() method."""
		
		jd_aa = (2451931.5, 2452015.5, 2452123.5, 2452192.5, 2443932.5)
		
		ijd = iter(jd_aa)
		for t in self.times:
			j = t.to_jd()
			self.assertAlmostEqual(j, next(ijd))
			
	def test_get_mean_sidereal_time(self):
		"""Test astro.get_mean_sidereal_time() function."""
		
		ACCURACY = 4
		
		gmst_aa = (\
			astro.hms( 8,  5, 39.1981).to_sec() / 3600.0, 
			astro.hms(13, 36, 49.8490).to_sec() / 3600.0, 
			astro.hms(20, 42, 37.8288).to_sec() / 3600.0,
			astro.hms( 1, 14, 40.1492).to_sec() / 3600.0,
			astro.hms(10, 28, 52.7540).to_sec() / 3600.0)
			
		igmst = iter(gmst_aa)
		for t in self.times:
			j = t.to_jd()
			st = astro.get_mean_sidereal_time(j)
			self.assertAlmostEqual(st, next(igmst), ACCURACY)
			
	def test_get_apparent_sidereal_time(self):
		"""Test astro.get_apparent_sidereal_time() function."""
		
		ACCURACY = 3
		
		gast_aa = (\
			astro.hms( 8,  5, 38.2531).to_sec() / 3600.0, 
			astro.hms(13, 36, 48.7661).to_sec() / 3600.0,
			astro.hms(20, 42, 36.8591).to_sec() / 3600.0,
			astro.hms( 1, 14, 39.0438).to_sec() / 3600.0,
			astro.hms(10, 28, 52.5980).to_sec() / 3600.0)
			
		igast = iter(gast_aa)
		for t in self.times:
			j = t.to_jd()
			st = astro.get_apparent_sidereal_time(j)
			self.assertAlmostEqual(st, next(igast), ACCURACY)
			
	def test_get_local_sidereal_time(self):
		"""Test astro.get_local_sidereal_time() function."""
		
		ACCURACY = 3
		
		lng = astro.deg_to_sec(self.geo.lng)
		
		gast_aa = (\
			astro.range_hours((astro.hms( 8,  5, 38.2531).to_sec() + lng) / 3600.0), 
			astro.range_hours((astro.hms(13, 36, 48.7661).to_sec() + lng) / 3600.0),
			astro.range_hours((astro.hms(20, 42, 36.8591).to_sec() + lng) / 3600.0),
			astro.range_hours((astro.hms( 1, 14, 39.0438).to_sec() + lng) / 3600.0),
			astro.range_hours((astro.hms(10, 28, 52.5980).to_sec() + lng) / 3600.0))
			
		igast = iter(gast_aa)
		for t in self.times:
			j = t.to_jd()
			st = astro.get_local_sidereal_time(self.geo.lng, j)
			self.assertAlmostEqual(st, next(igast), ACCURACY)
			
	def test_equ_posn_init(self):
		"""Test astro.equ_posn constructor."""
		
		e = astro.equ_posn(39.221, -24.543)
		self.assertAlmostEqual(e.ra, 39.221)
		self.assertAlmostEqual(e.dec, -24.543)
		
		ra = astro.deg_to_hms(39.221)
		dec = astro.deg_to_dms(-24.543)
		e = astro.equ_posn(ra, dec)
		self.assertAlmostEqual(e.ra, 39.221)
		self.assertAlmostEqual(e.dec, -24.543)
		
		self.assertRaises(ValueError, astro.equ_posn, 400, 0)
		self.assertRaises(ValueError, astro.equ_posn, -1, 0)
		self.assertRaises(ValueError, astro.equ_posn, 0, 100)
		self.assertRaises(ValueError, astro.equ_posn, 0, -100)
		
	def test_equ_posn_reduce(self):
		"""Test astro.equ_posn.__reduce__() method."""
		
		e = astro.equ_posn(39.221, -24.543)
		s = pickle.dumps(e)
		e = pickle.loads(s)
		self.assertAlmostEqual(e.ra, 39.221)
		self.assertAlmostEqual(e.dec, -24.543)
		
	def test_equ_posn_cmp(self):
		"""Test astro.equ_posn.__eq__() and astro.equ_posn.__ne__()
		methods."""
		
		e1 = astro.equ_posn(43, 12)
		e2 = astro.equ_posn(292, -21)
		e3 = astro.equ_posn(292, -21)
		self.assertTrue(operator.ne(e1, e2))
		self.assertFalse(operator.eq(e1, e2))
		self.assertTrue(operator.eq(e2, e3))
		self.assertFalse(operator.ne(e2, e3))
		
	def test_equ_posn_subscr(self):
		"""Test astro.equ_posn.__getitem__() and astro.equ_posn.__setitem__()
		methods."""
		
		e = astro.equ_posn(39.221, -24.543)
		self.assertAlmostEqual(e[0], 39.221)
		self.assertAlmostEqual(e[1], -24.543)
		
		e[0] = 132.043
		e[1] = 29.392
		self.assertAlmostEqual(e.ra, 132.043)
		self.assertAlmostEqual(e.dec, 29.392)
		
		self.assertRaises(ValueError, operator.itemgetter(2), e)
		self.assertRaises(ValueError, operator.itemgetter(-1), e)
        
	def test_equ_posn_format(self):
		"""Test astro.equ_posn.format() method."""
		
		e = astro.equ_posn(39.221, -24.543)
		(ra, dec) = e.format()
		self.assertAlmostEqual(ra.to_deg(), 39.221)
		self.assertAlmostEqual(dec.to_deg(), -24.543)
		
	def test_hrz_posn_init(self):
		"""Test astro.hrz_posn constructor."""
		
		h = astro.hrz_posn(39.221, 46.301)
		self.assertAlmostEqual(h.az, 39.221)
		self.assertAlmostEqual(h.alt, 46.301)
		
		self.assertRaises(ValueError, astro.hrz_posn, 400, 0)
		self.assertRaises(ValueError, astro.hrz_posn, -1, 0)
		self.assertRaises(ValueError, astro.hrz_posn, 0, 100)
		self.assertRaises(ValueError, astro.hrz_posn, 0, -100)
		
	def test_hrz_posn_reduce(self):
		"""Test astro.hrz_posn.__reduce__() method."""
		
		h = astro.hrz_posn(39.221, 46.301)
		s = pickle.dumps(h)
		h = pickle.loads(s)
		self.assertAlmostEqual(h.az, 39.221)
		self.assertAlmostEqual(h.alt, 46.301)
		
	def test_hrz_posn_cmp(self):
		"""Test astro.hrz_posn.__eq__() and astro.hrz_posn.__ne__() 
		methods."""
		
		h1 = astro.hrz_posn(300, 21)
		h2 = astro.hrz_posn(33, 43)
		h3 = astro.hrz_posn(33, 43)
		self.assertTrue(operator.ne(h1, h2))
		self.assertFalse(operator.eq(h1, h2))
		self.assertTrue(operator.eq(h2, h3))
		self.assertFalse(operator.ne(h2, h3))
		
	def test_hrz_posn_subscr(self):
		"""Test astro.hrz_posn.__getitem__() and astro.hrz_posn.__setitem__() 
		methods."""
		
		h = astro.hrz_posn(39.221, 46.301)
		self.assertAlmostEqual(h[0], 39.221)
		self.assertAlmostEqual(h[1], 46.301)
		h[0] = 319.293
		h[1] = -42.212
		self.assertAlmostEqual(h.az, 319.293)
		self.assertAlmostEqual(h.alt, -42.212)
		
		self.assertRaises(ValueError, operator.itemgetter(2), h)
		self.assertRaises(ValueError, operator.itemgetter(-1), h)
		
	def test_gal_posn_init(self):
		"""Test astro.gal_posn constructor."""
		
		g = astro.gal_posn(39.221, -24.543)
		self.assertAlmostEqual(g.l, 39.221)
		self.assertAlmostEqual(g.b, -24.543)
		
		l = astro.deg_to_dms(39.221)
		b = astro.deg_to_dms(-24.543)
		g = astro.gal_posn(l, b)
		self.assertAlmostEqual(g.l, 39.221)
		self.assertAlmostEqual(g.b, -24.543)
		
		self.assertRaises(ValueError, astro.gal_posn, 400, 0)
		self.assertRaises(ValueError, astro.gal_posn, -400, 0)
		self.assertRaises(ValueError, astro.gal_posn, 0, 100)
		self.assertRaises(ValueError, astro.gal_posn, 0, -100)
		
	def test_gal_posn_reduce(self):
		"""Test astro.gal_posn.__reduce__() method."""
		
		g = astro.gal_posn(39.221, -24.543)
		s = pickle.dumps(g)
		g = pickle.loads(s)
		self.assertAlmostEqual(g.l, 39.221)
		self.assertAlmostEqual(g.b, -24.543)
		
	def test_gal_posn_cmp(self):
		"""Test astro.gal_posn.__eq__() and astro.gal_posn.__ne__() 
		methods."""
		
		g1 = astro.gal_posn(232, 54)
		g2 = astro.gal_posn(123, -32)
		g3 = astro.gal_posn(123, -32)
		self.assertTrue(operator.ne(g1, g2))
		self.assertFalse(operator.eq(g1, g2))
		self.assertTrue(operator.eq(g2, g3))
		self.assertFalse(operator.ne(g2, g3))
		
	def test_gal_posn_subscr(self):
		"""Test astro.gal_posn.__getitem__() and astro.gal_posn.__setitem__()
		methods."""
		
		g = astro.gal_posn(39.221, -24.543)
		self.assertAlmostEqual(g[0], 39.221)
		self.assertAlmostEqual(g[1], -24.543)
		
		g[0] = 125.345
		g[1] = 43.232
		self.assertAlmostEqual(g.l, 125.345)
		self.assertAlmostEqual(g.b, 43.232)
		
		self.assertRaises(ValueError, operator.itemgetter(3), g)
		self.assertRaises(ValueError, operator.itemgetter(-1), g)
		
	def test_gal_posn_format(self):
		"""Test astro.gal_posn.format() method."""
		
		g = astro.gal_posn(39.221, -24.543)
		(l, b) = g.format()
		self.assertAlmostEqual(l.to_deg(), 39.221)
		self.assertAlmostEqual(b.to_deg(), -24.543)
		
	def test_geo_posn_init(self):
		"""Test astro.geo_posn constructor."""
		
		g = astro.geo_posn(39.221, -24.543, 2000.345)
		self.assertAlmostEqual(g.lng, 39.221)
		self.assertAlmostEqual(g.lat, -24.543)
		self.assertAlmostEqual(g.elv, 2000.345)
		
		lng = astro.deg_to_dms(39.221)
		lat = astro.deg_to_dms(-24.543)
		g = astro.geo_posn(lng, lat, 2000.345)
		self.assertAlmostEqual(g.lng, 39.221)
		self.assertAlmostEqual(g.lat, -24.543)
		self.assertAlmostEqual(g.elv, 2000.345)
		
		self.assertRaises(ValueError, astro.geo_posn, 400, 0)
		self.assertRaises(ValueError, astro.geo_posn, -400, 0)
		self.assertRaises(ValueError, astro.geo_posn, 0, 100)
		self.assertRaises(ValueError, astro.geo_posn, 0, -100)
		
	def test_geo_posn_reduce(self):
		"""Test astro.geo_posn.__reduce__() method."""
		
		g = astro.geo_posn(39.221, -24.543, 2000.345)
		s = pickle.dumps(g)
		g = pickle.loads(s)
		self.assertAlmostEqual(g.lng, 39.221)
		self.assertAlmostEqual(g.lat, -24.543)
		self.assertAlmostEqual(g.elv, 2000.345)
		
	def test_geo_posn_cmp(self):
		"""Test astro.geo_posn.__eq__() and astro.geo_posn.__ne__()
		methods."""
		
		g1 = astro.geo_posn(232, 53, 1322)
		g2 = astro.geo_posn(124, -24, 1243)
		g3 = astro.geo_posn(124, -24, 1243)
		self.assertTrue(operator.ne(g1, g2))
		self.assertFalse(operator.eq(g1, g2))
		self.assertTrue(operator.eq(g2, g3))
		self.assertFalse(operator.ne(g2, g3))
		
	def test_geo_posn_subscr(self):
		"""Test astro.geo_posn.__getitem__() and astro.geo_posn.__setitem__()
		methods."""
		
		g = astro.geo_posn(39.221, -24.543, 2000.345)
		self.assertAlmostEqual(g[0], 39.221)
		self.assertAlmostEqual(g[1], -24.543)
		self.assertAlmostEqual(g[2], 2000.345)
		
		g[0] = 200.221
		g[1] = 38.231
		g[2] = 1003.053
		self.assertAlmostEqual(g.lng, 200.221)
		self.assertAlmostEqual(g.lat, 38.231)
		self.assertAlmostEqual(g.elv, 1003.053)
		
		self.assertRaises(ValueError, operator.itemgetter(3), g)
		self.assertRaises(ValueError, operator.itemgetter(-1), g)
		
	def test_geo_posn_format(self):
		"""Test astro.geo_posn.format() method."""
		
		g = astro.geo_posn(39.221, -24.543, 2000.345)
		(lng, lat) = g.format()
		self.assertAlmostEqual(lng.to_deg(), 39.221)
		self.assertAlmostEqual(lat.to_deg(), -24.543)
		
	def test_rect_posn_init(self):
		"""Test astro.rect_posn constructor method."""
		
		r = astro.rect_posn(2.409, 9.324, 4.442)
		self.assertAlmostEqual(r.X, 2.409)
		self.assertAlmostEqual(r.Y, 9.324)
		self.assertAlmostEqual(r.Z, 4.442)
		
	def test_rect_posn_reduce(self):
		"""Test astro.rect_posn pickle.__reduce__() method."""
		
		r = astro.rect_posn(2.409, 9.324, 4.442)
		s = pickle.dumps(r)
		r = pickle.loads(s)
		self.assertAlmostEqual(r.X, 2.409)
		self.assertAlmostEqual(r.Y, 9.324)
		self.assertAlmostEqual(r.Z, 4.442)
		
	def test_rect_posn_cmp(self):
		"""Test astro.rect_posn.__eq__() and astro.rect_posn.__ne__()
		methods."""
		
		r1 = astro.rect_posn(4, 5, 2)
		r2 = astro.rect_posn(5, 3, 1)
		r3 = astro.rect_posn(5, 3, 1)
		self.assertTrue(operator.ne(r1, r2))
		self.assertFalse(operator.eq(r1, r2))
		self.assertTrue(operator.eq(r2, r3))
		self.assertFalse(operator.ne(r2, r3))
		
	def test_rect_posn_subscr(self):
		"""Test astro.rect_posn.__getitem__() and astro.rect_posn.__setitem__()
		methods."""
		
		r = astro.rect_posn(2.409, 9.324, 4.442)
		self.assertAlmostEqual(r[0], 2.409)
		self.assertAlmostEqual(r[1], 9.324)
		self.assertAlmostEqual(r[2], 4.442)
		r[0] = -2.343
		r[1] = -4.215
		r[2] = 0.432
		self.assertAlmostEqual(r.X, -2.343)
		self.assertAlmostEqual(r.Y, -4.215)
		self.assertAlmostEqual(r.Z, 0.432)
		
	def test_get_hrz_from_equ(self):
		"""Test astro.get_hrz_from_equ() function."""
		
		ACCURACY_AZ = 0
		ACCURACY_ALT = 0
		
		ra_aa = (\
			astro.hms( 8,  5, 38.2531).to_deg(), 
			astro.hms(13, 36, 48.7661).to_deg(),
			astro.hms(20, 42, 36.8591).to_deg(),
			astro.hms( 1, 14, 39.0438).to_deg(),
			astro.hms(10, 28, 52.5980).to_deg())
			
		obs = astro.lnlat_posn(0.0, 0.0)    
		
		ira = iter(ra_aa)
		for t in self.times:
			ra = next(ira)
			
			j = t.to_jd()
			equ = astro.equ_posn(ra, 30.0)
			hrz = astro.get_hrz_from_equ(equ, obs, j)
			az = hrz.az
			if az > 180:
					az -= 360
			self.assertAlmostEqual(math.radians(az), 0.0, ACCURACY_AZ)
			self.assertAlmostEqual(hrz.alt, 60.0, ACCURACY_ALT)               
			
	def test_get_equ_from_hrz(self):
		"""Test astro.get_equ_from_hrz() function."""
		
		ACCURACY_RA = 0
		ACCURACY_DEC = 0
		
		ra_aa = (\
			astro.hms( 8,  5, 38.2531).to_deg(), 
			astro.hms(13, 36, 48.7661).to_deg(),
			astro.hms(20, 42, 36.8591).to_deg(),
			astro.hms( 1, 14, 39.0438).to_deg(),
			astro.hms(10, 28, 52.5980).to_deg()) 
			
		hrz = astro.hrz_posn(0.0, 60.0)
		obs = astro.lnlat_posn(0.0, 0.0)
		
		ira = iter(ra_aa)
		for t in self.times:
			ra = next(ira)
			
			j = t.to_jd()
			equ = astro.get_equ_from_hrz(hrz, obs, j)
			self.assertAlmostEqual(math.radians(equ.ra), math.radians(ra), ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, 30.0, ACCURACY_DEC)
			
	def test_get_equ_prec(self):
		"""Test astro.get_equ_prec() function."""
		
		ACCURACY_RA = 4
		ACCURACY_DEC = 4
		
		equ_in = (\
			astro.equ_posn(  0.0,   0.0),
			astro.equ_posn(180.0,   0.0), 
			astro.equ_posn(  0.0,  80.0),
			astro.equ_posn(240.0, -55.0),
			astro.equ_posn(  0.0,   0.0))
			
		iequ = iter(equ_in)
		for t in self.times:
			j = t.to_jd() 
			equ0 = next(iequ)
			equ = astro.get_equ_prec(equ0, j)
			chk = astro.get_precession(astro.J2000_UTC_JD, equ0, j)
			self.assertAlmostEqual(equ.ra, chk.ra, ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, chk.dec, ACCURACY_DEC) 
			
	def test_get_equ_prec2(self):
		"""Test astro.get_equ_prec2() function."""
		
		ACCURACY_RA = 4
		ACCURACY_DEC = 4
		
		equ_in = (\
			astro.equ_posn(  0.0,   0.0),
			astro.equ_posn(180.0,   0.0), 
			astro.equ_posn(  0.0,  80.0),
			astro.equ_posn(240.0, -55.0),
			astro.equ_posn(  0.0,   0.0))
			
		iequ = iter(equ_in)
		for t in self.times:
			j = t.to_jd() 
			equ0 = next(iequ)
			equ = astro.get_equ_prec2(equ0, j, astro.J2000_UTC_JD)
			chk = astro.get_precession(j, equ0, astro.J2000_UTC_JD)
			self.assertAlmostEqual(equ.ra, chk.ra, ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, chk.dec, ACCURACY_DEC)
			
	def test_nutation_init(self):
		"""Test astro.nutation constructor."""
		
		n = astro.nutation(0.0234, -0.0421, 23.5656)
		self.assertAlmostEqual(n.longitude, 0.0234)
		self.assertAlmostEqual(n.obliquity, -0.0421)
		self.assertAlmostEqual(n.ecliptic, 23.5656)
		
		self.assertRaises(ValueError, astro.nutation, 400, 0, 0)
		self.assertRaises(ValueError, astro.nutation, -400, 0, 0)
		self.assertRaises(ValueError, astro.nutation, 0, 100, 0)
		self.assertRaises(ValueError, astro.nutation, 0, -100, 0)
		self.assertRaises(ValueError, astro.nutation, 0, 0, 100)
		self.assertRaises(ValueError, astro.nutation, 0, 0, -100)
		
	def test_get_nutation(self):
		"""Test astro.get_nutation() function."""
		
		ACCURACY_LNG = 4
		ACCURACY_OBL = 4
		ACCURACY_ECL = 4
		
		lng_aa = (\
			astro.dms(True,  0, 0, 15.453).to_deg(),
			astro.dms(True,  0, 0, 17.708).to_deg(),
			astro.dms(True,  0, 0, 15.856).to_deg(),
			astro.dms(True,  0, 0, 18.075).to_deg(),
			astro.dms(True,  0, 0,  2.557).to_deg())
			
		obl_aa = (\
			astro.dms(True,  0, 0, 2.569).to_deg(),
			astro.dms(True,  0, 0, 1.241).to_deg(),
			astro.dms(True,  0, 0, 0.797).to_deg(),
			astro.dms(False, 0, 0, 0.339).to_deg(),
			astro.dms(True,  0, 0, 8.606).to_deg())
			
		ecl_aa = (\
			astro.dms(False, 23, 26, 18.384).to_deg(),
			astro.dms(False, 23, 26, 19.604).to_deg(),
			astro.dms(False, 23, 26, 19.909).to_deg(),
			astro.dms(False, 23, 26, 20.957).to_deg(),
			astro.dms(False, 23, 26, 22.569).to_deg())
			
		ilng = iter(lng_aa)
		iobl = iter(obl_aa)
		iecl = iter(ecl_aa)
		for t in self.times:
			lng = next(ilng)
			obl = next(iobl)
			ecl = next(iecl)
			
			j = t.to_jd()
			nut = astro.get_nutation(j)
			self.assertAlmostEqual(nut.longitude, lng, ACCURACY_LNG)      
			self.assertAlmostEqual(nut.obliquity, obl, ACCURACY_OBL)
			self.assertAlmostEqual(nut.ecliptic, ecl, ACCURACY_ECL)
			
	def test_get_equ_nut(self):
		"""Test astro.get_equ_nut() function."""
		
		body = astro.equ_posn(41.547213, 49.348483)
		JD = 2462088.7
		
		ACCURACY = 3
		deltaRA = 15.842766509979356
		deltaDec = 6.2164102884679551
		
		body2 = astro.get_equ_nut(body, JD)
		self.assertAlmostEqual((body2.ra -body.ra )*3600, deltaRA,  ACCURACY)
		self.assertAlmostEqual((body2.dec-body.dec)*3600, deltaDec, ACCURACY)
		
	def get_equ_aber(self):
		"""Test astro.get_equ_aber() funciton."""
		
		body = astro.equ_posn(41.547213, 49.348483)
		JD = 2462088.7
		
		ACCURACY = 3
		deltaRA = 30.09128661282716
		deltaDec = 6.621256621590987
		
		body2 = astro.get_equ_aber(body, JD)
		self.assertAlmostEqual((body2.ra -body.ra )*3600, deltaRA,  ACCURACY)
		self.assertAlmostEqual((body2.dec-body.dec)*3600, deltaDec, ACCURACY)
		
	def test_rst_time_init(self):
		"""Test astro.rst_time constructor method."""
		
		rise = self.times[0].to_jd()
		transit = self.times[1].to_jd()
		set = self.times[2].to_jd()
		rst = astro.rst_time(rise, set, transit)
		self.assertAlmostEqual(rst.rise, rise)
		self.assertAlmostEqual(rst.transit, transit)
		self.assertAlmostEqual(rst.set, set)
		
	def test_get_object_rst(self):
		"""Test astro.get_object_rst() function."""
		
		body = astro.equ_posn(0, 33.0)
		
		ACCURACY = 2
		transit_aa = (\
			astro.date(2001,  1, 22, 21, 50, 48.04).to_jd(),
			astro.date(2001,  4, 16, 16, 20, 31.74).to_jd(),
			astro.date(2001,  8,  2,  9, 15, 56.73).to_jd(),
			astro.date(2001, 10, 10,  4, 44, 39.94).to_jd(),
			astro.date(1979,  2, 28, 19, 26, 50.45).to_jd())
			
		itr = iter(transit_aa) 
		for t in self.times:
			j = t.to_jd()
			trns = astro.get_object_rst(j, self.geo, body)
			self.assertAlmostEqual(astro.utc_to_tt(trns.transit), next(itr), ACCURACY)
			
		body = astro.equ_posn(0, 89.0)
		for t in self.times:
			j = t.to_jd()
			trns = astro.get_object_rst(j, self.geo, body)
			self.assertTrue(trns is None)
			
		body = astro.equ_posn(0, -89.0)
		for t in self.times:
			j = t.to_jd()
			trns = astro.get_object_rst(j, self.geo, body)
			self.assertTrue(trns is None)
			
	def test_get_solar_rst(self):
		"""Test astro.get_solar_rst() function."""
		
		ACCURACY = 2
		
		p = (-self.geo.lng / 360.0) + (1.002738 * astro.sec_to_jd(64)) 
		transit_aa = (\
			astro.date(2001,  1, 22, 12, 11, 37.57).to_jd() + p,
			astro.date(2001,  4, 16, 11, 59, 44.68).to_jd() + p,
			astro.date(2001,  8,  2, 12,  6, 14.09).to_jd() + p,
			astro.date(2001, 10, 10, 11, 46, 58.27).to_jd() + p,
			astro.date(1979,  2, 28, 12, 12, 29.32).to_jd() + p)
			
		itr = iter(transit_aa) 
		for t in self.times:
			j = t.to_jd()
			trns = astro.get_solar_rst(j, self.geo)
			self.assertAlmostEqual(astro.utc_to_tt(trns.transit), next(itr), ACCURACY)
			
	def test_get_solar_equ_coords(self):
		"""Test astro.get_solar_equ_coords() function."""
		
		ACCURACY_RA  = 2
		ACCURACY_DEC = 2
		
		ra_aa = (\
			astro.hms(20, 17,  7.63).to_deg(),
			astro.hms( 1, 36, 40.54).to_deg(),
			astro.hms( 8, 48, 53.15).to_deg(),
			astro.hms(13,  1, 45.13).to_deg(),
			astro.hms(22, 41, 38.54).to_deg())
			
		dec_aa = (\
			astro.dms(True,  19, 42, 30.7).to_deg(),
			astro.dms(False, 10,  3, 55.1).to_deg(),
			astro.dms(False, 17, 47, 59.1).to_deg(),
			astro.dms(True,   6, 34, 59.9).to_deg(),
			astro.dms(True,   8, 16, 14.4).to_deg())     
			
		ira = iter(ra_aa)
		idec = iter(dec_aa)
		for t in self.times:
			j = t.to_jd()
			equ = astro.get_solar_equ_coords(j)
			self.assertAlmostEqual(equ.ra, next(ira), ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, next(idec), ACCURACY_DEC)
			
	def test_get_jupiter_rst(self):
		"""Test astro.get_jupiter_rst() function."""
		
		ACCURACY = 2
		
		p = (-self.geo.lng / 360.0) + (1.002738 * astro.sec_to_jd(64))
		transit_aa = (\
			astro.date(2001,  1, 21, 19, 51, 59).to_jd() + p,
			astro.date(2001,  4, 16, 14, 56, 44).to_jd() + p,
			astro.date(2001,  8,  2,  9, 34, 53).to_jd() + p,
			astro.date(2001, 10, 10,  5, 48, 36).to_jd() + p,
			astro.date(1979,  2, 27, 21, 36, 57).to_jd() + p)
			
		itr = iter(transit_aa)
		for t in self.times:
			j = t.to_jd()
			trns = astro.get_jupiter_rst(j, self.geo)
			self.assertAlmostEqual(astro.utc_to_tt(trns.transit), next(itr), ACCURACY)
			
	def test_get_jupiter_equ_coord(self):
		"""Test astro.get_jupiter_equ_coords() function."""
		
		ACCURACY_RA  = 2
		ACCURACY_DEC = 2
		
		ra_aa = (\
			astro.hms( 3, 56, 56.236).to_deg(),
			astro.hms( 4, 35, 28.272).to_deg(),
			astro.hms( 6, 18, 42.909).to_deg(),
			astro.hms( 7,  4,  7.660).to_deg(),
			astro.hms( 8,  9, 40.552).to_deg()) 
			
		dec_aa = (\
			astro.dms(False, 19, 41, 39.73).to_deg(),
			astro.dms(False, 21, 34, 42.35).to_deg(),
			astro.dms(False, 23,  7,  6.00).to_deg(),
			astro.dms(False, 22, 27, 30.44).to_deg(),
			astro.dms(False, 20, 48, 29.55).to_deg())  
			
		ira = iter(ra_aa)
		idec = iter(dec_aa)
		for t in self.times:
			j = t.to_jd()
			equ = astro.get_jupiter_equ_coords(j)
			self.assertAlmostEqual(equ.ra, next(ira), ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, next(idec), ACCURACY_DEC)                                                     
			
	def test_get_lunar_rst(self):
		"""Test astro.get_lunar_rst() function."""
		
		ACCURACY = 1
		
		p = (-self.geo.lng / 360.0) + (1.002738 * astro.sec_to_jd(64))
		transit_aa = (\
			astro.date(2001,  1, 22, 10, 34, 47).to_jd() + p,
			astro.date(2001,  4, 16,  6, 46, 16).to_jd() + p,
			astro.date(2001,  8,  1, 22, 19, 42).to_jd() + p,
			astro.date(2001, 10, 10,  6,  3, 51).to_jd() + p,
			astro.date(1979,  2, 28, 13, 51, 32).to_jd() + p)
			
		itr = iter(transit_aa)
		for t in self.times:
			j = t.to_jd()
			trns = astro.get_lunar_rst(j, self.geo)
			self.assertAlmostEqual(astro.utc_to_tt(trns.transit), next(itr), ACCURACY)
			
	def test_get_lunar_equ_coords(self):
		"""Test astro.get_lunar_equ_coords() function."""
		
		ACCURACY_RA  = 1
		ACCURACY_DEC = 2
		
		ra_aa = (\
			astro.hms(18, 19, 24.98).to_deg(),
			astro.hms(20,  9, 55.87).to_deg(),
			astro.hms(19,  5, 43.33).to_deg(),
			astro.hms( 7,  3, 20.27).to_deg(),
			astro.hms(23, 50, 17.96).to_deg())
			
		dec_aa = (\
			astro.dms(True,  22, 20, 47.9).to_deg(),
			astro.dms(True,  21, 54, 24.7).to_deg(),
			astro.dms(True,  23, 23, 33.8).to_deg(),
			astro.dms(False, 23, 54,  8.3).to_deg(),
			astro.dms(True,   2,  0, 28.9).to_deg())    
			
		ira = iter(ra_aa)
		idec = iter(dec_aa)
		for t in self.times:
			j = t.to_jd()
			equ = astro.get_lunar_equ_coords(j)
			self.assertAlmostEqual(equ.ra, next(ira), ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, next(idec), ACCURACY_DEC) 
			
	def test_get_ecl_from_equ(self):
		"""Test astro.get_ecl_from_equ() function."""
		
		ACCURACY_LNG = 4
		ACCURACY_LAT = 4
		
		# Data set generated at http://ned.ipac.caltech.edu/forms/calculator.html
		lng_aa = (302.1000, 26.3000, 129.5000, 196.7000, 338.9000)
		lat_aa = (  0.0000, -1.0000,   3.0000,  -5.0000,   7.0000)
		
		ra_aa = (\
			astro.hms(20, 17, 26.98).to_deg(),
			astro.hms( 1, 39,  1.88).to_deg(),
			astro.hms( 8, 51, 10.38).to_deg(),
			astro.hms(12, 53, 44.69).to_deg(),
			astro.hms(22, 31, 32.27).to_deg())
		dec_aa = (\
			astro.dms(True,  19, 41, 50.94).to_deg(),
			astro.dms(False,  9, 13, 17.35).to_deg(),
			astro.dms(False, 20, 46, 10.16).to_deg(),
			astro.dms(True,  11, 10, 41.32).to_deg(),
			astro.dms(True,   1, 44, 23.34).to_deg())  
			
		ilng = iter(lng_aa)
		ilat = iter(lat_aa)
		ira = iter(ra_aa)
		idec = iter(dec_aa)
		for t in self.times:
			lng = next(ilng)
			lat = next(ilat)
			ra = next(ira)
			dec = next(idec)
			
			j = t.to_jd()
			equ = astro.equ_posn(ra, dec)
			ecl = astro.get_ecl_from_equ(equ, j)
			self.assertAlmostEqual(ecl.lng, lng, ACCURACY_LNG)
			self.assertAlmostEqual(ecl.lat, lat, ACCURACY_LAT)
			
	def test_get_equ_from_ecl(self):
		"""Test astro.get_equ_from_ecl() function."""
		
		ACCURACY_RA = 4
		ACCURACY_DEC = 4
		
		# Data set generated at http://ned.ipac.caltech.edu/forms/calculator.html
		lng_aa = (302.1000, 26.3000, 129.5000, 196.7000, 338.9000)
		lat_aa = (  0.0000, -1.0000,   3.0000,  -5.0000,   7.0000)
		
		ra_aa = (\
			astro.hms(20, 17, 26.98).to_deg(),
			astro.hms( 1, 39,  1.88).to_deg(),
			astro.hms( 8, 51, 10.38).to_deg(),
			astro.hms(12, 53, 44.69).to_deg(),
			astro.hms(22, 31, 32.27).to_deg())
		dec_aa = (\
			astro.dms(True,  19, 41, 50.94).to_deg(),
			astro.dms(False,  9, 13, 17.35).to_deg(),
			astro.dms(False, 20, 46, 10.16).to_deg(),
			astro.dms(True,  11, 10, 41.32).to_deg(),
			astro.dms(True,   1, 44, 23.34).to_deg()) 
			
		ilng = iter(lng_aa)
		ilat = iter(lat_aa)
		ira = iter(ra_aa)
		idec = iter(dec_aa)
		for t in self.times:
			lng = next(ilng)
			lat = next(ilat)
			ra = next(ira)
			dec = next(idec)
			
			j = t.to_jd()
			
			ecl = astro.ecl_posn(lng, lat)
			equ = astro.get_equ_from_ecl(ecl, j)
			self.assertAlmostEqual(equ.ra, ra, ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, dec, ACCURACY_DEC)
			
	def test_get_gal_from_equ(self):
		"""Test astro.get_gal_from_equ() function."""
		
		ACCURACY_L = 2
		ACCURACY_B = 2
		
		# Test set generated at http://fuse.pha.jhu.edu/support/tools/eqtogal.html
		ra_4c = (\
			astro.hms(13, 35, 31.41).to_deg(),
			astro.hms(15, 56,  3.01).to_deg(),
			astro.hms( 1, 35, 19.19).to_deg(),
			astro.hms(19, 45, 35.25).to_deg())  
		dec_4c = (\
			astro.dms(True,   6,  3, 38.0).to_deg(),
			astro.dms(False,  1, 49, 36.9).to_deg(),
			astro.dms(False, 22, 46,  3.9).to_deg(),
			astro.dms(False, 78, 37,  8.5).to_deg())
			
		l_4c = (323.3000, 11.5000, 136.7000, 110.9000)  
		b_4c = ( 54.7000, 38.5000, -38.6000,  24.1000)
		
		ira = iter(ra_4c)
		idec = iter(dec_4c)
		il = iter(l_4c)
		ib = iter(b_4c)
		for i in range(len(ra_4c)):
			ra = next(ira)
			dec = next(idec)
			l = next(il)
			b = next(ib)
			
			equ = astro.equ_posn(ra, dec)
			gal = astro.get_gal_from_equ(equ)
			self.assertAlmostEqual(gal.l, l, ACCURACY_L)
			self.assertAlmostEqual(gal.b, b, ACCURACY_B)
			
	def test_get_equ_from_gal(self):
		"""Test astro.get_equ_from_gal() function."""
		
		ACCURACY_RA = 2
		ACCURACY_DEC = 2
		
		# Test set generated at http://fuse.pha.jhu.edu/support/tools/eqtogal.html
		ra_4c = (\
			astro.hms(13, 35, 31.41).to_deg(),
			astro.hms(15, 56,  3.01).to_deg(),
			astro.hms( 1, 35, 19.19).to_deg(),
			astro.hms(19, 45, 35.25).to_deg())  
		dec_4c = (\
			astro.dms(True,   6,  3, 38.0).to_deg(),
			astro.dms(False,  1, 49, 36.9).to_deg(),
			astro.dms(False, 22, 46,  3.9).to_deg(),
			astro.dms(False, 78, 37,  8.5).to_deg())
			
		l_4c = (323.3000, 11.5000, 136.7000, 110.9000)  
		b_4c = ( 54.7000, 38.5000, -38.6000,  24.1000)
		
		ira = iter(ra_4c)
		idec = iter(dec_4c)
		il = iter(l_4c)
		ib = iter(b_4c)
		for i in range(len(ra_4c)):
			ra = next(ira)
			dec = next(idec)
			l = next(il)
			b = next(ib)
			
			gal = astro.gal_posn(l, b)
			equ = astro.get_equ_from_gal(gal)
			self.assertAlmostEqual(equ.ra, ra, ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, dec, ACCURACY_DEC)
			
	def test_get_gal_from_equ2000(self):
		"""Test astro.get_gal_from_equ2000() function."""
		
		ACCURACY_L = 2
		ACCURACY_B = 2
		
		# Test set generated at http://ned.ipac.caltech.edu/forms/calculator.html
		ra_4c = (\
			astro.hms(13, 35, 31.41).to_deg(),
			astro.hms(15, 56,  3.01).to_deg(),
			astro.hms( 1, 35, 19.19).to_deg(),
			astro.hms(19, 45, 35.25).to_deg())  
		dec_4c = (\
			astro.dms(True,   6,  3, 38.0).to_deg(),
			astro.dms(False,  1, 49, 36.9).to_deg(),
			astro.dms(False, 22, 46,  3.9).to_deg(),
			astro.dms(False, 78, 37,  8.5).to_deg())
			
		l_4c = (322.3585, 11.1982, 135.9796, 110.8142)  
		b_4c = ( 55.1396, 39.1055, -38.9896,  23.9482)
		
		ira = iter(ra_4c)
		idec = iter(dec_4c)
		il = iter(l_4c)
		ib = iter(b_4c)
		for i in range(len(ra_4c)):
			ra = next(ira)
			dec = next(idec)
			l = next(il)
			b = next(ib)
			
			equ = astro.equ_posn(ra, dec)
			gal = astro.get_gal_from_equ2000(equ)
			self.assertAlmostEqual(gal.l, l, ACCURACY_L)
			self.assertAlmostEqual(gal.b, b, ACCURACY_B)
			
	def test_get_equ2000_from_gal(self):
		"""Test astro.get_equ2000_from_gal() function."""
		
		ACCURACY_RA = 2
		ACCURACY_DEC = 2
		
		# Test set generated at http://ned.ipac.caltech.edu/forms/calculator.html
		ra_4c = (\
			astro.hms(13, 35, 31.41).to_deg(),
			astro.hms(15, 56,  3.01).to_deg(),
			astro.hms( 1, 35, 19.19).to_deg(),
			astro.hms(19, 45, 35.25).to_deg())  
		dec_4c = (\
			astro.dms(True,   6,  3, 38.0).to_deg(),
			astro.dms(False,  1, 49, 36.9).to_deg(),
			astro.dms(False, 22, 46,  3.9).to_deg(),
			astro.dms(False, 78, 37,  8.5).to_deg())
			
		l_4c = (322.3585, 11.1982, 135.9796, 110.8142)  
		b_4c = ( 55.1396, 39.1055, -38.9896,  23.9482)
		
		ira = iter(ra_4c)
		idec = iter(dec_4c)
		il = iter(l_4c)
		ib = iter(b_4c)
		for i in range(len(ra_4c)):
			ra = next(ira)
			dec = next(idec)
			l = next(il)
			b = next(ib)
			
			gal = astro.gal_posn(l, b)
			equ = astro.get_equ2000_from_gal(gal)
			self.assertAlmostEqual(equ.ra, ra, ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, dec, ACCURACY_DEC)
			
	def test_get_apparent_posn(self):
		"""Test astro.get_apparent_posn() function."""
		
		ACCURACY_RA = 4
		ACCURACY_DEC = 4
		
		ra_in  = (201.2500, 69.2517)  
		dec_in = (-43.0667,   29.6708)
		
		ra_out  = (201.26667, 69.27255)
		dec_out = (-43.06925, 29.67399)
		
		idec_in = iter(dec_in)
		ira_out = iter(ra_out)
		idec_out = iter(dec_out)
		j = self.times[0].to_jd()
		for ra in ra_in:
			equ = astro.equ_posn(ra, next(idec_in))
			equ = astro.get_apparent_posn(equ, j)
			self.assertAlmostEqual(equ.ra, next(ira_out), ACCURACY_RA)
			self.assertAlmostEqual(equ.dec, next(idec_out), ACCURACY_DEC)
			
	def test_B1950_to_J2000(self):
		"""Test astro.B1950_to_J2000() function."""
		
		ACCURACY_RA = 2
		ACCURACY_DEC = 2
		
		# Test set generated at http://ned.ipac.caltech.edu/forms/calculator.html
		pos_J2000 = astro.equ_posn(296.39687500, 78.61902778)
		pos_B1950 = astro.equ_posn(296.98456518, 78.49400140)
		
		equ = astro.B1950_to_J2000(pos_B1950)
		self.assertAlmostEqual(equ.ra,  pos_J2000.ra,  ACCURACY_RA)
		self.assertAlmostEqual(equ.dec, pos_J2000.dec, ACCURACY_DEC)
		
	def test_J2000_to_B1950(self):
		"""Test astro.J2000_to_B1950() function."""
		
		ACCURACY_RA = 2
		ACCURACY_DEC = 2
		
		# Test set generated at http://ned.ipac.caltech.edu/forms/calculator.html
		pos_J2000 = astro.equ_posn(296.39687500, 78.61902778)
		pos_B1950 = astro.equ_posn(296.98456518, 78.49400140)
		
		equ = astro.J2000_to_B1950(pos_J2000)
		self.assertAlmostEqual(equ.ra,  pos_B1950.ra,  ACCURACY_RA)
		self.assertAlmostEqual(equ.dec, pos_B1950.dec, ACCURACY_DEC)


class astro_test_suite(unittest.TestSuite):
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(astro_tests))
                   

if __name__ == '__main__':
	unittest.main()
