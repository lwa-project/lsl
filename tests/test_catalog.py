# -*- coding: utf-8 -*-

"""Unit test for lsl.catalog module."""

import unittest

from lsl import catalog


__revision__ = "$Revision$"
__version__  = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class catalog_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.catalog
	module."""
	
	def test_LWDA(self):
		"""Test catalog.LWA_Catalog constructor."""
		
		catalog.CatalogFactory.get_catalog('LWA')
		
	def test_PSR(self):
		"""Test catalog.PSR_Catalog contructor."""
		
		catalog.CatalogFactory.get_catalog('PSR')
		
	def test_PKS(self):
		"""Test catalog.PKS_Catalog constructor."""
		
		catalog.CatalogFactory.get_catalog('PKS')
		
	def test_PKS90(self):
		"""Test catalog.PKS90_Catalog constructor."""
		
		catalog.CatalogFactory.get_catalog('PKS90')
		
	def test_C3C(self):
		"""Test catalog.C3C_Catalog constructor."""
		
		catalog.CatalogFactory.get_catalog('3C')
		
	def test_C4C(self):
		"""Test catalog.C4C_Catalog constructor."""
		
		catalog.CatalogFactory.get_catalog('4C')
		
	def test_2FGL(self):
		"""Test catalog.F2FGL_Catalog constructor."""
		
		catalog.CatalogFactory.get_catalog('2FGL')
		
	def test_get_names(self):
		"""Test catalog.CatalogFactory.get_names() method."""
		
		names = catalog.CatalogFactory.get_names()
		self.assertTrue('LWA' in names)
		self.assertTrue('NONSENSE' not in names)
		
	def test_get_catalog(self):
		"""Test catalog.CatalogFactory.get_catalog() function."""
		
		cat = catalog.CatalogFactory.get_catalog('LWA')
		self.assertTrue(isinstance(cat, catalog.LWA_Catalog))
		
		self.assertRaises(ValueError, catalog.CatalogFactory.get_catalog, 'NONSENSE')
		
	def test_lookup(self):
		"""Test catalog.Catalog.lookup() method."""
		
		cat = catalog.CatalogFactory.get_catalog('LWA')
		
		source = cat.lookup('SgrA')
		self.assertEqual(source.name, 'SgrA')
		self.assertTrue('W24' in source.alias_list)
		self.assertAlmostEqual(source.position.j2000_equ.ra, 266.300, 3)
		self.assertAlmostEqual(source.position.j2000_equ.dec, -28.806, 3)
		
		source = cat.lookup('CasA')
		self.assertEquals(source.name, 'CasA')
		self.assertTrue('NRAO711' in source.alias_list)
		self.assertAlmostEqual(source.position.j2000_equ.ra, 350.850, 3)
		self.assertAlmostEqual(source.position.j2000_equ.dec, 58.815, 3)
		
		source = cat.lookup('CygA')
		self.assertEquals(source.name, 'CygA')
		self.assertTrue('3C405' in source.alias_list)
		self.assertAlmostEqual(source.position.j2000_equ.ra, 299.868, 3)
		self.assertAlmostEqual(source.position.j2000_equ.dec, 40.734, 3) 
		
		source = cat.lookup('NONSENSE')
		self.assertTrue(source is None) 


class catalog_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lwa_user.catalog
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(catalog_tests)) 


if __name__ == '__main__':
	unittest.main()
