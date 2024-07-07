"""
Unit test for the lsl.catalog module.
"""

import unittest

from lsl import catalog
import lsl.testing


__version__  = "0.2"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class catalog_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.catalog
    module."""
    
    def test_constructor(self):
        """Test catalog constructors."""
        
        for name in ('LWA', 'PSR', 'PKS', 'PKS90', '3C', '4C', '2FGL', '3FGL', '4FGL', '4FGL-DR4'):
            with self.subTest(name=name):
                catalog.CatalogFactory.get_catalog(name)
                
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
        self.assertEqual(source.name, 'CasA')
        self.assertTrue('NRAO711' in source.alias_list)
        self.assertAlmostEqual(source.position.j2000_equ.ra, 350.850, 3)
        self.assertAlmostEqual(source.position.j2000_equ.dec, 58.815, 3)
        
        source = cat.lookup('CygA')
        self.assertEqual(source.name, 'CygA')
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
