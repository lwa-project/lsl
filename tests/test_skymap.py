"""
Unit test for the lsl.skymap module.
"""

import unittest

from lsl import skymap, astro


__version__   = "0.2"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class skymap_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.skymap
    module."""
    
    def test_SkyMapGSM_init(self):
        """Test skymap.SkyMapGSM class constructor method.""" 
        
        skymap.SkyMapGSM()
        
    def test_SkyMapGSM_compute_total_power(self):
        """Test skymap.SkyMapGSM.compute_total_power() method."""
        
        s = skymap.SkyMapGSM()
        s.compute_total_power()
        
    def test_ProjectedSkyMap_init_GSM(self):
        """Test skymap.ProjectedSkyMap constructor method using SkyMapGSM."""
        
        s = skymap.SkyMapGSM()
        skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        
    def test_ProjectedSkyMap_compute_visible_power_GSM(self):
        """Test skymap.ProjectedSkyMap.compute_visible_power() method using SkyMapGSM."""
        
        s = skymap.SkyMapGSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys()) 
        p.compute_visible_power()
        
    def test_ProjectedSkyMap_get_direction_cosines_GSM(self):
        """Test skymap.ProjectedSkyMap.get_direction_cosines() method using SkyMapGSM."""
        
        s = skymap.SkyMapGSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        p.get_direction_cosines()
        
    def test_SkyMapLFSM_init(self):
        """Test skymap.SkyMapLFSM class constructor method.""" 
        
        skymap.SkyMapLFSM()
        
    def test_SkyMapLFSM_compute_total_power(self):
        """Test skymap.SkyMapLFSM.compute_total_power() method."""
        
        s = skymap.SkyMapLFSM()
        s.compute_total_power()

    def test_ProjectedSkyMap_init_LFSM(self):
        """Test skymap.ProjectedSkyMap constructor method using SkyMapLFSM."""
        
        s = skymap.SkyMapLFSM()
        skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())

    def test_ProjectedSkyMap_compute_visible_power_LFSM(self):
        """Test skymap.ProjectedSkyMap.compute_visible_power() method using SkyMapLFSM."""
        
        s = skymap.SkyMapLFSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys()) 
        p.compute_visible_power()

    def test_ProjectedSkyMap_get_direction_cosines_LFSM(self):
        """Test skymap.ProjectedSkyMap.get_direction_cosines() method using SkyMapLFSM."""
        
        s = skymap.SkyMapLFSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        p.get_direction_cosines()


class skymap_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lwa_user.skymap
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(skymap_tests))        
        
        
if __name__ == '__main__':
    unittest.main()
