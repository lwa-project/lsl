# -*- coding: utf-8 -*-

"""Unit test for lsl.skymap module."""

import unittest

from lsl import skymap, astro


__revision__  = "$Rev$"
__version__   = "0.2"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class skymap_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.skymap
    module."""
    
    def test_SkyMapGSM_init(self):
        """Test skymap.SkyMapGSM class constructor method.""" 
        
        skymap.SkyMapGSM()
        
    def test_SkyMapGSM_ComputeTotalPowerFromSky(self):
        """Test skymap.SkyMapGSM.ComputeTotalPowerFromSky() method."""
        
        s = skymap.SkyMapGSM()
        s.ComputeTotalPowerFromSky()
        
    def test_ProjectedSkyMap_init_GSM(self):
        """Test skymap.ProjectedSkyMap constructor method using SkyMapGSM."""
        
        s = skymap.SkyMapGSM()
        skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        
    def test_ProjectedSkyMap_ComputeTotalPowerFromVisibleSky_GSM(self):
        """Test skymap.ProjectedSkyMap.ComputeTotalPowerFromVisibleSky() method using SkyMapGSM."""
        
        s = skymap.SkyMapGSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys()) 
        p.ComputeTotalPowerFromVisibleSky()
        
    def test_ProjectedSkyMap_ComputeDirectionCosines_GSM(self):
        """Test skymap.ProjectedSkyMap.ComputeDirectionCosines() method using SkyMapGSM."""
        
        s = skymap.SkyMapGSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        p.ComputeDirectionCosines()
        
    def test_SkyMapLFSM_init(self):
        """Test skymap.SkyMapLFSM class constructor method.""" 
        
        skymap.SkyMapLFSM()
        
    def test_SkyMapLFSM_ComputeTotalPowerFromSky(self):
        """Test skymap.SkyMapLFSM.ComputeTotalPowerFromSky() method."""
        
        s = skymap.SkyMapLFSM()
        s.ComputeTotalPowerFromSky()

    def test_ProjectedSkyMap_init_LFSM(self):
        """Test skymap.ProjectedSkyMap constructor method using SkyMapLFSM."""
        
        s = skymap.SkyMapLFSM()
        skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())

    def test_ProjectedSkyMap_ComputeTotalPowerFromVisibleSky_LFSM(self):
        """Test skymap.ProjectedSkyMap.ComputeTotalPowerFromVisibleSky() method using SkyMapLFSM."""
        
        s = skymap.SkyMapLFSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys()) 
        p.ComputeTotalPowerFromVisibleSky()

    def test_ProjectedSkyMap_ComputeDirectionCosines_LFSM(self):
        """Test skymap.ProjectedSkyMap.ComputeDirectionCosines() method using SkyMapLFSM."""
        
        s = skymap.SkyMapLFSM()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        p.ComputeDirectionCosines()


class skymap_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lwa_user.skymap
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(skymap_tests))        
        
        
if __name__ == '__main__':
    unittest.main()
