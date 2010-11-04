# -*- coding: utf-8 -*-

"""
Unit test for lwa_user.skymap module.
"""

##########################################################
# $Id: test_skymap.py 88 2010-05-20 15:39:17Z dwood $
##########################################################


import unittest

from lsl import skymap, astro


__revision__  = "$Revision: 88 $"
__version__   = "dev"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class skymap_tests(unittest.TestCase):
    """
    A unittest.TestCase collection of unit tests for the lwa_user.skymap
    module.
    """
    
    def test_SkyMap_init(self):
        """
        Test skymap.SkyMap class constructor method.
        """ 
        
        skymap.SkyMap()
        
        
    def test_SkyMap_ComputeTotalPowerFromSky(self):
        """
        Test skymap.SkyMap.ComputeTotalPowerFromSky() method.
        """
        
        s = skymap.SkyMap()
        s.ComputeTotalPowerFromSky()
        
    def test_ProjectedSkyMap_init(self):
        """
        Test skymap.ProjectedSkyMap constructor method.
        """
        
        s = skymap.SkyMap()
        skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
     
       
    def test_ProjectedSkyMap_ComputeTotalPowerFromVisibleSky(self):
        """
        Test skymap.ProjectedSkyMap.ComputeTotalPowerFromVisibleSky() method.
        """
        
        s = skymap.SkyMap()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys()) 
        p.ComputeTotalPowerFromVisibleSky()
        
        
    def test_ProjectedSkyMap_ComputeDirectionCosines(self):
        """
        Test skymap.ProjectedSkyMap.ComputeDirectionCosines() method.
        """
        
        s = skymap.SkyMap()
        p = skymap.ProjectedSkyMap(s, 20, 30, astro.get_julian_from_sys())
        p.ComputeDirectionCosines()
        
                
        
class skymap_test_suite(unittest.TestSuite):
    """
    A unittest.TestSuite class which contains all of the lwa_user.skymap
    module unit tests.
    """
    
    def __init__(self):
    
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(skymap_tests))        
        
        
        
if __name__ == '__main__':

    unittest.main()
    
