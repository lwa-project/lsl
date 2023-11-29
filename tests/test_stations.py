"""
Unit test for the lsl.common.stations module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import ephem
import numpy
import pickle
import unittest
from datetime import datetime

from lsl.common.paths import DATA_BUILD
from lsl.common import stations, dp, mcs, sdf, metabundle, sdm
import lsl.testing


__version__  = "0.5"
__author__    = "Jayce Dowell"


class stations_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.stations
    module."""

    def test_station(self):
        """Test retrieving stations from the stations module."""

        for station in (stations.lwa1, stations.lwasv):
            with self.subTest(station=station.name):
                self.assertTrue(isinstance(station, stations.LWAStation))
                
    def test_observer(self):
        """Test the ephem.Observer portion of an LWAStation."""
        
        lwa1 = stations.lwa1
        jov = ephem.Jupiter()
        
        lwa1.date = '2013/7/10 22:07:07'
        lwa1.compute(jov)
        
        # RA/Dec
        self.assertAlmostEqual(jov.ra,  ephem.hours('6:14:41.01'), 6)
        self.assertAlmostEqual(jov.dec, ephem.degrees('23:11:49.1'), 6)
        
        #Az/Alt
        self.assertAlmostEqual(jov.az,  ephem.degrees('274:40:27.7'), 6)
        self.assertAlmostEqual(jov.alt, ephem.degrees( '37:24:10.5'), 6)
        
    def test_pickle(self):
        """Test pickling of LWAStation instances."""
        
        for station in (stations.lwa1, stations.lwasv):
            with self.subTest(station=station.name):
                # Pickle and re-load
                out  = pickle.dumps(station)
                stationPrime = pickle.loads(out)
                
                # Test similarity
                self.assertAlmostEqual(station.lat, stationPrime.lat)
                self.assertAlmostEqual(station.long, stationPrime.long)
                self.assertAlmostEqual(station.elev, stationPrime.elev)
                for i in range(len(station.antennas)):
                    self.assertEqual(station.antennas[i].id, stationPrime.antennas[i].id)
                    self.assertEqual(station.antennas[i].stand.id, stationPrime.antennas[i].stand.id)
                    self.assertEqual(station.antennas[i].digitizer, stationPrime.antennas[i].digitizer)
                self.assertEqual(station.interface.mcs, stationPrime.interface.mcs)
                self.assertEqual(station.interface.sdf, stationPrime.interface.sdf)
                
                # Check independence
                stationPrime.antennas[100].stand.id = 888
                self.assertTrue(station.antennas[100].stand.id != stationPrime.antennas[100].stand.id)
            
    def test_strings(self):
        """Test string representations in the stations module."""
        
        lwa1 = stations.lwa1
        str(lwa1)
        repr(lwa1)
        
        str(lwa1.antennas[0])
        repr(lwa1.antennas[0])
        
        str(lwa1.antennas[0].fee)
        repr(lwa1.antennas[0].fee)
        
        str(lwa1.antennas[0].stand)
        repr(lwa1.antennas[0].stand)
        
        str(lwa1.antennas[0].cable)
        repr(lwa1.antennas[0].cable)
        
        str(lwa1.antennas[0].arx)
        repr(lwa1.antennas[0].arx)
        
        str(lwa1.interface)
        repr(lwa1.interface)
        
    def test_ecef_conversion(self):
        """Test the stations.geo_to_ecef() function."""

        lat = 0.0
        lng = 0.0
        elev = 0.0
        x, y, z = stations.geo_to_ecef(lat, lng, elev)
        numpy.testing.assert_allclose((x,y,z), (6378137.0,0.0,0.0))
        
    def test_interfaces(self):
        """Test retrieving LSL interface information."""
        
        lwa1 = stations.lwa1
        self.assertEqual(lwa1.interface.backend, 'lsl.common.dp')
        self.assertEqual(lwa1.interface.mcs, 'lsl.common.mcs')
        self.assertEqual(lwa1.interface.sdf, 'lsl.common.sdf')
        self.assertEqual(lwa1.interface.metabundle, 'lsl.common.metabundle')
        self.assertEqual(lwa1.interface.sdm, 'lsl.common.sdm')
        
        lwasv = stations.lwasv
        self.assertEqual(lwasv.interface.backend, 'lsl.common.adp')
        self.assertEqual(lwasv.interface.mcs, 'lsl.common.mcsADP')
        self.assertEqual(lwasv.interface.sdf, 'lsl.common.sdfADP')
        self.assertEqual(lwasv.interface.metabundle, 'lsl.common.metabundleADP')
        self.assertEqual(lwasv.interface.sdm, 'lsl.common.sdmADP')
        
        lwana = stations.lwana
        self.assertEqual(lwana.interface.backend, None)
        self.assertEqual(lwana.interface.mcs, None)
        self.assertEqual(lwana.interface.sdf, None)
        self.assertEqual(lwana.interface.metabundle, None)
        self.assertEqual(lwana.interface.sdm, None)
        
    def test_interface_modules(self):
        """Test retrieving LSL interface modules."""
        
        lwa1 = stations.lwa1
        self.assertEqual(lwa1.interface.get_module('backend').__file__, dp.__file__)
        self.assertEqual(lwa1.interface.get_module('mcs').__file__, mcs.__file__)
        self.assertEqual(lwa1.interface.get_module('sdf').__file__, sdf.__file__)
        self.assertEqual(lwa1.interface.get_module('metabundle').__file__, metabundle.__file__)
        self.assertEqual(lwa1.interface.get_module('sdm').__file__, sdm.__file__)
        
        lwasv = stations.lwasv
        self.assertFalse(lwasv.interface.get_module('backend').__file__ == dp.__file__)
        self.assertFalse(lwasv.interface.get_module('mcs').__file__ == mcs.__file__)
        self.assertFalse(lwasv.interface.get_module('sdf').__file__ == sdf.__file__)
        self.assertFalse(lwasv.interface.get_module('metabundle').__file__ == metabundle.__file__)
        self.assertFalse(lwasv.interface.get_module('sdm').__file__ == sdm.__file__)
        
    def test_ssmif(self):
        """Test the SSMIF parser."""
        
        filenames = [os.path.join(DATA_BUILD, 'lwa1-ssmif.txt'),
                     os.path.join(DATA_BUILD, 'lwasv-ssmif.txt'),
                     os.path.join(DATA_BUILD, 'tests', 'ssmif.dat'),
                     os.path.join(DATA_BUILD, 'tests', 'ssmif-adp.dat')]
        sites = ['LWA1', 'LWA-SV', 'LWA1', 'LWA-SV']
        types = ['text', 'text', 'binary', 'binary']
        for filename,site,type in zip(filenames, sites, types):
            with self.subTest(station=site, type=type):
                out = stations.parse_ssmif(filename)
                
    def test_responses(self):
        """Test the various frequency responses."""
        
        for station in (stations.lwa1, stations.lwasv):
            with self.subTest(station=station.name):
                station[0].fee.response()
                station[0].cable.response()
                
            for filt in ('split', 'full', 'reduced', 'split@3MHz', 'full@3MHz'):
                with self.subTest(station=station.name, filter=filt):
                    station[0].arx.response(filt)


class stations_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.stations
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(stations_tests)) 


if __name__ == '__main__':
    unittest.main()
