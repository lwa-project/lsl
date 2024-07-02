"""
Unit test for the lsl.writer.virtual module.
"""

import os
import time
import ephem
import unittest
import numpy as np

from lsl.common import stations as lwa_common
from lsl.correlator import uvutils
from lsl.writer import virtual
from lsl.astro import unix_to_taimjd
import lsl.testing


__version__  = "0.1"
__author__   = "Jayce Dowell"


class virtual_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.virtual.VirtualWriter
    class."""
    
    def setUp(self):
        """Turn off all numpy warnings."""

        np.seterr(all='ignore')
        
    def _init_data(self):
        """Private function to generate a random set of data for writing a FITS
        IDI file.  The data is returned as a dictionary with keys:
         * freq - frequency array in Hz
         * site - lwa.common.stations object
         * stands - array of stand numbers
         * bl - list of baseline pairs in real stand numbers
         * vis - array of visibility data in baseline x freq format
        """

        # Frequency range
        freq = np.arange(0,512)*20e6/512 + 40e6
        # Site and stands
        site = lwa_common.lwa1
        antennas = site.antennas[0:40:2]
        
        # Set baselines and data
        blList = uvutils.get_baselines(antennas, include_auto=True, indicies=False)
        visData = np.random.rand(len(blList), len(freq))
        visData = visData.astype(np.complex64)

        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}

    def test_write_tables(self):
        """Test if the FITS IDI writer writes all of the tables."""

        testTime = time.time() - 2*86400.0
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = virtual.VirtualWriter(ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.set_observer('Dowell, Jayce', 'LD009', 'Test')
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        
        # Examine
        ds = fits.get_data_set(1)
        
    def test_multi_source(self):
        testTime = time.time() - 2*86400.0
        
        # Get some data
        data = self._init_data()
        
        # Create a source other than zenith to try
        source = ephem.FixedBody()
        source._ra = 0.0
        source._dec = np.pi/2
        source._epoch = ephem.J2000
        source.compute(data['site'])
        
        # Start the file
        fits = virtual.VirtualWriter(ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.set_observer('Dowell, Jayce', 'LD009', 'Test')
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.add_data_set(unix_to_taimjd(testTime+6.0), 6.0, data['bl'], data['vis'],
                          source=source)
        
        # Examine
        ds1 = fits.get_data_set(1)
        ds2 = fits.get_data_set(2)
        

class virtual_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.writer.virtual
    unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(virtual_tests))


if __name__ == '__main__':
    unittest.main()
