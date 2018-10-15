# -*- coding: utf-8 -*-

"""Unit test for the lsl.writer.measurementset module."""

import os
import time
import unittest
import tempfile
import numpy
import shutil

from lsl.writer import measurementset
from lsl.common import stations as lwa_common
from lsl.correlator import uvUtils

run_tests = False
try:
    import casacore
    run_tests = True
except ImportError:
    pass


__version__  = "0.1"
__revision__ = "$Rev$"
__author__   = "Jayce Dowell"


@unittest.skipUnless(run_tests, "requires the 'casacore' module")
class measurementset_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.measurementset.MS
    class."""
    
    testPath = None
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""
        
        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-measurementset-', suffix='.tmp')
        
    def __initData(self):
        """Private function to generate a random set of data for writing a UVFITS
        file.  The data is returned as a dictionary with keys:
         * freq - frequency array in Hz
         * site - lwa.common.stations object
         * stands - array of stand numbers
         * bl - list of baseline pairs in real stand numbers
         * vis - array of visibility data in baseline x freq format
        """

        # Frequency range
        freq = numpy.arange(0,512)*20e6/512 + 40e6
        # Site and stands
        site = lwa_common.lwa1
        antennas = site.antennas[0:40:2]
        
        # Set baselines and data
        blList = uvUtils.get_baselines(antennas, include_auto=True, indicies=False)
        visData = numpy.random.rand(len(blList), len(freq))
        visData = visData.astype(numpy.complex64)
        
        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}
        
    def test_write_tables(self):
        """Test if the MeasurementSet writer writes all of the tables."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'ms-test-W.ms')
        
        # Get some data
        data = self.__initData()
        
        # Start the table
        tbl = measurementset.MS(testFile, ref_time=testTime)
        tbl.set_stokes(['xx'])
        tbl.set_frequency(data['freq'])
        tbl.set_geometry(data['site'], data['antennas'])
        tbl.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        tbl.write()
        
        # Make sure everyone is there
        self.assertTrue(os.path.exists(testFile))
        for tbl in ('ANTENNA', 'DATA_DESCRIPTION', 'FEED', 'FIELD', 'FLAG_CMD', 'HISTORY', 
                    'OBSERVATION', 'POINTING', 'POLARIZATION', 'PROCESSOR', 'SOURCE', 
                    'SPECTRAL_WINDOW', 'STATE'):
            self.assertTrue(os.path.exists(os.path.join(testFile, tbl)))
            
    def test_main_table(self):
        """Test the primary data table."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'ms-test-UV.ms')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = measurementset.MS(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the table and examine
        ms = casacore.tables.table(testFile, ack=False)
        uvw  = ms.getcol('UVW')
        ant1 = ms.getcol('ANTENNA1')
        ant2 = ms.getcol('ANTENNA2')
        vis  = ms.getcol('DATA')
        
        ms2 = casacore.tables.table(os.path.join(testFile, 'ANTENNA'), ack=False)
        mapper = ms2.getcol('NAME')
        mapper = [int(m[3:], 10) for m in mapper]
        
        # Correct number of visibilities
        self.assertEqual(uvw.shape[0], data['vis'].shape[0])
        self.assertEqual(vis.shape[0], data['vis'].shape[0])
        
        # Correct number of uvw coordinates
        self.assertEqual(uvw.shape[1], 3)
        
        # Correct number of frequencies
        self.assertEqual(vis.shape[1], data['freq'].size)
            
        # Correct values
        for row in xrange(uvw.shape[0]):
            stand1 = ant1[row]
            stand2 = ant2[row]
            visData = vis[row,:,0]
           
            # Find out which visibility set in the random data corresponds to the 
            # current visibility
            i = 0
            for a1,a2 in data['bl']:
                if a1.stand.id == mapper[stand1] and a2.stand.id == mapper[stand2]:
                    break
                else:
                    i = i + 1
                    
            # Run the comparison
            for vd, sd in zip(visData, data['vis'][i,:]):
                self.assertAlmostEqual(vd, sd, 8)
            i = i + 1
            
        ms.close()
        ms2.close()
        
    def tearDown(self):
        """Remove the test path directory and its contents"""
        
        shutil.rmtree(self.testPath, ignore_errors=True)


class measurementset_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(measurementset_tests))


if __name__ == '__main__':
    unittest.main()
