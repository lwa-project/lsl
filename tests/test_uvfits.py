# -*- coding: utf-8 -*-

"""Unit test for the lsl.writer.uvfits modules."""

import os
import time
import unittest
import tempfile
import numpy
import shutil
from astropy.io import fits as astrofits

from lsl.common import stations as lwa_common
from lsl.correlator import uvutil
from lsl.writer import uvfits


__version__  = "0.2"
__revision__ = "$Rev$"
__author__   = "Jayce Dowell"


class uvfits_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.uvfits.uv
    class."""
    
    testPath = None
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""
        
        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-uvfits-', suffix='.tmp')
        
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
        blList = uvutil.get_baselines(antennas, include_auto=True, indicies=False)
        visData = numpy.random.rand(len(blList), len(freq))
        visData = visData.astype(numpy.complex64)
        
        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}

    def test_write_tables(self):
        """Test if the UVFITS writer writes all of the tables."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-W.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        # Check that all of the extensions are there
        extNames = [hdu.name for hdu in hdulist]
        for ext in ['AIPS AN', 'AIPS FQ', 'AIPS SU']:
            self.assertTrue(ext in extNames)
            
        hdulist.close()
        
    def test_array_geometry(self):
        """Test the 'AIPS AN' table, part 1."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-AG.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        ag = hdulist['AIPS AN'].data
        # Correct number of stands
        self.assertEqual(len(data['antennas']), len(ag.field('NOSTA')))
        
        # Correct stand names
        names = ['LWA%03i' % ant.stand.id for ant in data['antennas']]
        for name, anname in zip(names, ag.field('ANNAME')):
            self.assertEqual(name, anname)
            
        hdulist.close()
        
    def test_antenna(self):
        """Test the 'AIPS AN' table, part 2."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-AN.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        an = hdulist['AIPS AN'].data
        # Correct number of stands
        self.assertEqual(len(data['antennas']), len(an.field('NOSTA')))
        
        hdulist.close()
        
    def test_frequency(self):
        """Test the 'AIPS FQ' table."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-FQ.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        fq = hdulist['AIPS FQ'].data
        # Correct number of IFs
        self.assertEqual(len(fq.field('FRQSEL')), 1)
        # Correct channel width
        self.assertEqual(fq.field('CH WIDTH')[0], data['freq'][1]-data['freq'][0])
        # Correct bandwidth
        self.assertEqual(fq.field('TOTAL BANDWIDTH')[0], numpy.abs(data['freq'][-1]-data['freq'][0]).astype(numpy.float32), 4)
        
        hdulist.close()
        
    def test_source(self):
        """Test the 'AIPS SU' table."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-SU.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        su = hdulist['AIPS SU'].data
        
        hdulist.close()
        
    def test_bandpass(self):
        """Test the 'AIPS BP' table."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-BP.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        su = hdulist['AIPS BP'].data
        
        hdulist.close()
        
    def test_uvdata(self):
        """Test the primary data table."""
        
        testTime = time.time()
        testFile = os.path.join(self.testPath, 'uv-test-UV.fits')
        
        # Get some data
        data = self.__initData()
        
        # Start the file
        fits = uvfits.Uv(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(testTime, 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the file and examine
        hdulist = astrofits.open(testFile)
        uv = hdulist[0]
        
        # Load the mapper
        try:
            mp = hdulist['NOSTA_MAPPER'].data
            nosta = mp.field('NOSTA')
            noact = mp.field('NOACT')
        except KeyError:
            ag = hdulist['AIPS AN'].data
            nosta = ag.field('NOSTA')
            noact = ag.field('NOSTA')
        mapper = {}
        for s,a in zip(nosta, noact):
            mapper[s] = a
            
        # Correct number of visibilities
        self.assertEqual(len(uv.data), data['vis'].shape[0])
        
        # Correct number of frequencies
        for row in uv.data:
            vis = row['DATA'][0,0,:,:,:]
            self.assertEqual(len(vis), len(data['freq']))
            
        # Correct values
        for row in uv.data:
            bl = int(row['BASELINE'])
            vis = row['DATA'][0,0,:,:,:]
            
            # Unpack the baseline
            if bl >= 65536:
                a1 = int((bl - 65536) / 2048)
                a2 = int((bl - 65536) % 2048)
            else:
                a1 = int(bl / 256)
                a2 = int(bl % 256)
            
            # Convert mapped stands to real stands
            stand1 = mapper[a1]
            stand2 = mapper[a2]
            
            # Find out which visibility set in the random data corresponds to the 
            # current visibility
            i = 0
            for ant1,ant2 in data['bl']:
                if ant1.stand.id == stand1 and ant2.stand.id == stand2:
                    break
                else:
                    i = i + 1
                    
            # Extract the data and run the comparison
            visData = numpy.zeros(len(data['freq']), dtype=numpy.complex64)
            visData.real = vis[:,0,0]
            visData.imag = vis[:,0,1]
            for vd, sd in zip(visData, data['vis'][i,:]):
                self.assertAlmostEqual(vd, sd, 8)
            i = i + 1
            
        hdulist.close()
        
    def tearDown(self):
        """Remove the test path directory and its contents"""
        
        shutil.rmtree(self.testPath, ignore_errors=True)


class uvfits_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(uvfits_tests))


if __name__ == '__main__':
    unittest.main()
