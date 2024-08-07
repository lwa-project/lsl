"""
Unit test for the lsl.writer.fitsidi modules.
"""

import os
import time
import ephem
import unittest
import tempfile
import numpy as np
import shutil

from astropy.io import fits as astrofits
from astropy.coordinates import EarthLocation

from lsl.common import stations as lwa_common
from lsl.correlator import uvutils
from lsl.writer import fitsidi
from lsl.astro import unix_to_taimjd
import lsl.testing


__version__  = "0.2"
__author__   = "Jayce Dowell"


class fitsidi_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.fitsidi.idi
    class."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        np.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-fitsidi-', suffix='.tmp')

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
        blList = uvutils.get_baselines(antennas, include_auto=True)
        visData = np.random.rand(len(blList), len(freq))
        visData = visData.astype(np.complex64)

        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}

    def test_write_tables(self):
        """Test if the FITS IDI writer writes all of the tables."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-W.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.set_observer('Dowell, Jayce', 'LD009', 'Test')
        fits.add_header_keyword('EXAMPLE', 'example keyword')
        fits.add_comment('This is a comment')
        fits.add_history('This is history')
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        # Check that all of the extensions are there
        extNames = [hdu.name for hdu in hdulist]
        for ext in ['ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'BANDPASS', 'SOURCE', 'UV_DATA']:
            with self.subTest(table=ext):
                self.assertTrue(ext in extNames)
        # Check header values that we set
        self.assertEqual('Dowell, Jayce', hdulist[0].header['OBSERVER'])
        self.assertEqual('LD009', hdulist[0].header['PROJECT'])
        self.assertEqual('Test', hdulist[0].header['LWATYPE'])
        self.assertEqual('example keyword', hdulist[0].header['EXAMPLE'])
        # Check the comments and history
        self.assertTrue('This is a comment' in str(hdulist[0].header['COMMENT']).split('\n'))
        self.assertTrue('This is history' in str(hdulist[0].header['HISTORY']).split('\n'))
        
        hdulist.close()
        
    def test_writer_errors(self):
        """Test that common FITS IDI error conditions are caught."""
        
        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-ERR.fits')
        
        # Get some data
        data = self._init_data()
        
        for i in range(4):
            # Start the file
            fits = fitsidi.Idi(testFile, ref_time=testTime, overwrite=True)
            if i != 0:
                fits.set_stokes(['xx'])
            if i != 1:
                fits.set_frequency(data['freq'])
            if i != 2:
                fits.set_geometry(data['site'], data['antennas'])
            if i != 3:
                fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
            self.assertRaises(RuntimeError, fits.write)
            fits.close()
            
    def test_array_geometry(self):
        """Test the ARRAY_GEOMETRY table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-AG.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        
        # Correct station location
        ag = hdulist['ARRAY_GEOMETRY'].header
        site_ecef = data['site'].geocentric_location
        self.assertAlmostEqual(ag['ARRAYX'], site_ecef[0], 3)
        self.assertAlmostEqual(ag['ARRAYY'], site_ecef[1], 3)
        self.assertAlmostEqual(ag['ARRAYZ'], site_ecef[2], 3)
        
        # Correct number of stands
        ag = hdulist['ARRAY_GEOMETRY'].data
        self.assertEqual(len(data['antennas']), len(ag.field('NOSTA')))

        # Correct stand names
        names = ['LWA%03i' % ant.stand.id for ant in data['antennas']]
        for name, anname in zip(names, ag.field('ANNAME')):
            self.assertEqual(name, anname)
            
        # Correct stand locations
        for ant, ecef in zip(data['antennas'], ag.field('STABXYZ')):
            el = EarthLocation.from_geocentric(*ecef, unit='m')
            enz = data['site'].get_enz_offset(el)
            self.assertAlmostEqual(enz[0], ant.stand.x, 3)
            self.assertAlmostEqual(enz[1], ant.stand.y, 3)
            self.assertAlmostEqual(enz[2], ant.stand.z, 3)
            
        hdulist.close()

    def test_frequency(self):
        """Test the FREQUENCY table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-FQ.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        fq = hdulist['FREQUENCY'].data
        # Correct number of FREQIDs
        self.assertEqual(len(fq.field('FREQID')), 1)

        # Correct channel width
        self.assertAlmostEqual(fq.field('CH_WIDTH')[0], np.abs(data['freq'][1]-data['freq'][0]), 4)

        # Correct bandwidth
        self.assertAlmostEqual(fq.field('TOTAL_BANDWIDTH')[0], np.abs(data['freq'][-1]-data['freq'][0]).astype(np.float32), 4)

        # Correct sideband
        self.assertEqual(fq.field('SIDEBAND')[0], 1)

        hdulist.close()

    def test_antenna(self):
        """Test the ANTENNA table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-AN.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        an = hdulist['ANTENNA'].data
        # Correct number of stands
        self.assertEqual(len(data['antennas']), len(an.field('ANTENNA_NO')))

        # Correct FREQIDs
        for freqid in an.field('FREQID'):
            self.assertEqual(freqid, 1)

        hdulist.close()

    def test_bandpass(self):
        """Test the BANDPASS table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-BP.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        bp = hdulist['BANDPASS'].data
        # Correct number of entries
        self.assertEqual(len(data['antennas']), len(bp.field('ANTENNA_NO')))

        # Correct Source ID number
        for src in bp.field('SOURCE_ID'):
            self.assertEqual(src, 0)

        # Correct FREQIDs
        for freqid in bp.field('FREQID'):
            self.assertEqual(freqid, 1)

        hdulist.close()

    def test_source(self):
        """Test the SOURCE table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-SO.fits')
        
        # Get some data
        data = self._init_data()
        
        # Create a source other than zenith to try
        source = ephem.FixedBody()
        source._ra = 0.0
        source._dec = np.pi/2
        source._epoch = ephem.J2000
        source.compute(data['site'])
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.add_data_set(unix_to_taimjd(testTime+6.0), 6.0, data['bl'], data['vis'],
                          source=source)
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        so = hdulist['SOURCE'].data
        # Correct number of entries
        self.assertEqual(len(so.field('SOURCE_ID')), 2)

        # Correct Source ID number
        self.assertEqual(so.field('SOURCE_ID')[0], 1)
        self.assertEqual(so.field('SOURCE_ID')[1], 2)

        hdulist.close()

    def test_uvdata(self):
        """Test the UV_DATA table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-UV.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        uv = hdulist['UV_DATA'].data

        # Load the mapper
        try:
            mp = hdulist['NOSTA_MAPPER'].data
            nosta = mp.field('NOSTA')
            noact = mp.field('NOACT')
        except KeyError:
            ag = hdulist['ARRAY_GEOMETRY'].data
            nosta = ag.field('NOSTA')
            noact = ag.field('NOSTA')
        mapper = {}
        for s,a in zip(nosta, noact):
            mapper[s] = a

        # Correct number of visibilities
        self.assertEqual(len(uv.field('FLUX')), data['vis'].shape[0])
        
        # Correct number of frequencies
        for vis in uv.field('FLUX'):
            self.assertEqual(len(vis), 2*len(data['freq']))

        # Correct values
        for bl, vis in zip(uv.field('BASELINE'), uv.field('FLUX')):
            # Convert mapped stands to real stands
            stand1 = mapper[(bl >> 8) & 255]
            stand2 = mapper[bl & 255]

            # Find out which visibility set in the random data corresponds to the 
            # current visibility
            i = 0
            for ant1,ant2 in data['bl']:
                if ant1.stand.id == stand1 and ant2.stand.id == stand2:
                    break
                else:
                    i = i + 1
            
            # Extract the data and run the comparison
            visData = np.zeros(len(data['freq']), dtype=np.complex64)
            visData.real = vis[0::2]
            visData.imag = vis[1::2]
            np.testing.assert_allclose(visData, data['vis'][i,:])
            i = i + 1
        
        hdulist.close()
        
    def test_mapper(self):
        """Test the NOSTA_MAPPER table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-SM.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        extNames = [hdu.name for hdu in hdulist]
        maxStand = -1
        for ant in data['antennas']:
            if ant.stand.id > maxStand:
                maxStand = ant.stand.id
        if maxStand > 255:
            self.assertTrue('NOSTA_MAPPER' in extNames)

            # Make sure the mapper makes sense
            mp = hdulist['NOSTA_MAPPER'].data
            ag = hdulist['ARRAY_GEOMETRY'].data
            mNoSta = mp.field('NOSTA')
            aNoSta = ag.field('NOSTA')
            mNoAct = mp.field('NOACT')
            aAnNam = ag.field('ANNAME')
            for msta, mact, asta, anam in zip(mNoSta, mNoAct, aNoSta, aAnNam):
                self.assertEqual(msta, asta)
                self.assertEqual(mact, int(anam[3:]))

        hdulist.close()
        
    def test_multi_if(self):
        """Test writing more than one IF to a FITS IDI file."""
        
        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-MultiIF.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Idi(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_frequency(data['freq']+10e6)
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], 
                          np.concatenate([data['vis'], 10*data['vis']], axis=1))
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        # Check that all of the extensions are there
        extNames = [hdu.name for hdu in hdulist]
        for ext in ['ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'BANDPASS', 'SOURCE', 'UV_DATA']:
            self.assertTrue(ext in extNames)
            
        # Load the mapper
        try:
            mp = hdulist['NOSTA_MAPPER'].data
            nosta = mp.field('NOSTA')
            noact = mp.field('NOACT')
        except KeyError:
            ag = hdulist['ARRAY_GEOMETRY'].data
            nosta = ag.field('NOSTA')
            noact = ag.field('NOSTA')
        mapper = {}
        for s,a in zip(nosta, noact):
            mapper[s] = a
            
        # Correct number of visibilities
        uv = hdulist['UV_DATA'].data
        self.assertEqual(len(uv.field('FLUX')), data['vis'].shape[0])
        
        # Correct number of frequencies
        for vis in uv.field('FLUX'):
            self.assertEqual(len(vis), 2*2*len(data['freq']))

        # Correct values
        for bl, vis in zip(uv.field('BASELINE'), uv.field('FLUX')):
            # Convert mapped stands to real stands
            stand1 = mapper[(bl >> 8) & 255]
            stand2 = mapper[bl & 255]

            # Find out which visibility set in the random data corresponds to the 
            # current visibility
            i = 0
            for ant1,ant2 in data['bl']:
                if ant1.stand.id == stand1 and ant2.stand.id == stand2:
                    break
                else:
                    i = i + 1
            
            # Extract the data and run the comparison - IF 1
            visData = np.zeros(2*len(data['freq']), dtype=np.complex64)
            visData.real = vis[0::2]
            visData.imag = vis[1::2]
            np.testing.assert_allclose(visData[:len(data['freq'])], data['vis'][i,:])
            
            # Extract the data and run the comparison - IF 2
            visData = np.zeros(2*len(data['freq']), dtype=np.complex64)
            visData.real = vis[0::2]
            visData.imag = vis[1::2]
            np.testing.assert_allclose(visData[len(data['freq']):], 10*data['vis'][i,:])
                
        hdulist.close()
        
    def tearDown(self):
        """Remove the test path directory and its contents"""

        shutil.rmtree(self.testPath, ignore_errors=True)


class aipsidi_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.writer.fitsidi.aips
    class."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        np.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-aipsidi-', suffix='.tmp')

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
        blList = uvutils.get_baselines(antennas, include_auto=True)
        visData = np.random.rand(len(blList), len(freq))
        visData = visData.astype(np.complex64)

        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}

    def test_write_tables(self):
        """Test if the AIPS IDI writer writes all of the tables."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-W.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_comment('This is a comment')
        fits.add_history('This is history')
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        # Check that all of the extensions are there
        extNames = [hdu.name for hdu in hdulist]
        for ext in ['ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'BANDPASS', 'SOURCE', 'UV_DATA']:
            self.assertTrue(ext in extNames)
        # Check the comments and history
        self.assertTrue('This is a comment' in str(hdulist[0].header['COMMENT']).split('\n'))
        self.assertTrue('This is history' in str(hdulist[0].header['HISTORY']).split('\n'))
        
        hdulist.close()

    def test_array_geometry(self):
        """Test the AIPS IDI ARRAY_GEOMETRY table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-AG.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        
        # Correct station location
        ag = hdulist['ARRAY_GEOMETRY'].header
        site_ecef = data['site'].geocentric_location
        self.assertAlmostEqual(ag['ARRAYX'], site_ecef[0], 3)
        self.assertAlmostEqual(ag['ARRAYY'], site_ecef[1], 3)
        self.assertAlmostEqual(ag['ARRAYZ'], site_ecef[2], 3)
        
        # Correct number of stands
        ag = hdulist['ARRAY_GEOMETRY'].data
        self.assertEqual(len(data['antennas']), len(ag.field('NOSTA')))

        # Correct stand names
        names = ['L%03i' % ant.stand.id for ant in data['antennas']]
        for name, anname in zip(names, ag.field('ANNAME')):
            self.assertEqual(name, anname)
            
        # Correct stand locations
        for ant, ecef in zip(data['antennas'], ag.field('STABXYZ')):
            el = EarthLocation.from_geocentric(*ecef, unit='m')
            enz = data['site'].get_enz_offset(el)
            self.assertAlmostEqual(enz[0], ant.stand.x, 3)
            self.assertAlmostEqual(enz[1], ant.stand.y, 3)
            self.assertAlmostEqual(enz[2], ant.stand.z, 3)
            
        hdulist.close()

    def test_frequency(self):
        """Test the AIPS IDI FREQUENCY table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-FQ.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        fq = hdulist['FREQUENCY'].data
        # Correct number of FREQIDs
        self.assertEqual(len(fq.field('FREQID')), 1)

        # Correct channel width
        self.assertAlmostEqual(fq.field('CH_WIDTH')[0], np.abs(data['freq'][1]-data['freq'][0]), 4)

        # Correct bandwidth
        self.assertAlmostEqual(fq.field('TOTAL_BANDWIDTH')[0], np.abs(data['freq'][-1]-data['freq'][0]).astype(np.float32), 4)

        # Correct sideband
        self.assertEqual(fq.field('SIDEBAND')[0], 1)

        hdulist.close()

    def test_antenna(self):
        """Test the AIPS IDI ANTENNA table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-AN.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        an = hdulist['ANTENNA'].data
        # Correct number of stands
        self.assertEqual(len(data['antennas']), len(an.field('ANTENNA_NO')))

        # Correct FREQIDs
        for freqid in an.field('FREQID'):
            self.assertEqual(freqid, 1)

        hdulist.close()

    def test_bandpass(self):
        """Test the AIPS IDI BANDPASS table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-BP.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        bp = hdulist['BANDPASS'].data
        # Correct number of entries
        self.assertEqual(len(data['antennas']), len(bp.field('ANTENNA_NO')))

        # Correct Source ID number
        for src in bp.field('SOURCE_ID'):
            self.assertEqual(src, 0)

        # Correct FREQIDs
        for freqid in bp.field('FREQID'):
            self.assertEqual(freqid, 1)

        hdulist.close()

    def test_source(self):
        """Test the AIPS IDI SOURCE table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-SO.fits')
        
        # Get some data
        data = self._init_data()
        
        # Create a source other than zenith to try
        source = ephem.FixedBody()
        source._ra = 0.0
        source._dec = np.pi/2
        source._epoch = ephem.J2000
        source.compute(data['site'])
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.add_data_set(unix_to_taimjd(testTime+6.0), 6.0, data['bl'], data['vis'],
                          source=source)
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        so = hdulist['SOURCE'].data
        # Correct number of entries
        self.assertEqual(len(so.field('SOURCE_ID')), 2)

        # Correct Source ID number
        self.assertEqual(so.field('SOURCE_ID')[0], 1)
        self.assertEqual(so.field('SOURCE_ID')[1], 2)

        hdulist.close()

    def test_uvdata(self):
        """Test the AIPS IDI UV_DATA table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-UV.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        uv = hdulist['UV_DATA'].data

        # Load the mapper
        try:
            mp = hdulist['NOSTA_MAPPER'].data
            nosta = mp.field('NOSTA')
            noact = mp.field('NOACT')
        except KeyError:
            ag = hdulist['ARRAY_GEOMETRY'].data
            nosta = ag.field('NOSTA')
            noact = ag.field('NOSTA')
        mapper = {}
        for s,a in zip(nosta, noact):
            mapper[s] = a

        # Correct number of visibilities
        self.assertEqual(len(uv.field('FLUX')), data['vis'].shape[0])
        
        # Correct number of frequencies
        for vis in uv.field('FLUX'):
            self.assertEqual(len(vis), 2*len(data['freq']))

        # Correct values
        for bl, vis in zip(uv.field('BASELINE'), uv.field('FLUX')):
            # Convert mapped stands to real stands
            stand1 = mapper[(bl >> 8) & 255]
            stand2 = mapper[bl & 255]

            # Find out which visibility set in the random data corresponds to the 
            # current visibility
            i = 0
            for ant1,ant2 in data['bl']:
                if ant1.stand.id == stand1 and ant2.stand.id == stand2:
                    break
                else:
                    i = i + 1
            
            # Extract the data and run the comparison
            visData = np.zeros(len(data['freq']), dtype=np.complex64)
            visData.real = vis[0::2]
            visData.imag = vis[1::2]
            np.testing.assert_allclose(visData, data['vis'][i,:])
            i = i + 1
        
        hdulist.close()

    def test_mapper(self):
        """Test the AIPS IDI NOSTA_MAPPER table."""

        testTime = time.time() - 2*86400.0
        testFile = os.path.join(self.testPath, 'idi-test-SM.fits')
        
        # Get some data
        data = self._init_data()
        
        # Start the file
        fits = fitsidi.Aips(testFile, ref_time=testTime)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(unix_to_taimjd(testTime), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()

        # Open the file and examine
        hdulist = astrofits.open(testFile)
        extNames = [hdu.name for hdu in hdulist]
        maxStand = -1
        for ant in data['antennas']:
            if ant.stand.id > maxStand:
                maxStand = ant.stand.id
        if maxStand > 99:
            self.assertTrue('NOSTA_MAPPER' in extNames)

            # Make sure the mapper makes sense
            mp = hdulist['NOSTA_MAPPER'].data
            ag = hdulist['ARRAY_GEOMETRY'].data
            mNoSta = mp.field('NOSTA')
            aNoSta = ag.field('NOSTA')
            mNoAct = mp.field('NOACT')
            aAnNam = ag.field('ANNAME')
            for msta, mact, asta, anam in zip(mNoSta, mNoAct, aNoSta, aAnNam):
                self.assertEqual(msta, asta)
                self.assertEqual(mact, int(anam[1:]))

        hdulist.close()
        
    def tearDown(self):
        """Remove the test path directory and its contents"""

        shutil.rmtree(self.testPath, ignore_errors=True)


class fitsidi_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.reader units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(fitsidi_tests))
        self.addTests(loader.loadTestsFromTestCase(aipsidi_tests))


if __name__ == '__main__':
    unittest.main()
