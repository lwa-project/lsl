# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range
    
"""Unit test for lsl.imaging modules"""

import os
import copy
import time
import numpy
import shutil
import tempfile
import unittest

from lsl import astro
from lsl.common.paths import DATA_BUILD
from lsl.imaging import utils
from lsl.imaging import selfCal
from lsl.imaging.data import VisibilityData
from lsl.writer.fitsidi import Idi, NUMERIC_STOKES
from lsl.sim.vis import SOURCES as simSrcs
from lsl.common.stations import lwa1, parse_ssmif
from lsl.correlator import uvutil

run_ms_tests = False
try:
    import casacore
    from lsl.writer.measurementset import Ms
    run_ms_tests = True
except ImportError:
    pass


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


uvFile = os.path.join(DATA_BUILD, 'tests', 'uv-test.fits')
idiFile = os.path.join(DATA_BUILD, 'tests', 'idi-test.fits')
idiAltFile = os.path.join(DATA_BUILD, 'tests', 'idi-test-alt.fits')
idiSSMIFFile = os.path.join(DATA_BUILD, 'tests', 'idi-test-alt.txt')


class imaging_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.imaging
    modules."""
    
    testPath = None

    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        numpy.seterr(all='ignore')
        self.testPath = tempfile.mkdtemp(prefix='test-imaging-', suffix='.tmp')

    def test_CorrelatedDataIDI(self):
        """Test the utils.CorrelatedDataIDI class."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedDataIDI(idiFile)
        
        # Dates
        self.assertEqual(idi.date_obs.strftime("%Y-%m-%dT%H:%M:%S"), "2013-03-04T20:36:26")
        
        # Stand and baseline counts
        self.assertEqual(len(idi.stands), 5)
        self.assertEqual(idi.total_baseline_count, 5*(5+1)/2)
        
        # Basic functions (just to see that they run)
        junk = idi.get_antennaarray()
        junk = idi.get_observer()
        junk = idi.get_data_set(1)
        
        # Error checking
        self.assertRaises(IndexError, idi.get_data_set, 2)
        
        # 'with' statement support
        with utils.CorrelatedDataIDI(idiFile) as idi:
            ## Basic functions (just to see that they run)
            junk = idi.get_antennaarray()
            junk = idi.get_observer()
            junk = idi.get_data_set(1)
            
        # generator support
        idi = utils.CorrelatedDataIDI(idiFile)
        i = 0
        for ds in idi.data_set_sequence(include_auto=False):
            self.assertEqual(ds.XX.data.shape[0], 5*(5-1)/2)
            i += 1
        self.assertEqual(i, 1)
        
        # both at the same time
        with utils.CorrelatedDataIDI(idiFile) as idi:
            i = 0
            for ds in idi.data_set_sequence(include_auto=False):
                self.assertEqual(ds.XX.data.shape[0], 5*(5-1)/2)
                i += 1
            self.assertEqual(i, 1)
            
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
        site = lwa1
        antennas = site.antennas[0:40:2]
        
        # Set baselines and data
        blList = uvutil.get_baselines(antennas, include_auto=True, indicies=False)
        visData = numpy.random.rand(len(blList), len(freq))
        visData = visData.astype(numpy.complex64)
        
        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}
        
    def test_CorrelatedDataIDI_MultiIF(self):
        """Test the utils.CorrelatedDataIDI class on a file with multiple IFs."""
        
        # Get some data
        data = self.__initData()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join('idi-test-MultiIF.fits')
        
        # Start the file
        fits = Idi(testFile, ref_time=testTime, clobber=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_frequency(data['freq']+30e6)
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], 
                          numpy.concatenate([data['vis'], 10*data['vis']], axis=1))
        fits.write()
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(testFile)
        self.assertEqual(idi.freq.size, 2*data['freq'].size)
        ds = idi.get_data_set(1, include_auto=True)
        
        for i in xrange(len(data['bl'])):
            for j in xrange(data['freq'].size):
                self.assertAlmostEqual(ds.XX.data[i,j], data['vis'][i,j], 6)
                self.assertAlmostEqual(ds.XX.data[i,j+data['freq'].size], 10*data['vis'][i,j], 6)
                
    def test_CorrelatedDataIDI_Alt(self):
        """Test the utils.CorrelatedDataIDI class on a file with an unusual telescope."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedDataIDI(idiAltFile)
        
        # Dates
        self.assertEqual(idi.date_obs.strftime("%Y-%m-%dT%H:%M:%S"), "2013-03-04T20:36:26")
        
        # Stand and baseline counts
        self.assertEqual(len(idi.stands), 5)
        self.assertEqual(idi.total_baseline_count, 5*(5+1)/2)
        self.assertEqual(idi.integration_count, 1)
        
        # Basic functions (just to see that they run)
        junk = idi.get_antennaarray()
        junk = idi.get_observer()
        junk = idi.get_data_set(1)
        
        # Error checking
        self.assertRaises(IndexError, idi.get_data_set, 2)
        
    def test_CorrelatedDataIDI_AltArrayGeometry(self):
        """Test the utils.CorrelatedDataIDI class on determing array geometry."""
        
        # Open the FITS IDI files
        idi1 = utils.CorrelatedData(idiFile)
        idi2 = utils.CorrelatedData(idiAltFile)
        
        # Dates
        self.assertEqual(idi1.date_obs.strftime("%Y-%m-%dT%H:%M:%S"), idi2.date_obs.strftime("%Y-%m-%dT%H:%M:%S"))
        
        # Stand and baseline counts
        self.assertEqual(len(idi1.stands), len(idi2.stands))
        self.assertEqual(idi1.total_baseline_count, idi2.total_baseline_count)
        self.assertEqual(idi1.integration_count, idi2.integration_count)
        
        # Check stands
        for s1,s2 in zip(idi1.stands, idi2.stands):
            self.assertEqual(s1, s2)
            
        # Check stations
        station1 = parse_ssmif(idiSSMIFFile)
        station2 = idi2.station
        self.assertAlmostEqual(station1.lat, station2.lat, 3)
        self.assertAlmostEqual(station1.lon, station2.lon, 3)
        self.assertAlmostEqual(station1.elev, station2.elev, 1)
        
        # Check antennas
        ants1 = [a for a in station1.antennas if a.pol == 0]
        ants2 = station2.antennas
        for a1,a2 in zip(ants1, ants2):
            self.assertEqual(a1.id, a2.id)
            self.assertEqual(a1.stand.id, a2.stand.id)
            self.assertAlmostEqual(a1.stand.x, a2.stand.x, 2)
            self.assertAlmostEqual(a1.stand.y, a2.stand.y, 2)
            self.assertAlmostEqual(a1.stand.z, a2.stand.z, 2)
            
    def test_CorrelatedDataUV(self):
        """Test the utils.CorrelatedDataUV class."""
        
        # Open the UVFITS file
        uv = utils.CorrelatedDataUV(uvFile)
        
        # Dates
        self.assertEqual(uv.date_obs.strftime("%Y-%m-%dT%H:%M:%S"), "2013-03-04T20:36:26")
        
        # Stand and baseline counts
        self.assertEqual(len(uv.stands), 5)
        self.assertEqual(uv.total_baseline_count, 5*(5+1)/2)
        self.assertEqual(uv.integration_count, 1)
        
        # Basic functions (just to see that they run)
        junk = uv.get_antennaarray()
        junk = uv.get_observer()
        junk = uv.get_data_set(1)
        
        # Error checking
        self.assertRaises(IndexError, uv.get_data_set, 2)
        
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS(self):
        """Test the utils.CorrelatedDataMS class."""
        
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-W.ms')
        
        # Get some data
        data = self.__initData()
        
        # Start the table
        tbl = Ms(testFile, ref_time=testTime)
        tbl.set_stokes(['xx'])
        tbl.set_frequency(data['freq'])
        tbl.set_geometry(data['site'], data['antennas'])
        tbl.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], data['vis'])
        tbl.write()
        
        # Open the measurement set
        ms = utils.CorrelatedDataMS(testFile)
        
        # Basic functions (just to see that they run)
        junk = ms.get_antennaarray()
        junk = ms.get_observer()
        junk = ms.get_data_set(1)
        
        # Error checking
        self.assertRaises(IndexError, ms.get_data_set, 2)
        
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS_SingleIF(self):
        """Test the utils.CorrelatedDataMS class on a file with a single IF."""
        
        # Get some data
        data = self.__initData()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-SingleIF.ms')
        
        # Start the file
        fits = Ms(testFile, ref_time=testTime, clobber=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], data['vis'])
        fits.write()
        
        # Open the measurement set
        ms = utils.CorrelatedDataMS(testFile)
        self.assertEqual(ms.freq.size, data['freq'].size)
        ds = ms.get_data_set(1, include_auto=True)
        
        for i in xrange(len(data['bl'])):
            for j in xrange(data['freq'].size):
                self.assertAlmostEqual(ds.XX.data[i,j], data['vis'][i,j], 6)
                
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS_MultiIF(self):
        """Test the utils.CorrelatedDataMS class on a file with multiple IFs."""
        
        # Get some data
        data = self.__initData()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-MultiIF.ms')
        
        # Start the file
        fits = Ms(testFile, ref_time=testTime, clobber=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_frequency(data['freq']+30e6)
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], 
                          numpy.concatenate([data['vis'], 10*data['vis']], axis=1))
        fits.write()
        
        # Open the measurement set
        ms = utils.CorrelatedDataMS(testFile)
        self.assertEqual(ms.freq.size, 2*data['freq'].size)
        ds = ms.get_data_set(1, include_auto=True)
        
        for i in xrange(len(data['bl'])):
            for j in xrange(data['freq'].size):
                self.assertAlmostEqual(ds.XX.data[i,j], data['vis'][i,j], 6)
                self.assertAlmostEqual(ds.XX.data[i,j+data['freq'].size], 10*data['vis'][i,j], 6)
                
    def test_sort(self):
        """Test the utils.sort function."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1, sort=False)
        
        # Sort
        dss = ds.copy()
        dss.sort()
        for pol in ds.pols:
            p0 = getattr(ds,  pol)
            p1 = getattr(dss, pol)
            self.assertEqual(p0.nbaseline, p1.nbaseline)
            
    def test_sort_alt(self):
        """Test the utils.sort function - alternate FITS IDI file."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiAltFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1, sort=False)
        
        # Sort
        dss = ds.copy()
        dss.sort()
        for pol in ds.pols:
            p0 = getattr(ds,  pol)
            p1 = getattr(dss, pol)
            self.assertEqual(p0.nbaseline, p1.nbaseline)
            
    def test_sort_uvfits(self):
        """Test the utils.sort function - UVFITS file."""
        
        # Open the FITS IDI file
        uv = utils.CorrelatedData(uvFile)
        
        # Get some data to sort
        ds = uv.get_data_set(1, sort=False)
        
        # Sort
        dss = ds.copy()
        dss.sort()
        for pol in ds.pols:
            p0 = getattr(ds,  pol)
            p1 = getattr(dss, pol)
            self.assertEqual(p0.nbaseline, p1.nbaseline)
            
    def test_prune(self):
        """Test the utils.get_uv_range function."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        
        # Prune
        dsp1 = ds.get_uv_range(min_uv=10)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p1 = getattr(dsp1, pol)
            self.assertTrue(p1.nbaseline < p0.nbaseline)
            
        # Auto-prune
        dsp2 = idi.get_data_set(1, min_uv=10)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p2 = getattr(dsp2, pol)
            self.assertTrue(p2.nbaseline < p0.nbaseline)
            
        # Auto-prune that should result in no baselines
        dsp3 = idi.get_data_set(1, min_uv=100)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p3 = getattr(dsp3, pol)
            self.assertEqual(p3.nbaseline, 0)
            
    def test_subset(self):
        """Test the utils.get_antenna_subset function."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        
        # Indicies
        ## Include
        dss1 = ds.get_antenna_subset(include=(0,1,2), indicies=True)
        self.assertEqual(len(dss1.baselines), 3)
        for a1,a2 in dss1.baselines:
            self.assertTrue(a1 in (0,1,2))
            self.assertTrue(a2 in (0,1,2))
            
        ## Exclude
        dss2 = ds.get_antenna_subset(exclude=(4,), indicies=True)
        self.assertEqual(len(dss2.baselines), 6)
        for a1,a2 in dss2.baselines:
            self.assertTrue(a1 != 4)
            self.assertTrue(a2 != 4)
            
        # Stand numbers
        ## Include
        dss1 = ds.get_antenna_subset(include=(151,222,150), indicies=False)
        self.assertEqual(len(dss1.baselines), 3)
        for a1,a2 in dss1.baselines:
            self.assertTrue(a1 in [idi.stands.index(s) for s in (151,222,150)])
            self.assertTrue(a2 in [idi.stands.index(s) for s in (151,222,150)])
            
        ## Exclude
        dss2 = ds.get_antenna_subset(exclude=(173,), indicies=False)
        self.assertEqual(len(dss2.baselines), 6)
        for a1,a2 in dss2.baselines:
            self.assertTrue(a1 != idi.stands.index(173))
            self.assertTrue(a2 != idi.stands.index(173))
        
    def test_prune_alt(self):
        """Test the utils.get_uv_range function - alternate FITS IDI file."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiAltFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        
        # Prune
        dsp1 = ds.get_uv_range(min_uv=10)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p1 = getattr(dsp1, pol)
            self.assertTrue(p1.nbaseline < p0.nbaseline)
            
        # Auto-prune
        dsp2 = idi.get_data_set(1, min_uv=10)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p2 = getattr(dsp2, pol)
            self.assertTrue(p2.nbaseline < p0.nbaseline)
            
        # Auto-prune that should result in no baselines
        dsp3 = idi.get_data_set(1, min_uv=100)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p3 = getattr(dsp3, pol)
            self.assertEqual(p3.nbaseline, 0)
            
    def test_prune_uvfits(self):
        """Test the utils.get_uv_range function - UVFITS file."""
        
        # Open the FITS IDI file
        uv = utils.CorrelatedData(uvFile)
        
        # Get some data to sort
        ds = uv.get_data_set(1)
        
        # Prune
        dsp1 = ds.get_uv_range(min_uv=10)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p1 = getattr(dsp1, pol)
            self.assertTrue(p1.nbaseline < p0.nbaseline)
            
        # Auto-prune
        dsp2 = uv.get_data_set(1, min_uv=10)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p2 = getattr(dsp2, pol)
            self.assertTrue(p2.nbaseline < p0.nbaseline)
            
        # Auto-prune that should result in no baselines
        dsp3 = uv.get_data_set(1, min_uv=100)
        for pol in ds.pols:
            p0 = getattr(ds,   pol)
            p3 = getattr(dsp3, pol)
            self.assertEqual(p3.nbaseline, 0)
            
    def test_rephase(self):
        """Test the utils.rephase_data function."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiFile)
        
        # Get the AntennaArray instance
        aa = idi.get_antennaarray()
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        orig_bls = ds.baselines
        orig_dat = ds.XX.data.copy()
        orig_pc  = ds.phase_center
        
        # Rephase #1
        ds.rephase(new_phase_center=simSrcs['Sun'])
        for i in xrange(ds.nbaseline):
            self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
            self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
            
        # Rephase #2
        ds.rephase(new_phase_center=orig_pc)
        for i in xrange(ds.nbaseline):
            self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
            self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
            
            for j in xrange(ds.nchan):
                self.assertAlmostEqual(orig_dat[i][j], ds.XX.data[i][j], 2)
                
        # Bad rephase
        self.assertRaises(RuntimeError, ds.rephase, simSrcs['vir'])
        
    def test_rephase_alt(self):
        """Test the utils.rephase_data function - alternate FITS IDI file."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiAltFile)
        
        # Get the AntennaArray instance
        aa = idi.get_antennaarray()
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        orig_bls = ds.baselines
        orig_dat = ds.XX.data.copy()
        orig_pc  = ds.phase_center
        
        # Rephase #1
        ds.rephase(new_phase_center=simSrcs['Sun'])
        for i in xrange(ds.nbaseline):
            self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
            self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
            
        # Rephase #2
        ds.rephase(new_phase_center=orig_pc)
        for i in xrange(ds.nbaseline):
            self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
            self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
            
            for j in xrange(ds.nchan):
                self.assertAlmostEqual(orig_dat[i][j], ds.XX.data[i][j], 2)
                
        # Bad rephase
        self.assertRaises(RuntimeError, ds.rephase, simSrcs['vir'])
        
    def test_rephase_uvfits(self):
        """Test the utils.rephase_data function - UVFITS file."""
        
        # Open the UVFITS file
        uv = utils.CorrelatedData(uvFile)
        
        # Get the AntennaArray instance
        aa = uv.get_antennaarray()
        
        # Get some data to sort
        ds = uv.get_data_set(1)
        orig_bls = ds.baselines
        orig_dat = ds.XX.data.copy()
        orig_pc  = ds.phase_center
        
        # Rephase #1
        ds.rephase(new_phase_center=simSrcs['Sun'])
        for i in xrange(ds.nbaseline):
            self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
            self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
            
        # Rephase #2
        ds.rephase(new_phase_center=orig_pc)
        for i in xrange(ds.nbaseline):
            self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
            self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
            
            for j in xrange(ds.nchan):
                self.assertAlmostEqual(orig_dat[i][j], ds.XX.data[i][j], 2)
                
        # Bad rephase
        self.assertRaises(RuntimeError, ds.rephase, simSrcs['vir'])
        
    def test_gridding(self):
        """Test building a image from a visibility data set."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiFile)
        
        # Build the image
        ds = idi.get_data_set(1)
        junk = utils.build_gridded_image(ds, verbose=False)

        # Error checking
        self.assertRaises(RuntimeError, utils.build_gridded_image, ds, pol='XY')
        
        #
        # VisibilityData test
        #
        
        ds2 = VisibilityData()
        ds2.append( ds )
        junk = utils.build_gridded_image(ds, verbose=False)
        
    def test_gridding_alt(self):
        """Test building a image from a visibility data set - alternate FITS IDI file."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiAltFile)
        
        # Build the image
        ds = idi.get_data_set(1)
        junk = utils.build_gridded_image(ds, verbose=False)

        # Error checking
        self.assertRaises(RuntimeError, utils.build_gridded_image, ds, pol='XY')
        
    def test_gridding_uvfits(self):
        """Test building a image from a visibility data set - UVFITS file."""
        
        # Open the UVFITS file
        uv = utils.CorrelatedData(uvFile)
        
        # Build the image
        
        ds = uv.get_data_set(1)
        junk = utils.build_gridded_image(ds, verbose=False)
        
        # Error checking
        self.assertRaises(RuntimeError, utils.build_gridded_image, ds, pol='XY')
        
    def test_selfcal(self):
        """Test running a simple self calibration."""
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(idiFile)
        
        # Go for it!
        aa = idi.get_antennaarray()
        ds = idi.get_data_set(1)
        junk = selfCal.phase_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False)
        
        # Error checking
        self.assertRaises(RuntimeError, selfCal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=0  )
        self.assertRaises(RuntimeError, selfCal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=564)
        
    def test_selfcal_alt(self):
        """Test running a simple self calibration - alternate FITS IDI file."""
        
        # Open the alternate FITS IDI file
        idi = utils.CorrelatedData(idiAltFile)
        
        # Go for it!
        aa = idi.get_antennaarray()
        ds = idi.get_data_set(1)
        junk = selfCal.phase_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False)
        
        # Error checking
        self.assertRaises(RuntimeError, selfCal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=0  )
        self.assertRaises(RuntimeError, selfCal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=564)
        
    def test_selfcal_uvfits(self):
        """Test running a simple self calibration - UVFITS file."""
        
        # Open the alternate UVFITS file
        uv = utils.CorrelatedData(uvFile)
        
        # Go for it!
        aa = uv.get_antennaarray()
        ds = uv.get_data_set(1)
        junk = selfCal.phase_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False)
        
        # Error checking
        self.assertRaises(RuntimeError, selfCal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=0  )
        self.assertRaises(RuntimeError, selfCal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=564)
        
    def tearDown(self):
        """Remove the test path directory and its contents"""

        shutil.rmtree(self.testPath, ignore_errors=True)


class imaging_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.imaging units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(imaging_tests)) 


if __name__ == '__main__':
    unittest.main()
