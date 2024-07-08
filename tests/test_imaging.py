"""
Unit tests for the lsl.imaging modules.
"""

import os
import copy
import glob
import time
import numpy as np
import shutil
import tempfile
import unittest
import subprocess

from astropy.time import Time as AstroTime
from astropy.coordinates import SkyCoord, AltAz

from lsl import astro
from lsl.imaging import utils
from lsl.imaging import analysis
from lsl.imaging import deconv
from lsl.imaging import selfcal
from lsl.imaging import overlay
from lsl.imaging.data import VisibilityData, PolarizationDataSet
from lsl.writer.fitsidi import Idi, NUMERIC_STOKES
from lsl.sim import vis
from lsl.common.stations import lwa1, parse_ssmif
from lsl.correlator import uvutils
import lsl.testing

run_ms_tests = False
try:
    import casacore
    from lsl.writer.measurementset import Ms
    run_ms_tests = True
except ImportError:
    pass

run_plotting_tests = False
try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    run_plotting_tests = True
except ImportError:
    pass


__version__  = "0.3"
__author__    = "Jayce Dowell"


uvFile = os.path.join(os.path.dirname(__file__), 'data', 'uv-test.fits')
idiFile = os.path.join(os.path.dirname(__file__), 'data', 'idi-test.fits')
idiAltFile = os.path.join(os.path.dirname(__file__), 'data', 'idi-test-alt.fits')
idiSSMIFFile = os.path.join(os.path.dirname(__file__), 'data', 'idi-test-alt.txt')


class imaging_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.imaging
    modules."""
    
    def setUp(self):
        """Turn off all numpy warnings and create the temporary file directory."""

        np.seterr(all='ignore')
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
        
        idi.close()
        
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
        
        idi.close()
        
        # both at the same time
        with utils.CorrelatedDataIDI(idiFile) as idi:
            i = 0
            for ds in idi.data_set_sequence(include_auto=False):
                self.assertEqual(ds.XX.data.shape[0], 5*(5-1)/2)
                i += 1
            self.assertEqual(i, 1)
            
    def _init_data(self):
        """Private function to generate a random set of data for writing a UVFITS
        file.  The data is returned as a dictionary with keys:
         * freq - frequency array in Hz
         * site - lwa.common.stations object
         * stands - array of stand numbers
         * bl - list of baseline pairs in real stand numbers
         * vis - array of visibility data in baseline x freq format
        """

        # Frequency range
        freq = np.arange(0,512)*20e6/512 + 40e6
        # Site and stands
        site = lwa1
        antennas = site.antennas[0:40:2]
        
        # Set baselines and data
        blList = uvutils.get_baselines(antennas, include_auto=True, indicies=False)
        visData = np.random.rand(len(blList), len(freq))
        visData = visData.astype(np.complex64)
        
        return {'freq': freq, 'site': site, 'antennas': antennas, 'bl': blList, 'vis': visData}
        
    def test_CorrelatedDataIDI_MultiIF(self):
        """Test the utils.CorrelatedDataIDI class on a file with multiple IFs."""
        
        # Get some data
        data = self._init_data()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join('idi-test-MultiIF.fits')
        
        # Start the file
        fits = Idi(testFile, ref_time=testTime, overwrite=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_frequency(data['freq']+30e6)
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], 
                          np.concatenate([data['vis'], 10*data['vis']], axis=1))
        fits.write()
        fits.close()
        
        # Open the FITS IDI file
        idi = utils.CorrelatedData(testFile)
        self.assertEqual(idi.freq.size, 2*data['freq'].size)
        ds = idi.get_data_set(1, include_auto=True)
        
        np.testing.assert_allclose(ds.XX.data[:,:data['freq'].size], data['vis'])
        np.testing.assert_allclose(ds.XX.data[:,data['freq'].size:], 10*data['vis'])
        
        idi.close()
        
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
        
        idi.close()
        
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
        ants1 = idi1.station.antennas
        ants2 = station2.antennas
        for a1,a2 in zip(ants1, ants2):
            self.assertEqual(a1.id, a2.id)
            self.assertEqual(a1.stand.id, a2.stand.id)
            self.assertAlmostEqual(a1.stand.x, a2.stand.x, 2)
            self.assertAlmostEqual(a1.stand.y, a2.stand.y, 2)
            self.assertAlmostEqual(a1.stand.z, a2.stand.z, 2)
            
        idi1.close()
        idi2.close()
        
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
        
        uv.close()
        
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS(self):
        """Test the utils.CorrelatedDataMS class."""
        
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-W.ms')
        
        # Get some data
        data = self._init_data()
        
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
        
        tbl.close()
        
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS_SingleIF(self):
        """Test the utils.CorrelatedDataMS class on a file with a single IF."""
        
        # Get some data
        data = self._init_data()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-SingleIF.ms')
        
        # Start the file
        fits = Ms(testFile, ref_time=testTime, overwrite=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], data['vis'])
        fits.write()
        fits.close()
        
        # Open the measurement set
        ms = utils.CorrelatedDataMS(testFile)
        self.assertEqual(ms.freq.size, data['freq'].size)
        ds = ms.get_data_set(1, include_auto=True)
        
        np.testing.assert_allclose(ds.XX.data, data['vis'])
        
        ms.close()
        
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS_MultiIF(self):
        """Test the utils.CorrelatedDataMS class on a file with multiple IFs."""
        
        # Get some data
        data = self._init_data()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-MultiIF.ms')
        
        # Start the file
        fits = Ms(testFile, ref_time=testTime, overwrite=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_frequency(data['freq']+30e6)
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], 
                          np.concatenate([data['vis'], 10*data['vis']], axis=1))
        fits.write()
        fits.close()
        
        # Open the measurement set
        ms = utils.CorrelatedDataMS(testFile)
        self.assertEqual(ms.freq.size, 2*data['freq'].size)
        ds = ms.get_data_set(1, include_auto=True)
        
        np.testing.assert_allclose(ds.XX.data[:,:data['freq'].size], data['vis'])
        np.testing.assert_allclose(ds.XX.data[:,data['freq'].size:], 10*data['vis'])
        
        ms.close()
        
    @unittest.skipUnless(run_ms_tests, "requires the 'casacore' module")
    def test_CorrelatedDataMS_compressed(self):
        """Test the utils.CorrelatedDataMS class on a compressed file."""
        
        # Get some data
        data = self._init_data()
        
        # Filename and time
        testTime, testFile = time.time(), os.path.join(self.testPath, 'ms-test-MultiIF.ms')
        
        # Start the file
        fits = Ms(testFile, ref_time=testTime, overwrite=True)
        fits.set_stokes(['xx'])
        fits.set_frequency(data['freq'])
        fits.set_frequency(data['freq']+30e6)
        fits.set_geometry(data['site'], data['antennas'])
        fits.add_data_set(astro.utcjd_to_taimjd(astro.unix_to_utcjd(testTime)), 6.0, data['bl'], 
                          np.concatenate([data['vis'], 10*data['vis']], axis=1))
        fits.write()
        fits.close()
        
        # Compress
        compressedFile = os.path.splitext(testFile)[0]+'.tar.gz'
        cmd = ['tar', 'czf', compressedFile, '-C', self.testPath, os.path.basename(testFile)]
        subprocess.check_call(cmd)
        
        # Open the measurement set
        ms = utils.CorrelatedDataMS(compressedFile)
        self.assertEqual(ms.freq.size, 2*data['freq'].size)
        ds = ms.get_data_set(1, include_auto=True)
        
        np.testing.assert_allclose(ds.XX.data[:,:data['freq'].size], data['vis'])
        np.testing.assert_allclose(ds.XX.data[:,data['freq'].size:], 10*data['vis'])
        
        ms.close()
        
    def test_sort(self):
        """Test the utils.sort function."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                # Open the file
                idi = utils.CorrelatedData(filename)
                
                # Get some data to sort
                ds = idi.get_data_set(1, sort=False)
                
                # Sort
                dss = ds.copy()
                dss.sort()
                for pol in ds.pols:
                    p0 = getattr(ds,  pol)
                    p1 = getattr(dss, pol)
                    self.assertEqual(p0.nbaseline, p1.nbaseline)
                    
    def test_prune(self):
        """Test the utils.get_uv_range function."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                # Open the file
                idi = utils.CorrelatedData(filename)
                
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
                    
                idi.close()
                
    def test_subset(self):
        """Test the utils.get_antenna_subset function."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                # Open the file
                idi = utils.CorrelatedData(filename)
                
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
                    
                idi.close()
                
    def test_rephase(self):
        """Test the utils.rephase_data function."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                # Open the file
                idi = utils.CorrelatedData(filename)
                
                # Get the AntennaArray instance
                aa = idi.get_antennaarray()
                
                # Get some data to sort
                ds = idi.get_data_set(1)
                orig_bls = ds.baselines
                orig_dat = ds.XX.data.copy()
                orig_pc  = ds.phase_center
                
                # Rephase #1
                ds.rephase(new_phase_center=vis.SOURCES['Sun'])
                for i in range(ds.nbaseline):
                    self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
                    self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
                    
                # Rephase #2
                ds.rephase(new_phase_center=orig_pc)
                for i in range(ds.nbaseline):
                    self.assertEqual(orig_bls[i][0], ds.baselines[i][0])
                    self.assertEqual(orig_bls[i][1], ds.baselines[i][1])
                    
                    for j in range(ds.nchan):
                        self.assertAlmostEqual(orig_dat[i][j], ds.XX.data[i][j], 2)
                        
                # Bad rephase
                self.assertRaises(RuntimeError, ds.rephase, vis.SOURCES['vir'])
                
                idi.close()
                
    def test_convert_to_stokes(self):
        """Test the utils.convert_to_stokes function."""
        
        # Open the file
        idi = utils.CorrelatedData(idiFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        new_pol = PolarizationDataSet('YY', ds.XX.data, ds.XX.weight, ds.XX.mask)
        ds.append(new_pol)
        new_pol = PolarizationDataSet('XY', ds.XX.data, ds.XX.weight, ds.XX.mask)
        ds.append(new_pol)
        new_pol = PolarizationDataSet('YX', ds.XX.data, ds.XX.weight, ds.XX.mask)
        ds.append(new_pol)
        
        # Convert
        ds2 = utils.convert_to_stokes(ds)
        
        # Check
        self.assertTrue(getattr(ds2, 'I', None) is not None)
        self.assertTrue(getattr(ds2, 'Q', None) is not None)
        self.assertTrue(getattr(ds2, 'U', None) is not None)
        self.assertTrue(getattr(ds2, 'V', None) is not None)
        
        np.testing.assert_allclose(ds2.I.data, 2*ds.XX.data)
        np.testing.assert_allclose(ds2.Q.data, 0*ds.XX.data)
        np.testing.assert_allclose(ds2.U.data, 2*ds.XX.data)
        np.testing.assert_allclose(ds2.V.data, 0*ds.XX.data)
        
        idi.close()
        
    def test_convert_to_linear(self):
        """Test the utils.convert_to_linear function."""
        
        # Open the file
        idi = utils.CorrelatedData(idiFile)
        
        # Get some data to sort
        ds = idi.get_data_set(1)
        new_pol = PolarizationDataSet('YY', ds.XX.data, ds.XX.weight, ds.XX.mask)
        ds.append(new_pol)
        new_pol = PolarizationDataSet('XY', ds.XX.data, ds.XX.weight, ds.XX.mask)
        ds.append(new_pol)
        new_pol = PolarizationDataSet('YX', ds.XX.data, ds.XX.weight, ds.XX.mask)
        ds.append(new_pol)
        
        # Convert
        ds2 = utils.convert_to_stokes(ds)
        
        # Convert back
        ds3 = utils.convert_to_linear(ds2)
        
        # Check
        self.assertTrue(getattr(ds3, 'XX', None) is not None)
        self.assertTrue(getattr(ds3, 'YY', None) is not None)
        self.assertTrue(getattr(ds3, 'XY', None) is not None)
        self.assertTrue(getattr(ds3, 'YX', None) is not None)
        
        np.testing.assert_allclose(ds3.XX.data, ds.XX.data)
        np.testing.assert_allclose(ds3.YY.data, ds.XX.data)
        np.testing.assert_allclose(ds3.XY.data, ds.XX.data)
        np.testing.assert_allclose(ds3.YX.data, ds.XX.data)
        
        idi.close()
        
    def test_gridding(self):
        """Test building a image from a visibility data set."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                # Open the file
                idi = utils.CorrelatedData(filename)
                
                # Build the image
                ds = idi.get_data_set(1)
                junk = utils.build_gridded_image(ds, verbose=False)
                
                # Different weightings
                for weighting in ('natural', 'uniform', 'briggs'):
                    junk2 = junk.image(weighting=weighting)
                    
                # Different tapers
                for taper in ((0.0, 10.0), (0.1, 10.0)):
                    junk2 = junk.image(taper=taper)
                    
                # Error checking
                self.assertRaises(RuntimeError, utils.build_gridded_image, ds, pol='XY')
                
                #
                # VisibilityData test
                #
                
                ds2 = VisibilityData()
                ds2.append( ds )
                junk = utils.build_gridded_image(ds, verbose=False)
                
                idi.close()
                
    def test_image_coordinates(self):
        """Test getting per-pixel image coordinates."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                uv = utils.CorrelatedData(filename)
                
                # Build the image
                ds = uv.get_data_set(1)
                junk = utils.build_gridded_image(ds, verbose=False)
                
                # Get image properties
                junk.field_of_view
                junk.pixel_size
                
                radec = utils.get_image_radec(junk, uv.get_antennaarray())
                azalt = utils.get_image_azalt(junk, uv.get_antennaarray())
                
                uv.close()
        
    def test_selfcal(self):
        """Test running a simple self calibration."""
        
        for filename,type in zip((idiFile, idiAltFile, uvFile), ('FITS-IDI', 'Alt. FITS-IDI', 'UVFITS')):
            with self.subTest(filetype=type):
                # Open the file
                idi = utils.CorrelatedData(idiFile)
                
                # Go for it!
                aa = idi.get_antennaarray()
                ds = idi.get_data_set(1)
                junk = selfcal.phase_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False, amplitude=True, return_convergence=True)
                junk = selfcal.phase_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False, amplitude=True)
                junk = selfcal.phase_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False)
                junk = selfcal.delay_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False, amplitude=True, return_convergence=True)
                junk = selfcal.delay_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False, amplitude=True)
                junk = selfcal.delay_only(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False)
                junk = selfcal.delay_and_phase(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False, amplitude=True, return_convergence=True)
                junk = selfcal.delay_and_phase(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False, amplitude=True)
                junk = selfcal.delay_and_phase(aa, ds, ds, 173, 'XX', max_iter=1, verbose=False)
                
                # Error checking
                self.assertRaises(RuntimeError, selfcal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=0  )
                self.assertRaises(RuntimeError, selfcal.phase_only, aa, ds, ds, 173, 'YX', ref_ant=564)
                
                idi.close()
                
    def test_background(self):
        """Test the background estimation"""
        
        img = np.random.randn(256, 256)*0.5 + 10
        bkg = analysis.estimate_background(img)
        self.assertAlmostEqual(bkg.mean(), img.mean(), 0)
        
    def test_source_detection(self):
        """Test point source detection"""
        
        img = np.random.randn(256, 256)*0.5 + 10
        sx = ( 10, 56, 105)
        sy = (115, 35, 200)
        sf = ( 20, 30,  15)
        for i,j,f in zip(sx, sy, sf):
            for di in (-2, -1, 0, 1, 2):
                for dj in (-2, -1, 0, 1, 2):
                    s = np.exp(-(di**2+dj**2)/2.0/1.0**2)
                    img[i+di,j+dj] += f*s
        img = img - analysis.estimate_background(img)
        cx, cy, pf, sh, ro = analysis.find_point_sources(img, threshold=10, verbose=False)
        for x,y,f in zip(cx,cy,pf):
            self.assertTrue(int(round(x)) in sx)
            self.assertTrue(int(round(y)) in sy)
            
    def test_clean(self):
        """Test CLEAN"""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        with lsl.testing.SilentVerbose():
            # Build an image
            img = utils.build_gridded_image(out)
            
            # CLEAN
            deconv.clean(out, img, max_iter=5, verbose=False, plot=run_plotting_tests)
            
    def test_clean_sources(self):
        """Test CLEANing around specific sources"""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        with lsl.testing.SilentVerbose():
            # Build an image
            img = utils.build_gridded_image(out)
            
            # CLEAN
            deconv.clean_sources(out, img, vis.SOURCES, max_iter=5, verbose=False, plot=run_plotting_tests)
            
    def test_clean_leastsq(self):
        """Test CLEANing using least squares in the image plane"""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        with lsl.testing.SilentVerbose():
            # Build an image
            img = utils.build_gridded_image(out)
            
            # CLEAN
            deconv.lsq(out, img, max_iter=2, verbose=False, plot=run_plotting_tests)
            
    @unittest.skipUnless(run_plotting_tests, "requires the 'matplotlib' module")
    def test_plotting(self):
        """Test drawing an image."""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        # Build an image
        img = utils.build_gridded_image(out)
        
        # Plot
        fig = plt.figure()
        ax = fig.gca()
        utils.plot_gridded_image(ax, img)
        
    def test_radec_of(self):
        """Test finding the RA/dec of a topocentric position as viewed by an observer."""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # RA/dec -> az/alt
        el = lwa1.earth_location
        ot = AstroTime(lwa1.date, astro.DJD_OFFSET, format='jd', scale='utc')
        sc = SkyCoord('12h13m45.2s', '+15d10m13.4s', frame='fk5', equinox='J2000')
        tp = sc.transform_to(AltAz(location=el, obstime=ot))
        
        # Convert back
        eq = overlay._radec_of(aa, tp.az.deg, tp.alt.deg, degrees=True)
        
        # Compare with the original
        self.assertAlmostEqual(eq[0], sc.ra.deg, 6)
        self.assertAlmostEqual(eq[1], sc.dec.deg, 6)
        
    @unittest.skipUnless(run_plotting_tests, "requires the 'matplotlib' module")
    def test_plotting_horizon(self):
        """Test drawing the horizon on an image."""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        # Build an image
        img = utils.build_gridded_image(out)
        
        # Plot
        fig = plt.figure()
        ax = fig.gca()
        utils.plot_gridded_image(ax, img)
        overlay.horizon(ax, img)
        del fig
        
    @unittest.skipUnless(run_plotting_tests, "requires the 'matplotlib' module")
    def test_plotting_sources(self):
        """Test marking sources on an image."""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        # Build an image
        img = utils.build_gridded_image(out)
        
        # Plot
        fig = plt.figure()
        ax = fig.gca()
        utils.plot_gridded_image(ax, img)
        overlay.sources(ax, img, vis.SOURCES)
        del fig
        
    @unittest.skipUnless(run_plotting_tests, "requires the 'matplotlib' module")
    def test_plotting_graticules(self):
        """Test adding a graticule to an image."""
        
        # Setup
        antennas = lwa1.antennas[0:20]
        freqs = np.arange(30e6, 50e6, 1e6)
        aa = vis.build_sim_array(lwa1, antennas, freqs)
        
        # Build the data dictionary
        out = vis.build_sim_data(aa, vis.SOURCES, jd=2458962.16965)
        
        # Build an image
        img = utils.build_gridded_image(out)
        
        # Plot
        fig = plt.figure()
        ax = fig.gca()
        utils.plot_gridded_image(ax, img)
        with self.subTest(type='RA/Dec.'):
            overlay.graticule_radec(ax, img)
        with self.subTest(type='az/alt'):
            overlay.graticule_azalt(ax, img)
        del fig
        
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
