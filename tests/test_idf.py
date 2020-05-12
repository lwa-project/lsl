"""
Unit test for the lsl.common.idf module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import re
import copy
import pytz
import ephem
import tempfile
import unittest
from math import pi
from datetime import datetime, timedelta
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from astropy.coordinates import Angle as AstroAngle

from lsl.astro import MJD_OFFSET, DJD_OFFSET
from lsl.common.paths import DATA_BUILD
from lsl.common import idf
from lsl.common.stations import lwa1, lwasv


__version__  = "0.1"
__author__    = "Jayce Dowell"


drxFile = os.path.join(DATA_BUILD, 'tests', 'drx-idf.txt')
altFile = os.path.join(DATA_BUILD, 'tests', 'alt-idf.txt')
solFile = os.path.join(DATA_BUILD, 'tests', 'sol-idf.txt')
jovFile = os.path.join(DATA_BUILD, 'tests', 'jov-idf.txt')
sdfFile = os.path.join(DATA_BUILD, 'tests', 'drx-sdf.txt')


class _SilentVerbose(object):
    def __init__(self, stdout=True, stderr=False):
        self.stdout = stdout
        self.stderr = stderr
        
    def __enter__(self):
        if self.stdout:
            sys.stdout = StringIO()
        if self.stderr:
            sys.stderr = StringIO()
        return self
        
    def __exit__(self, exc_type, exc_value, exc_tb):
        if self.stdout:
            sys.stdout = sys.__stdout__
        if self.stderr:
            sys.stderr = sys.__stderr__


class idf_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.idf
    module."""
    
    testPath = None

    def setUp(self):
        """Create the temporary file directory."""

        self.testPath = tempfile.mkdtemp(prefix='test-idf-', suffix='.tmp')
        
    ### General ###
    
    def test_flat_projects(self):
        """Test single session/scans IDFs."""
        
        obs = idf.Observer('Test Observer', 99)
        targ = idf.DRX('Target', 'Target', '2019/1/1 00:00:00', '00:00:10', 0.0, 90.0, 40e6, 50e6, 7)
        sess = idf.Run('Test Session', 1, scans=targ)
        proj = idf.Project(obs, 'Test Project', 'COMTST', runs=sess)
        out = proj.render()
        
    def test_ucf_username(self):
        """Test setting the UCF username for auto-copy support."""
        
        obs = idf.Observer('Test Observer', 99)
        targ = idf.DRX('Target', 'Target', '2019/1/1 00:00:00', '00:00:10', 0.0, 90.0, 40e6, 50e6, 7)
        sess = idf.Run('Test Session', 1, scans=targ)
        sess.data_return_method = 'UCF'
        sess.ucf_username = 'test'
        proj = idf.Project(obs, 'Test Project', 'COMTST', runs=sess)
        out = proj.render()
        self.assertTrue(out.find('ucfuser:test') >= 0)
        
        obs = idf.Observer('Test Observer', 99)
        targ = idf.DRX('Target', 'Target', '2019/1/1 00:00:00', '00:00:10', 0.0, 90.0, 40e6, 50e6, 7)
        sess = idf.Run('Test Session', 1, scans=targ, comments='This is a comment')
        sess.data_return_method = 'UCF'
        sess.ucf_username = 'test/dir1'
        proj = idf.Project(obs, 'Test Project', 'COMTST', runs=sess)
        out = proj.render()
        self.assertTrue(out.find('ucfuser:test/dir1') >= 0)
        
    ### DRX - TRK_RADEC ###
    
    def test_drx_parse(self):
        """Test reading in a TRK_RADEC IDF file."""
        
        project = idf.parse_idf(drxFile)
        
        # Basic file structure
        self.assertEqual(len(project.runs), 1)
        self.assertEqual(len(project.runs[0].scans), 1)
        
        # Correlator setup
        self.assertEqual(project.runs[0].correlator_channels, 256)
        self.assertAlmostEqual(project.runs[0].correlator_inttime, 1.0, 6)
        self.assertEqual(project.runs[0].correlator_basis, 'linear')
        
        # Observational setup - 1
        self.assertEqual(project.runs[0].scans[0].mode, 'TRK_RADEC')
        self.assertEqual(project.runs[0].scans[0].mjd,  58490)
        self.assertEqual(project.runs[0].scans[0].mpm,  74580000)
        self.assertEqual(project.runs[0].scans[0].dur,  7200000)
        self.assertEqual(project.runs[0].scans[0].freq1, 766958446)
        self.assertEqual(project.runs[0].scans[0].freq2, 1643482384)
        self.assertEqual(project.runs[0].scans[0].filter,   6)
        self.assertAlmostEqual(project.runs[0].scans[0].ra, 19.991210200, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].dec, 40.733916000, 6)
        
        # Ordering
        scan = copy.deepcopy(project.runs[0].scans[0])
        scan.start = "UTC 2019/1/7 22:43:05"
        project.runs[0].append(scan)
        
        self.assertTrue(project.runs[0].scans[0] < project.runs[0].scans[1])
        self.assertFalse(project.runs[0].scans[0] > project.runs[0].scans[1])
        self.assertTrue(project.runs[0].scans[0] != project.runs[0].scans[1])
        
    def test_drx_update(self):
        """Test updating TRK_RADEC values."""
        
        project = idf.parse_idf(drxFile)
        project.runs[0].scans[0].start = "MST 2011 Feb 23 17:00:15"
        project.runs[0].scans[0].duration = timedelta(seconds=15)
        project.runs[0].scans[0].frequency1 = 75e6
        project.runs[0].scans[0].frequency2 = 76e6
        project.runs[0].scans[0].ra = AstroAngle('5:30:00', unit='hourangle')
        project.runs[0].scans[0].dec = ephem.degrees('+22:30:00')
        
        self.assertEqual(project.runs[0].scans[0].mjd,  55616)
        self.assertEqual(project.runs[0].scans[0].mpm,  15000)
        self.assertEqual(project.runs[0].scans[0].dur,  15000)
        self.assertEqual(project.runs[0].scans[0].freq1, 1643482384)
        self.assertEqual(project.runs[0].scans[0].freq2, 1665395482)
        self.assertAlmostEqual(project.runs[0].scans[0].ra, 5.5, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].dec, 22.5, 6)
        
        dt0, dt1 = idf.get_scan_start_stop(project.runs[0].scans[0])
        self.assertEqual(dt0.year, 2011)
        self.assertEqual(dt0.month, 2)
        self.assertEqual(dt0.day, 24)
        self.assertEqual(dt0.hour, 0)
        self.assertEqual(dt0.minute, 0)
        self.assertEqual(dt0.second, 15)
        self.assertEqual(dt0.microsecond, 0)
        self.assertEqual(dt1.year, 2011)
        self.assertEqual(dt1.month, 2)
        self.assertEqual(dt1.day, 24)
        self.assertEqual(dt1.hour, 0)
        self.assertEqual(dt1.minute, 0)
        self.assertEqual(dt1.second, 30)
        self.assertEqual(dt1.microsecond, 0)
        
        for scan in project.runs[0].scans:
            scan.mjd += 1
            scan.mpm += 1000
        self.assertEqual(project.runs[0].scans[0].mjd,  55617)
        self.assertEqual(project.runs[0].scans[0].mpm,  16000)
        self.assertEqual(project.runs[0].scans[0].start, 'UTC 2011/02/25 00:00:16.000000')
        
        project.runs[0].scans[0].duration = 16.0
        self.assertEqual(project.runs[0].scans[0].dur,  16000)
        
        project.runs[0].scans[0].duration = '16.1'
        self.assertEqual(project.runs[0].scans[0].dur,  16100)
        
        project.runs[0].scans[0].duration = '0:01:01.501'
        self.assertEqual(project.runs[0].scans[0].dur,  61501)
        
    def test_drx_write(self):
        """Test writing a TRK_RADEC IDF file."""
        
        project = idf.parse_idf(drxFile)
        with _SilentVerbose() as sv:
            out = project.render(verbose=True)
        
    def test_drx_proper_motion(self):
        """Test proper motion handling in a TRK_RADEC IDF file."""
        
        project = idf.parse_idf(drxFile)
        project.runs[0].scans[0].pm = [3182.7, 592.1]
        
        self.assertAlmostEqual(project.runs[0].scans[0].pm[0], 3182.7, 1)
        self.assertAlmostEqual(project.runs[0].scans[0].pm[1], 592.1, 1)
        
        ## TODO: Coordinate test?
        sdfs = project.generate_sdfs()
        for sdf in sdfs:
            for o in range(len(project.runs[0].scans)):
                bdy = project.runs[0].scans[o].fixed_body
                bdy.compute(project.runs[0].scans[o].mjd + MJD_OFFSET - DJD_OFFSET + project.runs[0].scans[o].mjd/1000.0/86400.0)
                self.assertAlmostEqual(bdy.a_ra, sdf.sessions[0].observations[o].ra*pi/12.0, 5)
                self.assertAlmostEqual(bdy.a_dec, sdf.sessions[0].observations[o].dec*pi/180.0, 5)
                
    def test_drx_errors(self):
        """Test various TRK_RADEC IDF errors."""
        
        project = idf.parse_idf(drxFile)
        
        # Bad project
        old_id = project.id
        project.id = 'ThisIsReallyLong'
        self.assertFalse(project.validate())
        
        # Bad session
        project.id = old_id
        old_id = project.runs[0].id
        project.runs[0].id = 10001
        self.assertFalse(project.validate())
        
        # Bad interferometer
        project.runs[0].id = old_id
        with self.assertRaises(ValueError):
            project.runs[0].stations = [lwa1,]
            
        # Bad correlator channel count
        project.runs[0].stations = [lwa1,lwasv]
        with self.assertRaises(ValueError):
            project.runs[0].correlator_channels = 129
            
        # Bad correlator integration time
        project.runs[0].correlator_channels = 128
        with self.assertRaises(ValueError):
            project.runs[0].correlator_inttime = 1e-6
            
        # Bad correlator output polarization basis
        project.runs[0].correlator_inttime = 1.0
        with self.assertRaises(ValueError):
            project.runs[0].correlator_basis = 'cats'
            
        # Bad intent
        project.runs[0].correlator_basis = 'linear'
        with self.assertRaises(ValueError):
            project.runs[0].scans[0].intent = 'cats'
            
        # Good filter
        project.runs[0].scans[0].intent = 'Target'
        project.runs[0].scans[0].filter = 7
        project.runs[0].scans[0].update()
        self.assertTrue(project.validate())
        
        # Bad filter
        project.runs[0].scans[0].filter = 8
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.runs[0].scans[0].filter = 6
        project.runs[0].scans[0].frequency1 = 10.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        project.runs[0].scans[0].filter = 6
        project.runs[0].scans[0].frequency1 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        project.runs[0].scans[0].frequency1 = 38.0e6
        project.runs[0].scans[0].frequency2 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.runs[0].scans[0].frequency2 = 38.0e6
        project.runs[0].scans[0].duration = '96:00:00.000'
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad pointing
        project.runs[0].scans[0].duration = '00:00:01.000'
        project.runs[0].scans[0].dec = -72.0
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - TRK_RADEC + Alternate phase centers###
    
    def test_drx_alt_parse(self):
        """Test reading in a TRK_RADEC IDF file with other phase centers."""
        
        project = idf.parse_idf(altFile)
        
        # Basic file structure
        self.assertEqual(len(project.runs), 1)
        self.assertEqual(len(project.runs[0].scans), 1)
        
        # Correlator setup
        self.assertEqual(project.runs[0].correlator_channels, 256)
        self.assertAlmostEqual(project.runs[0].correlator_inttime, 1.0, 6)
        self.assertEqual(project.runs[0].correlator_basis, 'linear')
        
        # Observational setup - 1
        self.assertEqual(project.runs[0].scans[0].mode, 'TRK_RADEC')
        self.assertEqual(project.runs[0].scans[0].mjd,  58490)
        self.assertEqual(project.runs[0].scans[0].mpm,  74580000)
        self.assertEqual(project.runs[0].scans[0].dur,  7200000)
        self.assertEqual(project.runs[0].scans[0].freq1, 766958446)
        self.assertEqual(project.runs[0].scans[0].freq2, 1643482384)
        self.assertEqual(project.runs[0].scans[0].filter,   6)
        self.assertAlmostEqual(project.runs[0].scans[0].ra, 19.991210200, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].dec, 40.733916000, 6)
        
        # Phase center - 1
        self.assertEqual(len(project.runs[0].scans[0].alt_phase_centers), 2)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[0].ra, 19.991310200, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[0].dec, 40.733916000, 6)
        
        # Phase center - 2
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[1].ra, 19.991210200, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[1].dec, 40.734016000, 6)
        
    def test_drx_alt_update(self):
        """Test updating TRK_RADEC values with other phase centers."""
        
        project = idf.parse_idf(altFile)
        project.runs[0].scans[0].start = "MST 2011 Feb 23 17:00:15"
        project.runs[0].scans[0].duration = timedelta(seconds=15)
        project.runs[0].scans[0].frequency1 = 75e6
        project.runs[0].scans[0].frequency2 = 76e6
        project.runs[0].scans[0].ra = ephem.hours('5:30:00')
        project.runs[0].scans[0].dec = ephem.degrees('+22:30:00')
        project.runs[0].scans[0].alt_phase_centers[0].ra = ephem.hours('5:35:00')
        project.runs[0].scans[0].alt_phase_centers[1].ra = ephem.hours('5:25:00')
        
        self.assertEqual(project.runs[0].scans[0].mjd,  55616)
        self.assertEqual(project.runs[0].scans[0].mpm,  15000)
        self.assertEqual(project.runs[0].scans[0].dur,  15000)
        self.assertEqual(project.runs[0].scans[0].freq1, 1643482384)
        self.assertEqual(project.runs[0].scans[0].freq2, 1665395482)
        self.assertAlmostEqual(project.runs[0].scans[0].ra, 5.5, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].dec, 22.5, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[0].ra, 5.583333, 6)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[1].ra, 5.416667, 6)
        
    def test_drx_alt_write(self):
        """Test writing a TRK_RADEC IDF file with other phase centers."""
        
        project = idf.parse_idf(altFile)
        out = project.render()
        
    def test_drx_alt_proper_motion(self):
        """Test proper motion handling in a TRK_RADEC IDF file with other phase centers."""
        
        project = idf.parse_idf(altFile)
        project.runs[0].scans[0].alt_phase_centers[0].pm = [3182.7, 592.1]
        
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[0].pm[0], 3182.7, 1)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[0].pm[1], 592.1, 1)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[1].pm[0], 0.0, 1)
        self.assertAlmostEqual(project.runs[0].scans[0].alt_phase_centers[1].pm[1], 0.0, 1)
        
        ## TODO: Coordinate test?
        sdfs = project.generate_sdfs()
        for sdf in sdfs:
            for o in range(len(project.runs[0].scans)):
                sdf_phase_centers = sdf.project_office.observations[0][o]
                for i,phase_center in enumerate(project.runs[0].scans[o].alt_phase_centers):
                    bdy = phase_center.fixed_body
                    bdy.compute(project.runs[0].scans[o].mjd + MJD_OFFSET - DJD_OFFSET + project.runs[0].scans[o].mjd/1000.0/86400.0)
                    
                    ra = re.search("altra%i:(?P<ra>\d+(.\d*)?)" % (i+1,), sdf_phase_centers)
                    ra = float(ra.group('ra'))
                    dec = re.search("altdec%i:(?P<dec>[-+]?\d+(.\d*)?)" % (i+1,), sdf_phase_centers)
                    dec = float(dec.group('dec'))
                    self.assertAlmostEqual(bdy.a_ra, ra*pi/12.0, 5)
                    self.assertAlmostEqual(bdy.a_dec, dec*pi/180.0, 5)
                    
    def test_drx_alt_errors(self):
        """Test various TRK_RADEC IDF errors with other phase centers."""
        
        project = idf.parse_idf(altFile)
        
        # Bad interferometer
        with self.assertRaises(ValueError):
            project.runs[0].stations = [lwa1,]
            
        # Bad correlator channel count
        with self.assertRaises(ValueError):
            project.runs[0].correlator_channels = 129
            
        # Bad correlator integration time
        with self.assertRaises(ValueError):
            project.runs[0].correlator_inttime = 1e-6
            
        # Bad correlator output polarization basis
        with self.assertRaises(ValueError):
           project.runs[0].correlator_basis = 'cats'
            
        # Bad intent
        with self.assertRaises(ValueError):
            project.runs[0].scans[0].intent = 'cats'
            
        # Bad filter
        project.runs[0].scans[0].intent = 'Target'
        project.runs[0].scans[0].filter = 8
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.runs[0].scans[0].filter = 6
        project.runs[0].scans[0].frequency1 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        project.runs[0].scans[0].frequency1 = 38.0e6
        project.runs[0].scans[0].frequency2 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.runs[0].scans[0].frequency2 = 38.0e6
        project.runs[0].scans[0].duration = '96:00:00.000'
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad pointing
        project.runs[0].scans[0].duration = '00:00:01.000'
        project.runs[0].scans[0].dec = -72.0
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad alternate phase center
        project.runs[0].scans[0].dec = 40.733916000
        project.runs[0].scans[0].alt_phase_centers[0].dec = 45.0
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad alternate phase center intent
        project.runs[0].scans[0].alt_phase_centers[0].dec = 40.733916000
        with self.assertRaises(ValueError):
            project.runs[0].scans[0].alt_phase_centers[0].intent = 'cats'
            
        # Too many phase centers
        project.runs[0].scans[0].alt_phase_centers[0].intent = 'PhaseCal'
        for i in range(40):
            project.runs[0].scans[0].add_alt_phase_center('test', 'Target', 19.991210200, 40.733916000)
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - TRK_SOL ###
    
    def test_sol_parse(self):
        """Test reading in a TRK_SOL IDF file."""
        
        project = idf.parse_idf(solFile)
        
        # Basic file structure
        self.assertEqual(len(project.runs), 1)
        self.assertEqual(len(project.runs[0].scans), 1)
        
        # Correlator setup
        self.assertEqual(project.runs[0].correlator_channels, 512)
        self.assertAlmostEqual(project.runs[0].correlator_inttime, 0.5, 6)
        self.assertEqual(project.runs[0].correlator_basis, 'stokes')
        
        # Observational setup - 1
        self.assertEqual(project.runs[0].scans[0].mode, 'TRK_SOL')
        self.assertEqual(project.runs[0].scans[0].mjd,  58490)
        self.assertEqual(project.runs[0].scans[0].mpm,  74580000)
        self.assertEqual(project.runs[0].scans[0].dur,  7200000)
        self.assertEqual(project.runs[0].scans[0].freq1, 766958446)
        self.assertEqual(project.runs[0].scans[0].freq2, 1643482384)
        self.assertEqual(project.runs[0].scans[0].filter,   6)
        
    def test_sol_update(self):
        """Test updating TRK_SOL values."""
        
        project = idf.parse_idf(solFile)
        project.runs[0].scans[0].start = "MST 2011 Feb 23 17:00:15"
        project.runs[0].scans[0].duration = timedelta(seconds=15)
        project.runs[0].scans[0].frequency1 = 75e6
        project.runs[0].scans[0].frequency2 = 76e6
        
        self.assertEqual(project.runs[0].scans[0].mjd,  55616)
        self.assertEqual(project.runs[0].scans[0].mpm,  15000)
        self.assertEqual(project.runs[0].scans[0].dur,  15000)
        self.assertEqual(project.runs[0].scans[0].freq1, 1643482384)
        self.assertEqual(project.runs[0].scans[0].freq2, 1665395482)
        
    def test_sol_write(self):
        """Test writing a TRK_SOL IDF file."""
        
        project = idf.parse_idf(solFile)
        out = project.render()
        
    def test_sol_errors(self):
        """Test various TRK_SOL IDF errors."""
        
        project = idf.parse_idf(solFile)
        
        # Bad interferometer
        with self.assertRaises(ValueError):
            project.runs[0].stations = [lwa1,]
            
        # Bad correlator channel count
        with self.assertRaises(ValueError):
            project.runs[0].correlator_channels = 129
            
        # Bad correlator integration time
        with self.assertRaises(ValueError):
            project.runs[0].correlator_inttime = 1e-6
            
        # Bad correlator output polarization basis
        with self.assertRaises(ValueError):
            project.runs[0].correlator_basis = 'cats'
            
        # Bad intent
        with self.assertRaises(ValueError):
            project.runs[0].scans[0].intent = 'cats'
            
        # Bad filter
        project.runs[0].scans[0].intent = 'Target'
        project.runs[0].scans[0].filter = 8
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.runs[0].scans[0].filter = 6
        project.runs[0].scans[0].frequency1 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        project.runs[0].scans[0].frequency1 = 38.0e6
        project.runs[0].scans[0].frequency2 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.runs[0].scans[0].frequency2 = 38.0e6
        project.runs[0].scans[0].duration = '96:00:00.000'
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - TRK_JOV ###
    
    def test_jov_parse(self):
        """Test reading in a TRK_JOV IDF file."""
        
        project = idf.parse_idf(jovFile)
        
        # Basic file structure
        self.assertEqual(len(project.runs), 1)
        self.assertEqual(len(project.runs[0].scans), 1)
        
        # Correlator setup
        self.assertEqual(project.runs[0].correlator_channels, 512)
        self.assertAlmostEqual(project.runs[0].correlator_inttime, 0.1, 6)
        self.assertEqual(project.runs[0].correlator_basis, 'circular')
        
        # Observational setup - 1
        self.assertEqual(project.runs[0].scans[0].mode, 'TRK_JOV')
        self.assertEqual(project.runs[0].scans[0].mjd,  58490)
        self.assertEqual(project.runs[0].scans[0].mpm,  74580000)
        self.assertEqual(project.runs[0].scans[0].dur,  3600000)
        self.assertEqual(project.runs[0].scans[0].freq1, 766958446)
        self.assertEqual(project.runs[0].scans[0].freq2, 1643482384)
        self.assertEqual(project.runs[0].scans[0].filter,   6)
        
    def test_jov_update(self):
        """Test updating TRK_JOV values."""
        
        project = idf.parse_idf(jovFile)
        project.runs[0].scans[0].start = "MST 2011 Feb 23 17:00:15"
        project.runs[0].scans[0].duration = timedelta(seconds=15)
        project.runs[0].scans[0].frequency1 = 75e6
        project.runs[0].scans[0].frequency2 = 76e6
        
        self.assertEqual(project.runs[0].scans[0].mjd,  55616)
        self.assertEqual(project.runs[0].scans[0].mpm,  15000)
        self.assertEqual(project.runs[0].scans[0].dur,  15000)
        self.assertEqual(project.runs[0].scans[0].freq1, 1643482384)
        self.assertEqual(project.runs[0].scans[0].freq2, 1665395482)
        
    def test_jov_write(self):
        """Test writing a TRK_JOV IDF file."""
        
        project = idf.parse_idf(jovFile)
        out = project.render()
        
    def test_jov_errors(self):
        """Test various TRK_JOV IDF errors."""
        
        project = idf.parse_idf(jovFile)
        
        # Bad interferometer
        with self.assertRaises(ValueError):
            project.runs[0].stations = [lwa1,]
            
        # Bad correlator channel count
        with self.assertRaises(ValueError):
            project.runs[0].correlator_channels = 129
            
        # Bad correlator integration time
        with self.assertRaises(ValueError):
            project.runs[0].correlator_inttime = 1e-6
            
        # Bad correlator output polarization basis
        with self.assertRaises(ValueError):
            project.runs[0].correlator_basis = 'cats'
            
        # Bad intent
        with self.assertRaises(ValueError):
            project.runs[0].scans[0].intent = 'cats'
            
        # Bad filter
        project.runs[0].scans[0].intent = 'Target'
        project.runs[0].scans[0].filter = 8
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.runs[0].scans[0].filter = 6
        project.runs[0].scans[0].frequency1 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        project.runs[0].scans[0].frequency1 = 38.0e6
        project.runs[0].scans[0].frequency2 = 90.0e6
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.runs[0].scans[0].frequency2 = 38.0e6
        project.runs[0].scans[0].duration = '96:00:00.000'
        project.runs[0].scans[0].update()
        self.assertFalse(project.validate())
        
    ### Misc. ###
    
    def test_generate_sdfs(self):
        """Test generated IDFs from the IDF."""
        
        project = idf.parse_idf(drxFile)
        sdfs = project.generate_sdfs()
        
        project = idf.parse_idf(solFile)
        sdfs = project.generate_sdfs()
        
        project = idf.parse_idf(jovFile)
        sdfs = project.generate_sdfs()
        
    def test_auto_update(self):
        """Test project auto-update on render."""
        
        # Part 1 - frequency and duration
        project = idf.parse_idf(drxFile)
        project.runs[0].scans[0].frequency1 = 75e6
        project.runs[0].scans[0].duration = '00:01:31.000'
        
        fh = open(os.path.join(self.testPath, 'idf.txt'), 'w')
        fh.write(project.render())
        fh.close()
        
        project = idf.parse_idf(os.path.join(self.testPath, 'idf.txt'))
        self.assertEqual(project.runs[0].scans[0].freq1, 1643482384)
        self.assertEqual(project.runs[0].scans[0].dur, 91000)
        
        # Part 2 - frequency and duration (timedelta)
        project = idf.parse_idf(drxFile)
        project.runs[0].scans[0].frequency1 = 75e6
        project.runs[0].scans[0].duration = timedelta(minutes=1, seconds=31, microseconds=1000)
        
        fh = open(os.path.join(self.testPath, 'idf.txt'), 'w')
        fh.write(project.render())
        fh.close()
        
        project = idf.parse_idf(os.path.join(self.testPath, 'idf.txt'))
        self.assertEqual(project.runs[0].scans[0].freq1, 1643482384)
        self.assertEqual(project.runs[0].scans[0].dur, 91001)
        
        # Part 3 - frequency and start time
        project = idf.parse_idf(drxFile)
        project.runs[0].scans[0].frequency2 = 75e6
        project.runs[0].scans[0].start = "MST 2011 Feb 23 14:00:15"
        
        fh = open(os.path.join(self.testPath, 'idf.txt'), 'w')		
        fh.write(project.render())
        fh.close()
        
        project = idf.parse_idf(os.path.join(self.testPath, 'idf.txt'))
        self.assertEqual(project.runs[0].scans[0].freq2, 1643482384)
        self.assertEqual(project.runs[0].scans[0].mjd,  55615)
        self.assertEqual(project.runs[0].scans[0].mpm,  75615000)
        
        # Part 4 - frequency and start time (timedelta)
        project = idf.parse_idf(drxFile)
        _MST = pytz.timezone('US/Mountain')
        project.runs[0].scans[0].frequency2 = 75e6
        project.runs[0].scans[0].start = _MST.localize(datetime(2011, 2, 23, 14, 00, 30, 1000))
        
        fh = open(os.path.join(self.testPath, 'idf.txt'), 'w')		
        fh.write(project.render())
        fh.close()
        
        project = idf.parse_idf(os.path.join(self.testPath, 'idf.txt'))
        self.assertEqual(project.runs[0].scans[0].freq2, 1643482384)
        self.assertEqual(project.runs[0].scans[0].mjd,  55615)
        self.assertEqual(project.runs[0].scans[0].mpm,  75630001)
        
    def test_set_stations(self):
        """Test the set stations functionlity."""
        
        project = idf.parse_idf(drxFile)
        project.runs[0].stations = [lwasv, lwa1]
        self.assertTrue(project.validate())
        
    def test_is_valid(self):
        """Test whether or not is_valid works."""
        
        self.assertTrue(idf.is_valid(drxFile))
        self.assertTrue(idf.is_valid(solFile))
        self.assertTrue(idf.is_valid(jovFile))
        
    def test_is_not_valid(self):
        """Test whether or not is_valid works on LWA1 and IDF files."""
        
        self.assertFalse(idf.is_valid(sdfFile))
        
    def test_username(self):
        """Test setting auto-copy parameters."""
        
        project = idf.parse_idf(drxFile)
        project.runs[0].data_return_method = 'UCF'
        project.runs[0].ucf_username = 'jdowell'
        out = project.render()
        
        self.assertTrue(out.find('Requested data return method is UCF') > 0)
        self.assertTrue(out.find('ucfuser:jdowell') > 0)
        
        project.writeto(os.path.join(self.testPath, 'idf.txt'))
        
        project = idf.parse_idf(os.path.join(self.testPath, 'idf.txt'))
        out = project.render()
        
        self.assertTrue(out.find('Requested data return method is UCF') > 0)
        self.assertTrue(out.find('ucfuser:jdowell') > 0)
        
    def tearDown(self):
        """Remove the test path directory and its contents"""

        tempFiles = os.listdir(self.testPath)
        for tempFile in tempFiles:
            os.unlink(os.path.join(self.testPath, tempFile))
        os.rmdir(self.testPath)
        self.testPath = None


class idf_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.idf units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(idf_tests)) 


if __name__ == '__main__':
    unittest.main()
