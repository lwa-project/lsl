# -*- coding: utf-8 -*-

"""
Unit test for the lsl.common.sdf module.
"""

import os
import pytz
import ephem
import tempfile
import unittest
from datetime import datetime, timedelta
try:
    import cStringIO as StringIO
except ImprotError:
    import StringIO

from lsl.common.paths import dataBuild as dataPath
from lsl.common import sdfADP
from lsl.common.stations import lwa1, lwasv


__revision__ = "$Rev$"
__version__  = "0.1"
__author__    = "Jayce Dowell"


tbwFile = os.path.join(dataPath, 'tests', 'tbw-sdf.txt')
tbnFile = os.path.join(dataPath, 'tests', 'tbn-sdf.txt')
drxFile = os.path.join(dataPath, 'tests', 'drx-sdf.txt')
solFile = os.path.join(dataPath, 'tests', 'sol-sdf.txt')
jovFile = os.path.join(dataPath, 'tests', 'jov-sdf.txt')
stpFile = os.path.join(dataPath, 'tests', 'stp-sdf.txt')
spcFile = os.path.join(dataPath, 'tests', 'spc-sdf.txt')
tbfFile = os.path.join(dataPath, 'tests', 'tbf-sdf.txt')


class sdf_adp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.sdf
    module."""
    
    testPath = None

    def setUp(self):
        """Create the temporary file directory."""

        self.testPath = tempfile.mkdtemp(prefix='test-sdf-', suffix='.tmp')
        
    ### General ###
    
    def test_time(self):
        """Test the sdfADP.parseTime() function."""
        
        _UTC = pytz.utc
        _EST = pytz.timezone('US/Eastern')
        
        # Different realizations of the same thing
        s1 = "EST 2011-01-01 12:13:14.567"
        s2 = "EST 2011 01 01 12:13:14.567"
        s3 = "EST 2011 Jan 01 12:13:14.567"
        s4 = _EST.localize(datetime(2011, 1, 1, 12, 13, 14, 567000))
        s5 = _EST.localize(datetime(2011, 1, 1, 12, 13, 14, 567123))
        
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s2))
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s3))
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s4))
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s5))
        self.assertEqual(sdfADP.parseTime(s2), sdfADP.parseTime(s3))
        self.assertEqual(sdfADP.parseTime(s2), sdfADP.parseTime(s4))
        self.assertEqual(sdfADP.parseTime(s2), sdfADP.parseTime(s5))
        self.assertEqual(sdfADP.parseTime(s3), sdfADP.parseTime(s4))
        self.assertEqual(sdfADP.parseTime(s3), sdfADP.parseTime(s5))
        self.assertEqual(sdfADP.parseTime(s4), sdfADP.parseTime(s5))
        
        # Month name and month number agreement
        for n,m in enumerate(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            s1 = "UTC 2011-%s-14 12:13:14.000" % m
            s2 = "UTC 2011-%02i-14 12:13:14.000" % (n+1)
            s3 = _UTC.localize(datetime(2011, n+1, 14, 12, 13, 14, 0))
            self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s2))
            self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s3))
            
        # Time zone agreement - UTC
        s1 = "2011-01-01 12:13:14.567"
        s2 = "2011 01 01 12:13:14.567"
        s3 = _UTC.localize(datetime(2011, 1, 1, 12, 13, 14, 567000))
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s2))
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s3))
        self.assertEqual(sdfADP.parseTime(s2), sdfADP.parseTime(s3))
        
        # Time zone agreement - local
        for o,z in enumerate(['EST', 'CST', 'MST', 'PST']):
            h = 12
            o = -5 - o
            s1 = "%s 2011 01 01 %02i:13:14.567" % ('UTC', h)
            s2 = "%s 2011 01 01 %02i:13:14.567" % (z, h+o)
            self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s2))
            
        # Something strange
        s1 = "CET 2013-01-08 19:42:00.000"
        s2 = "2013-01-08 18:42:00.000+00:00"
        s3 = "2013-01-08 11:42:00-0700"
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s2))
        self.assertEqual(sdfADP.parseTime(s1), sdfADP.parseTime(s3))
        self.assertEqual(sdfADP.parseTime(s2), sdfADP.parseTime(s3))
        
        # Details
        s1 = "2011 01 02 03:04:05.678"
        out = sdfADP.parseTime(s1)
        ## Date
        self.assertEqual(out.year, 2011)
        self.assertEqual(out.month, 1)
        self.assertEqual(out.day, 2)
        ## Time
        self.assertEqual(out.hour, 3)
        self.assertEqual(out.minute, 4)
        self.assertEqual(out.second, 5)
        self.assertEqual(out.microsecond, 678000)
        
        # LST at LWA1
        s1 = "LST 2013-01-08 19:42:00.000"
        s2 = "UTC 2013-01-08 19:35:28.994"
        self.assertEqual(sdfADP.parseTime(s1, station=lwasv), sdfADP.parseTime(s2))
        
    ### TBW ###
    
    def test_tbw_parse(self):
        """Test reading in a TBW SDF file."""
        
        self.assertRaises(RuntimeError, sdfADP.parseSDF, tbwFile)
        
    ### TBN ###
    
    def test_tbn_parse(self):
        """Test reading in a TBN SDF file."""
        
        project = sdfADP.parseSDF(tbnFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 2)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode, 'TBN')
        self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[0].mpm,      0)
        self.assertEqual(project.sessions[0].observations[0].dur,  10000)
        self.assertEqual(project.sessions[0].observations[0].freq1, 438261968)
        self.assertEqual(project.sessions[0].observations[0].filter,   7)
        
        # Observational setup - 2
        self.assertEqual(project.sessions[0].observations[1].mode, 'TBN')
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
        self.assertEqual(project.sessions[0].observations[1].dur,  10000)
        self.assertEqual(project.sessions[0].observations[1].freq1, 832697741)
        self.assertEqual(project.sessions[0].observations[1].filter,   7)
        
    def test_tbn_update(self):
        """Test updating TBN values."""
        
        project = sdfADP.parseSDF(tbnFile)
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:00:15")
        project.sessions[0].observations[1].setDuration(timedelta(seconds=15))
        project.sessions[0].observations[1].setFrequency1(75e6)
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        self.assertEqual(project.sessions[0].observations[1].dur,  15000)
        self.assertEqual(project.sessions[0].observations[1].freq1, 1643482384)
        
    def test_tbn_write(self):
        """Test writing a TBN SDF file."""
        
        project = sdfADP.parseSDF(tbnFile)
        out = project.render()
        
    def test_tbn_errors(self):
        """Test various TBN SDF errors."""
        
        project = sdfADP.parseSDF(tbnFile)
        
        # Bad filter
        project.sessions[0].observations[0].filter = 10
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.sessions[0].observations[0].filter = 7
        project.sessions[0].observations[0].frequency1 = 95.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.sessions[0].observations[0].frequency1 = 38.0e6
        project.sessions[0].observations[0].duration = '96:00:00.000'
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - TRK_RADEC ###
    
    def test_drx_parse(self):
        """Test reading in a TRK_RADEC SDF file."""
        
        project = sdfADP.parseSDF(drxFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 2)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode, 'TRK_RADEC')
        self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[0].mpm,      0)
        self.assertEqual(project.sessions[0].observations[0].dur,  10000)
        self.assertEqual(project.sessions[0].observations[0].freq1,  438261968)
        self.assertEqual(project.sessions[0].observations[0].freq2, 1928352663)
        self.assertEqual(project.sessions[0].observations[0].filter,   7)
        self.assertAlmostEqual(project.sessions[0].observations[0].ra, 5.6, 6)
        self.assertAlmostEqual(project.sessions[0].observations[0].dec, 22.0, 6)
        
        # Observational setup - 2
        self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_RADEC')
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
        self.assertEqual(project.sessions[0].observations[1].dur,  10000)
        self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
        self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
        self.assertEqual(project.sessions[0].observations[1].filter,   7)
        self.assertAlmostEqual(project.sessions[0].observations[1].ra, 5.6, 6)
        self.assertAlmostEqual(project.sessions[0].observations[1].dec, 22.0, 6)
        
    def test_drx_update(self):
        """Test updating TRK_RADEC values."""
        
        project = sdfADP.parseSDF(drxFile)
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:00:15")
        project.sessions[0].observations[1].setDuration(timedelta(seconds=15))
        project.sessions[0].observations[1].setFrequency1(75e6)
        project.sessions[0].observations[1].setFrequency2(76e6)
        project.sessions[0].observations[1].setRA(ephem.hours('5:30:00'))
        project.sessions[0].observations[1].setDec(ephem.degrees('+22:30:00'))
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        self.assertEqual(project.sessions[0].observations[1].dur,  15000)
        self.assertEqual(project.sessions[0].observations[1].freq1, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].freq2, 1665395482)
        self.assertAlmostEqual(project.sessions[0].observations[1].ra, 5.5, 6)
        self.assertAlmostEqual(project.sessions[0].observations[1].dec, 22.5, 6)
        
    def test_drx_write(self):
        """Test writing a TRK_RADEC SDF file."""
        
        project = sdfADP.parseSDF(drxFile)
        # Fix for LWA-SV only going up to filter code 6
        for obs in project.sessions[0].observations:
            obs.setFilter(6)
        out = project.render()
        
    def test_drx_errors(self):
        """Test various TRK_RADEC SDF errors."""
        
        project = sdfADP.parseSDF(drxFile)
        
        # Bad beam
        project.sessions[0].drxBeam = 6
        self.assertFalse(project.validate())
        
        # No beam
        project.sessions[0].drxBeam = -1
        self.assertFalse(project.validate())
        
        # Bad filter
        project.sessions[0].observations[0].filter = 7
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.sessions[0].observations[0].filter = 6
        project.sessions[0].observations[0].frequency1 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        project.sessions[0].observations[0].frequency1 = 38.0e6
        project.sessions[0].observations[0].frequency2 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.sessions[0].observations[0].frequency2 = 38.0e6
        project.sessions[0].observations[0].duration = '96:00:00.000'
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        # Bad pointing
        project.sessions[0].observations[0].duration = '00:00:01.000'
        project.sessions[0].observations[0].dec = -72.0
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - TRK_SOL ###
    
    def test_sol_parse(self):
        """Test reading in a TRK_SOL SDF file."""
        
        project = sdfADP.parseSDF(solFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 2)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode, 'TRK_SOL')
        self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[0].mpm,      0)
        self.assertEqual(project.sessions[0].observations[0].dur,  10000)
        self.assertEqual(project.sessions[0].observations[0].freq1,  438261968)
        self.assertEqual(project.sessions[0].observations[0].freq2, 1928352663)
        self.assertEqual(project.sessions[0].observations[0].filter,   7)
        
        # Observational setup - 2
        self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_SOL')
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
        self.assertEqual(project.sessions[0].observations[1].dur,  10000)
        self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
        self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
        self.assertEqual(project.sessions[0].observations[1].filter,   7)
        
    def test_sol_update(self):
        """Test updating TRK_SOL values."""
        
        project = sdfADP.parseSDF(solFile)
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:00:15")
        project.sessions[0].observations[1].setDuration(timedelta(seconds=15))
        project.sessions[0].observations[1].setFrequency1(75e6)
        project.sessions[0].observations[1].setFrequency2(76e6)
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        self.assertEqual(project.sessions[0].observations[1].dur,  15000)
        self.assertEqual(project.sessions[0].observations[1].freq1, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].freq2, 1665395482)
        
    def test_sol_write(self):
        """Test writing a TRK_SOL SDF file."""
        
        project = sdfADP.parseSDF(solFile)
        # Fix for LWA-SV only going up to filter code 6
        for obs in project.sessions[0].observations:
            obs.setFilter(6)
        out = project.render()
        
    def test_sol_errors(self):
        """Test various TRK_SOL SDF errors."""
        
        project = sdfADP.parseSDF(solFile)
        
        # Bad beam
        project.sessions[0].drxBeam = 6
        self.assertFalse(project.validate())
        
        # No beam
        project.sessions[0].drxBeam = -1
        self.assertFalse(project.validate())
        
        # Bad filter
        project.sessions[0].observations[0].filter = 7
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.sessions[0].observations[0].filter = 6
        project.sessions[0].observations[0].frequency1 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        project.sessions[0].observations[0].frequency1 = 38.0e6
        project.sessions[0].observations[0].frequency2 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.sessions[0].observations[0].frequency2 = 38.0e6
        project.sessions[0].observations[0].duration = '96:00:00.000'
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - TRK_JOV ###
    
    def test_jov_parse(self):
        """Test reading in a TRK_JOV SDF file."""
        
        project = sdfADP.parseSDF(jovFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 2)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode, 'TRK_JOV')
        self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[0].mpm,      0)
        self.assertEqual(project.sessions[0].observations[0].dur,  10000)
        self.assertEqual(project.sessions[0].observations[0].freq1,  438261968)
        self.assertEqual(project.sessions[0].observations[0].freq2, 1928352663)
        self.assertEqual(project.sessions[0].observations[0].filter,   7)
        
        # Observational setup - 2
        self.assertEqual(project.sessions[0].observations[1].mode, 'TRK_JOV')
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  10000)
        self.assertEqual(project.sessions[0].observations[1].dur,  10000)
        self.assertEqual(project.sessions[0].observations[1].freq1,  832697741)
        self.assertEqual(project.sessions[0].observations[1].freq2, 1621569285)
        self.assertEqual(project.sessions[0].observations[1].filter,   7)
        
    def test_jov_update(self):
        """Test updating TRK_JOV values."""
        
        project = sdfADP.parseSDF(jovFile)
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:00:15")
        project.sessions[0].observations[1].setDuration(timedelta(seconds=15))
        project.sessions[0].observations[1].setFrequency1(75e6)
        project.sessions[0].observations[1].setFrequency2(76e6)
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        self.assertEqual(project.sessions[0].observations[1].dur,  15000)
        self.assertEqual(project.sessions[0].observations[1].freq1, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].freq2, 1665395482)
        
    def test_jov_write(self):
        """Test writing a TRK_JOV SDF file."""
        
        project = sdfADP.parseSDF(jovFile)
        # Fix for LWA-SV only going up to filter code 6
        for obs in project.sessions[0].observations:
            obs.setFilter(6)
        out = project.render()
        
    def test_jov_errors(self):
        """Test various TRK_JOV SDF errors."""
        
        project = sdfADP.parseSDF(jovFile)
        
        # Bad beam
        project.sessions[0].drxBeam = 6
        self.assertFalse(project.validate())
        
        # No beam
        project.sessions[0].drxBeam = -1
        self.assertFalse(project.validate())
        
        # Bad filter
        project.sessions[0].observations[0].filter = 7
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.sessions[0].observations[0].filter = 6
        project.sessions[0].observations[0].frequency1 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        project.sessions[0].observations[0].frequency1 = 38.0e6
        project.sessions[0].observations[0].frequency2 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.sessions[0].observations[0].frequency2 = 38.0e6
        project.sessions[0].observations[0].duration = '96:00:00.000'
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - STEPPED ###
    
    def test_stp_parse(self):
        """Test reading in a STEPPED SDF file."""
        
        project = sdfADP.parseSDF(stpFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 2)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode, 'STEPPED')
        self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[0].mpm, 440000)
        self.assertEqual(project.sessions[0].observations[0].dur, 300000)
        self.assertEqual(project.sessions[0].observations[0].filter,   7)
        self.assertEqual(project.sessions[0].observations[0].obsFEE[0], [1,1])
        self.assertEqual(project.sessions[0].observations[0].aspFlt[0], 2)
        self.assertEqual(project.sessions[0].observations[0].aspAT1[0], 10)
        self.assertEqual(project.sessions[0].observations[0].aspAT2[0], 12)
        self.assertEqual(project.sessions[0].observations[0].aspATS[0], 14)
        
        # Steps - 1
        self.assertEqual(len(project.sessions[0].observations[0].steps), 4)
        for i in xrange(4):
            self.assertEqual(project.sessions[0].observations[0].steps[i].RADec, project.sessions[0].observations[0].RADec)
            self.assertEqual(project.sessions[0].observations[0].steps[i].freq1,  832697741)
            self.assertEqual(project.sessions[0].observations[0].steps[i].freq2, 1621569285)
        self.assertAlmostEqual(project.sessions[0].observations[0].steps[0].c1, 90.0, 6)
        self.assertAlmostEqual(project.sessions[0].observations[0].steps[0].c2, 45.0, 6)
        self.assertEqual(project.sessions[0].observations[0].steps[0].dur, 60000)
        self.assertAlmostEqual(project.sessions[0].observations[0].steps[-1].c1, 0.0, 6)
        self.assertAlmostEqual(project.sessions[0].observations[0].steps[-1].c2, 1.0, 6)
        self.assertEqual(project.sessions[0].observations[0].steps[-1].dur, 120000)
        
        # Observational setup - 2
        self.assertEqual(project.sessions[0].observations[1].mode, 'STEPPED')
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm, 800000)
        self.assertEqual(project.sessions[0].observations[1].dur, 180000)
        self.assertEqual(project.sessions[0].observations[1].filter,   7)
        self.assertEqual(project.sessions[0].observations[1].obsFEE[0], [1,0])
        self.assertEqual(project.sessions[0].observations[1].aspFlt[0], 1)
        self.assertEqual(project.sessions[0].observations[1].aspAT1[0], 11)
        self.assertEqual(project.sessions[0].observations[1].aspAT2[0], 13)
        self.assertEqual(project.sessions[0].observations[1].aspATS[0], 15)
        
        # Steps - 2
        self.assertEqual(len(project.sessions[0].observations[1].steps), 2)
        for i in xrange(2):
            self.assertEqual(project.sessions[0].observations[1].steps[i].RADec, project.sessions[0].observations[1].RADec)
            self.assertEqual(project.sessions[0].observations[1].steps[i].freq1,  832697741)
            self.assertEqual(project.sessions[0].observations[1].steps[i].freq2, 1621569285)
        self.assertAlmostEqual(project.sessions[0].observations[1].steps[0].c1, 0.0, 6)
        self.assertAlmostEqual(project.sessions[0].observations[1].steps[0].c2, 90.0, 6)
        self.assertEqual(project.sessions[0].observations[1].steps[0].dur, 60000)
        self.assertAlmostEqual(project.sessions[0].observations[1].steps[-1].c1, 12.0, 6)
        self.assertAlmostEqual(project.sessions[0].observations[1].steps[-1].c2, 80.0, 6)
        self.assertEqual(project.sessions[0].observations[1].steps[-1].dur, 120000)
        
    def test_stp_update(self):
        """Test updating a STEPPED SDF file."""
        
        project = sdfADP.parseSDF(stpFile)
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:00:15")
        for step in project.sessions[0].observations[1].steps:
            step.setDuration(timedelta(seconds=15))
            step.setFrequency1(75e6)
            step.setFrequency2(76e6)
            step.setC1(ephem.hours('10:30:00'))
            step.setC2(ephem.degrees('89:30:00'))
        project.sessions[0].observations[1].update()
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        self.assertEqual(project.sessions[0].observations[1].dur,  30000)
        for step in project.sessions[0].observations[1].steps:
            self.assertEqual(step.dur, 15000)
            self.assertEqual(step.freq1, 1643482384)
            self.assertEqual(step.freq2, 1665395482)
            self.assertEqual(step.c1, 10.5)
            self.assertEqual(step.c2, 89.5)
            
        project = sdfADP.parseSDF(stpFile)
        project.sessions[0].observations[1].RADec = False
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:00:15")
        for step in project.sessions[0].observations[1].steps:
            step.RADec = False
            step.setDuration(timedelta(seconds=15))
            step.setFrequency1(75e6)
            step.setFrequency2(76e6)
            step.setC1(ephem.hours('10:30:00'))
            step.setC2(ephem.degrees('89:30:00'))
        project.sessions[0].observations[1].update()
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        self.assertEqual(project.sessions[0].observations[1].dur,  30000)
        for step in project.sessions[0].observations[1].steps:
            self.assertFalse(step.RADec)
            self.assertEqual(step.dur, 15000)
            self.assertEqual(step.freq1, 1643482384)
            self.assertEqual(step.freq2, 1665395482)
            self.assertEqual(step.c1, 10.5*15)
            self.assertEqual(step.c2, 89.5)
        
    def test_stp_write(self):
        """Test writing a STEPPED SDF file."""
        
        project = sdfADP.parseSDF(stpFile)
        # Fix for LWA-SV only going up to filter code 6
        for obs in project.sessions[0].observations:
            obs.setFilter(6)
        out = project.render()
        
    def test_stp_errors(self):
        """Test various STEPPED SDF errors."""
        
        project = sdfADP.parseSDF(stpFile)
        
        # Bad beam
        project.sessions[0].drxBeam = 6
        self.assertFalse(project.validate())
        
        # No beam
        project.sessions[0].drxBeam = -1
        self.assertFalse(project.validate())
        
        # Bad filter
        project.sessions[0].observations[0].filter = 7
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.sessions[0].observations[0].filter = 6
        project.sessions[0].observations[0].steps[0].frequency1 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        project.sessions[0].observations[0].steps[0].frequency1 = 38.0e6
        project.sessions[0].observations[0].steps[1].frequency2 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        # Bad duration
        project.sessions[0].observations[0].steps[1].frequency2 = 38.0e6
        project.sessions[0].observations[0].steps[2].duration = '96:00:00.000'
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
    ### DRX - STEPPED with delays and gains ###
    
    def test_spc_parse(self):
        """Test reading in a STEPPED Delay and Gain SDF file."""
        
        project = sdfADP.parseSDF(spcFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 1)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode, 'STEPPED')
        self.assertEqual(project.sessions[0].observations[0].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[0].mpm, 440000)
        self.assertEqual(project.sessions[0].observations[0].dur,  60000)
        self.assertEqual(project.sessions[0].observations[0].filter,   7)
        
        # Steps - 1
        self.assertEqual(len(project.sessions[0].observations[0].steps), 1)
        self.assertEqual(project.sessions[0].observations[0].steps[0].RADec, project.sessions[0].observations[0].RADec)
        self.assertAlmostEqual(project.sessions[0].observations[0].steps[0].c1, 90.0, 6)
        self.assertAlmostEqual(project.sessions[0].observations[0].steps[0].c2, 45.0, 6)
        self.assertEqual(project.sessions[0].observations[0].steps[0].freq1,  832697741)
        self.assertEqual(project.sessions[0].observations[0].steps[0].freq2, 1621569285)
        self.assertEqual(project.sessions[0].observations[0].steps[0].dur, 60000)
        
        # Delays - 1
        for i in xrange(256):
            self.assertEqual(project.sessions[0].observations[0].steps[0].delays[i], 0)
            
        # Gains - 1
        for i in xrange(256):
            self.assertEqual(project.sessions[0].observations[0].steps[0].gains[i][0][0], 1)
            self.assertEqual(project.sessions[0].observations[0].steps[0].gains[i][0][1], 0)
            self.assertEqual(project.sessions[0].observations[0].steps[0].gains[i][1][0], 0)
            self.assertEqual(project.sessions[0].observations[0].steps[0].gains[i][1][1], 1)
        
    ### TBF ###
    
    def test_tbf_parse(self):
        """Test reading in a TBF SDF file."""
        
        project = sdfADP.parseSDF(tbfFile)
        
        # Basic file structure
        self.assertEqual(len(project.sessions), 1)
        self.assertEqual(len(project.sessions[0].observations), 5)
        
        # Observational setup - 1
        self.assertEqual(project.sessions[0].observations[0].mode,      'TBF')
        self.assertEqual(project.sessions[0].observations[0].mjd,       57865)
        self.assertEqual(project.sessions[0].observations[0].mpm,    46800000)
        self.assertEqual(project.sessions[0].observations[0].freq1, 876523938)
        self.assertEqual(project.sessions[0].observations[0].filter,        6)
        
        # Observational setup - 2
        self.assertEqual(project.sessions[0].observations[1].mode,       'TBF')
        self.assertEqual(project.sessions[0].observations[1].mjd,        57865)
        self.assertEqual(project.sessions[0].observations[1].mpm,     47100000)
        self.assertEqual(project.sessions[0].observations[1].freq1, 1095654922)
        self.assertEqual(project.sessions[0].observations[1].filter,         6)
        
    def test_tbf_update(self):
        """Test updating TRK_SOL values."""
        
        project = sdfADP.parseSDF(tbfFile)
        project.sessions[0].observations[1].setStart("MST 2011 Feb 23 17:10:15")
        
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  615000)
        
    def test_tbf_write(self):
        """Test writing a TBF SDF file."""
        
        project = sdfADP.parseSDF(tbfFile)
        out = project.render()
        
    def test_tbf_errors(self):
        """Test various TBF SDF errors."""
        
        project = sdfADP.parseSDF(tbfFile)
        
        # Bad number of TBF samples
        project.sessions[0].observations[0].samples = 6*196000000
        self.assertFalse(project.validate())
        
        # Bad filter
        project.sessions[0].observations[0].filter = 7
        self.assertFalse(project.validate())
        
        # Bad frequency
        project.sessions[0].observations[0].filter = 6
        project.sessions[0].observations[0].frequency1 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
        project.sessions[0].observations[0].frequency1 = 38.0e6
        project.sessions[0].observations[0].frequency2 = 90.0e6
        project.sessions[0].observations[0].update()
        self.assertFalse(project.validate())
        
    ### Misc. ###
    
    def test_auto_update(self):
        """Test project auto-update on render."""
        
        # Part 1 - frequency and duration
        project = sdfADP.parseSDF(drxFile)
        project.sessions[0].observations[1].frequency1 = 75e6
        project.sessions[0].observations[1].duration = '00:01:31.000'

        fh = open(os.path.join(self.testPath, 'sdf.txt'), 'w')		
        fh.write(project.render())
        fh.close()
        
        project = sdfADP.parseSDF(os.path.join(self.testPath, 'sdf.txt'))
        self.assertEqual(project.sessions[0].observations[1].freq1, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].dur, 91000)
        
        # Part 2 - frequency and duration (timedelta)
        project = sdfADP.parseSDF(drxFile)
        project.sessions[0].observations[1].frequency1 = 75e6
        project.sessions[0].observations[1].duration = timedelta(minutes=1, seconds=31, microseconds=1000)

        fh = open(os.path.join(self.testPath, 'sdf.txt'), 'w')		
        fh.write(project.render())
        fh.close()
        
        project = sdfADP.parseSDF(os.path.join(self.testPath, 'sdf.txt'))
        self.assertEqual(project.sessions[0].observations[1].freq1, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].dur, 91001)
        
        # Part 3 - frequency and start time
        project = sdfADP.parseSDF(drxFile)
        project.sessions[0].observations[1].frequency2 = 75e6
        project.sessions[0].observations[1].start = "MST 2011 Feb 23 17:00:15"

        fh = open(os.path.join(self.testPath, 'sdf.txt'), 'w')		
        fh.write(project.render())
        fh.close()
        
        project = sdfADP.parseSDF(os.path.join(self.testPath, 'sdf.txt'))
        self.assertEqual(project.sessions[0].observations[1].freq2, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  15000)
        
        # Part 4 - frequency and start time (timedelta)
        project = sdfADP.parseSDF(drxFile)
        _MST = pytz.timezone('US/Mountain')
        project.sessions[0].observations[1].frequency2 = 75e6
        project.sessions[0].observations[1].start = _MST.localize(datetime(2011, 2, 23, 17, 00, 30, 1000))

        fh = open(os.path.join(self.testPath, 'sdf.txt'), 'w')		
        fh.write(project.render())
        fh.close()
        
        project = sdfADP.parseSDF(os.path.join(self.testPath, 'sdf.txt'))
        self.assertEqual(project.sessions[0].observations[1].freq2, 1643482384)
        self.assertEqual(project.sessions[0].observations[1].mjd,  55616)
        self.assertEqual(project.sessions[0].observations[1].mpm,  30001)
        
    def test_set_station(self):
        """Test the set stations functionlity."""
        
        project = sdfADP.parseSDF(drxFile)
        project.sessions[0].setStation(lwasv)
        # Fix for LWA-SV only going up to filter code 6
        for obs in project.sessions[0].observations:
            obs.setFilter(6)
        self.assertTrue(project.validate())
        
        self.assertRaises(RuntimeError, project.sessions[0].setStation, lwa1)
        
    def test_is_valid(self):
        """Test whether or not isValid works."""
        
        self.assertTrue(sdfADP.isValid(tbnFile))
        self.assertFalse(sdfADP.isValid(drxFile))
        self.assertFalse(sdfADP.isValid(solFile))
        self.assertFalse(sdfADP.isValid(jovFile))
        self.assertFalse(sdfADP.isValid(stpFile))
        self.assertFalse(sdfADP.isValid(spcFile))
        self.assertTrue(sdfADP.isValid(tbfFile))
        
    def test_is_not_valid(self):
        """Test whether or not isValid works on LWA1 files."""
        
        self.assertFalse(sdfADP.isValid(tbwFile))
        
    def test_username(self):
        """Test setting auto-copy parameters."""
        
        project = sdfADP.parseSDF(drxFile)
        project.sessions[0].setDataReturnMethod('UCF')
        project.sessions[0].setUCFUsername('jdowell')
        # Fix for LWA-SV only going up to filter code 6
        for obs in project.sessions[0].observations:
            obs.setFilter(6)
        out = project.render()
        
        self.assertTrue(out.find('Requested data return method is UCF') > 0)
        self.assertTrue(out.find('ucfuser:jdowell') > 0)
        
        fh = open(os.path.join(self.testPath, 'sdf.txt'), 'w')		
        fh.write(out)
        fh.close()
        
        project = sdfADP.parseSDF(os.path.join(self.testPath, 'sdf.txt'))
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


class sdf_adp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.sdf units 
    tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(sdf_adp_tests)) 


if __name__ == '__main__':
    unittest.main()
