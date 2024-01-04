"""
Unit test for the lsl.common.metabundleADP module.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import unittest

from lsl.common import metabundleADP


__version__  = "0.3"
__author__    = "Jayce Dowell"

mdbFile = os.path.join(os.path.dirname(__file__), 'data', 'metadata.tgz')
mdbFileOld0 = os.path.join(os.path.dirname(__file__), 'data', 'metadata-old-0.tgz')
mdbFileOld1 = os.path.join(os.path.dirname(__file__), 'data', 'metadata-old-1.tgz')
mdbFileADP = os.path.join(os.path.dirname(__file__), 'data', 'metadata-adp.tgz')
mdbFileNDP = os.path.join(os.path.dirname(__file__), 'data', 'metadata-ndp.tgz')
mdbFileGDB = os.path.join(os.path.dirname(__file__), 'data', 'metadata-gdb.tgz')
mdbFileGDBOld0 = os.path.join(os.path.dirname(__file__), 'data', 'metadata-gdb-old-0.tgz')


class metabundle_tests_adp(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.metabundle
    module."""
    
    def test_ss(self):
        """Test the session specification utilties."""
        
        ses = metabundleADP.get_session_spec(mdbFileADP)
        obs = metabundleADP.get_observation_spec(mdbFileADP)
        
        # Check session start time
        self.assertEqual(ses['mjd'], 57774)
        self.assertEqual(ses['mpm'], 29970000)
        
        # Check the duration
        self.assertEqual(ses['dur'], obs[0]['dur'] + 10000)
        
        # Check the number of observations
        self.assertEqual(ses['nobs'], len(obs))
    
    def test_os(self):
        """Test the observation specification utilities."""
        
        obs1 = metabundleADP.get_observation_spec(mdbFileADP)
        obs2 = metabundleADP.get_observation_spec(mdbFileADP, obs_id=1)
        
        # Check if the right observation is returned
        self.assertEqual(obs1[0], obs2)
        
        # Check the mode
        self.assertEqual(obs2['mode'], 1)
        
        # Check the time
        self.assertEqual(obs2['mjd'], 57774)
        self.assertEqual(obs2['mpm'], 29975000)
        
    def test_cs(self):
        """Test the command script utilities."""
        
        cmnds = metabundleADP.get_command_script(mdbFileADP)
        
        # Check number of command
        self.assertEqual(len(cmnds), 50)
        
        # Check the first and last commands
        self.assertEqual(cmnds[ 0]['command_id'], 'NUL')
        self.assertEqual(cmnds[-2]['command_id'], 'STP')
        self.assertEqual(cmnds[-1]['command_id'], 'ESN')
        
        # Check the counds of DP BAM commands
        nBAM = 0
        for cmnd in cmnds:
            if cmnd['command_id'] == 'BAM':
                nBAM += 1
        self.assertEqual(nBAM, 40)
        
    def test_sm(self):
        """Test the session metadata utilties."""
        
        sm = metabundleADP.get_session_metadata(mdbFileADP)
        
        # Make sure all of the observations are done
        self.assertEqual(len(sm.keys()), 1)
        
    def test_sdf(self):
        """Test building a SDF from a tarball."""
        
        sdf = metabundleADP.get_sdf(mdbFileADP)
        
    def test_station(self):
        """Test building a station from a tarball."""
        
        station = metabundleADP.get_station(mdbFileADP)
        
    def test_sdm(self):
        """Test the station dynamic MIB utilties."""
        
        sm = metabundleADP.get_sdm(mdbFileADP)
        
    def test_metadata(self):
        """Test the observation metadata utility."""
        
        fileInfo = metabundleADP.get_session_metadata(mdbFileADP)
        self.assertEqual(len(fileInfo.keys()), 1)
        
        # File tag
        self.assertEqual(fileInfo[1]['tag'], '057774_000770030')
        
        # DRSU barcode
        self.assertEqual(fileInfo[1]['barcode'], 'S10TCC13S0016')
        
    def test_aspconfig(self):
        """Test retrieving the ASP configuration."""
        
        # Beginning config.
        aspConfig = metabundleADP.get_asp_configuration_summary(mdbFileADP, which='beginning')
        self.assertEqual(aspConfig['asp_filter'],      0)
        self.assertEqual(aspConfig['asp_atten_1'],     6)
        self.assertEqual(aspConfig['asp_atten_2'],     5)
        self.assertEqual(aspConfig['asp_atten_split'],15)
        
        # End config.
        aspConfig = metabundleADP.get_asp_configuration_summary(mdbFileADP, which='End')
        self.assertEqual(aspConfig['asp_filter'],      0)
        self.assertEqual(aspConfig['asp_atten_1'],     6)
        self.assertEqual(aspConfig['asp_atten_2'],     5)
        self.assertEqual(aspConfig['asp_atten_split'],15)
        
        # Unknown code
        self.assertRaises(ValueError, metabundleADP.get_asp_configuration_summary, mdbFileADP, 'middle')
        
    def test_aspconfig_gdbm(self):
        """Test retrieving the ASP configuration from a GDBM MIB."""
        
        # Beginning config.
        aspConfig = metabundleADP.get_asp_configuration_summary(mdbFileGDB, which='beginning')
        self.assertEqual(aspConfig['asp_filter'],      0)
        self.assertEqual(aspConfig['asp_atten_1'],     6)
        self.assertEqual(aspConfig['asp_atten_2'],     5)
        self.assertEqual(aspConfig['asp_atten_split'],15)
        
        # End config.
        aspConfig = metabundleADP.get_asp_configuration_summary(mdbFileGDB, which='End')
        self.assertEqual(aspConfig['asp_filter'],      0)
        self.assertEqual(aspConfig['asp_atten_1'],     6)
        self.assertEqual(aspConfig['asp_atten_2'],     5)
        self.assertEqual(aspConfig['asp_atten_split'],15)
        
        # Unknown code
        self.assertRaises(ValueError, metabundleADP.get_asp_configuration_summary, mdbFileGDB, 'middle')
        
    def test_is_valid(self):
        """Test whether or not is_valid works."""
        
        self.assertTrue(metabundleADP.is_valid(mdbFileADP))
        self.assertTrue(metabundleADP.is_valid(mdbFileGDB))
        self.assertTrue(metabundleADP.is_valid(mdbFileGDBOld0))
        
    def test_is_not_valid(self):
        """Test whether or not is_valid works on LWA1 files."""
        
        self.assertFalse(metabundleADP.is_valid(mdbFile))
        self.assertFalse(metabundleADP.is_valid(mdbFileNDP))
        self.assertFalse(metabundleADP.is_valid(mdbFileOld0))
        self.assertFalse(metabundleADP.is_valid(mdbFileOld1))


class metabundle_adp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.metabundleADP
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(metabundle_tests_adp))
        
if __name__ == '__main__':
    unittest.main()
