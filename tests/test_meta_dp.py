"""
Unit test for the DP portion of the lsl.common.metabundle module.
"""

import os
import unittest

from lsl.common.mcs import CommandID, ObservingMode
from lsl.common import metabundle, metabundleDP, stations

run_gdbm_tests = False
try:
    from dbm import gnu
    run_gdbm_tests = True
except ImportError:
    pass


__version__  = "0.7"
__author__    = "Jayce Dowell"

mdbFile = os.path.join(os.path.dirname(__file__), 'data', 'metadata.tgz')
mdbFileOld0 = os.path.join(os.path.dirname(__file__), 'data', 'metadata-old-0.tgz')
mdbFileOld1 = os.path.join(os.path.dirname(__file__), 'data', 'metadata-old-1.tgz')
mdbFileADP = os.path.join(os.path.dirname(__file__), 'data', 'metadata-adp.tgz')
mdbFileNDP = os.path.join(os.path.dirname(__file__), 'data', 'metadata-ndp.tgz')
mdbFileGDB = os.path.join(os.path.dirname(__file__), 'data', 'metadata-gdb.tgz')
mdbFileGDBOld0 = os.path.join(os.path.dirname(__file__), 'data', 'metadata-gdb-old-0.tgz')

class metabundle_dp_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.metabundle
    module."""
    
    def test_ss(self):
        """Test the session specification utilties."""
        
        ses = metabundle.get_session_spec(mdbFile)
        obs = metabundle.get_observation_spec(mdbFile)
        
        # Check session start time
        self.assertEqual(ses['mjd'], 56742)
        self.assertEqual(ses['mpm'], 4914000)
        
        # Check the duration
        self.assertEqual(ses['dur'], obs[0]['dur'] + 10000)
        
        # Check the number of observations
        self.assertEqual(ses['nobs'], len(obs))
    
    def test_os(self):
        """Test the observation specification utilities."""
        
        obs1 = metabundle.get_observation_spec(mdbFile)
        obs2 = metabundle.get_observation_spec(mdbFile, obs_id=1)
        
        # Check if the right observation is returned
        self.assertEqual(obs1[0], obs2)
        
        # Check the mode
        self.assertEqual(obs2['mode'], ObservingMode.TRK_RADEC)
        
        # Check the time
        self.assertEqual(obs2['mjd'], 56742)
        self.assertEqual(obs2['mpm'], 4919000)
        
    def test_cs(self):
        """Test the command script utilities."""
        
        cmnds = metabundle.get_command_script(mdbFile)
        
        # Check number of command
        self.assertEqual(len(cmnds), 150)
        
        # Check the first and last commands
        self.assertEqual(cmnds[ 0]['command_id'], CommandID.NUL)
        self.assertEqual(cmnds[-2]['command_id'], CommandID.OBE)
        self.assertEqual(cmnds[-1]['command_id'], CommandID.ESN)
        
        # Check the counds of DP BAM commands
        nBAM = 0
        for cmnd in cmnds:
            if cmnd['command_id'] == CommandID.BAM:
                nBAM += 1
        self.assertEqual(nBAM, 143)
        
    def test_sm(self):
        """Test the session metadata utilties."""
        
        sm = metabundle.get_session_metadata(mdbFile)
        
        # Make sure all of the observations are done
        self.assertEqual(len(sm.keys()), 1)
        
    def test_sdf(self):
        """Test building a SDF from a tarball."""
        
        sdf = metabundle.get_sdf(mdbFile)
        
    def test_beamformer_min_delay(self):
        """Test reading the beamformer minimum delay info."""
        
        md = metabundle.get_beamformer_min_delay(mdbFile)
        
    def test_station(self):
        """Test building a station from a tarball."""
        
        station = metabundle.get_station(mdbFile, apply_sdm=False)
        station = metabundle.get_station(mdbFile)
        
    def test_sdm(self):
        """Test the station dynamic MIB utilties."""
        
        sm = metabundle.get_sdm(mdbFile)
        
    def test_sdm_dynamic_update(self):
        """Test applying a station dynamic MIB to a LWAStation object."""
        
        station = stations.lwa1
        sm = metabundle.get_sdm(mdbFile)
        newAnts = sm.update_antennas(station.antennas)
        
    def test_metadata(self):
        """Test the observation metadata utility."""
        
        fileInfo = metabundle.get_session_metadata(mdbFile)
        self.assertEqual(len(fileInfo.keys()), 1)
        
        # File tag
        self.assertEqual(fileInfo[1]['tag'], '056742_000440674')
        
        # DRSU barcode
        self.assertEqual(fileInfo[1]['barcode'], 'S15TCV23S0001')
        
    @unittest.skipUnless(run_gdbm_tests, "requires the 'dbm.gnu' module")
    def test_aspconfig(self):
        """Test retrieving the ASP configuration."""
        
        # Beginning config.
        aspConfig = metabundle.get_asp_configuration_summary(mdbFile, which='beginning')
        self.assertEqual(aspConfig['asp_filter'],      1)
        self.assertEqual(aspConfig['asp_atten_1'],    13)
        self.assertEqual(aspConfig['asp_atten_2'],    13)
        self.assertEqual(aspConfig['asp_atten_split'],15)
        
        # End config.
        aspConfig = metabundle.get_asp_configuration_summary(mdbFile, which='End')
        self.assertEqual(aspConfig['asp_filter'],      1)
        self.assertEqual(aspConfig['asp_atten_1'],    13)
        self.assertEqual(aspConfig['asp_atten_2'],    13)
        self.assertEqual(aspConfig['asp_atten_split'],15)
        
        # Unknown code
        self.assertRaises(ValueError, metabundle.get_asp_configuration_summary, mdbFile, 'middle')
        
    def test_is_valid(self):
        """Test whether or not is_valid works."""
        
        for filename in (mdbFile, mdbFileADP, mdbFileNDP):
            self.assertTrue(metabundle.is_valid(mdbFile))
            
    def test_get_style(self):
        """Test whether or not get_style works."""
        
        self.assertEqual(metabundle.get_style(mdbFile),    'lsl.common.metabundleDP')
        self.assertEqual(metabundle.get_style(mdbFileADP), 'lsl.common.metabundleADP')
        self.assertEqual(metabundle.get_style(mdbFileNDP), 'lsl.common.metabundleNDP')
        
    def test_is_valid_dp(self):
        """Test whether or not the DP-specific is_valid works."""
        
        self.assertTrue(metabundleDP.is_valid(mdbFile))
        
    def test_is_not_valid_dp(self):
        """Test whether or not the DP-specific is_valid works on LWA-SV and LWA-NA files."""
        
        self.assertFalse(metabundleDP.is_valid(mdbFileADP))
        self.assertFalse(metabundleDP.is_valid(mdbFileNDP))
        self.assertFalse(metabundleDP.is_valid(mdbFileGDB))
        self.assertFalse(metabundleDP.is_valid(mdbFileGDBOld0))


class metabundle_dp_tests_old_0(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.metabundle
    module based on the tarball format supported in LSL 0.5.x."""
    
    def test_ss(self):
        """Test the session specification utilties."""
        
        ses = metabundle.get_session_spec(mdbFileOld0)
        obs = metabundle.get_observation_spec(mdbFileOld0)
        
        # Check session start time
        self.assertEqual(ses['mjd'], 56013)
        self.assertEqual(ses['mpm'], 25855000)
        
        # Check the duration
        self.assertEqual(ses['dur'], obs[0]['dur'] + 10000)
        
        # Check the number of observations
        self.assertEqual(ses['nobs'], len(obs))
        
    def test_os(self):
        """Test the observation specification utilities."""
        
        obs1 = metabundle.get_observation_spec(mdbFileOld0)
        obs2 = metabundle.get_observation_spec(mdbFileOld0, obs_id=1)
        
        # Check if the right observation is returned
        self.assertEqual(obs1[0], obs2)
        
        # Check the mode
        self.assertEqual(obs2['mode'], ObservingMode.TRK_RADEC)
        
        # Check the time
        self.assertEqual(obs2['mjd'], 56013)
        self.assertEqual(obs2['mpm'], 25860000)
        
    def test_cs(self):
        """Test the command script utilities."""
        
        cmnds = metabundle.get_command_script(mdbFileOld0)
        
        # Check number of command
        self.assertEqual(len(cmnds), 491)
        
        # Check the first and last commands
        self.assertEqual(cmnds[ 0]['command_id'], CommandID.NUL)
        self.assertEqual(cmnds[-2]['command_id'], CommandID.OBE)
        self.assertEqual(cmnds[-1]['command_id'], CommandID.ESN)
        
        # Check the counds of DP BAM commands
        nBAM = 0
        for cmnd in cmnds:
            if cmnd['command_id'] == CommandID.BAM:
                nBAM += 1
        self.assertEqual(nBAM, 484)
        
    def test_sm(self):
        """Test the session metadata utilties."""
        
        sm = metabundle.get_session_metadata(mdbFileOld0)
        
        # Make sure all of the observations are done
        self.assertEqual(len(sm.keys()), 1)
        
    def test_sdf(self):
        """Test building a SDF from a tarball."""
        
        sdf = metabundle.get_sdf(mdbFileOld0)
        
    def test_sdm(self):
        """Test the station dynamic MIB utilties."""
        
        sm = metabundle.get_sdm(mdbFileOld0)
        
    def test_is_valid(self):
        """Test whether or not is_valid works."""
        
        self.assertTrue(metabundleDP.is_valid(mdbFileOld0))


class metabundle_dp_tests_old_1(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.metabundle
    module."""
    
    def test_ss(self):
        """Test the session specification utilties."""
        
        ses = metabundle.get_session_spec(mdbFileOld1)
        obs = metabundle.get_observation_spec(mdbFileOld1)
        
        # Check session start time
        self.assertEqual(ses['mjd'], 56492)
        self.assertEqual(ses['mpm'], 68995000)
        
        # Check the duration
        self.assertEqual(ses['dur'], obs[0]['dur'] + 10000)
        
        # Check the number of observations
        self.assertEqual(ses['nobs'], len(obs))
        
    def test_os(self):
        """Test the observation specification utilities."""
        
        obs1 = metabundle.get_observation_spec(mdbFileOld1)
        obs2 = metabundle.get_observation_spec(mdbFileOld1, obs_id=1)
        
        # Check if the right observation is returned
        self.assertEqual(obs1[0], obs2)
        
        # Check the mode
        self.assertEqual(obs2['mode'], ObservingMode.TRK_RADEC)
        
        # Check the time
        self.assertEqual(obs2['mjd'], 56492)
        self.assertEqual(obs2['mpm'], 69000000)
        
    def test_cs(self):
        """Test the command script utilities."""
        
        cmnds = metabundle.get_command_script(mdbFileOld1)
        
        # Check number of command
        self.assertEqual(len(cmnds), 8)
        
        # Check the first and last commands
        self.assertEqual(cmnds[ 0]['command_id'], CommandID.NUL)
        self.assertEqual(cmnds[-2]['command_id'], CommandID.OBE)
        self.assertEqual(cmnds[-1]['command_id'], CommandID.ESN)
        
        # Check the counds of DP BAM commands
        nBAM = 0
        for cmnd in cmnds:
            if cmnd['command_id'] == CommandID.BAM:
                nBAM += 1
        self.assertEqual(nBAM, 1)
        
    def test_sm(self):
        """Test the session metadata utilties."""
        
        sm = metabundle.get_session_metadata(mdbFileOld1)
        
        # Make sure all of the observations are done
        self.assertEqual(len(sm.keys()), 1)
        
    def test_sdf(self):
        """Test building a SDF from a tarball."""
        
        sdf = metabundle.get_sdf(mdbFileOld1)
        
    def test_sdm(self):
        """Test the station dynamic MIB utilties."""
        
        sm = metabundle.get_sdm(mdbFileOld1)
        
    def test_metadata(self):
        """Test the observation metadata utility."""
        
        fileInfo = metabundle.get_session_metadata(mdbFileOld1)
        self.assertEqual(len(fileInfo.keys()), 1)
        
        # File tag
        self.assertEqual(fileInfo[1]['tag'], '056492_000000094')
        
        # DRSU barcode
        self.assertEqual(fileInfo[1]['barcode'], 'S10TCC13S0007')
        
    @unittest.skipUnless(run_gdbm_tests, "requires the 'dbm.gnu' module")
    def test_aspconfig(self):
        """Test retrieving the ASP configuration."""
        
        # Beginning config.
        aspConfig = metabundle.get_asp_configuration_summary(mdbFileOld1, which='beginning')
        self.assertEqual(aspConfig['asp_filter'],      3)
        self.assertEqual(aspConfig['asp_atten_1'],     0)
        self.assertEqual(aspConfig['asp_atten_2'],     0)
        self.assertEqual(aspConfig['asp_atten_split'], 0)
        
        # End config.
        aspConfig = metabundle.get_asp_configuration_summary(mdbFileOld1, which='End')
        self.assertEqual(aspConfig['asp_filter'],      1)
        self.assertEqual(aspConfig['asp_atten_1'],    13)
        self.assertEqual(aspConfig['asp_atten_2'],    13)
        self.assertEqual(aspConfig['asp_atten_split'], 0)
        
        # Unknown code
        self.assertRaises(ValueError, metabundle.get_asp_configuration_summary, mdbFileOld1, 'middle')
        
        # Not a summary
        aspConfig = metabundle.get_asp_configuration(mdbFileOld1, which='End')
        self.assertEqual(aspConfig['asp_filter'][0],      1)
        self.assertEqual(aspConfig['asp_atten_1'][0],    13)
        self.assertEqual(aspConfig['asp_atten_2'][0],    13)
        self.assertEqual(aspConfig['asp_atten_split'][0], 0)
        
    def test_is_valid(self):
        """Test whether or not is_valid works."""
        
        self.assertTrue(metabundleDP.is_valid(mdbFileOld1))


class metabundle_dp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.metabundle
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(metabundle_dp_tests))        
        self.addTests(loader.loadTestsFromTestCase(metabundle_dp_tests_old_0))
        self.addTests(loader.loadTestsFromTestCase(metabundle_dp_tests_old_1))
        
if __name__ == '__main__':
    unittest.main()
