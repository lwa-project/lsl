"""
Unit test for the NDP portion of the lsl.common.metabundle module.
"""

import os
import unittest

from lsl.common.mcs import CommandID, ObservingMode
from lsl.common import metabundle, metabundle, stations
from lsl.common.paths import DATA_BUILD

run_gdbm_tests = False
try:
    from dbm import gnu
    run_gdbm_tests = True
except ImportError:
    pass


__version__  = "0.2"
__author__    = "Jayce Dowell"

mdbFile = os.path.join(os.path.dirname(__file__), 'data', 'metadata.tgz')


class metabundle_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.metabundle
    module."""
    
    def test_ss(self):
        """Test the session specification utilties."""
        
        ses = metabundle.get_session_spec(mdbFile)
        obs = metabundle.get_observation_spec(mdbFile)
        
        # Check session start time
        self.assertEqual(ses['mjd'], 60902)
        self.assertEqual(ses['mpm'], 71499000-5000)
        
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
        self.assertEqual(obs2['mjd'], 60902)
        self.assertEqual(obs2['mpm'], 71499000)
        
    def test_cs(self):
        """Test the command script utilities."""
        
        cmnds = metabundle.get_command_script(mdbFile)
        
        # Check number of command
        self.assertEqual(len(cmnds), 10)
        
        # Check the first and last commands
        self.assertEqual(cmnds[ 0]['command_id'], CommandID.NUL)
        self.assertEqual(cmnds[-2]['command_id'], CommandID.STP)
        self.assertEqual(cmnds[-1]['command_id'], CommandID.ESN)
        
        # Check the counds of DP BAM commands
        nBAM = 0
        for cmnd in cmnds:
            if cmnd['command_id'] == CommandID.BAM:
                nBAM += 1
        self.assertEqual(nBAM, 1)
        
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
        
        station = stations.lwana
        sm = metabundle.get_sdm(mdbFile)
        newAnts = sm.update_antennas(station.antennas)
        
    def test_metadata(self):
        """Test the observation metadata utility."""
        
        fileInfo = metabundle.get_session_metadata(mdbFile)
        self.assertEqual(len(fileInfo.keys()), 1)
        
        # File tag
        self.assertEqual(fileInfo[1]['tag'], '060902_000067516')
        
        # DRSU barcode
        self.assertEqual(fileInfo[1]['barcode'], 'UNK')
        
    @unittest.skipUnless(run_gdbm_tests, "requires the 'dbm.gnu' module")
    def test_aspconfig(self):
        """Test retrieving the ASP configuration."""
        
        # Beginning config.
        aspConfig = metabundle.get_asp_configuration_summary(mdbFile, which='beginning')
        self.assertEqual(aspConfig['asp_filter'],      3)
        self.assertEqual(aspConfig['asp_atten_1'],     0)
        self.assertEqual(aspConfig['asp_atten_2'],     0)
        self.assertEqual(aspConfig['asp_atten_split'], 0)
        
        # End config.
        aspConfig = metabundle.get_asp_configuration_summary(mdbFile, which='End')
        self.assertEqual(aspConfig['asp_filter'],      3)
        self.assertEqual(aspConfig['asp_atten_1'],     0)
        self.assertEqual(aspConfig['asp_atten_2'],     0)
        self.assertEqual(aspConfig['asp_atten_split'], 0)
        
        # Unknown code
        self.assertRaises(ValueError, metabundle.get_asp_configuration_summary, mdbFile, 'middle')
        
        # Not a summary
        aspConfig = metabundle.get_asp_configuration(mdbFile, which='End')
        self.assertEqual(aspConfig['asp_filter'][0],      3)
        self.assertEqual(aspConfig['asp_atten_1'][0],     0)
        self.assertEqual(aspConfig['asp_atten_2'][0],     0)
        self.assertEqual(aspConfig['asp_atten_split'][0], 0)
        
    def test_is_valid(self):
        """Test whether or not is_valid works."""
        
        self.assertTrue(metabundle.is_valid(mdbFile))


class metabundle_ndp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.metabundle
    module unit tests."""
    
    def __init__(self):
        unittest.TestSuite.__init__(self)
        
        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(metabundle_tests))
        
if __name__ == '__main__':
    unittest.main()
