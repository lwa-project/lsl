# -*- coding: utf-8 -*-

"""Unit test for lsl.common.metabundle"""

import os
import unittest

from lsl.common import metabundle
from lsl.common.paths import dataBuild as dataPath


__revision__ = "$ Revision: 2 $"
__version__  = "0.1"
__author__    = "Jayce Dowell"

mdbFile = os.path.join(dataPath, 'tests', 'metadata.tgz')


class metabundle_tests(unittest.TestCase):
	"""A unittest.TestCase collection of unit tests for the lsl.common.metabundle
	module."""
	
	def test_ss(self):
		"""Test the session specification utilties."""
		
		ses = metabundle.getSessionSpec(mdbFile)
		obs = metabundle.getObservationSpec(mdbFile)
		
		# Check session start time
		self.assertEqual(ses['MJD'], 56013)
		self.assertEqual(ses['MPM'], 25855000)
		
		# Check the duration
		self.assertEqual(ses['Dur'], obs[0]['Dur'] + 10000)
		
		# Check the number of observations
		self.assertEqual(ses['nObs'], len(obs))
	
	def test_os(self):
		"""Test the observation specification utilities."""
		
		obs1 = metabundle.getObservationSpec(mdbFile)
		obs2 = metabundle.getObservationSpec(mdbFile, selectObs=1)
		
		# Check if the right observation is returned
		self.assertEqual(obs1[0], obs2)
		
		# Check the mode
		self.assertEqual(obs2['Mode'], 1)
		
		# Check the time
		self.assertEqual(obs2['MJD'], 56013)
		self.assertEqual(obs2['MPM'], 25860000)
		
	
	def test_cs(self):
		"""Test the command script utilities."""
		
		cmnds = metabundle.getCommandScript(mdbFile)
		
		# Check number of command
		self.assertEqual(len(cmnds), 491)
		
		# Check the first and last commands
		self.assertEqual(cmnds[ 0]['commandID'], 'NUL')
		self.assertEqual(cmnds[-2]['commandID'], 'OBE')
		self.assertEqual(cmnds[-1]['commandID'], 'ESN')
		
		# Check the counds of DP BAM commands
		nBAM = 0
		for cmnd in cmnds:
			if cmnd['commandID'] == 'BAM':
				nBAM += 1
		self.assertEqual(nBAM, 484)
		
	def test_sm(self):
		"""Test the session metadata utilties."""
		
		sm = metabundle.getSessionMetaData(mdbFile)
		
		# Make sure all of the observations are done
		self.assertEqual(len(sm.keys()), 1)
		
	def test_sdm(self):
		"""Test the station dynamic MIB utilties."""
		
		sm = metabundle.getSDM(mdbFile)

    
class metabundle_test_suite(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the lsl.common.metabundle
	module unit tests."""
	
	def __init__(self):
		unittest.TestSuite.__init__(self)
		
		loader = unittest.TestLoader()
		self.addTests(loader.loadTestsFromTestCase(metabundle_tests))        
        
        
if __name__ == '__main__':
	unittest.main()
