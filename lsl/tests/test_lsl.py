# -*- coding: utf-8 -*-

"""Unit test suite for the lsl package."""

import unittest

from lsl.tests import test_paths
from lsl.tests import test_astro
from lsl.tests import test_skymap
from lsl.tests import test_mathutil
from lsl.tests import test_nec_util
from lsl.tests import test_catalog
from lsl.tests import test_stations
from lsl.tests import test_robust
from lsl.tests import test_reader
from lsl.tests import test_uvUtils
from lsl.tests import test_fx
from lsl.tests import test_visUtils
from lsl.tests import test_fakedata
from lsl.tests import test_simdp
from lsl.tests import test_simvis
from lsl.tests import test_geodesy
from lsl.tests import test_fitsidi
from lsl.tests import test_tsfits
from lsl.tests import test_vdif
from lsl.tests import test_beamformer


__revision__  = "$Revision: 95 $"
__version__   = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


class lsl_tests(unittest.TestSuite):
	"""A unittest.TestSuite class which contains all of the package unit tests."""

	def __init__(self):
		"""Setup the lsl package unit test suite."""

		unittest.TestSuite.__init__(self)
		
		self.addTest(test_paths.paths_test_suite())
		self.addTest(test_astro.astro_test_suite())
		self.addTest(test_skymap.skymap_test_suite())
		self.addTest(test_mathutil.mathutil_test_suite())
		self.addTest(test_nec_util.nec_util_test_suite())
		self.addTest(test_catalog.catalog_test_suite())
		self.addTest(test_stations.stations_test_suite())
		self.addTest(test_robust.robust_test_suite())
		self.addTest(test_reader.reader_test_suite())
		self.addTest(test_uvUtils.uvUtils_test_suite())
		self.addTest(test_fx.fx_test_suite())
		self.addTest(test_visUtils.visUtils_test_suite())
		self.addTest(test_fakedata.fakedata_test_suite())
		self.addTest(test_simdp.simdp_test_suite())
		self.addTest(test_simvis.simvis_test_suite())
		self.addTest(test_geodesy.geodesy_test_suite())
		self.addTest(test_fitsidi.fitsidi_test_suite())
		self.addTest(test_tsfits.tsfits_test_suite())
		self.addTest(test_vdif.vdif_test_suite())
		self.addTest(test_beamformer.beamformer_test_suite())


def main(opts=None, args=None):
	"""Function to call all of the lsl tests."""

	if opts is not None:
		if opts.verbose:
			level = 2
		else:
			level = 1
	else:
		level = 2
			
	suite = lsl_tests()
	runner = unittest.TextTestRunner(verbosity = level)
	runner.run(suite)


if __name__  == '__main__':
	import optparse
	
	parser = optparse.OptionParser(usage = "python %prog [options]", description = __doc__)
	parser.add_option("-v", "--verbose", action = "store_true", dest = "verbose", default = False,
		help = "extra print output")
	(opts, args) = parser.parse_args()

	main(opts, args)
