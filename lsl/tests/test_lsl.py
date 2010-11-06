# -*- coding: utf-8 -*-

"""Unit test suite for the lsl package."""

import unittest

from lsl.tests import test_astro
from lsl.tests import test_skymap
from lsl.tests import test_mathutil
from lsl.tests import test_nec_util
from lsl.tests import test_catalog


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
		self.addTest(test_reader.reader_test_suite())
		self.addTest(test_uvUtils.uvUtils_test_suite())


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
