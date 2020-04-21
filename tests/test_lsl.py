"""
Unit test suite for the LSL package.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import unittest

from lsl.misc import telemetry
telemetry.ignore()

import test_paths
import test_astro
import test_skymap
import test_transform
import test_mathutil
import test_necutil
import test_catalog
import test_dp
import test_adp
import test_stations
import test_robust
import test_kurtosis
import test_stattests
import test_reader
import test_reader_adp
import test_buffer
import test_ldp
import test_ldp_adp
import test_uvutil
import test_fx
import test_filterbank
import test_fakedata
import test_simdp
import test_simvis
import test_fitsidi
import test_uvfits
import test_sdfits
import test_vdif
import test_measurementset
import test_beamformer
import test_dedispersion
import test_imaging
import test_progress
import test_busy
import test_mcs
import test_mcs_adp
import test_sdf
import test_sdf_adp
import test_meta
import test_meta_adp
import test_parser
import test_scripts
import test_notebooks
import test_idf


__revision__  = "$Rev$"
__version__   = "0.5"
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
        self.addTest(test_transform.transform_test_suite())
        self.addTest(test_mathutil.mathutil_test_suite())
        self.addTest(test_necutil.necutil_test_suite())
        self.addTest(test_catalog.catalog_test_suite())
        self.addTest(test_dp.dp_test_suite())
        self.addTest(test_adp.adp_test_suite())
        self.addTest(test_stations.stations_test_suite())
        self.addTest(test_robust.robust_test_suite())
        self.addTest(test_kurtosis.kurtosis_test_suite())
        self.addTest(test_stattests.stattests_test_suite())
        self.addTest(test_reader.reader_test_suite())
        self.addTest(test_reader_adp.reader_adp_test_suite())
        self.addTest(test_buffer.buffer_test_suite())
        self.addTest(test_ldp.ldp_test_suite())
        self.addTest(test_ldp_adp.ldp_adp_test_suite())
        self.addTest(test_uvutil.uvutil_test_suite())
        self.addTest(test_fx.fx_test_suite())
        self.addTest(test_filterbank.filterbank_test_suite())
        self.addTest(test_fakedata.fakedata_test_suite())
        self.addTest(test_simdp.simdp_test_suite())
        self.addTest(test_simvis.simvis_test_suite())
        self.addTest(test_fitsidi.fitsidi_test_suite())
        self.addTest(test_uvfits.uvfits_test_suite())
        self.addTest(test_sdfits.sdfits_test_suite())
        self.addTest(test_vdif.vdif_test_suite())
        self.addTest(test_measurementset.measurementset_test_suite())
        self.addTest(test_beamformer.beamformer_test_suite())
        self.addTest(test_dedispersion.dedispersion_test_suite())
        self.addTest(test_imaging.imaging_test_suite())
        self.addTest(test_progress.progress_test_suite())
        self.addTest(test_busy.busy_test_suite())
        self.addTest(test_mcs.mcs_test_suite())
        self.addTest(test_mcs_adp.mcs_adp_test_suite())
        self.addTest(test_sdf.sdf_test_suite())
        self.addTest(test_sdf_adp.sdf_adp_test_suite())
        self.addTest(test_meta.metabundle_test_suite())
        self.addTest(test_meta_adp.metabundle_adp_test_suite())
        self.addTest(test_parser.parser_test_suite())
        self.addTest(test_scripts.scripts_test_suite())
        self.addTest(test_notebooks.notebooks_test_suite())
        self.addTest(test_idf.idf_test_suite())


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
