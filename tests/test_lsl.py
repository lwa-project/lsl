# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import absolute_import

"""Unit test suite for the lsl package."""

import unittest

from . import test_paths
from . import test_astro
from . import test_skymap
from . import test_transform
from . import test_mathutil
from . import test_necutil
from . import test_catalog
from . import test_dp
from . import test_stations
from . import test_robust
from . import test_stattests
from . import test_reader
from . import test_reader_adp
from . import test_buffer
from . import test_ldp
from . import test_ldp_adp
from . import test_uvUtils
from . import test_fx
from . import test_filterbank
from . import test_fakedata
from . import test_simdp
from . import test_simvis
from . import test_fitsidi
from . import test_uvfits
from . import test_sdfits
from . import test_vdif
from . import test_measurementset
from . import test_beamformer
from . import test_imaging
from . import test_progress
from . import test_busy
from . import test_mcs
from . import test_mcs_adp
from . import test_sdf
from . import test_sdf_adp
from . import test_meta
from . import test_meta_adp


__revision__  = "$Rev$"
__version__   = "0.4"
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
        self.addTest(test_stations.stations_test_suite())
        self.addTest(test_robust.robust_test_suite())
        self.addTest(test_stattests.stattests_test_suite())
        self.addTest(test_reader.reader_test_suite())
        self.addTest(test_reader_adp.reader_adp_test_suite())
        self.addTest(test_buffer.buffer_test_suite())
        self.addTest(test_ldp.ldp_test_suite())
        self.addTest(test_ldp_adp.ldp_adp_test_suite())
        self.addTest(test_uvUtils.uvUtils_test_suite())
        self.addTest(test_fx.fx_test_suite())
        self.addTest(test_filterbank.filterbank_test_suite())
        self.addTest(test_fakedata.fakedata_test_suite())
        self.addTest(test_simdp.simdp_test_suite())
        self.addTest(test_simvis.simvis_test_suite())
        self.addTest(test_fitsidi.fitsidi_test_suite())
        self.addTest(test_uvfits.uvfits_test_suite())
        self.addTest(test_sdfits.sdfits_test_suite())
        self.addTest(test_vdif.vdif_test_suite())
        self.addTest(test_measurementset.measurement_test_suite())
        self.addTest(test_beamformer.beamformer_test_suite())
        self.addTest(test_imaging.imaging_test_suite())
        self.addTest(test_progress.progress_test_suite())
        self.addTest(test_busy.busy_test_suite())
        self.addTest(test_mcs.mcs_test_suite())
        self.addTest(test_mcs_adp.mcs_adp_test_suite())
        self.addTest(test_sdf.sdf_test_suite())
        self.addTest(test_sdf_adp.sdf_adp_test_suite())
        self.addTest(test_meta.metabundle_test_suite())
        self.addTest(test_meta_adp.metabundle_adp_test_suite())


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
