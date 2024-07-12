""
Unit test for the lsl.common.ndp module.
"""

import warnings
import unittest
import numpy as np

from lsl.common import ndp
from lsl.common import stations


__version__  = "0.1"
__author__    = "Jayce Dowell"


class ndp_bandpass_tests(unittest.TestCase):
    """A unittest.TestCase collection of unit tests for the lsl.common.ndp
    module function for the bandpass."""

    def setUp(self):
        """Turn off all numpy and python warnings."""

        np.seterr(all='ignore')
        warnings.simplefilter('ignore')

    def test_drx_bandpass(self):
        """Test that the DRX bandpass generator actually runs."""

        fnc = ndp.drx_filter(sample_rate=19.6e6, npts=256)
        junk = fnc(1e3)


class ndp_test_suite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the lsl.common.ndp
    module unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(ndp_bandpass_tests)) 


if __name__ == '__main__':
    unittest.main()

