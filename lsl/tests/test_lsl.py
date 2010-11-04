# -*- coding: utf-8 -*-

"""
Unit test suite for lwa_user package.
"""


###########################################################
# $Id: test_lwa_user.py 93 2010-05-21 18:29:27Z dwood $
###########################################################


import unittest

from lsl.tests import test_astro
from lsl.tests import test_skymap
from lsl.tests import test_mathutil
from lsl.tests import test_nec_util
from lsl.tests import test_catalog


__revision__  = "$Revision: 93 $"
__version__   = "dev"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"



class lwa_user_tests(unittest.TestSuite):
    """
    A unittest.TestSuite class which contains all of the package unit tests.
    """

    def __init__(self):
        """
        Setup the lwa_user package unit test suite.
        """

        unittest.TestSuite.__init__(self)
        
        self.addTest(test_astro.astro_test_suite())
        self.addTest(test_skymap.skymap_test_suite())
        self.addTest(test_mathutil.mathutil_test_suite())
        self.addTest(test_nec_util.nec_util_test_suite())
        self.addTest(test_catalog.catalog_test_suite())


if __name__  == '__main__':

    import optparse
    
    parser = optparse.OptionParser(usage = "python %prog [options]", description = __doc__)
    parser.add_option("-v", "--verbose", action = "store_true", dest = "verbose", default = False,
        help = "extra print output")
    (opts, args) = parser.parse_args()
    
    if opts.verbose:
        level = 2
    else:
        level = 1    
          
    suite = lwa_user_tests()
    runner = unittest.TextTestRunner(verbosity = level)
    runner.run(suite)
    
    
       
        
        
        
