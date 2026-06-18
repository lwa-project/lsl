import sys
import logging
import numpy as np
import unittest
import contextlib
from io import StringIO


__version__ = '0.3'
__all__ = ['assert_allclose', 'assert_spatially_close', 'SilentVerbose',
           'quiet_lsl_logging']


def assert_allclose(actual, desired, rtol=1e-01, atol=1e-6, err_msg='', verbose=True):
    """
    Similar to numpy.testing.assert_allclose but compares using an absolute
    tolerance equal to a fraction of the mean magnitude. This ignores large
    relative errors on values with magnitudes much smaller than the mean.
    
    From:
      https://github.com/ledatelescope/bifrost/blob/master/test/test_fft.py
    """
    
    absmean = np.abs(desired).mean()
    np.testing.assert_allclose(actual, desired, rtol=rtol, atol=atol * absmean,
                               err_msg=err_msg, verbose=verbose)


def assert_spatially_close(actual, desired, degrees=False, decimal=7, err_msg='', verbose=True):
    """
    Similar to numpy.testing.assert_almost_equal but compares ra/dec or az/alt
    pairs to see if they are spatially coincident.
    """
    
    lon0, lat0 = actual
    lon1, lat1 = desired
    if degrees:
        lon0 *= np.pi/180
        lat0 *= np.pi/180
        lon1 *= np.pi/180
        lat1 *= np.pi/180
        
    ang = np.sin((lat1-lat0)/2.)**2
    ang += np.cos(lat0)*np.cos(lat1)*np.sin((lon1-lon0)/2.)**2
    ang = min([1.0, ang])
    ang = 2*np.arcsin(ang)
    
    if degrees:
        ang *= 180/np.pi
        
    np.testing.assert_almost_equal(ang, 0.0, decimal=decimal,
                                   err_msg=err_msg, verbose=verbose)


class SilentVerbose(object):
    """
    Class that works as a context manager to capture the output to stdout/stderr
    and re-direct them to StringIO instances.
    """
    
    def __init__(self, stdout=True, stderr=False):
        self.stdout = stdout
        self.stderr = stderr
        
    def __enter__(self):
        if self.stdout:
            self._orig_stdout = sys.stdout
            sys.stdout = StringIO()
        if self.stderr:
            self._orig_stderr = sys.stderr
            sys.stderr = StringIO()
        return self
        
    def __exit__(self, exc_type, exc_value, exc_tb):
        if self.stdout:
            buffer = sys.stdout
            sys.stdout = self._orig_stdout
            buffer.flush()
            buffer.close()
        if self.stderr:
            buffer = sys.stderr
            sys.stderr = self._orig_stderr
            buffer.flush()
            buffer.close()


def quiet_lsl_logging(kls=None, *, level=logging.WARNING):
    """
    Class decorator for unittest.TestCase subclasses that raises the LSL logger
    threshold to `level` (default WARNING) for the duration of each test so
    routine INFO/DEBUG output (e.g. the per-build beam model provenance from
    lsl.sim.vis) does not flood the test output.  The previous level is restored
    after each test via TestCase.addCleanup.  Any existing setUp is preserved.
    """
    
    def decorate(kls):
        orig_setUp = kls.setUp
        
        def setUp(self):
            from lsl.logger import get_log_level, set_log_level
            
            original = get_log_level()
            set_log_level(level)
            self.addCleanup(set_log_level, original)
            orig_setUp(self)
            
        kls.setUp = setUp
        return kls
        
    return decorate if kls is None else decorate(kls)
