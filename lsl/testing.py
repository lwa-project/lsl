import sys
import numpy
import unittest
import contextlib
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


__version__ = '0.2'
__all__ = ['assert_allclose', 'assert_spatially_close', 'SilentVerbose']


# Check if unittest has a 'subTest' method.  If it does not add a dummy one
# to help out Python2.7.
if not hasattr(unittest.TestCase, "subTest"):
    @contextlib.contextmanager
    def subTest(self, msg=None, **params):
        yield
        return
        
    unittest.TestCase.subTest = subTest


def assert_allclose(actual, desired, rtol=1e-01, atol=1e-6, err_msg='', verbose=True):
    """
    Similar to numpy.testing.assert_allclose but compares using an absolute
    tolerance equal to a fraction of the mean magnitude. This ignores large
    relative errors on values with magnitudes much smaller than the mean.
    
    From:
      https://github.com/ledatelescope/bifrost/blob/master/test/test_fft.py
    """
    
    absmean = numpy.abs(desired).mean()
    numpy.testing.assert_allclose(actual, desired, rtol=rtol, atol=atol * absmean,
                                  err_msg=err_msg, verbose=verbose)


def assert_spatially_close(actual, desired, degrees=False, decimal=7, err_msg='', verbose=True):
    """
    Similar to numpy.testing.assert_almost_equal but compares ra/dec or az/alt
    pairs to see if they are spatially coincident.
    """
    
    lon0, lat0 = actual
    lon1, lat1 = desired
    if degrees:
        lon0 *= numpy.pi/180
        lat0 *= numpy.pi/180
        lon1 *= numpy.pi/180
        lat1 *= numpy.pi/180
        
    ang = numpy.sin((lat1-lat0)/2.)**2
    ang += numpy.cos(lat0)*numpy.cos(lat1)*numpy.sin((lon1-lon0)/2.)**2
    ang = min([1.0, ang])
    ang = 2*numpy.arcsin(ang)
    
    if degrees:
        ang *= 180/numpy.pi
        
    numpy.testing.assert_almost_equal(ang, 0.0, decimal=decimal,
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
