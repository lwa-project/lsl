"""
Common utilities shared by the various metabundle modules.
"""

import re
import shutil
import tarfile
import tempfile

from functools import lru_cache


__version__ = '0.1'
__all__ = ['FILENAME_RE', '_open_tarball', '_get_members', 'managed_mkdtemp']


# Regular expression for figuring out filenames
FILENAME_RE = re.compile(r'(?P<projectID>[a-zA-Z0-9]{1,8})_(?P<sessionID>\d+)(_(?P<obsID>\d+)(_(?P<obsOutcome>\d+))?)?.*\..*')


@lru_cache(maxsize=1)
def _open_tarball(tarname):
    return tarfile.open(tarname, mode='r:*')


@lru_cache(maxsize=1)
def _get_members(tarname):
    tf = _open_tarball(tarname)
    return tf.getmembers()


class managed_mkdtemp(object):
    """
    Wrapper class around tempfile.mkdtemp to enable 'with' statements with 
    automatic cleanup.
    """
    
    def __init__(self, suffix='', prefix='tmp', dir=None):
        self._dir = tempfile.mkdtemp(suffix, prefix, dir)
        
    def __enter__(self):
        return self._dir
        
    def __exit__(self, type, value, tb):
        shutil.rmtree(self._dir, ignore_errors=True)
        
    @property
    def name(self):
        return self._dir
