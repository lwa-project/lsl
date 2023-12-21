
import os
import sys
import time
import shutil
import warnings
import contextlib
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

from lsl.common.paths import DATA as DATA_PATH
from lsl.common.progress import DownloadBar
from lsl.misc.file_cache import FileCache, MemoryCache
from lsl.common.color import colorfy

from lsl.config import LSL_CONFIG
DOWN_CONFIG = LSL_CONFIG.view('download')

__version__ = '0.1'
__all__ = ['DataAccess']


class _DataAccess(object):
    """
    Class to streamline LSL data access whether on disk as part of the LSL
    install or remotely hosted.
    """
    
    # Base URL for remote data downloads
    _BASE_URL = 'https://fornax.phys.unm.edu/lwa/data/lsl'
    
    # Backup URL for remote data downloads
    _BACKUP_URL = 'https://leo.phys.unm.edu/data/lsl'
    
    # Minimum age before checking for an update
    _MIN_AGE_SEC = 14*86400
    
    def __init__(self):
        # Create the cache directory
        try:
            self._CACHE_DIR = FileCache(os.path.join(os.path.expanduser('~'), '.lsl', 'data_cache'), max_size=0)
        except OSError:
            self._CACHE_DIR = MemoryCache(max_size=0)
            warnings.warn(colorfy("{{%yellow Cannot create or write to on-disk data cache, using in-memory data cache}}"), RuntimeWarning)
            
    def _local_copy_mtime(self, relative_url):
        local_path = os.path.join(DATA_PATH, relative_url)
        
        mtime = 0
        if os.path.exists(local_path):
            mtime = os.path.getmtime(local_path)
            
        return mtime
        
    def _local_copy_worker(self, relative_url, filename):
        """
        Look for a local copy of the file in lsl.common.paths.DATA and save it to a
        file.
        """
        
        local_path = os.path.join(DATA_PATH, relative_url)
        
        received = 0
        if os.path.exists(local_path):
            with open(local_path, 'rb') as uh:
                with self._CACHE_DIR.open(filename, 'wb') as fh:
                    while True:
                        data = uh.read(DOWN_CONFIG.get('block_size'))
                        if len(data) == 0:
                            break
                        received += len(data)
                        
                        fh.write(data)
                        
        # Did we get anything or, at least, enough of something like it looks like 
        # a real file?
        if received < 3:
            ## Fail
            self._CACHE_DIR.remove(filename)
            return False
            
        return True
        
    def _download_mtime(self, relative_url, use_backup=False):
        url = self._BASE_URL+'/'+relative_url
        if use_backup:
            url = self._BACKUP_URL+'/'+relative_url
            
        mtime = 0
        try:
            with urlopen(url, timeout=DOWN_CONFIG.get('timeout')):
                mtime = uh.headers['etag']
        except socket.timeout:
            pass
            
        return mtime
        
    def _download_worker(self, relative_url, filename, use_backup=False):
        """
        Download the URL and save it to a file.
        """
        
        is_interactive = sys.__stdin__.isatty()
        
        # Attempt to download the data
        url = self._BASE_URL+'/'+relative_url
        if use_backup:
            url = self._BACKUP_URL+'/'+relative_url
        print("Downloading %s" % url)
        try:
            uh = urlopen(url, timeout=DOWN_CONFIG.get('timeout'))
            remote_size = 1
            try:
                remote_size = int(uh.headers["Content-Length"])
            except AttributeError:
                pass
            try:
                meta = uh.info()
                remote_size = int(meta.getheaders("Content-Length")[0])
            except AttributeError:
                pass
            pbar = DownloadBar(max=remote_size)
            received = 0
            with self._CACHE_DIR.open(filename, 'wb') as fh:
                while True:
                    data = uh.read(DOWN_CONFIG.get('block_size'))
                    if len(data) == 0:
                        break
                    received += len(data)
                    pbar.inc(len(data))
                    
                    fh.write(data)
                    
                if is_interactive:
                    sys.stdout.write(pbar.show()+'\r')
                    sys.stdout.flush()
            uh.close()
            if is_interactive:
                sys.stdout.write(pbar.show()+'\n')
                sys.stdout.flush()
        except IOError as e:
            warnings.warn(colorfy("{{%%yellow Error downloading file from %s: %s}}" % (url, str(e))), RuntimeWarning)
            received = 0
        except socket.timeout:
            received = 0
            
        # Did we get anything or, at least, enough of something like it looks like 
        # a real file?
        if received < 3:
            ## Fail
            self._CACHE_DIR.remove(filename)
            return False
            
        return True
        
    def fetch_data_file(self, filename):
        """
        Make sure that a LSL data file exists.  If it does not copy the file
        from the base LSL installation or download it.
        """
        
        if filename not in self._CACHE_DIR:
            # No file, go get one.
            status = self._local_copy_worker(filename, filename)
            if not status:
                status = self._download_worker(filename, filename, use_backup=False)
            if not status:
                status = self._download_worker(filename, filename, use_backup=True)
            if not status:
                raise RuntimeError("Failed to download '%s'" % filename)
                
        else:
            # There is a file.  Make sure that it is up to date.
            cache_mtime = self._CACHE_DIR.getmtime(filename)
            if time.time() - cache_mtime > self._MIN_AGE_SEC:
                ## Yep, looks like it could be in need of a refresh.  See
                ## what our options are.
                local_mtime = self._local_copy_mtime(filename)
                remote_mtime = self._download_mtime(filename, backup=False)
                if remote_mtime == 0:
                    remote_mtime = self._download_mtime(filename, backup=True)
                    
                if max([local_mtime, remote_mtime]) > cache_mtime:
                    ## Found something newer.  Delete the current version and
                    ## copy/download again.
                    self._CACHE_DIR.remove(filename)
                    self.fetch_data_file(filename)
                    
    @contextlib.contextmanager
    def open(self, filename, mode='r'):
        """
        Open a file with the given mode and return the file handle.
        """
        
        # Make sure we have a file
        self.fetch_data_file(filename)
        
        # Open the file
        with self._CACHE_DIR.open(filename, mode=mode) as fh:
            yield fh


#: DataAccess instance
DataAccess = _DataAccess()
