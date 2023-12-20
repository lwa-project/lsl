
import os
import sys
import shutil
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
    # Base URL for remote data downloads
    _BASE_URL = 'https://fornax.phys.unm.edu/lwa/data/lsl'
    
    def __init__(self):
        # Create the cache directory
        try:
            self._CACHE_DIR = FileCache(os.path.join(os.path.expanduser('~'), '.lsl', 'data_cache'), max_size=0)
        except OSError:
            self._CACHE_DIR = MemoryCache(max_size=0)
            warnings.warn(colorfy("{{%yellow Cannot create or write to on-disk data cache, using in-memory data cache}}"), RuntimeWarning)
            
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
        
    def _download_worker(self, relative_url, filename):
        """
        Download the URL and save it to a file.
        """
        
        is_interactive = sys.__stdin__.isatty()
        
        # Attempt to download the data
        url = self._BASE_URL+'/'+relative_url
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
        
    def establish_data_file(self, filename):
        if filename not in self._CACHE_DIR:
            status = self._local_copy_worker(filename, filename)
            if not status:
                status = self._download_worker(filename, filename)
            if not status:
                raise RuntimeError("Failed to download '%s'" % filename)
                
    @contextlib.contextmanager
    def open(self, filename, mode='r'):
        self.establish_data_file(filename)
        
        fh = self._CACHE_DIR.open(filename, mode=mode)
        try:
            yield fh
        finally:
            fh.close()


DataAccess = _DataAccess()
