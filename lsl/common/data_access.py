
import os
import sys
import time
import shutil
import socket
import calendar
import warnings
import contextlib
from datetime import datetime
from urllib.request import Request, urlopen
from urllib.error import HTTPError

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
    
    def __init__(self):
        # Create the cache directory
        try:
            self._data_cache = FileCache(os.path.join(os.path.expanduser('~'), '.lsl', 'data_cache'), max_size=0)
        except OSError:
            self._data_cache = MemoryCache(max_size=0)
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
        
        metaname = filename+'_meta'
        
        local_path = os.path.join(DATA_PATH, relative_url)
        
        received = 0
        if os.path.exists(local_path):
            with open(local_path, 'rb') as uh:
                with self._data_cache.open(filename, 'wb') as fh:
                    while True:
                        data = uh.read(DOWN_CONFIG.get('block_size'))
                        if len(data) == 0:
                            break
                        received += len(data)
                        
                        fh.write(data)
                        
            with self._data_cache.open(metaname, 'w') as fh:
                fh.write(f"created: {time.time():.0f}\n")
                fh.write(f"source: {local_path}\n")
                fh.write(f"source size: {os.path.getsize(local_path)} B\n")
                fh.write(f"source last modified: {os.path.getmtime(local_path):.0f}\n")
                
        # Did we get anything or, at least, enough of something like it looks like 
        # a real file?
        if received < 3:
            ## Fail
            self._data_cache.remove(filename)
            self._data_cache.remove(metaname)
            return False
            
        return True
        
    def _download_mtime(self, relative_url, use_backup=False):
        url = self._BASE_URL+'/'+relative_url
        if use_backup:
            url = self._BACKUP_URL+'/'+relative_url
            
        mtime = 0
        try:
            with urlopen(url, timeout=DOWN_CONFIG.get('timeout')) as uh:
                mtime = uh.headers['Last-Modified']    
                
                mtime = datetime.strptime(mtime, "%a, %d %b %Y %H:%M:%S GMT")
                mtime = calendar.timegm(mtime.timetuple())
        except IOError as e:
            warnings.warn(colorfy("{{%%yellow Error finding modification time of file from %s: %s}}" % (url, str(e))), RuntimeWarning)
        except (socket.timeout, TimeoutError):
            pass
        except HTTPError:
            pass
            
        return mtime
        
    def _download_worker(self, relative_url, filename, use_backup=False):
        """
        Download the URL and save it to a file.
        """
        
        is_interactive = sys.__stdin__.isatty()
        
        metaname = filename+'_meta'
        
        # Attempt to download the data
        url = self._BASE_URL+'/'+relative_url
        if use_backup:
            url = self._BACKUP_URL+'/'+relative_url
        print(f"Downloading {url}")
        try:
            mtime = 0.0
            remote_size = 1
            req = Request(url)
            if os.path.splitext(url)[1] not in ('.Z', '.gz'):
                req.add_header('Accept-encoding', 'gzip')
            with urlopen(req, timeout=DOWN_CONFIG.get('timeout')) as uh:
                remote_size = int(uh.headers["Content-Length"])
                mtime = uh.headers['Last-Modified']
                
                mtime = datetime.strptime(mtime, "%a, %d %b %Y %H:%M:%S GMT")
                mtime = calendar.timegm(mtime.timetuple())
                
                pbar = DownloadBar(max=remote_size)
                received = 0
                with self._data_cache.open(filename, 'wb') as fh:
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
                        
                if is_interactive:
                    sys.stdout.write(pbar.show()+'\n')
                    sys.stdout.flush()
        except IOError as e:
            warnings.warn(colorfy("{{%%yellow Error downloading file from %s: %s}}" % (url, str(e))), RuntimeWarning)
            received = 0
        except (socket.timeout, TimeoutError):
            received = 0
        except HTTPError:
            received = 0
            
        with self._data_cache.open(metaname, 'w') as fh:
            fh.write(f"created: {time.time():.0f}\n")
            fh.write(f"source: {url}\n")
            fh.write(f"source size: {remote_size} B\n")
            fh.write(f"source last modified: {mtime:.0f}\n")
            
        # Did we get anything or, at least, enough of something like it looks like 
        # a real file?
        if received < 3:
            ## Fail
            self._data_cache.remove(filename)
            self._data_cache.remove(metaname)
            return False
            
        return True
        
    def fetch_data_file(self, filename):
        """
        Make sure that a LSL data file exists.  If it does not copy the file
        from the base LSL installation or download it.
        """
        
        if filename not in self._data_cache:
            # No file, go get one.
            status = self._local_copy_worker(filename, filename)
            if not status:
                status = self._download_worker(filename, filename, use_backup=False)
            if not status:
                status = self._download_worker(filename, filename, use_backup=True)
            if not status:
                raise RuntimeError(f"Failed to download '{filename}'")
                
        else:
            # There is a file.  Make sure that it is up to date.
            cache_mtime = self._data_cache.getmtime(filename)
            
            if time.time() - cache_mtime > DOWN_CONFIG.get('refresh_age')*86400:
                ## Yep, looks like it could be in need of a refresh.  See
                ## what our options are.
                local_mtime = self._local_copy_mtime(filename)
                remote_mtime = self._download_mtime(filename, use_backup=False)
                if remote_mtime == 0:
                    remote_mtime = self._download_mtime(filename, use_backup=True)
                    
                if max([local_mtime, remote_mtime]) > cache_mtime:
                    ## Found something newer.  Delete the current version and
                    ## copy/download again.
                    self._data_cache.remove(filename)
                    self.fetch_data_file(filename)
                    
    @contextlib.contextmanager
    def open(self, filename, mode='r'):
        """
        Open a file with the given mode and return the file handle.
        """
        
        # Make sure we have a file
        self.fetch_data_file(filename)
        
        # Open the file
        with self._data_cache.open(filename, mode=mode) as fh:
            yield fh
            
    def remove(self, filename):
        """
        Remove a file from the data cache.
        """
        
        self._data_cache.remove(filename)
        
    def getsize(self, filename):
        """
        Return the size of a file in the cache.
        """
        
        return self._data_cache.getsize(filename)
        
    def getmtime(self, filename):
        """
        Return the last modification time of a file in the cache.
        """
        
        return self._data_cache.getmtime(filename)
        
    def stat(self, filename):
        """
        Return the os.stat_result value for a file in the cache.
        """
        
        return self._data_cache.stat(filename)
    


#: DataAccess instance for accessing LSL software data
DataAccess = _DataAccess()
