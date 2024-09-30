"""
Module that provides classes for file cache management.

.. versionadded:: 2.1.3
"""

import os
import glob
import time
import contextlib
from threading import Lock
from io import StringIO, BytesIO
from collections import namedtuple
    
from lsl.misc.file_lock import FileLock, MemoryLock


__version__ = '0.1'
__all__ = ['FileCache', 'MemoryCache']


class FileCache(object):
    """
    Class to hold a file cache on disk and enforce a maximum size for that cache
    if needed.
    """
    
    def __init__(self, cache_dir, max_size=0):
        """
        Create a new FileCache instance using the specified cache directory and
        maximum cache size in MB.  If the maximum size is 0, no size limites are
        enforced.
        """
        
        self.cache_dir = os.path.abspath(cache_dir)
        self.max_size = max_size
        self._lock = FileLock(os.path.join(self.cache_dir, 'cache_access'))
        
        # Create the cache directory as needed
        needed = []
        while not os.path.exists(self.cache_dir):
            self.cache_dir, dir = self.cache_dir.rsplit(os.path.sep, 1)
            needed.append(dir)
        needed.reverse()
        
        for dir in needed:
            self.cache_dir = os.path.join(self.cache_dir, dir)
            os.mkdir(self.cache_dir)
            
        # Make sure we can write here
        with self._lock:
            if not self._lock.locked():
                raise RuntimeError(f"Could not acquire {self.cache_dir} lock.  Is another process holding the lock?")
            with open(os.path.join(self.cache_dir, 'write.test'), 'w') as fh:
                fh.write('test')
            os.unlink(os.path.join(self.cache_dir, 'write.test'))
            
    def __contains__(self, filename):
        return os.path.exists(os.path.join(self.cache_dir, filename))
        
    def _cache_management(self):
        try:
            max_size = self.max_size()
        except TypeError:
            max_size = self.max_size
            
        if max_size > 0:
            filenames = glob.glob(os.path.join(self.cache_dir, '*'))
            filenames = list(filter(lambda x:os.path.basename(x) != 'cache_access.lock', filenames))
            sizes = [os.path.getsize(f)/1024.**2 for f in filenames]
            mtimes = [os.path.getmtime(f) for f in filenames]
            
            # Delete until we mee the cache size limit or until there is only one
            # file left (the most recent one)
            while sum(sizes) > max_size and len(filenames) > 1:
                oldest = min(mtimes)
                idx = mtimes.index(oldest)
                try:
                    os.unlink(filenames[idx])
                    del filenames[idx]
                    del sizes[idx]
                    del mtimes[idx]
                except OSError:
                    pass
                    
    @contextlib.contextmanager
    def open(self, filename, mode='r'):
        """
        Open a file from the cache with the given mode and return the file
        handle.  If the file is opened in a writing mode ('w' or 'a') then the
        cache is checked against the 'max_size' value and cleaned before the 
        handle is returned.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        if mode.startswith('w') or mode.startswith('a'):
            parent_dirs = os.path.dirname(filename)
            if parent_dirs != '':
                parent_dirs = os.path.join(self.cache_dir, parent_dirs)
                if not os.path.exists(parent_dirs):
                    os.makedirs(parent_dirs, exist_ok=True)
                    
        fh = open(os.path.join(self.cache_dir, filename), mode)
        if mode.startswith('w') or mode.startswith('a'):
            self._cache_management()
            
        try:
            yield fh
        finally:
            fh.close()
            self._lock.release()
            
    def remove(self, filename):
        """
        Remove the specified filename from the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            os.unlink(os.path.join(self.cache_dir, filename))
        except OSError:
            pass
        finally:
            self._lock.release()
            
    def getsize(self, filename):
        """
        Return the size of a file in the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            size = os.path.getsize(os.path.join(self.cache_dir, filename))
        except OSError as e:
            raise e
        finally:
            self._lock.release()
            
        return size
        
    def getmtime(self, filename):
        """
        Return the last modificaton time of a file in the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            mtime = os.path.getmtime(os.path.join(self.cache_dir, filename))
        except OSError as e:
            raise e
        finally:
            self._lock.release()
            
        return mtime
        
    def stat(self, filename):
        """
        Return the os.stat_result value for a file in the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            fstat = os.stat(os.path.join(self.cache_dir, filename))
        except OSError as e:
            raise e
        finally:
            self._lock.release()
            
        return fstat


#: os.stat_result - like namedtuple for stat-ing a MemoryFile
mf_stat_result = namedtuple('mf_stat_result', ['st_mode', 'st_ino', 'st_dev', 'st_nlink',
                                               'st_uid', 'st_gid', 'st_size',
                                               'st_atime', 'st_mtime', 'st_ctime'])


class MemoryFile(object):
    """
    Wrapper around StringIO/BytesIO for creating an in-memory file-like object.
    Unlike StringIO/BytesIO calling 'close' does not destroy the internal buffer
    making it suitable as a file replacement.
    
    In order to ensure consistent buffer contents a lock is used to make sure
    that only one thread can have the buffer open and accessible at any time.
    """
    
    def __init__(self, name):
        self.name = name
        self.mode = ''
        
        try:
            kls = BytesIO
        except NameError:
            kls = StringIO
        self._buffer = kls()
        self._lock = Lock()
        self._closed = True
        self._is_binary = False
        
        self._ctime = 0.0
        self._atime = 0.0
        self._mtime = 0.0
        
    def __iter__(self):
        contents = self.readline()
        if len(contents) > 0:
            yield contents
            
    @property
    def closed(self):
        """
        Whether or not the buffer is currently closed.
        """
        
        return self._closed
        
    @property
    def size(self):
        """
        Current size of the buffer in bytes.
        """
        
        mark = self._buffer.tell()
        last = self._buffer.seek(0, 2)
        self._buffer.seek(mark, 0)
        return last
        
    @property
    def stat(self):
        """
        The os.stat_result-like (really a mf_stat_result) value for the buffer.
        """
        
        return mf_stat_result(st_mode=0, st_ino=0, st_dev=0, st_nlink=1,
                              st_uid=0, st_gid=0, st_size=self.size,
                              st_atime=int(self._atime), st_mtime=int(self._mtime), st_ctime=int(self._ctime))
    
    def open(self, mode='r'):
        """
        Prepare the buffer and lock it for access. 
        """
        
        if not self._closed:
            raise IOError(f"MemoryFile:{self.name} is already open")
            
        self._lock.acquire(True)
        
        self.mode = mode
        self._atime = time.time()
        if mode.startswith('w'):
            self._ctime = time.time()
        self._closed = False
        self._is_binary = (mode.find('b') != -1)
        
        return self
            
    def seek(self, pos, whence=0):
        """
        Change the position in the buffer and return the new position.  The
        interpretation of what this new position means is determined by the
        'whence' keyword:
         * 0 - Start of stream (the default).  pos should be >= 0;
         * 1 - Current position - pos must be 0;
         * 2 - End of stream - pos must be 0.
        """
        
        if self._closed:
            raise IOError(f"MemoryFile:{self.name} is closed")
            
        self._buffer.seek(pos, whence)
        
    def tell(self):
        """
        Return the current position in the buffer.
        """
        
        if self._closed:
            raise IOError(f"MemoryFile:{self.name} is closed")
            
        return self._buffer.tell()
        
    def truncate(self, pos=None):
        """
        Truncate the size of the buffer to 'pos'.  If 'pos' is None then the
        current position in the buffer is used.  Returns the new absolute
        position in the buffer.
        """
        
        return self._buffer.truncate(pos)
        
    def read(self, n=-1):
        """
        Read at most 'n' bytes from the buffer and return them in as an object
        that is consistent with how the buffer was opened.  If n is <0 then the
        entire remaining contents of the buffer are returned.
        """
        
        if self._closed:
            raise IOError(f"MemoryFile:{self.name} is closed")
            
        contents = self._buffer.read(n)
        if not self._is_binary:
            contents = contents.decode()
        return contents
        
    def readline(self, size=-1):
        """
        Read from the buffer until the next new line character and return the
        line as an object that is conistent with how the buffer was opened.
        """
        
        if self._closed:
            raise IOError(f"MemoryFile:{self.name} is closed")
            
        contents = self._buffer.readline(size)
        if not self._is_binary:
            contents = contents.decode()
        return contents
        
    def write(self, s):
        """
        Write the given data to the buffer and return the number of bytes
        written.
        """
        
        if self._closed:
            raise IOError(f"MemoryFile:{self.name} is closed")
            
        if not self._is_binary:
            try:
                s = s.encode()
            except AttributeError:
                # Already encoded as binary
                pass
        return self._buffer.write(s)
        
    def flush(self):
        """
        Flush the write buffers, if applicable.
        """
        
        if self._closed:
            raise IOError(f"MemoryFile:{self.name} is closed")
            
        self._buffer.flush()
        if self.mode.startswith('w') or self.mode.startswith('a'):
            self._mtime = time.time()
            
    def close(self):
        """
        Close the buffer and release the access lock.
        """
        
        self._buffer.flush()
        self._buffer.seek(0)
        self._closed = True
        if self._lock.locked():
            self._lock.release()


# Lock for accessing the MemoryCache
_MEMORY_CACHE_LOCK = MemoryLock('cache_access')


class MemoryCache(object):
    """
    Class to hold a file cache in memory using StringIO and enforce a maximum
    size for that cache if needed.
    """
    
    def __init__(self, max_size=0):
        """
        Create a new MemoryCache instance using the specified maximum cache
        size.  If the maximum size is 0 in MB, no size limites are enforced.
        """
        
        self.max_size = max_size
        self._cache = {}
        self._lock = _MEMORY_CACHE_LOCK
        
    def __contains__(self, filename):
        return filename in self._cache
            
    def _cache_management(self):
        try:
            max_size = self.max_size()
        except TypeError:
            max_size = self.max_size
            
        if max_size > 0:
            filenames = list(self._cache.keys())
            sizes = [self._cache[filename].size/1024.**2 for filename in filenames]
            mtimes = [self._cache[filename]._mtime for filename in filenames]
            
            # Delete until we mee the cache size limit or until there is only one
            # file left (the most recent one)
            while sum(sizes) > max_size and len(filenames) > 1:
                oldest = min(mtimes)
                idx = mtimes.index(oldest)
                
                filename = filenames[idx]
                del self._cache[filename]
                
                del filenames[idx]
                del sizes[idx]
                del mtimes[idx]
                
    @contextlib.contextmanager
    def open(self, filename, mode='r'):
        """
        Return a StringIO instance corresponding to the provided filename.
        If the call is made with a writing mode ('w' or 'a') then the cache
        is checked against the 'max_size' value and cleaned before the 
        instance is returned.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            fh = self._cache[filename]
        except KeyError:
            self._cache[filename] = MemoryFile(filename)
            fh = self._cache[filename]
        if mode.startswith('w'):
            fh.truncate(0)
        if mode.startswith('w') or mode.startswith('a'):
            self._cache_management()
            
        try:
            yield fh.open(mode=mode)
        finally:
            fh.close()
            self._lock.release()
            
    def remove(self, filename):
        """
        Remove the specified filename from the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            del self._cache[filename]
        except KeyError:
            pass
        finally:
            self._lock.release()
            
    def getsize(self, filename):
        """
        Return the size of a file in the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            size = self._cache[filename].size
        except KeyError:
            raise OSError()
        finally:
            self._lock.release()
            
        return size
        
    def getmtime(self, filename):
        """
        Return the last modification time of a file in the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            mtime = self._cache[filename].stat['st_mtime']
        except KeyError:
            raise OSError()
        finally:
            self._lock.release()
            
        return mtime
        
    def stat(self, filename):
        """
        Return the mf_stat_result value for a file in the cache.
        """
        
        self._lock.acquire()
        assert(self._lock.locked())
        
        try:
            fstat = self._cache[filename].stat
        except KeyError:
            raise OSError()
        finally:
            self._lock.release()
            
        return fstat
