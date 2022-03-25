"""
Basic file-based lock for filesystem operations.

.. versionadded:: 2.1.3
"""

import os
import time
import errno
import fcntl
import warnings
from threading import current_thread, RLock


__version__ = '0.1'
__all__ = ['FileLock', 'MemoryLock']


class FileLock(object):
    """
    Create a file-based lock that can be used to control access to a particular
    file.   For files that reside on read-only filesystems or have permission
    problems this class initially generates a warning and then silently allows
    all locks/unlocks.
    """
    
    def __init__(self, filename):
        self._lockname = filename+'.lock'
        self._locked = False
        self._our_lock = False
        
    def __del__(self):
        if self._locked:
            self.release()
            
    def __enter__(self):
        self.acquire()
        return self
        
    def __exit__(self, type, value, tb):
        self.release()
        
    def locked(self):
        return self._locked
        
    def acquire(self, blocking=True, timeout=120):
        t0 = time.time()
        
        ident = current_thread().ident
        while not self._locked:
            try:
                if os.path.exists(self._lockname):
                    with open(self._lockname, 'r') as fh:
                        try:
                            owner_ident = int(fh.read(), 10)
                        except ValueError:
                            owner_ident = 0
                else:
                    owner_ident = 0
                    
                if ident != owner_ident:
                    if os.path.exists(self._lockname):
                        err = IOError()
                        err.errno = errno.EAGAIN
                        raise err
                        
                    with open(self._lockname, 'a+') as fh:
                        fcntl.flock(fh, fcntl.LOCK_EX|fcntl.LOCK_NB)
                        fh.truncate(0)
                        fh.write("%i" % ident)
                        fh.flush()
                    self._our_lock = True
                else:
                    self._our_lock = False
                self._locked = True
                
            except IOError as e:
                if e.errno != errno.EAGAIN:
                    raise
                if blocking:
                    if time.time()-t0 > timeout:
                        break
                else:
                    break
                time.sleep(0.01)
                
        return self._locked
        
    def release(self):
        if self._our_lock:
            os.unlink(self._lockname)
            self._our_lock = False
        self._locked = False
        return True


class MemoryLock(object):
    def __init__(self, filename):
        self._lock = RLock()
        self._locked = False
        
    def __enter__(self):
        self.acquire()
        return self
        
    def __exit__(self, type, value, tb):
        self.release()
        
    def locked(self):
        return self._locked
        
    def acquire(self, blocking=True, timeout=120):
        try:
            self._locked = self._lock.acquire(blocking, timeout)
        except TypeError:
            self._locked = self._lock.acquire(blocking)
        return self._locked
        
    def release(self):
        self._lock.release()
        self._locked = False
