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
        self._our_lock = None
        
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
        emit_waiting_warning = 0
        
        pid = os.getpid()
        ident = current_thread().ident
        while not self._locked:
            try:
                # Open the lock file and try to read the thread ID of who owns
                # it
                fh = open(self._lockname, 'a+')
                try:
                    owner_info = fh.read().split()
                    owner_pid, owner_ident = int(owner_info[0], 10), int(owner_info[1], 10)
                except (IndexError, ValueError):
                    owner_pid = owner_ident = 0
                    
                if pid != owner_pid or ident != owner_ident:
                    ## If we don't own it, try to claim an exclusive lock on it.
                    fcntl.flock(fh, fcntl.LOCK_EX|fcntl.LOCK_NB)
                    
                    ## Verify that the lock file and what we just opened are the
                    ## the same file.  If not then give up.
                    fh_stat = os.fstat(fh.fileno())
                    lock_stat  = os.stat(fh.name)
                    if fh_stat.st_ino != lock_stat.st_ino:
                        err = IOError()
                        err.errno = errno.EAGAIN
                        raise err
                        
                    ## Write our thread ID to the file and save the file handle
                    ## to _our_lock so that we know that we need to clean things
                    ## up when we are done.
                    fh.truncate(0)
                    fh.write(f"{pid} {ident}")
                    fh.flush()
                    self._our_lock = fh
                else:
                    ## If we do already own it then there is nothing else to do.
                    fh.close()
                    self._our_lock = None
                self._locked = True
                
            except IOError as e:
                fh.close()
                if e.errno != errno.EAGAIN:
                    raise
                if blocking:
                    tElapsed = time.time() - t0
                    if int(tElapsed/15) != emit_waiting_warning:
                        warnings.warn(f"Waiting {tElapsed:.0f} s to acquire lock on '{self._lockname}'")
                        emit_waiting_warning = int(tElapsed/15)
                        
                    if tElapsed > timeout:
                        break
                else:
                    break
                time.sleep(0.01)
               
        return self._locked
        
    def release(self):
        if self._our_lock is not None:
            # To cleanup we truncate the file and then release the lock.
            self._our_lock.truncate(0)
            self._our_lock.flush()
            fcntl.flock(self._our_lock, fcntl.LOCK_UN)
            self._our_lock.close()
            self._our_lock = None
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
