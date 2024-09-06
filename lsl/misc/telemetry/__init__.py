"""
Basic telemetry client for LSL to help establish usage patterns

.. versionadded:: 2.0.0
"""

import os
import sys
import time
import uuid
import atexit
import socket
import inspect
import warnings
from urllib.request import urlopen
from urllib.parse import urlencode
from threading import RLock
from functools import wraps

from lsl.version import version as lsl_version

from lsl.config import LSL_CONFIG
TELE_CONFIG = LSL_CONFIG.view('telemetry')

from typing import Callable, Dict, List


__version__ = '0.3'
__all__ = ['is_active', 'enable', 'disable', 'ignore',
           'track_script', 'track_module',
           'track_function', 'track_function_timed',
           'track_method', 'track_method_timed']


try:
    # Create the cache directory
    if not os.path.exists(os.path.join(os.path.expanduser('~'), '.lsl')):
        os.mkdir(os.path.join(os.path.expanduser('~'), '.lsl'))
    _CACHE_DIR = os.path.join(os.path.expanduser('~'), '.lsl', 'telemetry_cache')
    if not os.path.exists(_CACHE_DIR):
        os.mkdir(_CACHE_DIR)

    # Load the install ID key, creating it if it doesn't exist
    _INSTALL_KEY = os.path.join(_CACHE_DIR, 'install.key')
    if not os.path.exists(_INSTALL_KEY):
        with open(_INSTALL_KEY, 'w') as fh:
            fh.write(str(uuid.uuid4()))

    with open(_INSTALL_KEY, 'r') as fh:
        _INSTALL_KEY = fh.read().rstrip()
        
    _IS_READONLY = False
except OSError as e:
    _INSTALL_KEY = ''
    _IS_READONLY = True
    
    warnings.warn("Could not create telemetry cache, telemetry will be disabled for this session: %s" % str(e),
                  RuntimeWarning)


class _TelemetryClient(object):
    """
    LSL telemetry client to help understand usage of the LSL library.
    """
    _lock = RLock()
    
    def __init__(self, key: str, version: str=lsl_version):
        # Setup
        self.key = key
        self.version = version
        self.py_version = "%i.%i" % (sys.version_info.major, sys.version_info.minor)
        self.max_entries = TELE_CONFIG.get('max_entries')
        self.timeout = TELE_CONFIG.get('timeout')
        
        # Session reference
        self._session_start = time.time()
        
        # Telemetry cache
        self._cache: Dict[str,List] = {}
        self._cache_count = 0
        
        if not _IS_READONLY:
            # Reporting lockout
            self.active = TELE_CONFIG.get('enabled')
            
            # Register the "send" method to be called by atexit... at exit
            if not _IS_READONLY:
                atexit.register(self.send, True)
                
        else:
            self.active = False
            
    def track(self, name: str, timing: float=0.0) -> bool:
        """
        Add an entry to the telemetry cache with optional timing information.
        """
        
        if name[:3] != 'lsl' or not self.active:
            return False
            
        with self._lock:
            try:
                self._cache[name][0] += 1
                self._cache[name][1] += (1 if timing > 0 else 0)
                self._cache[name][2] += timing
            except KeyError:
                self._cache[name] = [1, 0, timing]
                self._cache[name][1] += (1 if timing > 0 else 0)
                self._cache_count += 1
                
            if self._cache_count >= self.max_entries:
                self.send()
                
        return True
                
    def send(self, final: bool=False) -> bool:
        """
        Send the current cache of telemetry data back to the maintainers for 
        analysis.
        """
        
        success = False
        with self._lock:
            if self.active and self._cache_count > 0:
                try:
                    tNow = time.time()
                    payload = ';'.join(["%s;%i;%i;%.6f" % (name,
                                                           self._cache[name][0],
                                                           self._cache[name][1],
                                                           self._cache[name][2]) for name in self._cache])
                    payload = urlencode({'timestamp'   : int(tNow),
                                         'key'         : self.key, 
                                         'version'     : self.version,
                                         'py_version'  : self.py_version,
                                         'session_time': "%.6f" % ((tNow-self._session_start) if final else 0.0,),
                                         'payload'     : payload})
                    with urlopen('https://fornax.phys.unm.edu/telemetry/log.php', payload.encode(), timeout=self.timeout) as uh:
                        status = uh.read()
                    if status == '':
                        self.clear()
                        success = True
                    uh.close()
                except Exception as e:
                    warnings.warn("Failed to send telemetry data: %s" % str(e))
            else:
                self.clear()
                
        return success
                
    def clear(self):
        """
        Clear the current telemetry cache.
        """
        
        with self._lock:
            self._cache.clear()
            self._cache_count = 0
            
    @property
    def is_active(self) -> bool:
        """
        Whether or not the cache is active and sending data back.
        """
        
        return self.active
        
    def enable(self):
        """
        Enable saving data to the telemetry cache.
        """
        
        TELE_CONFIG.set('enabled', True)
        if not _IS_READONLY:
            self.active = TELE_CONFIG.get('enabled')
            
    def disable(self):
        """
        Disable saving data to the telemetry cache in a persistent way.
        """
        
        TELE_CONFIG.set('enabled', False)
        self.active = TELE_CONFIG.get('enabled')
            
    def ignore(self):
        """
        Disable saving data to the telemetry cache in a temporary way.
        """
        
        TELE_CONFIG.set_temp('enabled', False)
        self.active = TELE_CONFIG.get('enabled')


# Create an instance of the telemetry client to use.
_telemetry_client = _TelemetryClient(_INSTALL_KEY)


# Telemetry control
def is_active() -> bool:
    """
    Return a boolean of whether or not the LSL telemetry client is active.
    """
    
    global _telemetry_client
    return _telemetry_client.is_active


def enable():
    """
    Enable logging of usage data via the LSL telemetry client.
    """
    
    global _telemetry_client
    _telemetry_client.enable()


def disable():
    """
    Disable logging of usage data via the LSL telemetry client.
    
    .. note::
        This function disables in a global way that persists across
        invocations.
    """
    
    global _telemetry_client
    _telemetry_client.disable()


def ignore():
    """
    Temporarily disable logging of usage data via the LSL telemetry client.
    
    .. note::
        This function only stops logging after it is called and until enable()
        is called again.  For a persistent call use disable().
    """
    
    global _telemetry_client
    _telemetry_client.ignore()


def track_script():
    """
    Record the use of a LSL script.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    name = os.path.basename(caller.f_globals['__file__'])
    _telemetry_client.track('lsl.scripts.'+name)


def track_module():
    """
    Record the import of an LSL module.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    _telemetry_client.track(caller.f_globals['__name__'])


def track_function(user_function: Callable) -> Callable:
    """
    Record the use of a function in LSL without execution time information.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back  # type: ignore
    mod = caller.f_globals['__name__']  # type: ignore
    fnc = user_function.__name__
    name = mod+'.'+fnc+'()'
    
    @wraps(user_function)
    def wrapper(*args, **kwds):
        global _telemetry_client
        result =  user_function(*args, **kwds)
        
        _telemetry_client.track(name)
        return result
        
    return wrapper


def track_function_timed(user_function):
    """
    Record the use of a function in LSL with execution time information.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    mod = caller.f_globals['__name__']
    fnc = user_function.__name__
    name = mod+'.'+fnc+'()'
    
    @wraps(user_function)
    def wrapper(*args, **kwds):
        global _telemetry_client
        t0 = time.time()
        result = user_function(*args, **kwds)
        t1 = time.time()
        
        _telemetry_client.track(name, t1-t0)
        return result
        
    return wrapper


def track_method(user_method):
    """
    Record the use of a method in LSL with execution time information.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    mod = caller.f_globals['__name__']
    cls = None
    fnc = user_method.__name__
    name = mod+'.'+'%s'+'.'+fnc+'()'
    
    @wraps(user_method)
    def wrapper(*args, **kwds):
        global _telemetry_client
        result =  user_method(*args, **kwds)
        
        cls = type(args[0]).__name__
        _telemetry_client.track(name % cls)
        return result
        
    return wrapper


def track_method_timed(user_method):
    """
    Record the use of a method in LSL with execution time information.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    mod = caller.f_globals['__name__']
    cls = None
    fnc = user_method.__name__
    name = mod+'.'+'%s'+'.'+fnc+'()'
    
    @wraps(user_method)
    def wrapper(*args, **kwds):
        global _telemetry_client
        t0 = time.time()
        result =  user_method(*args, **kwds)
        t1 = time.time()
        
        cls = type(args[0]).__name__
        _telemetry_client.track(name % cls, t1-t0)
        return result
        
    return wrapper
