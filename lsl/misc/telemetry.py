# -*- coding: utf-8 -*-

"""
Basic telemetry client for LSL to help establish usage patterns

.. versionadded:: 1.3.0
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import time
import uuid
import atexit
import socket
import inspect
import warnings
try:
    from urllib2 import urlopen
    from urllib import urlencode
except ImportError:
    from urllib.request import urlopen
    from urllib.parse import urlencode
from threading import RLock
from functools import update_wrapper

from lsl.version import version as lsl_version


__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['is_active', 'enable', 'disable', 'ignore',
           'track_script', 'track_module',
           'track_function', 'track_function_timed',
           'track_method', 'track_method_timed']


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


class _TelemetryClient(object):
    """
    LSL telemetry client to help understand usage of the LSL library.
    """
    _lock = RLock()
    
    _lockout_file = os.path.join(_CACHE_DIR, 'lockout.key')
    
    def __init__(self, key, version=lsl_version, max_entries=50, timeout=1.0):
        # Setup
        self.key = key
        self.version = version
        self.max_entries = max_entries
        self.timeout = timeout
        
        # Session reference
        self._session_start = time.time()
        
        # Telemetry cache
        self._cache = {}
        self._cache_count = 0
        
        # Reporting lockout
        self.active = True
        if os.path.exists(self._lockout_file):
            self.active = False
            
        # Register the "send" method to be called by atexit... at exit
        atexit.register(self.send, True)
        
    def track(self, name, timing=0.0):
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
                
    def send(self, final=False):
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
                                         'session_time': "%.6f" % ((tNow-self._session_start) if final else 0.0,),
                                         'payload'     : payload})
                    try:
                        payload = payload.encode()
                    except AttributeError:
                        pass
                    uh = urlopen('https://fornax.phys.unm.edu/telemetry/log.php', payload, 
                                 timeout=self.timeout)
                    status = uh.read()
                    if status == '':
                        self.clear()
                        success = True
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
    def is_active(self):
        """
        Whether or not the cache is active and sending data back.
        """
        
        return self.active
        
    def enable(self):
        """
        Enable saving data to the telemetry cache.
        """
        
        try:
            self.active = True
            os.unlink(self._lockout_file)
        except OSError:
            pass
            
    def disable(self):
        """
        Disable saving data to the telemetry cache in a persistent way.
        """
        
        self.active = False
        with open(self._lockout_file, 'w') as fh:
            fh.write('disable')
            
    def ignore(self):
        """
        Disable saving data to the telemetry cache in a temporary way.
        """
        
        self.active = False


# Create an instance of the telemetry client to use.
_telemetry_client = _TelemetryClient(_INSTALL_KEY)


# Telemetry control
def is_active():
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


def track_function(user_function):
    """
    Record the use of a function in LSL without execution time information.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    mod = caller.f_globals['__name__']
    fnc = user_function.__name__
    name = mod+'.'+fnc+'()'
    
    def wrapper(*args, **kwds):
        global _telemetry_client
        result =  user_function(*args, **kwds)
        
        _telemetry_client.track(name)
        return result
        
    wrapper.__wrapped__ = user_function
    return update_wrapper(wrapper, user_function)


def track_function_timed(user_function):
    """
    Record the use of a function in LSL with execution time information.
    """
    
    global _telemetry_client
    
    caller = inspect.currentframe().f_back
    mod = caller.f_globals['__name__']
    fnc = user_function.__name__
    name = mod+'.'+fnc+'()'
    
    def wrapper(*args, **kwds):
        global _telemetry_client
        t0 = time.time()
        result = user_function(*args, **kwds)
        t1 = time.time()
        
        _telemetry_client.track(name, t1-t0)
        return result
        
    wrapper.__wrapped__ = user_function
    return update_wrapper(wrapper, user_function)


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
    
    def wrapper(*args, **kwds):
        global _telemetry_client
        result =  user_method(*args, **kwds)
        
        cls = type(args[0]).__name__
        _telemetry_client.track(name % cls)
        return result
        
    wrapper.__wrapped__ = user_method
    return update_wrapper(wrapper, user_method)


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
    
    def wrapper(*args, **kwds):
        global _telemetry_client
        t0 = time.time()
        result =  user_method(*args, **kwds)
        t1 = time.time()
        
        cls = type(args[0]).__name__
        _telemetry_client.track(name % cls, t1-t0)
        return result
        
    wrapper.__wrapped__ = user_method
    return update_wrapper(wrapper, user_method)
