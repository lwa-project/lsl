"""
Module to make a blinking ASCII busy indicator.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import

import sys
import time
import threading

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.1'
__all__ = ['BusyIndicator',]


class BusyIndicator(object):
    """
    Object to make a ASCII busy indicator for use with various long-
    run tasks that don't have a way of calculating how long they will
    take.

    Example Usage:
        >>> from busy import BusyIndicator
        >>> bi = BusyIndicator()
        >>> bi.start()
        >>> longRunningTask()
        >>> bi.stop()
    """
    
    def __init__(self, message='Busy', interval=0.5, dots=3):
        """
        Initialize the BusyIndicator class with various parameters:
        * message: message to display
        * interval: interval in seconds between message displays
        * dots: number of dots to cycle through at the end of the
            message
        """
            
        self.message = message
        self.interval = interval
        self.dots = dots
        
        self.thread = None
        self.alive = threading.Event()
        self._i = 0
        
    def __enter__(self):
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_value, exc_tb):
        success = True
        if exc_type is not None:
            success = False
        self.stop(success=success)
        
    def start(self):
        """
        Start the indicator running.
        """
        
        if self.thread is not None:
            self.stop()
            
        self.thread = threading.Thread(target=self._run, name='indicator')
        self.thread.setDaemon(1)
        self.alive.set()
        self.thread.start()
        
    def stop(self, success=True):
        """
        Stop the indicator and display a 'Done'  or 'Failed' message depending on
        whether or not the 'success' keyword is True.
        
        .. note::
            This can take up to one BusyIndicator.interval to complete.
        """
        
        if self.thread is not None:
            sys.stdout.write('%s%s%s%s\n' % (self.message, 
                                             '.'*self._i, 
                                             'Done' if success else 'Failed', 
                                             ' '*self.dots))
            
            self.alive.clear()
            self.thread.join()
            self.thread = None
            self._i = 0
            
    def _run(self):
        """
        Internal function used by the thread to make/change the displayed text.
        """
        
        while self.alive.isSet():
            sys.stdout.write('%s%s%s\r' % (self.message, '.'*self._i, ' '*(self.dots-self._i)))
            sys.stdout.flush()
            
            self._i += 1
            self._i %= (self.dots+1)
            time.sleep(self.interval)
