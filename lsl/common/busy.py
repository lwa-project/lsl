"""
Module to make a blinking ASCII busy indicator.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import

import sys
import time
import threading

from lsl.common.color import colorfy

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['BusyIndicator', 'BusyIndicatorPlus']


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
    
    def __init__(self, message='Busy', interval=0.5, dots=3, color=None):
        """
        Initialize the BusyIndicator class with various parameters:
         * message: message to display
         * interval: interval in seconds between message displays
         * dots: number of dots to cycle through at the end of the
                 message
         * color: color to use for the indicator (default: None)
        """
            
        self.message = message
        self.interval = interval
        self.dots = dots
        self.color = color
        
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
            if self.color is None:
                out = "%s%s%s%s\n" % (self.message, 
                                      '.'*self._i, 
                                      'Done' if success else 'Failed', 
                                      ' '*self.dots)
            else:
                out = colorfy("%s{{%%%s %s}}%s%s\n" % (self.message, 
                                                       self.color,
                                                       '.'*self._i, 
                                                       'Done' if success else 'Failed', 
                                                       ' '*self.dots))
            sys.stdout.write(out)
            
            self.alive.clear()
            self.thread.join()
            self.thread = None
            self._i = 0
            
    def _run(self):
        """
        Internal function used by the thread to make/change the displayed text.
        """
        
        while self.alive.isSet():
            if self.color is None:
                out = "%s%s%s\r" % (self.message,
                                    '.'*self._i,
                                    ' '*(self.dots-self._i))
            else:
                out = colorfy("%s{{%%%s %s}}%s\r" % (self.message,
                                                     self.color,
                                                     '.'*self._i,
                                                     ' '*(self.dots-self._i)))
            sys.stdout.write(out)
            sys.stdout.flush()
            
            self._i += 1
            self._i %= (self.dots+1)
            time.sleep(self.interval)


class BusyIndicatorPlus(BusyIndicator):
    _styles = ('boomerang', 'pingpong', 'flow')
    
    def __init__(self, message='Busy', interval=0.1, width=10, style='flow', color=None):
        BusyIndicator.__init__(self, message, interval, 0, color)
        
        if style not in self._styles:
            raise ValueError("Unknown BusyIndicatorPlus style '%s'" % style)
        self.style = style
        self.width = width
        self._dir = 1
        
    def _render_boomerang(self, active=True):
        out = [' ',]*self.width
        out[self._i] = ('|', '/', '-', '\\')[self._i % 4] 
        out = ''.join(out)
        if self.color is not None:
            out = colorfy("{{%%%s %s}}" % (self.color, out))
            
        self._i += self._dir
        if self._i < 0:
            self._i += 1
            self._dir = 1
        if self._i == self.width:
            self._i -= 2
            self._dir = -1
        return out
        
    def _render_pingpong(self, active=True):
        sym = 'o' if active else '.'
        out = ''
        if self._i == 0:
            out += ")%s" % sym
            out += ' '*(self.width-3)
            out += '|'
            self._dir = 1
        elif self._i == self.width-3:
            out += '|'
            out += ' '*(self.width-3)
            out += "%s(" % sym
            self._dir = -1
        else:
            out += '|'
            out += ' '*(self._i)
            out += sym
            out += ' '*(self.width-self._i-3)
            out += '|'
        if self.color is not None:
            out = colorfy(out.replace(sym, "{{%%%s %s}}" % (self.color, sym)))
            
        self._i += self._dir
        return out
        
    def _render_flow(self, active=True):
        out = [' ',]*self.width
        for i in (-2,-1,0):
            out[self._i+i] = '>'
        out = ''.join(out)
        if self.color is not None:
            out = colorfy("{{%%%s %s}}" % (self.color, out))
            
        self._i += 1
        self._i %= self.width
        return out
        
    def stop(self, success=True):
        out = getattr(self, "_render_%s" % self.style)(active=False)
        out = "%s%s%s\n" % (self.message, out, 'Done  ' if success else 'Failed')
        sys.stdout.write(out)
        
        self.alive.clear()
        self.thread.join()
        self.thread = None
        self._i = 0
        self._dir = 1
        
    def _pprint(self):
        t = time.time() - self.t0
        m = int(t)//60
        s = t % 60
        if m == 0:
            out = "%is" % s
        else:
            out = "%im%02is" % (m, s)
        out = "%6s" % out
        return out
        
    def _run(self):
        """
        Internal function used by the thread to make/change the displayed text.
        """
        
        self.t0 = time.time()
        while self.alive.isSet():
            out = getattr(self, "_render_%s" % self.style)(active=True)
            out = "%s%s%s\r" % (self.message, out, self._pprint())
            sys.stdout.write(out)
            sys.stdout.flush()
            
            time.sleep(self.interval)

