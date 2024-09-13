"""
Module to make an ASCII progress bar.
"""

import sys
import copy
import time

from lsl.common.color import colorfy

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.4'
__all__ = ['ProgressBar', 'ProgressBarPlus', 'DownloadBar']


class ProgressBar(object):
    """
    Object to make a ASCII progress bar for use with various long-
    run tasks.
    
    Example Usage:
     >>> import sys
     >>> from progess import ProgressBar
     >>> pb = ProgressBar()
     >>> pb.inc()
     >>> sys.stdout.write(pb.show())
     >>> sys.stdout.flush()
    """
    
    def __init__(self, max=100, span=70, sym='=', print_percent=True, color=None):
        """
        Initialize the ProgressBar class with various parameters:
         * max: maximum count for the progress bar (default: 100)
         * span: width in characters to make the bar (default: 70)
         * sym: character to use in the progress bar (default: '=')
         * print_percent: whether or not to print the percentage in addition to
                          the bar or not (default: True)
         * color: color to use for the progress bar (default: None)
        """
        
        self.amount = 0
        self.max = max
        self.span = span
        self.sym = sym
        self.rotations = ['-', '\\', '|', '/', self.sym]
        self.print_percent = print_percent
        self.color = color
        
    def inc(self, amount=1):
        """
        Increment the progress bar's internal counter by some amount.  The
        default is one.
        """
        
        self.__iadd__(amount)
        
    def dec(self, amount=1):
        """
        Decrement the progress bar's internal counter by some amount.  The
        default is one.
        """
            
        self.__isub__(amount)
        
    def show(self):
        """
        Build a string representation of the progress bar and return it.
        """
        
        if self.print_percent:
            # If we want the percentage also displayed, trim a little 
            # more from the progress bar's wdith
            barSpan = self.span - 9
            nMarks = float(self.amount)/self.max * barSpan
            nMarksFull = int(nMarks)
            if nMarksFull < barSpan:
                partial = nMarks - nMarksFull
                lastMark = self.rotations[int(partial*len(self.rotations))]
            else:
                lastMark = ''
            bar = self.sym * nMarksFull
            bar = bar + lastMark
            bar = bar+(' ' * (barSpan-(nMarksFull+len(lastMark))))
            nte = "%5.1f%%" % (float(self.amount)/self.max*100)
            
            if self.color is None:
                out = "[%s] %s" % (bar, nte)
            else:
                out = colorfy("[{{%%%s %s}}] %s" % (self.color, bar, nte))
        else:
            # Progress bar only
            barSpan = self.span - 2
            nMarks = float(self.amount)/self.max * barSpan
            nMarksFull = int(nMarks)
            if nMarksFull < barSpan:
                partial = nMarks - nMarksFull
                lastMark = self.rotations[int(partial*len(self.rotations))]
            else:
                lastMark = ''
            bar = self.sym * nMarksFull
            bar = bar + lastMark
            bar = bar+(' ' * (barSpan-(nMarksFull+len(lastMark))))
            
            if self.color is None:
                out = "[%s]" % bar
            else:
                out = colorfy("[{{%%%s %s}}]" % (self.color, bar))
            
        return out
        
    def __add__(self, amount):
        """
        Increment the internal counter by a certain amount, return a new
        ProgressBar object.
        """
        
        newBar = copy.deepcopy(self)
        newBar += amount
        return newBar
        
    def __iadd__(self, amount):
        """
        Increment the internal counter by a certain amount.
        """
        
        self.amount += amount
        return self
        
    def __sub__(self, amount):
        """
        Decrement the internal counter by a certain amount, return a new
        ProgressBar object.
        """
        
        newBar = copy.deepcopy(self)
        if newBar.amount >= amount:
            newBar -= amount
        return newBar
        
    def __isub__(self, amount):
        """
        Decrement the internal counter by a certain amount.
        """
        
        if self.amount >= amount:
            self.amount -= amount
        return self
        
    def __str__(self):
        """
        Alternative to self.show().
        """
        
        return self.show()


class ProgressBarPlus(ProgressBar):
    """
    Extended version of the ProgressBar class that has a crude time 
    estimator.
    
    Example Usage:
     >>> import sys
     >>> from progess import ProgressBarPlus
     >>> pb = ProgressBarPlus()
     >>> pb.inc()
     >>> sys.stdout.write(pb.show())
     >>> sys.stdout.flush()
        
    .. note::
        The timing feature is only active when the inc()/dec() functions are called.
        
    .. versionadded:: 0.6.4
    """
    
    t0 = None
    t1 = None
    
    def startTimer(self):
        """
        Initialize the timer.
        """
        
        self.t0 = time.time()
        
    def __iadd__(self, amount):
        """
        Increment the internal counter by a certain amount.
        """
        
        self.amount += amount
        
        if self.t0 is None:
            self.t0 = time.time()
        self.t1 = time.time()
        
        return self
        
    def __isub__(self, amount):
        """
        Decrement the internal counter by a certain amount.
        """
        
        if self.amount >= amount:
            self.amount -= amount
            
        if self.t0 is None:
            self.t0 = time.time()
        self.t1 = time.time()
        
        return self
        
    @staticmethod
    def _pprint(value):
        """
        Nice units for printing.
        """
        
        m, s = value / 60, value % 60
        return '%4im%02is' % (m, s)
        
    def show(self):
        """
        Build a string representation of the progress bar and return it.
        """
        
        if self.t0 is None:
            # Have we started?
            cte = '----m--s'
        elif self.amount == 0:
            # Have we gone far enough to get a "good" estimate?
            cte = '----m--s'
        elif self.amount >= self.max:
            # Are we done?
            cte = self.t1 - self.t0
            cte = self._pprint(cte)
        elif self.t1 - self.t0 < 0.2:
            # Have we running long enough to get a "good" estimate?
            cte = '----m--s'
        else:
            cte = (self.max - self.amount) * (self.t1 - self.t0)/self.amount
            cte = self._pprint(cte)
            
        if self.print_percent:
            # If we want the percentage also displayed, trim a little 
            # more from the progress bar's wdith
            barSpan = self.span - 9
            nMarks = float(self.amount)/self.max * barSpan
            nMarksFull = int(nMarks)
            if nMarksFull < barSpan:
                partial = nMarks - nMarksFull
                lastMark = self.rotations[int(partial*len(self.rotations))]
            else:
                lastMark = ''
            bar = self.sym * nMarksFull
            bar = bar + lastMark
            bar = bar+(' ' * (barSpan-(nMarksFull+len(lastMark))))
            nte = "%5.1f%%" % (float(self.amount)/self.max*100)
            
            if self.color is None:
                out = "[%s] %s %s" % (bar, nte, cte)
            else:
                out = colorfy("[{{%%%s %s}}] %s %s" % (self.color, bar, nte, cte))
        else:
            # Progress bar only
            barSpan = self.span - 2
            nMarks = float(self.amount)/self.max * barSpan
            nMarksFull = min([int(nMarks), barSpan])
            if nMarksFull < barSpan:
                partial = nMarks - nMarksFull
                lastMark = self.rotations[int(partial*len(self.rotations))]
            else:
                lastMark = ''
            bar = self.sym * min([nMarksFull, barSpan])
            bar = bar + lastMark
            bar = bar+(' ' * (barSpan-(nMarksFull+len(lastMark))))
            
            if self.color is None:
                out = "[%s] %s" % (bar, cte)
            else:
                out = colorfy("[{{%%%s %s}}] %s" % (self.color, bar, cte))
            
        return out


class DownloadBar(ProgressBarPlus):
    """
    Modified version of the ProgressBarPlus class that has a crude bandwidth
    estimator.  At the end of the download the total file size (`max`) is 
    displayed instead of the average rate.
    
    Example Usage:
     >>> import sys
     >>> from progess import DownloadBar
     >>> db = DownloadBar()
     >>> db.inc()
     >>> sys.stdout.write(db.show())
     >>> sys.stdout.flush()
        
    .. note::
        The timing feature is only active when the inc()/dec() functions are called.
        
    .. versionadded:: 2.0.2
    """
    
    @staticmethod
    def _pprint(value):
        """
        Nice units for printing.
        """
        
        units = 'B/s'
        if value > 0.9*1024**3:
            value = value/1024.**3
            units = 'GB/s'
        elif value > 0.9*1024**2:
            value = value/1024.**2
            units = 'MB/s'
        elif value > 0.9*1024:
            value = value/1024.
            units = 'kB/s'
        return "%5.1f%4s" % (value, units)
        
    def show(self):
        """
        Build a string representation of the download bar and return it.
        """
        
        if self.t0 is None:
            # Have we started?
            cte = '----- B/s'
        elif self.amount == 0:
            # Have we gone far enough to get a "good" estimate?
            cte = '----- B/s'
        elif self.amount >= self.max:
            # Are we done?
            cte = self._pprint(self.amount)[:-2]
        elif self.t1 - self.t0 < 0.01:
            # Have we running long enough to get a "good" estimate?
            cte = '----- B/s'
        else:
            cte = self.amount / (self.t1 - self.t0)
            cte = self._pprint(cte)
            
        if self.print_percent:
            # If we want the percentage also displayed, trim a little 
            # more from the progress bar's wdith
            barSpan = self.span - 10
            nMarks = float(self.amount)/self.max * barSpan
            nMarksFull = min([int(nMarks), barSpan])
            if nMarksFull < barSpan:
                partial = nMarks - nMarksFull
                lastMark = self.rotations[int(partial*len(self.rotations))]
            else:
                lastMark = ''
            bar = self.sym * nMarksFull
            bar = bar + lastMark
            bar = bar+(' ' * (barSpan-(nMarksFull+len(lastMark))))
            nte = "%5.1f%%" % (float(self.amount)/self.max*100)
            if self.amount > self.max:
                nte = "-----%"
                
            if self.color is None:
                out = "[%s] %s %s" % (bar, nte, cte)
            else:
                out = colorfy("[{{%%%s %s}}] %s %s" % (self.color, bar, nte, cte))
        else:
            # Progress bar only
            barSpan = self.span - 3
            nMarks = float(self.amount)/self.max * barSpan
            nMarksFull = int(nMarks)
            if nMarksFull < barSpan:
                partial = nMarks - nMarksFull
                lastMark = self.rotations[int(partial*len(self.rotations))]
            else:
                lastMark = ''
            bar = self.sym * nMarksFull
            bar = bar + lastMark
            bar = bar+(' ' * (barSpan-(nMarksFull+len(lastMark))))
            
            if self.color is None:
                out = "[%s] %s" % (bar, cte)
            else:
                out = colorfy("[{{%%%s %s}}] %s" % (self.color, bar, cte))
            
        return out
    
