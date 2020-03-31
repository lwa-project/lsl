from __future__ import print_function, division

from astropy.time import Time

__all__ = ['B1950', 'J2000', 'Date']


_DJD_OFFSET = 2415020.0


class _FloatableTime(Time):
    def __float__(self):
        return self.jd - _DJD_OFFSET
        
    def __neg__(self):
        return -float(self)
        
    def __eq__(self, other):
        return float(self) == other
        
    def __gt__(self, other):
        return float(self) > float(other)
        
    def __ge__(self, other):
        return float(self) >= float(other)
        
    def __lt__(self, other):
        return float(self) < float(other)
        
    def __le__(self, other):
        return float(self) <= float(other)
        
    def tuple(self):
        dt = self.datetime
        return (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)


B1950 = _FloatableTime(1950.0, format='byear', scale='utc')
J2000 = _FloatableTime(2000.0, format='jyear', scale='utc')


def Date(value):
    date = None
    if isinstance(value, Time):
        date = value
    if isinstance(value, (int, float)):
        date = _FloatableTime(value+_DJD_OFFSET, format='jd', scale='utc')
    elif isinstance(value, str):
        value = value.replace('/', '-')
        date = _FloatableTime(value, format='iso', scale='utc')
    if date is None:
        raise ValueError("Cannot parse '%s'" % str(value))
    return date
