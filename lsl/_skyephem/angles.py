"""
Python module providing ephem-like hours, degrees, and separation functions as
well as the hour and degree constants.
"""

import numpy
from functools import total_ordering

from skyfield.units import Angle as SkyAngle
from skyfield.starlib import Star
from astropy.coordinates import Angle as AstroAngle

from lsl._skyephem.config import PYEPHEM_REPR


__all__ = ['hour', 'degree', 'Angle', 'hours', 'degrees', 'separation']


# Fixed quantity to radian values
hour = 1 / 24.
degree = numpy.pi / 180.


@total_ordering
class Angle(SkyAngle):
    """
    Base class for representing angles in a way that behaves like ephem.angle.
    """
    
    def __repr__(self):
        if PYEPHEM_REPR:
            return str(float(self))
        else:
            return SkyAngle.__repr__(self)
            
    def __str__(self):
        if PYEPHEM_REPR:
            output = SkyAngle.__str__(self)
            for u in ('deg ', "' ", 'h ', 'm '):
                output = output.replace(u, ':')
            output = output.replace('"', '')
            output = output.replace('s', '')
            return output
        else:
            return SkyAngle.__str__(self)
            
    def __float__(self):
        return float(self.radians)
        
    def __round__(self, ndigits=0):
        return round(float(self), ndigits)
        
    def __abs__(self):
        return abs(float(self))
        
    def __add__(self, other):
        if isinstance(other, (int, float)):
            return float(self) + other
        elif isinstance(other, SkyAngle):
            return Angle(radians=float(self) + float(other), preference=self.preference)
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return float(self) - other
        elif isinstance(other, SkyAngle):
            return Angle(radians=float(self) - float(other), preference=self.preference)
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return float(self) * other
        elif isinstance(other, SkyAngle):
            return float(self) * float(other)
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return float(self) / other
        elif isinstance(other, SkyAngle):
            return float(self) / float(other)
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __floordiv__(self, other):
        return Angle.__truediv__(self, other)
        
    def __neg__(self):
        return -float(self)
        
    def __eq__(self, other):
        if isinstance(other, (int, float)):
            return float(self) == float(other)
        elif isinstance(other, SkyAngle):
            return float(self) == float(other)
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
        
    def __lt__(self, other):
        if isinstance(other, (int, float)):
            return float(self) < float(other)
        elif isinstance(other, SkyAngle):
            return float(self) < float(other)
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)


def hours(value, wrap=True):
    """
    Build an Angle instance measured in hours.
    """
    
    ang = None
    if isinstance(value, SkyAngle):
        ang = Angle(value, preference='hours')
    elif isinstance(value, (int, float)):
        ang = Angle(radians=value, preference='hours')
    elif isinstance(value, str):
        value = AstroAngle(value, 'hourangle')
        ang = Angle(radians=value.rad, preference='hours')
    if ang is None:
        raise ValueError("Cannot parse '%s'" % str(value))
        
    if wrap:
        value = ang.radians % (2*numpy.pi)
        ang = Angle(radians=value, preference='hours')
    return ang


def degrees(value, wrap360=False, wrap180=False):
    """
    Build an Angle instance measured in degrees.
    """
    
    ang = None
    if isinstance(value, SkyAngle):
        ang = Angle(value, preference='degrees')
    elif isinstance(value, (int, float)):
        ang = Angle(radians=value, preference='degrees')
    elif isinstance(value, str):
        value = AstroAngle(value, 'deg')
        ang = Angle(radians=value.rad, preference='degrees')
    if ang is None:
        raise ValueError("Cannot parse '%s'" % str(value))
        
    if wrap180:
        value = ang.radians % (2*numpy.pi)
        if value > numpy.pi:
            value -= 2*numpy.pi
        ang = Angle(radians=value, preference='degrees')
    elif wrap360:
        value = ang.radians % (2*numpy.pi)
        ang = Angle(radians=value, preference='degrees')
    return ang


def separation(pos1, pos2):
    """
    Return the angular separation between two objects or positions as an
    Angle instance.
    """
    
    c1 = None
    if isinstance(pos1, (tuple, list)):
        ra, dec = pos1
        ra, dec = hours(ra), degrees(dec)
        c1 = Star(ra=ra, dec=dec)
    elif isinstance(pos1, Star):
        c1 = pos1
    if c1 is None:
        raise TypeError("Unexpected type for pos1: %s" % type(pos1))
        
    c2 = None
    if isinstance(pos2, (tuple, list)):
        ra, dec = pos2
        ra, dec = hours(ra), degrees(dec)
        c2 = Star(ra=ra, dec=dec)
    elif isinstance(pos2, Star):
        c2 = pos2
    if c2 is None:
        raise TypeError("Unexpected type for pos2: %s" % type(pos1))
        
    distance = c1.separation_from(c2)
    return degrees(distance)
    
