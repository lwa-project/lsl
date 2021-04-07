"""
Python module providing ephem-like hours, degrees, and separation functions as
well as the hour and degree constants.
"""

import numpy

from astropy.coordinates import Angle as AstroAngle, SkyCoord
from astropy import units as u

from functools import total_ordering

from lsl._astroephem.config import PYEPHEM_REPR


__all__ = ['hour', 'degree', 'Angle', 'hours', 'degrees', 'separation']


# Fixed quantity to radian values
hour = 1 / 24.
degree = numpy.pi / 180.


@total_ordering
class Angle(AstroAngle):
    """
    Base class for representing angles in a way that behaves like ephem.angle.
    """
    
    def __repr__(self):
        if PYEPHEM_REPR:
            return str(float(self))
        else:
            return AstroAngle.__repr__(self)
            
    def __str__(self):
        if PYEPHEM_REPR:
            output = AstroAngle.__str__(self)
            for u in ('d', 'm', 'h'):
                output = output.replace(u, ':')
            output = output.replace('s', '')
            return output
        else:
            return AstroAngle.__str__(self)
            
    def __float__(self):
        return float(self.to('radian').value)
        
    def __round__(self, ndigits=0):
        return round(float(self), ndigits)
        
    def __add__(self, other):
        if isinstance(other, (int, float)):
            return float(self) + other
        else:
            return AstroAngle.__add__(self, other)
            
    def __sub__(self, other):
        if isinstance(other, (int, float)):
            return float(self) - other
        else:
            return AstroAngle.__sub__(self, other)
            
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return float(self) * other
        else:
            return AstroAngle.__mul__(self, other)
            
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return float(self) / other
        else:
            return AstroAngle.__truediv__(self, other)
            
    def __floordiv__(self, other):
        return Angle.__truediv__(self, other)
        
    def __neg__(self):
        return -float(self)
        
    def __eq__(self, other):
        if isinstance(other, (int, float)):
            return float(self) == float(other)
        else:
            return AstroAngle.__eq__(self, other)
        
    def __lt__(self, other):
        if isinstance(other, (int, float)):
            return float(self) < float(other)
        else:
            return AstroAngle.__lt__(self, other)


def hours(value, wrap=False):
    """
    Build an Angle instance measured in hours.
    """
    
    ang = None
    if isinstance(value, AstroAngle):
        ang = Angle(value).to(u.hourangle)
    elif isinstance(value, (int, float)):
        ang = Angle(value, u.rad).to(u.hourangle)
    elif isinstance(value, str):
        ang = Angle(value, u.hourangle)
    if ang is None:
        raise ValueError("Cannot parse '%s'" % str(value))
        
    if wrap:
        ang.wrap_at(24*u.hourangle, inplace=True)
    return ang


def degrees(value, wrap360=False, wrap180=False):
    """
    Build an Angle instance measured in degrees.
    """
    
    ang = None
    if isinstance(value, AstroAngle):
        ang = Angle(value).to(u.deg)
    elif isinstance(value, (int, float)):
        ang = Angle(value, u.rad).to(u.deg)
    elif isinstance(value, str):
        ang = Angle(value, u.deg)
    if ang is None:
        raise ValueError("Cannot parse '%s'" % str(value))
        
    if wrap360:
        ang.wrap_at(360*u.deg, inplace=True)
    elif wrap180:
        ang.wrap_at(180*u.deg, inplace=True)
    return ang


def separation(pos1, pos2):
    """
    Return the angular separation between two objects or positions as an
    Angle instance.
    """
    
    if isinstance(pos1, SkyCoord):
        c1 = pos1
    else:
        pos1_1, pos1_2 = pos1
        if not isinstance(pos1_1, u.quantity.Quantity):
            pos1_1 = pos1_1*u.radian
        if not isinstance(pos1_2, u.quantity.Quantity):
            pos1_2 = pos1_2*u.radian
        c1 = SkyCoord(pos1_1, pos1_2, frame='icrs')
    if isinstance(pos2, SkyCoord):
        c2 = pos2
    else:
        pos2_1, pos2_2 = pos2
        if not isinstance(pos2_1, u.quantity.Quantity):
            pos2_1 = pos2_1*u.radian
        if not isinstance(pos2_2, u.quantity.Quantity):
            pos2_2 = pos2_2*u.radian
        c2 = SkyCoord(pos2_1, pos2_2, frame='icrs')
    return Angle(c1.separation(c2))
