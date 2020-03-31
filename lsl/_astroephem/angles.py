import numpy

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u


__all__ = ['hours', 'degrees', 'separation']


class _FloatableAngle(Angle):
    def __float__(self):
        return self.to('radian').value


def hours(value):
    ang = None
    if isinstance(value, Angle):
        ang = value
    elif isinstance(value, (int, float)):
        ang = _FloatableAngle(value, u.rad).to(u.hourangle)
    elif isinstance(value, str):
        ang = _FloatableAngle(value, u.hourangle)
    if ang is None:
        raise ValueError("Cannot parse '%s'" % str(value))
    return ang


def degrees(value):
    ang = None
    if isinstance(value, Angle):
        ang = value
    elif isinstance(value, (int, float)):
        ang = _FloatableAngle(value, u.rad).to(u.deg)
    elif isinstance(value, str):
        ang = _FloatableAngle(value, u.deg)
    if ang is None:
        raise ValueError("Cannot parse '%s'" % str(value))
    return ang


def separation(pos1, pos2):
    if isinstance(pos1, SkyCoord):
        c1 = pos1
    else:
        c1 = SkyCoord(*pos0, frame='icrs')
    if isinstance(pos2, SkyCoord):
        c2 = pos2
    else:
        c2 = SkyCoord(*pos1, frame='icrs')
    return c1.separation(c2)

