"""
Python module providing ephem-like conversions between coordinate systems.
"""

from __future__ import print_function, division

import numpy
from scipy.optimize import minimize

from skyfield import api
from skyfield.starlib import Star
from skyfield.positionlib import Astrometric, position_of_radec
ts = api.load.timescale()

from lsl._skyephem.angles import hours, degrees
from lsl._skyephem.dates import Date, B1950, J2000
from lsl._skyephem.cache import load_planetary_ephemeris
from lsl._skyephem.bodies import Body, FixedBody


__all__ = ['Coordinate', 'Equatorial', 'Ecliptic', 'Galactic', 'separation']


_solar_system = load_planetary_ephemeris()


def _equ_from_ecl(lon, lat, epoch):
    best = 1e9
    for r in numpy.linspace(0, 2*numpy.pi, 25):
        r = hours(r)
        for d in numpy.linspace(-numpy.pi/2, numpy.pi/2, 13):
            d = degrees(d)
            s = Star(r, d, epoch=epoch)
            p = _solar_system['earth'].at(epoch).observe(s).apparent().ecliptic_latlon()
            o = (lon.radians-p[1].radians)**2 + (lat.radians-p[0].radians)**2
            if o < best:
                best = o
                ra, dec = r, d
                
    def _search(x):
         p = _solar_system['earth'].at(epoch).observe(Star(hours(x[0]), degrees(x[1]), epoch=epoch)).apparent().ecliptic_latlon()
         return (lon.radians-p[1].radians)**2 + (lat.radians-p[0].radians)**2
    
    x0 = [ra.radians, dec.radians]
    sol = minimize(_search, x0, tol=(0.1/3600*numpy.pi/180)**2)
    return (hours(sol.x[0]), degrees(sol.x[1]))


def _equ_from_gal(lon, lat, epoch):
    best = 1e9
    for r in numpy.linspace(0, 2*numpy.pi, 49):
        r = hours(r)
        for d in numpy.linspace(-numpy.pi/2, numpy.pi/2, 25):
            d = degrees(d)
            s = Star(r, d, epoch=epoch)
            p = _solar_system['earth'].at(J2000).observe(s).apparent().galactic_latlon()
            o = (lon.radians-p[1].radians)**2 + (lat.radians-p[0].radians)**2
            if o < best:
                best = o
                ra, dec = r, d
                
    def _search(x):
         p = _solar_system['earth'].at(J2000).observe(Star(hours(x[0]), degrees(x[1]), epoch=epoch)).apparent().galactic_latlon()
         return (lon.radians-p[1].radians)**2 + (lat.radians-p[0].radians)**2

    x0 = [ra.radians, dec.radians]
    sol = minimize(_search, x0, tol=(0.01/3600*numpy.pi/180)**2)
    return (hours(sol.x[0]), degrees(sol.x[1]))


class Coordinate(object):
    """
    Base class to facilitate the coordinate transforms in a way that behaves
    like ephem.Coordinate.
    """
    
    def __init__(self, c1_or_class, c2=None, epoch=J2000):
        epoch = Date(epoch)
        
        if isinstance(c1_or_class, Coordinate):
            ra, dec, _ = c1_or_class._apparent_position(epoch).radec('date')
        elif isinstance(c1_or_class, Body):
            ra, dec = c1_or_class.a_ra, c1_or_class.a_dec
            epoch = c1_or_class.a_epoch
        else:
            ra, dec = hours(c1_or_class), degrees(c2)
        self._coord = Star(ra, dec, epoch=epoch)
        self._epoch = epoch
        
    def _apparent_position(self, epoch=None):
        if epoch is None:
            epoch = self._epoch
        return _solar_system['earth'].at(epoch).observe(self._coord).apparent()
        
    @classmethod
    def from_radec(cls, ra, dec):
        """
        Set the coordinates for the instance using the provided ICRS right
        ascension and declination.
        """
        
        return cls(ra, dec)
        
    def to_radec(self):
        """
        Return the ICRS right assencsion and declination for the coordinate.
        """
        
        ra, dec, _ = self._apparent_position().radec()
        return (hours(ra), degrees(dec))
               

class Equatorial(Coordinate):
    """
    Class for representing an Equatorial position on the sky.  The default frame
    is ICRS.  Setting an epoch of J2000 sets the frame to FK5 and setting it to
    B1950 set the frame to FK4.
    """
    
    @property
    def ra(self):
        """
        Right ascension in the current frame.
        """
        
        return self.get()[0]
        
    @property
    def dec(self):
        """
        Declination in the current frame.
        """
        
        return self.get()[1]
        
    def get(self):
        """
        Get the right ascension and declination as two-element tuple of Angle
        instances.
        """
        
        ra, dec, _ = self._apparent_position().radec()
        return (hours(ra), degrees(dec))
        
    def set(self, ra, dec):
        """
        Set the right ascension and declination coordinates for the instance.
        """
        
        ra, dec = hours(c1_or_class), degrees(c2)
        self._coord = Star(ra, dec, epoch=self._epoch)
        

class Ecliptic(Coordinate):
    """
    Class for representing an Ecliptic position on the sky in the geocentric true
    ecliptic frame.
    """
    
    def __init__(self, c1_or_class, c2=None, epoch=J2000):
        epoch = Date(epoch)
        
        try:
            ra, dec = c1_or_class.to_radec()
            epoch = c1_or_class._epoch
        except AttributeError:
            ra, dec = _equ_from_ecl(degrees(c1_or_class), degrees(c2), epoch)
        Coordinate.__init__(self, ra, dec, epoch=epoch)
        
    @property
    def lon(self):
        """
        Ecliptic longitude.
        """
        
        return self.get()[0]
        
    @property
    def lat(self):
        """
        Ecliptic latitude.
        """
        
        return self.get()[1]
        
    def get(self):
        """
        Get the longitude and latitude as two-element tuple of Angle instances.
        """
        
        lat, lon, _ = self._apparent_position().ecliptic_latlon('date')
        lat, lon = degrees(lat), degrees(lon)
        return lon, lat
        
    def set(self, lon, lat):
        """
        Set the longitude and latitude coordinates for the instance.
        """
        
        lon, lat = degrees(lon), degrees(lat)
        ra, dec = _equ_from_ecl(lon, lat, self._epoch)
        self.set(ra, dec)
        

class Galactic(Coordinate):
    """
    Class for representing a Galactic position on the sky in the geocentric true
    ecliptic frame.
    """
    
    def __init__(self, c1_or_class, c2=None, epoch=J2000):
        epoch = Date(epoch)
        
        try:
            ra, dec = c1_or_class.to_radec()
            epoch = c1_or_class._epoch
        except AttributeError:
            ra, dec = _equ_from_gal(degrees(c1_or_class), degrees(c2), epoch)
        Coordinate.__init__(self, ra, dec, epoch=epoch)
        
    @property
    def lon(self):
        """
        Galactic longitude.
        """
        
        return self.get()[0]
        
    @property
    def lat(self):
        """
        Galactic latitude.
        """
        
        return self.get()[1]
        
    def get(self):
        """
        Get the longitude and latitude as two-element tuple of Angle instances.
        """
        
        lat, lon, _ = self._apparent_position(J2000).galactic_latlon()
        lat, lon = degrees(lat), degrees(lon)
        return lon, lat
        
    def set(self, lon, lat):
        """
        Set the longitude and latitude coordinates for the instance.
        """
        
        lon, lat = degrees(lon), degrees(lat)
        ra, dec = _equ_from_gal(lon, lat, self._epoch)
        self.set(ra, dec)


def separation(pos1, pos2):
    """
    Return the angular separation between two objects or positions as an
    Angle instance.
    """
    
    c1 = None
    e1 = None
    if isinstance(pos1, (tuple, list)):
        ra, dec = pos1
        ra, dec = hours(ra), degrees(dec)
        c1 = FixedBody()
        c1._ra = ra
        c1._dec = dec
        c1._epoch = J2000
        e1 = c1._epoch
    elif isinstance(pos1, Body):
        ra, dec = pos1.g_ra, pos1.g_dec
        c1 = FixedBody()
        c1._ra = ra
        c1._dec = dec
        c1._epoch = J2000
        e1 = c1._epoch
    elif isinstance(pos1, Coordinate):
        ra, dec = pos1.to_radec()
        c1 = FixedBody()
        c1._ra = ra
        c1._dec = dec
        c1._epoch = pos1._epoch
        e1 = c1._epoch
    if c1 is None:
        raise TypeError("Unexpected type for pos1: %s" % type(pos1))
        
    c2 = None
    e2 = None
    if isinstance(pos2, (tuple, list)):
        ra, dec = pos2
        ra, dec = hours(ra), degrees(dec)
        c2 = FixedBody()
        c2._ra = ra
        c2._dec = dec
        c2._epoch = J2000
        e2 = c2._epoch
    elif isinstance(pos2, Body):
        ra, dec = pos2.g_ra, pos2.g_dec
        c2 = FixedBody()
        c2._ra = ra
        c2._dec = dec
        c2._epoch = J2000
        e2 = c2._epoch
    elif isinstance(pos2, Coordinate):
        ra, dec = pos2.to_radec()
        c2 = FixedBody()
        c2._ra = ra
        c2._dec = dec
        c2._epoch = pos2._epoch
        e2 = c2._epoch
    if c2 is None:
        raise TypeError("Unexpected type for pos2: %s" % type(pos1))
        
    e = e1 if e2 == J2000 else e2
    c1 = _solar_system['earth'].at(e).observe(c1._body)
    c2 = _solar_system['earth'].at(e).observe(c2._body)
    distance = c1.separation_from(c2)
    return degrees(distance)
