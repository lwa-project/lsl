"""
Python module providing ephem-like conversions between coordinate systems.
"""

from __future__ import print_function, division

from skyfield import api
from skyfield.positionlib import Astrometric, position_of_radec
ts = api.load.timescale()

from lsl._skyephem.angles import hours, degrees
from lsl._skyephem.dates import Date, B1950, J2000

__all__ = ['Coordinate', 'Equatorial', 'Ecliptic', 'Galactic']


class Coordinate(object):
    """
    Base class to facilitate the coordinate transforms in a way that behaves
    like ephem.Coordinate.
    """
    
    def __init__(self, c1_or_class, c2=None, epoch=None):
        if epoch is None:
            epoch = Time.now()
        self._epoch = Date(epoch)
        
        try:
            ra, dec = c1_or_class.to_radec()
        except AttributeError:
            ra, dec = hours(c1_or_class), degrees(c2)
        self._coord = position_of_radec(ra.hours, dec.degrees, epoch=self._epoch)
            
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
        
        ra, dec, _ = self._coord.radec()
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
        
        ra, dec, _ = self._coord.radec()
        ra, dec = hours(ra), degrees(dec)
        return ra, dec
        
    def set(self, ra, dec):
        """
        Set the coordinates for the instance.  Exactly what is set is dependent
        on the frame used for the instance.
        """
        
        ra, dec = hours(c1_or_class), degrees(c2)
        self._coord = position_of_radec(ra, dec, epoch=self._epoch)
        

class Ecliptic(Coordinate):
    """
    Class for representing an Ecliptic position on the sky in the geocentric true
    ecliptic frame.
    """
    
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
        
        lat, lon, _ = self._coord.ecliptic_latlon()
        lat, lon = degrees(lat), degrees(lon)
        return lon, lat
        

class Galactic(Coordinate):
    """
    Class for representing a Galactic position on the sky in the geocentric true
    ecliptic frame.
    """
    
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
        
        lat, lon, _ = self._coord.galactic_latlon()
        lat, lon = degrees(lat), degrees(lon)
        return lon, lat
