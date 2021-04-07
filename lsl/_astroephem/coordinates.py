"""
Python module providing ephem-like conversions between coordinate systems.
"""

from __future__ import print_function, division

from astropy.time import Time
from astropy.coordinates import ICRS, ITRS, CIRS, FK4, FK5, GeocentricTrueEcliptic, Galactic as AstroGalactic
import astropy.units as u

from lsl._astroephem.angles import hours, degrees
from lsl._astroephem.dates import Date, B1950, J2000

__all__ = ['Equatorial', 'Ecliptic', 'Galactic']


class EphemCoordinate(object):
    """
    Base class to facilitate the coordinate transforms.
    """
    
    __astropy_base__ = ICRS
    
    def __init__(self, c1_or_class, c2=None, epoch=None, hours=False):
        if epoch is None:
            epoch = Time.now()
        self._epoch = Date(epoch)
        self._hours = hours
        
        if isinstance(c1_or_class, EphemCoordinate):
            try:
                new_frame = self.__astropy_base__(equinox=self._epoch,
                                                  obstime=self._epoch)
            except TypeError:
                try:
                    new_frame = self.__astropy_base__(obstime=self._epoch)
                except TypeError:
                    new_frame = self.__astropy_base__()
            self._coord = c1_or_class._coord.transform_to(new_frame)
        else:
            assert(c2 is not None)
            self.set(c1_or_class, c2)
            
    def get(self):
        """
        Get the coordinate pair as two-element tuple of Angle instances.
        Exactly what is returned is dependent on the frame used for the
        instance.
        """
        
        try:
            c1, c2 = degrees(self._coord.lon), degrees(self._coord.lat)
        except AttributeError:
            try:
                c1, c2 = degrees(self._coord.l), degrees(self._coord.b)
            except AttributeError:
                c1 = hours(self._coord.ra) if self._hours else degrees(self._coord.ra)
                c2 = degrees(self._coord.dec)
        return c1, c2
        
    def to_radec(self):
        """
        Return the ICRS right assencsion and declination for the coordinate.
        """
        
        radec = self._coord.transform_to(ICRS())
        return (hours(radec.ra), degrees(radec.dec))
        
    def set(self, c1, c2):
        """
        Set the coordinates for the instance.  Exactly what is set is dependent
        on the frame used for the instance.
        """
        
        c1 = hours(c1) if self._hours else degrees(c1)
        c2 = degrees(c2)
        try:
            self._coord = self.__astropy_base__(c1, c2, equinox=self._epoch,
                                                        obstime=self._epoch)
        except TypeError:
            try:
                self._coord = self.__astropy_base__(c1, c2, obstime=self._epoch)
            except TypeError:
                self._coord = self.__astropy_base__(c1, c2)
                
    def from_radec(self, ra, dec):
        """
        Set the coordinates for the instance using the provided ICRS right
        ascension and declination.
        """
        
        radec = ICRS(hours(ra), degres(dec))
        try:
            new_frame = self.__astropy_base__(equinox=self._epoch,
                                              obstime=self._epoch)
        except TypeError:
            try:
                new_frame = self.__astropy_base__(obstime=self._epoch)
            except TypeError:
                new_frame = self.__astropy_base__()
        self._coord = c1_or_class.transform_to(new_frame)
               

class Equatorial(EphemCoordinate):
    """
    Class for representing an Equatorial position on the sky.  The default frame
    is ICRS.  Setting an epoch of J2000 sets the frame to FK5 and setting it to
    B1950 set the frame to FK4.
    """
    
    __astropy_base__ = ICRS
    
    def __init__(self, ra_or_class, dec=None, epoch=None):
        if epoch == B1950:
            self.__astropy_base__ = FK4
        elif epoch == J2000:
            self.__astropy_base__ = FK5
            
        EphemCoordinate.__init__(self, ra_or_class, dec, epoch=epoch, hours=True)
        
    @property
    def ra(self):
        """
        Right ascension in the current frame.
        """
        
        return hours(self._coord.ra)
        
    @property
    def dec(self):
        """
        Declination in the current frame.
        """
        
        return degrees(self._coord.dec)


class Ecliptic(EphemCoordinate):
    """
    Class for representing an Ecliptic position on the sky in the geocentric true
    ecliptic frame.
    """
    
    __astropy_base__ = GeocentricTrueEcliptic
    
    def __init__(self, lon_or_class, lat=None, epoch=None):
        EphemCoordinate.__init__(self, lon_or_class, lat, epoch=epoch, hours=False)
        
    @property
    def lon(self):
        """
        Ecliptic longitude.
        """
        
        return degrees(self._coord.lon)
        
    @property
    def lat(self):
        """
        Ecliptic latitude.
        """
        
        return degrees(self._coord.lat)


class Galactic(EphemCoordinate):
    """
    Class for representing a Galactic position on the sky in the geocentric true
    ecliptic frame.
    """
    
    __astropy_base__ = AstroGalactic
    
    def __init__(self, lon_or_class, lat=None, epoch=None):
        EphemCoordinate.__init__(self, lon_or_class, lat, epoch=epoch, hours=False)
        
    @property
    def lon(self):
        """
        Galactic longitude.
        """
        
        return degrees(self._coord.l)
        
    @property
    def lat(self):
        """
        Galactic latitude.
        """
        
        return degrees(self._coord.b)
