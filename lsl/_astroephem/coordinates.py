from __future__ import print_function, division

from astropy.time import Time
from astropy.coordinates import ICRS, GeocentricTrueEcliptic, Galactic as AstroGalactic
import astropy.units as u

from angles import hours, degrees
from dates import Date

__all__ = ['Equatorial', 'Ecliptic', 'Galactic']


class _Coordinate(object):
    __astropy_base__ = ICRS
    
    def __init__(self, c1_or_class, c2=None, epoch=None, hours=False):
        if epoch is None:
            epoch = Time.now()
        self._epoch = epoch
        self._hours = hours
        
        if isinstance(c1_or_class, _Coordinate):
            try:
                new_frame = self.__astropy_base__(obstime=self._epoch)
            except TypeError:
                new_frame = self.__astropy_base__()
            self._coord = c1_or_class._coord.transform_to(new_frame)
        else:
            assert(c2 is not None)
            self.set(c1_or_class, c2)
            
    def get(self):
        try:
            c1, c2 = self._coord.lon, self._coord.lat
        except AttributeError:
            try:
                c1, c2 = self._coord.l, self._coord.b
            except AttributeError:
                c1, c2 = self._coord.ra, self._coord.dec
        return c1, c2
      
    def to_radec(self):
        radec = self._coord.transform_to(ICRS())
        return (radec.ra, radec.dec)
        
    def set(self, c1, c2):
        c1 = hours(c1) if self._hours else degrees(c1)
        c2 = degrees(c2)
        try:
            self._coord = self.__astropy_base__(c1, c2, obstime=self._epoch)
        except TypeError:
            self._coord = self.__astropy_base__(c1, c2)
            
    def from_radec(self, ra, dec):
        radec = ICRS(hours(ra), degres(dec))
        try:
            new_frame = self.__astropy_base__(obstime=self._epoch)
        except TypeError:
            new_frame = self.__astropy_base__()
        self._coord = c1_or_class.transform_to(new_frame)
               

class Equatorial(_Coordinate):
    __astropy_base__ = ICRS
    
    def __init__(self, ra_or_class, dec=None, epoch=None):
        _Coordinate.__init__(self, ra_or_class, dec, epoch=epoch, hours=True)
        
    @property
    def ra(self):
        return self._coord.ra.to('radian').value
        
    @property
    def dec(self):
        return self._coord.dec.to('radian').value


class Ecliptic(_Coordinate):
    __astropy_base__ = GeocentricTrueEcliptic
    
    def __init__(self, lon_or_class, lat=None, epoch=None):
        _Coordinate.__init__(self, lon_or_class, lat, epoch=epoch, hours=False)
        
    @property
    def lon(self):
        return self._coord.lon.to('radian').value
        
    @property
    def lat(self):
        return self._coord.lat.to('radian').value


class Galactic(_Coordinate):
    __astropy_base__ = AstroGalactic
    
    def __init__(self, lon_or_class, lat=None, epoch=None):
        _Coordinate.__init__(self, lon_or_class, lat, epoch=epoch, hours=False)
        
    @property
    def lon(self):
        return self._coord.l.to('radian').value
        
    @property
    def lat(self):
        return self._coord.b.to('radian').value
