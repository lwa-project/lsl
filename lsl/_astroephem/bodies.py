from __future__ import print_function, division

import numpy

from astropy.time import Time
from astropy.coordinates import SkyCoord, ITRS, AltAz, get_body, get_moon
from astropy.coordinates import CartesianRepresentation
import astropy.units as u

from angles import hours, degrees
from dates import Date
from observers import Observer


__all__ = ['FixedBody', 'readdb', 'Sun', 'Mercury', 'Venus', 'Moon', 'Mars', 
           'Jupiter', 'Saturn', 'Uranus', 'Neptune']


class _Body(object):
    """
    Base class that provides the compute() method for FixedBody and Planet.
    """
    
    def __init__(self, skycoord_or_body=SkyCoord(0*u.hourangle, 0*u.deg, frame='icrs')):
        self._sc = skycoord_or_body
        
    def compute(self, date_or_observer=None):
        if date_or_observer is None:
            date_or_observer = Time.now()
        elif isinstance(date_or_observer, str):
            date_or_observer = Date(date_or_observer)
        elif isinstance(date_or_observer, float):
            date_or_observer = Date(date_or_observer)
            
        if isinstance(date_or_observer, Time):
            self.a_epoch = date_or_observer
            gast = self.a_epoch.sidereal_time('apparent', longitude=0)
            ogc = CartesianRepresentation(0.0, 0.0, 0.0)
            
            try:
                _ac = self._sc.apply_space_motion(new_obstime=date_or_observer)
            except (AttributeError, ValueError):
                _ac = self._sc
            _gc = _ac.transform_to(ITRS(obstime=self.a_epoch))
            _lc = _ac.transform_to(ITRS(ogc,
                                        obstime=self.a_epoch))
        else:
            self.a_epoch = date_or_observer.date
            gast = self.a_epoch.sidereal_time('apparent', longitude=0)
            ogc = CartesianRepresentation(*date_or_observer._el.to_geocentric())
            
            try:
                _ac = self._sc.apply_space_motion(new_obstime=date_or_observer)
            except (AttributeError, ValueError):
                _ac = self._sc
            _gc = _ac.transform_to(ITRS(obstime=self.a_epoch))
            _lc = _ac.transform_to(ITRS(ogc,
                                        obstime=self.a_epoch))
            _tc = _lc.transform_to(AltAz(obstime=date_or_observer.date,
                                         location=date_or_observer._el))
        self.a_ra  = hours(_ac.ra, wrap=True)
        self.a_dec = degrees(_ac.dec, wrap180=True)
        self.g_ra  = hours(gast + _gc.spherical.lon, wrap=True)
        self.g_dec = degrees(_gc.spherical.lat, wrap180=True)
        try:
            self.ra  = hours(gast + _lc.spherical.lon, wrap=True)
            self.dec = degrees(_lc.spherical.lat, wrap180=True)
        except NameError:
            pass
        try:
            self.az = degrees(_tc.az)
            self.alt = degrees(_tc.alt)
        except NameError:
            pass
            
        self._rise_transit_set(date_or_observer)
            
    def _rise_transit_set(self, observer):
        """
        From:
            https://github.com/brandon-rhodes/pyephem/blob/master/libastro-3.7.7/riset.c
        """
        
        if not isinstance(observer, Observer):
            return False
            
        try:
            self.g_ra
        except AttributeError:
            return False
            
        # Load the values we need in radians
        ra = self.g_ra.to('radian').value
        dec = self.g_dec.to('radian').value
        lat = observer.lat.to('radian').value
        
        # Flip things around if we are south of the equator
        southern = (lat < 0)
        if southern:
            lat = -lat
            dec = -dec
            
        # Figure out how high it gets, if it is circumpolar, and if it never rises
        z = numpy.pi/2
        zmin = abs(dec - lat)
        zmax = numpy.pi - abs(dec + lat)
        if zmax <= z + 1e-9:
            circumpolar = True
        else:
            circumpolar = False
        if zmin >= z - 1e-9:
            neverup = True
        else:
            neverup = False
        transit_alt = zmin
            
        # Find the hour angle that it rises at
        cos_h = (numpy.cos(z)-numpy.sin(lat)*numpy.sin(dec))/(numpy.cos(lat)*numpy.cos(dec));
        if cos_h >= 1:
            rising_ha =  0.0
        elif cos_h <= -1:
            rising_ha = numpy.pi
        else:
            rising_ha = numpy.arccos(cos_h)
            
        # Find the setting azimuth
        xaz = numpy.sin(dec)*numpy.cos(lat)-numpy.cos(dec)*numpy.cos(rising_ha)*numpy.sin(lat);
        yaz = -numpy.cos(dec)*numpy.sin(rising_ha);
        if xaz == 0.:
	        if yaz > 0:
	            setting_az = numpy.pi/2
	        else:
	            setting_az = -numpy.pi/2
        else:
	        setting_az = numpy.arctan2(yaz, xaz)
        if southern:
	        setting_az = numpy.pi - setting_az
        rising_az = 2*numpy.pi - setting_az
	    
        # Find the rising and setting local sidereal times
        rising_lst = ra - rising_ha
        setting_lst = ra + rising_ha

        # Save it all
        ## Direct access
        self.circumpolar = circumpolar
        self.neverup = neverup
        self.transit_alt = degrees(transit_alt)
        self.rising_az = degrees(rising_az, wrap360=True)
        self.setting_az = degrees(setting_az, wrap360=True)
        ## Used later
        self._rising_lst = hours(rising_lst, wrap=True)
        self._setting_lst = hours(setting_lst, wrap=True)


class FixedBody(_Body):
    """
    A celestial body, that can compute() its sky position.
    """
    
    name = ''
    __ra = 0.0*u.hourangle
    __dec = 0.0*u.deg
    __epoch = "J2000.0"
    __pmra = 0.0*u.mas/u.yr
    __pmdec = 0.0*u.mas/u.yr
    
    def __repr__(self):
        return "<%s.%s %s at 0x%x>" % (type(self).__module__, type(self).__name__, '"%s"' % self.name if self.name else 'None', id(self))
        
    def _update(self):
        self._sc = SkyCoord(self.__ra, self.__dec,
                            frame='icrs', obstime=self.__epoch)
    
    @property
    def _ra(self):
        return self.__ra
    @_ra.setter
    def _ra(self, value):
        self.__ra = hours(value)
        self._update()
        
    @property
    def _dec(self):
        return self.__dec
    @_dec.setter
    def _dec(self, value):
        self.__dec = degrees(value)
        self._update()
        
    @property
    def _epoch(self):
        return self.__epoch
    @_epoch.setter
    def _epoch(self, value):
        self.__epoch = Date(value)
        self._update()
        
    @property
    def _pmra(self):
        return self.__pmra.value
    @_pmra.setter
    def _pmra(self, value):
        self.__pmra = value*u.mas/u.yr
        self._update()
    @property
    def _pmdec(self):
        return self.__pmdec.value
    @_pmdec.setter
    def _pmdec(self, value):
        self.__pmdec = value*u.mas/u.yr
        self._update() 


def readdb(line):
    fields = line.split(',')
    name, type, ra, dec, _ = fields
    if type != 'f|J':
        raise ValueError("Can only deal with database entries that are 'f|J'")
    bdy = FixedBody()
    bdy.name = name
    bdy._ra = ra
    bdy._dec = dec
    return bdy


class Planet(_Body):
    """
    A solar system body, that can compute() its sky position.
    """
    
    def __init__(self, func):
        self.func = func
        
    def __repr__(self):
        return "<%s.%s %s at 0x%x>" % (type(self).__module__, type(self).__name__, '"%s"' % self.name if self.name else 'None', id(self))
        
    def compute(self, date_or_observer=None):
        if date_or_observer is None:
            date_or_observer = Time.now()
        elif isinstance(date_or_observer, str):
            date_or_observer = Date(date_or_observer)
        elif isinstance(date_or_observer, float):
            date_or_observer = Date(date_or_observer)
            
        if isinstance(date_or_observer, Time):
            self._sc = self.func(date_or_observer)
        else:
            self._sc = self.func(date_or_observer.date)
            
        _Body.compute(self, date_or_observer=date_or_observer)            


class Sun(Planet):
    """
    A Planet instance representing the Sun.
    """
    
    name = 'Sun'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('sun', x))


class Mercury(Planet):
    """
    A Planet instance representing Mercury.
    """
    
    name = 'Mercury'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('mercury', x))


class Venus(Planet):
    """
    A Planet instance representing Venus.
    """
    
    name = 'Venus'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('venus', x))


class Moon(Planet):
    """
    A Planet instance representing Earth's moon"""
    
    name = 'Moon'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_moon(x))


class Mars(Planet):
    """
    A Planet instance representing Mars.
    """
    
    name = 'Mars'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('mars', x))


class Jupiter(Planet):
    """
    A Planet instance representing Jupiter.
    """
    
    name = 'Jupiter'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('jupiter', x))


class Saturn(Planet):
    """
    A Planet instance representing Saturn.
    """
    
    name = 'Saturn'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('saturn', x))


class Uranus(Planet):
    """
    A Planet instance representing Uranus.
    """
    
    name = 'Uranus'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('uranus', x))


class Neptune(Planet):
    """
    A Planet instance representing Neptune.
    """
    
    name = 'Neptune'
    
    def __init__(self):
        Planet.__init__(self, lambda x: get_body('neptune', x))
