from __future__ import print_function, division

from astropy.time import Time
from astropy.coordinates import SkyCoord, GCRS, AltAz, get_body, get_moon
from astropy.coordinates import CartesianRepresentation
import astropy.units as u
import numpy
import math

from angles import hours, degrees
from dates import Date
from observers import Observer


__all__ = ['FixedBody', 'readdb', 'Sun', 'Mercury', 'Venus', 'Moon', 'Mars', 
           'Jupiter', 'Saturn', 'Uranus', 'Neptune']


class _Body(object):
    def __init__(self, skycoord_or_body=SkyCoord(0*u.hourangle, 0*u.deg, frame='icrs')):
        self._sc = skycoord_or_body
        
    def compute(self, date_or_observer=None):
        if date_or_observer is None:
            date_or_observer = Time.now()
        elif isinstance(date_or_observer, float):
            date_or_observer = Date(date_or_observer)
            
        if isinstance(date_or_observer, Time):
            self.a_epoch = date_or_observer
            try:
                _ac = self._sc.apply_space_motion(new_obstime=date_or_observer)
            except (AttributeError, ValueError):
                _ac = self._sc
            _gc = _ac.transform_to(GCRS(obstime=self.a_epoch))
        else:
            self.a_epoch = date_or_observer.date
            try:
                _ac = self._sc.apply_space_motion(new_obstime=date_or_observer)
            except (AttributeError, ValueError):
                _ac = self._sc
            _gc = _ac.transform_to(GCRS(obstime=self.a_epoch))
            cr = CartesianRepresentation(*date_or_observer._el.to_geocentric())
            _lc = _ac.transform_to(GCRS(obstime=self.a_epoch,
                                        obsgeoloc=cr))
            _tc = _lc.transform_to(AltAz(obstime=date_or_observer.date,
                                         location=date_or_observer._el))
        self.a_ra  = _ac.ra.to('radian').value
        self.a_dec = _ac.dec.to('radian').value
        self.g_ra  = _gc.ra.to('radian').value
        self.g_dec = _gc.dec.to('radian').value
        try:
            self.ra  = _lc.ra.to('radian').value
            self.dec = _lc.dec.to('radian').value
        except NameError:
            pass
        try:
            self.az = _tc.az.to('radian').value
            self.alt = _tc.alt.to('radian').value
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
        ra = self.g_ra
        dec = self.g_dec
        lat = observer.lat
        
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
            rising_ha = numpy.py;
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
        self.rising_az = degrees(rising_az % (2*numpy.pi))
        self.setting_az = degrees(setting_az % (2*numpy.pi))
        ## Used later
        self._rising_lst = hours(rising_lst % (2*numpy.pi))
        self._setting_lst = hours(setting_lst % (2*numpy.pi))


class FixedBody(_Body):
    name = ''
    __ra = 0.0*u.hourangle
    __dec = 0.0*u.deg
    __epoch = "J2000.0"
    __pmra = 0.0*u.mas/u.yr
    __pmdec = 0.0*u.mas/u.yr
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


class Sun(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('sun', Time.now()))


class Mercury(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('sun', Time.now()))


class Venus(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('venus', Time.now()))


class Moon(_Body):
    def __init__(self):
        _Body.__init__(self, get_moon(Time.now()))


class Mars(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('mars', Time.now()))


class Jupiter(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('jupiter', Time.now()))


class Saturn(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('saturn', Time.now()))


class Uranus(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('uranus', Time.now()))


class Neptune(_Body):
    def __init__(self):
        _Body.__init__(self, get_body('neptune', Time.now()))
