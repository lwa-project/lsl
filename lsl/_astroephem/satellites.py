import numpy
from sgp4 import api as sgp4_api

from astropy.time import Time
from astropy.coordinates import ITRS, AltAz, ICRS, CIRS, TEME, CartesianDifferential, CartesianRepresentation
from astropy import units as u

from angles import hours, degrees
from dates import Date

__all__ = ['EphemEarthSatellte', 'readtle']


class EphemEarthSatellte(object):
    """
    An Earth orbiting satellite, that can compute() its sky position.
    """
    
    def __init__(self, name, sgp4_api_tle):
        self.name = name.strip().rstrip()
        self._tle = sgp4_api_tle
        
    def compute(self, date_or_observer=None):
        """
        Given an EphemTime or Observer instance or None, compute the location of the
        satellite.
        
        For an EphemTime, this computes:
         * the sub-body location - sublat, sublon, and elevation
         * the geocentric position of the body - g_ra and g_dec
         
        For an Observer, this computes everything as in the case of EphemTime
        and also includes:
         * the astrometric position of the body - a_ra and a_dec
         * the apparent position of the body - ra and dec
         * the topocentric position of the body - az and alt
         * the range to the body - range
         * the line-of-sight velocity of the body - range_velocity
         
        If None is provided, the current time is used.
        """
        
        if date_or_observer is None:
            date_or_observer = Time.now()
        elif isinstance(date_or_observer, str):
            date_or_observer = Date(date_or_observer)
        elif isinstance(date_or_observer, float):
            date_or_observer = Date(date_or_observer)
            
        if isinstance(date_or_observer, Time):
            self.a_epoch = date_or_observer
        else:
            self.a_epoch = date_or_observer.date
        error_code, teme_p, teme_v = self._tle.sgp4(self.a_epoch.jd1, self.a_epoch.jd2)
        if error_code != 0:
            raise RuntimeError(sgp4_api.SGP4_ERRORS[error_code])
            
        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        _ac = TEME(teme_p.with_differentials(teme_v), obstime=self.a_epoch)
        
        if isinstance(date_or_observer, Time):
            gast = self.a_epoch.sidereal_time('apparent', longitude=0*u.deg)
            ogc = CartesianRepresentation(0.0*u.m, 0.0*u.m, 0.0*u.m)
            
            _gc = _ac.transform_to(ITRS(obstime=self.a_epoch))
            _lc = _ac.transform_to(ITRS(ogc,
                                        obstime=self.a_epoch))
        else:
            gast = self.a_epoch.sidereal_time('apparent', longitude=0*u.deg)
            ogc = CartesianRepresentation(*date_or_observer._el.to_geocentric())
            
            _gc = _ac.transform_to(ITRS(obstime=self.a_epoch))
            _lc = _ac.transform_to(ITRS(ogc,
                                        obstime=self.a_epoch))
            _tc = _ac.transform_to(AltAz(obstime=date_or_observer.date,
                                         location=date_or_observer._el))
            
        sat = _gc.earth_location.geodetic
        self.sublat = degrees(sat.lat)
        self.sublon = degrees(sat.lon)
        self.elevation = sat.height.to(u.m).value
        
        self.g_ra  = hours(gast + _gc.spherical.lon, wrap=True)
        self.g_dec = degrees(_gc.spherical.lat, wrap180=True)
        try:
            self.az = degrees(_tc.az)
            self.alt = degrees(_tc.alt)
            
            self.a_ra, self.a_dec = date_or_observer.radec_of(self.az, self.alt)
            _ec = ICRS(self.a_ra, self.a_dec).transform_to(CIRS(obstime=self.a_epoch))
            self.ra = hours(_ec.ra, wrap=True)
            self.dec = degrees(_ec.dec, wrap180=True)
            
            self.range = _tc.distance.to(u.m).value
            self.range_velocity = _tc.radial_velocity.to(u.m/u.s).value
        except NameError:
            pass


def readtle(name, line1, line2):
    """
    Read TLE-format satellite elements.
    """
    
    tle = sgp4_api.Satrec.twoline2rv(line1, line2)
    return EphemEarthSatellte(name, tle)
