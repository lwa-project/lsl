import numpy

from skyfield import api
from skyfield.toposlib import Topos
from skyfield.timelib import Time as SkyTime
ts = api.load.timescale()

from lsl._skyephem.angles import hours, degrees
from lsl._skyephem.dates import Date
from lsl._skyephem.bodies import Body, prepare_date_or_observer
from lsl._skyephem.cache import load_planetary_ephemeris


__all__ = ['EarthSatellite', 'readtle']


_solar_system = load_planetary_ephemeris()


class EarthSatellite(Body):
    """
    An Earth orbiting satellite, that can compute() its sky position.
    """
    
    def __init__(self, name, skyfield_earthsatellite):
        self.name = name.strip().rstrip()
        self._body = skyfield_earthsatellite
        
    @prepare_date_or_observer
    def compute(self, date_or_observer=None):
        """
        Given a Date or Observer instance or None, compute the location of the
        satellite.
        
        For a Date, this computes:
         * the sub-body location - sublat, sublon, and elevation
         
        For an Observer, this computes everything as in the case of a Date
        and also includes:
         * the astrometric position of the body - a_ra and a_dec
         * the apparent position of the body - ra and dec
         * the topocentric position of the body - az and alt
         * the range to the body - range
         * the line-of-sight velocity of the body - range_velocity
         
        If None is provided, the current time is used.
        """
        
        if isinstance(date_or_observer, SkyTime):
            t = date_or_observer
            obs = None
        else:
            t = date_or_observer.date
            obs = _solar_system['earth'] + date_or_observer._wgs84
            
        geo = self._body.at(t)
        sub = geo.subpoint()
        self.sublat = degrees(sub.latitude)
        self.sublon = degrees(sub.longitude)
        self.elevation = sub.elevation.m
        
        if obs is not None:
            diff = self._body - date_or_observer._wgs84
            topo = diff.at(t)
            
            alt, az, dist = topo.altaz()
            self.az = degrees(az)
            self.alt = degrees(alt)
            self.range = dist.m
            
            vel = topo.velocity.km_per_s * 1000 # m/s
            dot = topo.position.m / numpy.sqrt(((topo.position.m)**2).sum())
            self.range_velocity = (vel*dot).sum()
            
            ra, dec, _ = topo.radec()
            self.a_ra = hours(ra)
            self.a_dec = degrees(dec)
            
            ra, dec, _ = topo.radec(epoch=t)
            self.ra = hours(ra)
            self.dec = degrees(dec)


def readtle(name, line1, line2):
    """
    Read TLE-format satellite elements.
    """
    
    tle = api.EarthSatellite(line1, line2, name, ts)
    return EarthSatellite(name, tle)
