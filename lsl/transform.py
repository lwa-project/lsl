"""
Time and position transform objects.
"""

from collections.abc import Sequence as SequenceABC
import copy
import datetime
import math
import abc
from functools import total_ordering

from astropy.time import Time as AstroTime

from lsl import astro
from lsl.common.ndp import fS

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.4'
__all__ = ['Time', 'SkyPosition', 'CelestialPosition', 'PlanetaryPosition', 
           'GeographicalPosition', 'PointingDirection']
__author__ = "Unknown"
__maintainer__ = "Jayce Dowell"


@total_ordering
class Time(object):
    """
    Holds an absolute time value and can manipulate the value in 
    various representations.
    
    After creation, a particular time format and time system may
    be accessed using the appropriate instance member.  If marked
    with '(S)', the Time value may also be updated by setting the
    member to a new value.
    
    utc_jd (S)        - UTC standard Julian day
    utc_mjd (S)       - UTC modified Julian day
    utc_timet (S)     - UTC UNIX timet seconds
    utc_py_date (S)   - UTC python datetime.datetime object
    utc_dp (S)        - UTC DP samples at 196 MHz
    utc_mcs           - UTC MCS MJD/MPM pair
    utc_str           - UTC ISO8601 calendar string format
    
    tai_jd (S)        - TAI standard Julian day
    tai_mjd (S)       - TAI modified Julian day
    tai_timet (S)     - TAI UNIX timet seconds
    """

    # time format types

    FORMAT_PY_DATE    = 'PY_DATE'
    FORMAT_STR        = 'STR'
    FORMAT_JD         = 'JD'
    FORMAT_MJD        = 'MJD'
    FORMAT_DP         = 'DP'
    FORMAT_MCS        = 'MCS'
    FORMAT_TIMET      = 'TIMET'
    FORMAT_ASTROPY    = 'ASTROPY'
    
    known_formats = (FORMAT_PY_DATE, FORMAT_STR, FORMAT_JD, 
    FORMAT_MJD, FORMAT_DP, FORMAT_MCS, FORMAT_TIMET, FORMAT_ASTROPY)
    
    # time system types
    TIMESYS_UTC       = 'UTC'
    TIMESYS_TAI       = 'TAI'
    
    known_timesys = (TIMESYS_UTC, TIMESYS_TAI)
    
    @classmethod
    def from_system(klass):
        """
        Factory method to create a Time instance from current system clock value.
        """
        
        return klass(astro.get_julian_from_sys(), klass.FORMAT_JD, klass.TIMESYS_UTC)
    
    def __init__(self, value, format = FORMAT_MJD, timesys = TIMESYS_UTC):
        """
        Create a Time instance, using 'value' as the initial time.
        
        'format' describes the type of 'value'
            Time.FORMAT_PY_DATE - python datetime.datetime class calendar format
            Time.FORMAT_STR     - ISO 8601 (YYYY-MM-DD hh:mm:ss.s) calendar format 
                                  string 
            Time.FORMAT_JD      - standard julian day as float
            Time.FORMAT_MJD     - modified julian day as float
            Time FORMAT_DP      - samples of the DP 196 MHz clock as integer
            Time FORMAT_MCS     - MCS MJD, MPM pair as integers
            Time.FORMAT_TIMET   - UNIX timet seconds as integer
            Time.FORMAT_ASTROPY - astropy.time.Time instance
        
        'timesys' defines the time system
            Time.TIMESYS_UTC  - UTC time system
            Time.TIMESYS_TAI  - TAI time system
            
        For Time.FORMAT_ASTROPY the 'timesys' is ignored
        """
        
        # check parameters
        if format not in self.known_formats:
            raise ValueError(f"unknown format {format}")
            
        if timesys not in self.known_timesys:
            raise ValueError(f"unknown timesys {timesys}")
        
        # parse init value base on format type
        # time value is held internally as UTC JD float
        if format == self.FORMAT_PY_DATE:
            if timesys == self.TIMESYS_UTC:
                self.utc_py_date = value
            else:
                raise ValueError(f"{self.FORMAT_PY_DATE} not supported for TIMESYS_TAI")
                
        elif format == self.FORMAT_STR:
            if timesys == self.TIMESYS_UTC:
                self.utc_str = value
            else:
                raise ValueError(f"{self.FORMAT_STR} not supported for TIMESYS_TAI")
                
        elif format == self.FORMAT_JD:
            if timesys == self.TIMESYS_UTC:
                self.utc_jd = value
            else:
                self.tai_jd = value
                
        elif format == self.FORMAT_MJD:
            if timesys == self.TIMESYS_UTC:
                self.utc_mjd = value
            else:
                self.tai_mjd = value
                
        elif format == self.FORMAT_DP:
            if timesys == self.TIMESYS_UTC:
                self.utc_dp = value
            else:
                raise ValueError(f"{self.FORMAT_DP} not supported for TIMESYS_TAI")
                
        elif format == self.FORMAT_MCS:
            if timesys == self.TIMESYS_UTC:
                self.utc_mcs = value
            else:
                raise ValueError(f"{self.FORMAT_MCS} not supported for TIMESYS_TAI")
                
        elif format == self.FORMAT_TIMET:
            if timesys == self.TIMESYS_UTC:
                self.utc_timet = value
            else:
                self.tai_timet = value
            
        elif format == self.FORMAT_ASTROPY:
            self.astropy = value
            
    def __cmp__(self, other):
        """
        Determine ordering for two Time instances.
        """
        
        if self._time < other._time:
            return -1
        elif self._time > other._time:
            return 1
        else:
            return 0
            
    def __eq__(self, other):
        return True if self.__cmp__(other) == 0 else False
        
    def __lt__(self, other):
        return True if self.__cmp__(other) < 0 else False
        
    def __hash__(self):
        """
        Return a hash key for the Time instance.
        """
        
        return hash(self._time)
    
    def __repr__(self):
        """
        Return low-level string representation of Time instance.
        """
        
        return "%s.%s(%s)" % (type(self).__module__, type(self).__name__, repr(self.utc_mjd))
        
    def __str__(self):
        """
        Return high-level string representation of Time instance.
        """
        
        return self.utc_str
        
    @property
    def utc_jd(self):
        """
        Time value formatted as UTC standard julian day (float).
        """
        
        return self._time.utc.jd
        
    @utc_jd.setter
    def utc_jd(self, value):

        if not isinstance(value, (int, float)):
            raise TypeError("value must be type int or float")
            
        self._time = AstroTime(value, format='jd', scale='utc')
        
    @property
    def utc_mjd(self):
        """
        Time value formatted as UTC modified julian day (float).
        """
        
        return self._time.utc.mjd
        
    @utc_mjd.setter
    def utc_mjd(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("value must be type int or float")
            
        self._time = AstroTime(value, format='mjd', scale='utc')
        
    @property
    def utc_dp(self):
        return int(self._time.utc.unix * fS)
        
    @utc_dp.setter
    def utc_dp(self, value):
        if not isinstance(value, int):
            raise TypeError("value must be type int")
            
        self._time = AstroTime(value // fS, value % fS, format='unix', scale='utc')
        
    @property
    def utc_mcs(self):
        mjd = int(self.utc_mjd)
        mpm = int(round((self.utc_mjd - mjd) * 86400 * 1000))
        return (mjd, mpm)
        
    @utc_mcs.setter
    def utc_mcs(self, value):
        if isinstance(value, (tuple, list)):
            if len(value) != 2:
                raise ValueError("value must be a two-element tuple or list")
            if not isinstance(value[0], (int, float)):
                raise TypeError("value[0] must be type int or float")
            if not isinstance(value[1], (int, float)):
                raise TypeError("value[1] must be type int or float")
            mjd, mpm = value
        else:
            raise TypeError("Expected a two-element tuple of int/float values")
        self._time = AstroTime(float(mjd) + float(mpm)/86400/1000, format='mjd', scale='utc')
        
    @property
    def utc_py_date(self):
        """
        Time value formatted as UTC calendar datetime.datetime object.
        """
        
        return self._time.utc.datetime
        lnDate = astro.get_date(self._time)
        (usec, sec) = math.modf(lnDate.seconds)
        pyDate = datetime.datetime(lnDate.years, lnDate.months, lnDate.days,
            lnDate.hours, lnDate.minutes, int(sec), int(usec * 1e6))
        return pyDate
        
    @utc_py_date.setter
    def utc_py_date(self, value):
        if not isinstance(value, datetime.datetime):
            raise ValueError("value must be type datetime.datetime")
            
        self._time = AstroTime(value, format='datetime', scale='utc')
        
    @property
    def utc_str(self):
        """
        Time value formatted as UTC ISO 8601 calendar string.
        """
        
        return self._time.utc.iso
        
    @utc_str.setter
    def utc_str(self, value):
        if not isinstance(value, str):
            raise TypeError("value must be type str")
            
        self._time = AstroTime(value, format='iso', scale='utc')
    
    @property
    def tai_jd(self):
        """
        Time value formatted as TAI standard julian day (float).
        """
        
        return self._time.tai.jd
        
    @tai_jd.setter
    def tai_jd(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("value must be type int or float")
            
        self._time = AstroTime(value, format='jd', scale='tai')
        
    @property
    def tai_mjd(self):
        """
        Time value formatted as TAI modified julian day (float).
        """
        
        return self._time.tai.mjd
        
    @tai_mjd.setter
    def tai_mjd(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("value must be type int or float")
            
        self._time = AstroTime(value, format='mjd', scale='tai')
        
    @property
    def utc_timet(self):
        """
        Time value formatted as UTC UNIX timet seconds.
        """
        
        return self._time.utc.unix
        
    @utc_timet.setter
    def utc_timet(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("value must be type int or float")
            
        self._time = AstroTime(value, format='unix', scale='utc')
        
    @property
    def tai_timet(self):
        """
        Time value formatted as TAI UNIX timet seconds.
        """
        
        return self._time.tai.unix_tai
        
    @tai_timet.setter
    def tai_timet(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("value must be type int or float")
            
        self._time = AstroTime(value, format='unix_tai', scale='tai')
        
    @property
    def astropy(self):
        """
        Time value as a astropy.time.Time instance.
        """
        
        return self._time
        
    @astropy.setter
    def astropy(self, value):
        if not isinstance(value, AstroTime):
            raise TypeError("value must be type astropy.time.Time")
            
        self._time = value


class SkyPosition(object):
    """
    Base abstract class for representing positions on the sky.
    """
    
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def apparent_equ(self, time_):
        """
        Return position formatted as apparent equatorial coordinates.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return value is object of type astro.equ_posn.
        """
        
        raise NotImplementedError()
        
    @abc.abstractmethod
    def apparent_ecl(self, time_):
        """
        Return position formatted as apparent ecliptic coordinates.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return alue is object of type astro.ecl_posn.
        """
        
        raise NotImplementedError()


class CelestialPosition(SkyPosition):
    """
    Holds a celestial object position value and can manipulate
    the value in various representations.
    
    After creation, the celestial coordinates may be accessed in different
    formats and epoch by using the appropriate instance member.  If marked
    with '(S)', the Time value may also be updated by setting the member to 
    a new value. 
    
    j2000_equ (S)     - J2000 equatorial coordinates
    j2000_gal (S)     - J2000 galactic coordinates
    j2000_ecl (S)     - J2000 ecliptic coordinates
    
    b1950_equ (S)     - B1950 equatorial coordinates
    
    The instance methods apparent_equ() and apparent_ecl() may be used to
    get the apparent equatorial or ecliptic coordinates for a particular
    time.
    """
    
    # format types
    FORMAT_EQU      = 'EQU'
    FORMAT_GAL      = 'GAL'
    FORMAT_ECL      = 'ECL'
    
    known_formats = (FORMAT_EQU, FORMAT_ECL, FORMAT_GAL)
    
    # epoch types
    EPOCH_J2000     = 'J2000'
    EPOCH_B1950     = 'B1950'
    
    known_epochs = (EPOCH_J2000, EPOCH_B1950)
    
    def __init__(self, value, format = FORMAT_EQU, epoch = EPOCH_J2000, name = ''):
        """
        Create a CelestialPosition object, using 'value' as the initial coordinates. 
        
        'format' describes the type of 'value'
        * CelestialPosition.FORMAT_EQU  - equatorial coordinates RA and DEC
        * CelestialPosition.FORMAT_GAL  - galactic coordinates l and b
        * CelestialPosition.FORMAT_ECL  - ecliptic longitude and latitude
        
        'ephoch' defines the coordinate system equinox epoch
        * CelestialPosition.EPOCH_J2000   - J2000 epoch
        * CelestialPosition.EPOCH_B1950   - B1950 epoch 
        
        For format type EQU, the 'value' parameter may be either an instance of 
        the astro.equ_posn class or a sequence of length 2, in which case 
        RA=value[0] and DEC=value[1].
        
        For format type GAL, the 'value' parameter may be either an instance of
        the astro.gal_posn class or a sequence of length 2, in which case
        l=value[0] and b=value[1].
        
        For format type ECL, the 'value' parameter may be either an instance of
        the astro.ecl_posn class or a sequence of length 2, in which case
        long=value[0], lat=value[1]. 
        """
        
        # check parameters
        if format not in self.known_formats:
            raise ValueError(f"unknown format {format}")
            
        if epoch not in self.known_epochs:
            raise ValueError(f"unknown epoch {epoch}")
            
        self.name = name
        
        # parse init value based on format type
        # position is held internally as J2000 astro.equ_posn
        
        if format == self.FORMAT_EQU:
            if epoch == self.EPOCH_J2000:
                self.j2000_equ = value
            elif epoch == self.EPOCH_B1950:
                self.b1950_equ = value
                
        elif format == self.FORMAT_GAL:
            if epoch == self.EPOCH_J2000:
                self.j2000_gal = value
            else:
                raise ValueError(f"epoch {epoch} not supported for GAL format")
        elif format == self.FORMAT_ECL:
            if epoch == self.EPOCH_J2000:
                self.j2000_ecl = value
            else:
                raise ValueError(f"epoch {epoch} not supported for ECL format")
                
    def __repr__(self):
        return "%s.%s(%s, name=%s)" % (type(self).__module__, type(self).__name__, repr(self._posn), repr(self.name))
        
    @property
    def j2000_equ(self):
        """
        Position formatted as J2000 equatorial coordinates.
        Value is object of type astro.equ_posn.
        """
        
        return self._posn
        
    @j2000_equ.setter
    def j2000_equ(self, value):
        if not isinstance(value, (astro.equ_posn, SequenceABC)):
            raise TypeError("value must be type astro.equ_posn or sequence of length 2")
        if isinstance(value, SequenceABC):
            if len(value) != 2:
                raise TypeError("value sequence must be length 2")
            value = astro.equ_posn(*value)
            
        self._posn = copy.copy(value)
        
    def apparent_equ(self, time_):
        """
        Return position formatted as apparent equatorial coordinates.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return value is object of type astro.equ_posn.
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be of type transform.Time")
            
        return astro.get_apparent_posn(self._posn, time_.utc_jd)
        
    @property
    def b1950_equ(self):
        """
        Position formatted as B1950 equatorial coordinates.
        Value is object of type astro.equ_posn.
        """
        
        return astro.J2000_to_B1950(self._posn)
        
    @b1950_equ.setter
    def b1950_equ(self, value):
        if not isinstance(value, (astro.equ_posn, SequenceABC)):
            raise TypeError("value must be type astro.equ_posn or sequence of length 2")
        if isinstance(value, SequenceABC):
            if len(value) != 2:
                raise TypeError("value sequence must be length 2")
            value = astro.equ_posn(*value)
        
        self._posn = astro.B1950_to_J2000(value)
        
    @property
    def j2000_ecl(self):
        """
        Position formatted as J2000 ecliptic coordinates.
        Value is object of type astro.ecl_posn.
        """
        
        return astro.get_ecl_from_equ(self._posn, astro.J2000_UTC_JD)
        
    @j2000_ecl.setter
    def j2000_ecl(self, value):
        if not isinstance(value, (astro.ecl_posn, SequenceABC)):
            raise TypeError("value must be type astro.ecl_posn or sequence of length 2")
        if isinstance(value, SequenceABC):
            if len(value) != 2:
                raise TypeError("value sequence must be length 2")
            value = astro.ecl_posn(*value)
            
        self._posn = astro.get_equ_from_ecl(value, astro.J2000_UTC_JD)
        
    def apparent_ecl(self, time_):
        """
        Return position formatted as apparent ecliptic coordinates.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return alue is object of type astro.ecl_posn.
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be of type transform.Time")
            
        equ = astro.get_apparent_posn(self._posn, time_.utc_jd)
        return astro.get_ecl_from_equ(equ, time_.utc_jd)
        
    @property
    def j2000_gal(self):
        """
        Position formatted as J2000 galactic coordinates.
        Value is object of type astro.gal_posn.
        """
        
        return astro.get_gal_from_equ(self._posn)
        
    @j2000_gal.setter
    def j2000_gal(self, value):
        if not isinstance(value, (astro.gal_posn, SequenceABC)):
            raise TypeError("value must be type astro.gal_posn or sequence of length 2")
        if isinstance(value, SequenceABC):
            if len(value) != 2:
                raise TypeError("value sequence must be length 2")
            value = astro.gal_posn(*value)
            
        self._posn = astro.get_equ_from_gal(value)


class PlanetaryPosition(SkyPosition):
    """
    Holds a solar, lunar, or planetary position value and can manipulate
    the value in various representations.
    
    The instance methods apparent_equ() and apparent_ecl() may be used to
    get the apparent equatorial or ecliptic coordinates for a particular
    time.
    """
    
    # planetary names
    NAME_SUN      = 'Sun'
    NAME_MOON     = 'Moon'
    NAME_VENUS    = 'Venus'
    NAME_MARS     = 'Mars'
    NAME_JUPITER  = 'Jupiter'
    NAME_SATURN   = 'Saturn'
    
    known_names = (NAME_SUN, NAME_MOON, NAME_VENUS, NAME_MARS, NAME_JUPITER, NAME_SATURN)
    
    def __init__(self, name):
        """
        Create a PlanetaryPosition instance.
        """
        
        # check parameters
        if name not in self.known_names:
            raise ValueError(f"unknown name {name}")
            
        self.name = name
        
        # cache the appropriate ephemeris function
        if name == self.NAME_SUN:
            self._posn_func = astro.get_solar_equ_coords
        elif name == self.NAME_MOON:
            self._posn_func = astro.get_lunar_equ_coords
        elif name == self.NAME_VENUS:
            self._posn_func = astro.get_venus_equ_coords
        elif name == self.NAME_MARS:
            self._posn_func = astro.get_mars_equ_coords
        elif name == self.NAME_JUPITER:
            self._posn_func = astro.get_jupiter_equ_coords
        elif name == self.NAME_SATURN:
            self._posn_func = astro.get_saturn_equ_coords
            
    def __repr__(self):
        return "%s.%s(%s)" % (type(self).__module__, type(self).__name__, repr(self.name))
        
    def apparent_equ(self, time_):
        """
        Position formatted as apparent equatorial coordinates.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return value is object of type astro.equ_posn.
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be of type transform.Time")
        
        return self._posn_func(time_.utc_jd)
        
    def apparent_ecl(self, time_):
        """
        Position formatted as apparent ecliptic coordinates.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return value is object of type astro.ecl_posn.
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be of type transform.Time")
        
        equ = self._posn_func(time_.utc_jd)
        return astro.get_ecl_from_equ(equ, time_.utc_jd)


class GeographicalPosition(object):
    """
    Holds a geographical position valuee and can manipulate the value in
    various representations.
    
    After creation, the celestial coordinates may be accessed in different
    formats and epoch by using the appropriate instance member.  If marked
    with '(S)', the Time value may also be updated by setting the member to 
    a new value.
    
    geo (S)     - longitude and latitude degrees, elevation meters
    ecef (S)    - Earth centered rectilinear coordinates
    """
    
    # format types
    FORMAT_GEO  = 'GEO'
    FORMAT_ECEF = 'ECEF'
    
    known_formats = (FORMAT_GEO, FORMAT_ECEF)
    
    def __init__(self, value, format = FORMAT_GEO, name = ''):
        """
        Create a new GeographicalPosition instance, using 'value' as the
        initial coordinates.
        
        'format' describes the type of 'value':
        * GeographicalPoistion.FORMAT_GEO   - longitude,latitude, and elevation
        * GeographicalPoistion.FORMAT_ECEF  - ECEF rectangular geodetic
        """
        
        # check parameters
        if format not in self.known_formats:
            raise ValueError(f"unknown format {format}") 
            
        self.name = name
        
        # parse init value based on format type
        # value is held internally as type astro.geo_posn
        if format == self.FORMAT_GEO:
            self.geo = value
            
        elif format == self.FORMAT_ECEF:
            self.ecef = value
            
    def __repr__(self):
        return "%s.%s(%s, name=%s)" % (type(self).__module__, type(self).__name__, repr(self._posn), repr(self.name))
        
    @property
    def geo(self):
        """
        Position formatted as geodedic longitude, latitude, elevation.
        Value is object of type astro.geo_posn.
        """
        
        return self._posn
        
    @geo.setter
    def geo(self, value):
        if not isinstance(value, (astro.geo_posn, SequenceABC)):
            raise TypeError("value must be type astro.geo_posn or sequence of length 2/3")
        if isinstance(value, SequenceABC):
            if (len(value) != 2) and (len(value) != 3):
                raise TypeError("value sequence must be length 2/3")
            value = astro.geo_posn(*value)
            
        self._posn = copy.copy(value)
        
    @property
    def ecef(self):
        """
        Position formatted as ECEF rectagular coordinates.
        Value is object of type astro.rect_posn.
        """
        
        return astro.get_rect_from_geo(self._posn)
        
    @ecef.setter
    def ecef(self, value):
        if not isinstance(value, (astro.rect_posn, SequenceABC)):
            raise TypeError("value must be type astro.rect_posn or sequence of length 3")
        if isinstance(value, SequenceABC):
            if len(value) != 3:
                raise TypeError("value sequence must be length 3")
            value = astro.rect_posn(*value)
            
        self._posn = astro.get_geo_from_rect(value)
        
    def sidereal(self, time_):
        """
        Return the apparent sidereal time for this location.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Returns sidereal time as a float in range [0.0, 24.0).
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be type transform.Time")
            
        return astro.get_local_sidereal_time(self._posn.lng, time_.utc_jd)


class PointingDirection(object):
    """
    Represent a target source object pointing direction for a particular ground
    site location.
    
    Each PointingDirection instance pairs a celestial position with a ground
    observation site position.  Each instance contains a 'source' member reference
    of type CelestialPosition or PlanetaryPosition and a 'site' member reference
    of type GeographicalPosition.
    
    The instance methods hrz() and dir_cos() may be called to get the
    pointing direction as horizontal coordinates or direction cosine components
    for a particular time of interest.  The rst() method may be called to get
    the rise, set, and transit ephemeris times.
    """
    
    def __init__(self, source, site):
        """
        Create a new pointing direction instance.
        
        source  - instance of type CelestialPosition or PlanetaryPosition 
                represeting target source location
        site    - instance of GeographicalPosition representing observer 
                location
        """
        
        # save parameters
        self.source = source
        self.site = site
        
    def __repr__(self):
        return "%s.%s(source=%s, site=%s)" % (type(self).__module__, type(self).__name__, repr(self.source), repr(self.site))
        
    def __setattr__(self, name, value):
        # make surce 'source' and 'site' member are correct type
        
        if (name == 'source') and (not isinstance(value, SkyPosition)):
            raise TypeError("'source' must be type SkyPosition")
            
        elif (name == 'site') and (not isinstance(value, GeographicalPosition)):
            raise TypeError("'site' must be type GeographicalPosition")
            
        object.__setattr__(self, name, value)
        
    def hrz(self, time_):
        """
        Return the pointing direction in horizontal coordinates as type 
        astro.hrz_posn.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be type transform.Time")
        
        return astro.get_hrz_from_equ(self.source.apparent_equ(time_), 
            self.site.geo, time_.utc_jd)
            
    def dir_cos(self, time_):
        """
        Return the pointing direction as three direction cosine components.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return value is a tuple (l, m, n) of direction cosine components.
        """
        return astro.dir_cos(self.hrz(time_))
        
    def rst(self, time_):
        """
        Return the rise, set, and transit ephemeris times.
        The 'time_' parameter should be set to a Time instance providing
        the time of interest.
        Return value is an object of type astro.rst_time, or None if
        the source is circumpolar.
        """
        
        if not isinstance(time_, Time):
            raise TypeError("time_ must be type transform.Time")
            
        return astro.get_object_rst(time_.utc_jd, self.site.geo, self.source.apparent_equ(time_))
