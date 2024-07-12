"""
Astronomical utility functions and classes based on libnova library.
"""

import os
import sys
import time
import math
import numpy as np
from urllib.request import urlopen
from urllib.parse import quote_plus
from calendar import timegm
from datetime import datetime
from xml.etree import ElementTree
from functools import total_ordering

from scipy.interpolate import interp1d

from astropy import units as astrounits
from astropy.time import Time as AstroTime
from astropy.coordinates import EarthLocation, SkyCoord, ICRS, FK4, FK5, Galactic, GeocentricTrueEcliptic, PrecessedGeocentric, CartesianRepresentation, AltAz, solar_system_ephemeris, get_body

from lsl.misc import telemetry
telemetry.track_module()

from typing import Optional, Tuple, Union


__version__   = '0.6'
__all__ = ['dms', 'hms', 'date', 'zonedate', 'rst_time', 'hrz_posn', 'equ_posn', 
           'gal_posn', 'rect_posn', 'ecl_posn', 'get_gmtoff', 'date_to_zonedate',
           'zonedate_to_date', 'rad_to_deg', 'deg_to_rad', 'dms_to_rad',
           'dms_to_deg', 'deg_to_dms', 'rad_to_dms', 'hms_to_deg', 'hms_to_rad',
           'deg_to_hms', 'rad_to_hms', 'add_secs_hms', 'add_hms', 'hrz_to_nswe',
           'range_degrees', 'get_julian_day', 'get_julian_local_date', 'get_date',
           'get_day_of_week', 'get_julian_from_sys', 'get_date_from_sys',
           'get_julian_from_timet', 'get_timet_from_julian', 'get_hrz_from_equ',
           'get_equ_from_hrz', 'get_ecl_from_rect', 'get_equ_from_ecl',
           'get_ecl_from_equ', 'get_equ_from_gal', 'get_gal_from_equ',
           'get_apparent_sidereal_time', 'get_mean_sidereal_time',
           'get_angular_separation', 'get_rel_posn_angle', 'get_apparent_posn',
           'get_equ_prec', 'get_equ_prec2', 'get_object_rst',
           'get_solar_equ_coords', 'get_solar_rst', 'get_jupiter_equ_coords',
           'get_jupiter_rst', 'get_saturn_equ_coords', 'get_saturn_rst',
           'get_lunar_equ_coords', 'get_lunar_rst', 'get_venus_equ_coords',
           'get_venus_rst', 'get_mars_equ_coords', 'get_mars_rst', 'sec_to_jd', 
           'jd_to_sec',  'range_hours', 'jd_to_mjd', 'mjd_to_jd', 'leap_secs', 
           'utc_to_tai', 'tai_to_utc', 'taimjd_to_utcjd', 'utcjd_to_taimjd', 
           'unix_to_utcjd', 'unix_to_taimjd', 'utcjd_to_unix', 'taimjd_to_unix', 
           'tai_to_tt', 'tt_to_tai', 'utc_to_tt', 'tt_to_utc', 'tt_to_tdb', 
           'get_tai_from_sys', 'hms_to_sec', 'deg_to_sec',
           'get_local_sidereal_time', 'geo_posn', 'dir_cos', 'get_rect_from_equ',
           'get_equ_from_rect', 'get_geo_from_rect', 'get_rect_from_geo',
           'get_precession', 'B1950_to_J2000', 'J2000_to_B1950', 'resolve_name']
__author__    = 'D. L. Wood'
__maintainer__ = 'Jayce Dowell'


@total_ordering
class dms(object):
    """
    Represents angles in degrees, minutes, seconds.
    
    Public members:
      neg     - True if measured west of GM, south of EQ
                False if measured east of GM, north of EQ.
      degrees - Angle degrees (integer).
      minutes - Angle minutes (integer).
      seconds - Angle seconds (float).
    """
    
    def __init__(self, neg: bool=False, degrees: int=0, minutes: int=0, seconds: float=0.0):
        """
        Create a dms object.
        
        Param: neg      - True if measured west of GM, south of EQ
                          False if measured east of GM, north of EQ.
        Param: degrees  - Angle degrees (integer [0, 359)).
        Param: minutes  - Angle minutes (integer [0, 59]).
        Param: seconds  - Angle seconds (float [0.0, 60.0)).
        """

        if neg is not None:
            self.neg = neg
            
        if degrees is not None:
            if degrees < 0 or degrees > 359:
                raise ValueError(f"degrees parameter range is [0, 359], is set to {degrees}")
            self.degrees = degrees
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError(f"minutes parameter range is [0, 59], is set to {minutes}")
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError(f"seconds parameter range is [0.0, 60.0), is set to {seconds:0.3f}")
            self.seconds = seconds
            
    def __neg__(self):
        """
        Return negative of object.
        """
        
        if self.neg:
            neg = False
        else:
            neg = True
            
        return dms(neg, self.degrees, self.minutes, self.seconds)
        
    def __str__(self):
        """
        dms object print/str method.
        """
            
        if self.neg:
            sign = '-'
        else:
            sign = '+' 
            
        return f"{sign}{self.degrees:02d} {self.minutes:02d} {self.seconds:05.2f}"
        
    def __repr__(self):
        """
        dms object repr string method
        """
            
        return "%s.%s(%s,%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.neg), repr(self.degrees), repr(self.minutes), repr(self.seconds))
        
    def __reduce__(self):
        """
        dms object pickle reduce method.
        """
            
        return (dms, (self.neg, self.degrees, self.minutes, self.seconds))
        
    def __cmp__(self, other):
        """
        dms comparison tests
        """
            
        if isinstance(other, (float, int)):
            o = other
        elif isinstance(other, (dms, hms)):
            o = other.to_deg()
        else:
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
                
        if self.to_deg() < o:
            return -1
        elif self.to_deg() > o:
            return 1
        else:
            return 0
            
    def __eq__(self, other):
        return True if self.__cmp__(other) == 0 else False
        
    def __lt__(self, other):
        return True if self.__cmp__(other) < 0 else False
        
    def to_deg(self) -> float:
        """
        Convert angles degrees, minutes, seconds to float degrees.
        Returns angle in degrees (float).
        """
        
        return dms_to_deg(self)
        
    def to_rad(self) -> float:
        """
        Convert angles degrees, minutes, seconds to float radians.
        Returns angle in radians (float).
        """
        
        return dms_to_rad(self)
        
    def to_hms(self) -> "hms":
        """
        Convert angles degrees, minutes seconds to hours, minutes, seconds.
        Returns: object of type hms representing angle.
        """
        
        return deg_to_hms(self.to_deg())


@total_ordering
class hms(object):
    """
    Represents times/angles in hours, minutes, seconds.
    
    Public members:
      hours - Angle/time hours (integer).
      minutes - Angle/time minutes (integer).
      seconds - Angle/time seconds (float).
    """
    
    def __init__(self, hours: int=0, minutes: int=0, seconds: float=0.0):
        """
        Create a hms object.
        
        Param: hours    - Angle/time hours (integer [0, 23]).
        Param: minutes  - Angle/time minutes (integer [0, 59]).
        Param: seconds  - Angle/time seconds (float [0.0, 60.0)).
        """
        
        if hours is not None:
            if hours < 0 or hours > 23:
                raise ValueError(f"hours paramerer range is [0, 23], is set to {hours}")
            self.hours = hours
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError(f"minutes paramerer range is [0, 59], is set to {minutes}")
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError(f"seconds paramerer range is [0.0, 60.0), is set to {seconds:0.3f}")
            self.seconds = seconds
            
    def __iadd__(self, other):
        """
        Add hms objects.
        """
        
        return add_hms(other, self) 
        
    def __repr__(self):
        """
        hms object repr string method
        """
        
        return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.hours), repr(self.minutes), repr(self.seconds))           
        
    def __str__(self):	
        """
        hms object print/str method.
        """
        
        return f"{self.hours:02d} {self.minutes:02d} {self.seconds:05.2f}"
        
    def __reduce__(self):
        """
        hms object pickle reduce method.
        """
        
        return (hms, (self.hours, self.minutes, self.seconds))
        
    def __cmp__(self, other):
        """
        hms comparison tests.
        """
        
        if isinstance(other, (float, int)):
            o = other
        elif isinstance(other, (hms, dms)):
            o = other.to_deg()
        else:
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        if self.to_deg() < o:
            return -1
        elif self.to_deg() > o:
            return 1
        else:
            return 0
            
    def __eq__(self, other):
        return True if self.__cmp__(other) == 0 else False
        
    def __lt__(self, other):
        return True if self.__cmp__(other) < 0 else False
        
    def to_deg(self) -> float:
        """
        Convert angles hours, minutes, seconds to float degrees.
        Returns angle in degrees (float).
        """
        
        return hms_to_deg(self)
        
    def to_rad(self) -> float:
        """
        Convert angles hours, minutes, second to float radians.
        Returns angle in radians (float).
        """
        
        return hms_to_rad(self)
        
    def to_dms(self) -> dms:
        """
        Convert angle hours, minutes, seconds to degrees, minutes, seconds.
        Returns: object of type dms representing angle.
        """
        
        return deg_to_dms(self.to_deg())
        
    def to_sec(self) -> float:
        """
        Convert angle hours, minutes, seconds to seconds.
        Returns: time/angle as seconds.
        """
        
        return hms_to_sec(self)


@total_ordering
class date(object):
    """
    Represents UT time in calendar units.
    
    Public members:
      years   - Date years (integer).
      months  - Date months (integer).
      days    - Date days (integer).
      hours   - Date hours (integer).
      minutes - Date minutes (integer).
      seconds - Date seconds (float).
    """
    
    def __init__(self, years: int=2000, months: int=1, days: int=1, hours: int=0, minutes: int=0, seconds: float=0.0):
        """
        Create a date object.
        
        Param: years    - Date years (integer).
        Param: months   - Date months (integer [1, 12]).
        Param: days     - Date days (integer [1, 31]).
        Param: hours    - Date hours (integer [0, 23]).
        Param: minutes  - Date minutes (integer [0, 59]).
        Param: seconds  - Date seconds (float [0.0, 60.0)).
        """
        
        if years is not None:
            self.years = years
            
        if months is not None:
            if months < 1 or months > 12:
                raise ValueError(f"months paramerer range is [1, 12], is set to {months}")
            self.months = months
            
        if days is not None:
            if days < 1 or days > 31:
                raise ValueError(f"days paramerer range is [1, 31], is set to {days}")
            self.days = days
            
        if hours is not None:
            if hours < 0 or hours > 23:
                raise ValueError(f"hours paramerer range is [0, 23], is set to {hours}")
            self.hours = hours
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError(f"minutes paramerer range is [0, 59], is set to {minutes}")
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError(f"seconds paramerer range is [0.0, 60.0), is set to {seconds:0.3f}")
            self.seconds = seconds
            
    def __str__(self):
        """
        date object print/str method.
        """
        
        return f"{self.years:04d}-{self.months:02d}-{self.days:02d} {self.hours:02d}:{self.minutes:02d}:{self.seconds:06.3f}"
        
    def __repr__(self):
        """
        date object repr string method
        """
        
        return "%s.%s(%s,%s,%s,%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.years), repr(self.months), repr(self.days), repr(self.hours), repr(self.minutes), repr(self.seconds))
        
    def __reduce__(self):
        """
        date object pickle reduce method.
        """
        
        return (date, (self.years, self.months, self.days, self.hours, self.minutes, self.seconds))
        
    def __cmp__(self, other):
        """
        date comparison tests.
        """
        
        if not isinstance(other, (date, zonedate)):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        if self.to_jd() < other.to_jd():
            return -1
        elif self.to_jd() > other.to_jd():
            return 1
        else:
            return 0
            
    def __eq__(self, other):
        return True if self.__cmp__(other) == 0 else False
        
    def __lt__(self, other):
        return True if self.__cmp__(other) < 0 else False
        
    def to_zone(self, gmtoff: Optional[int]= None) -> "zonedate":
        """
        Convert UTC calendar time to local calendar time.
        
        Param: gmtoff - local seconds offset from UTC (integer -43200 to 43200).
                    if set to None, the time module offset value is used for
                    current location 
        
        Returns object of type zonedate representing local time.
        """
        
        if gmtoff is None:
            gmtoff = get_gmtoff() 
            
        return date_to_zonedate(self, gmtoff)             
    
    def to_jd(self) -> float:
        """
        Convert calendar time to Julian day.
        Returns UTC time in Julian days (float).
        """
        
        return get_julian_day(self)
        
    def load(self, dateStr: str, timeStr: str):
        """
        Load date object from formatted strings.
        
        Param: dateStr - A string of format YYYY-MM-DD giving date.
        Param: timeStr - A string of format HH:MM:SS.S giving time.
        """
        
        # parse date string
        try:
            (yearStr, monthStr, dayStr) = dateStr.split('-')
        except ValueError:
            raise ValueError(f"date incorrectly formated: {dateStr}") 
            
        try:
            year = int(yearStr)
        except ValueError:
            raise ValueError(f"year incorrectly formated: {yearStr}")
            
        try:
            month = int(monthStr)
        except ValueError:
            raise ValueError(f"month incorrectly formated: {monthStr}")
        if month < 1 or month > 12:
            raise ValueError(f"months paramerer range is [1, 12], is set to {month}")
            
        try:
            day = int(dayStr)
        except ValueError:
            raise ValueError(f"day incorrectly formated: {dayStr}")
        if day < 1 or day > 31:
            raise ValueError(f"days paramerer range is [1, 31], is set to {day}")
            
        # parse time sting 
        try:    
            (hourStr, minuteStr, secondStr) = timeStr.split(':') 
        except ValueError:
            raise ValueError(f"time incorrectly formated: {timeStr}")         
            
        try:
            hour = int(hourStr)
        except ValueError:
            raise ValueError(f"hour incorrectly formated: {hourStr}")
        if hour < 0 or hour > 23:
            raise ValueError(f"hours paramerer range is [0, 23], is set to {hour}")
            
        try:
            minute = int(minuteStr)
        except ValueError:
            raise ValueError(f"minutes incorrectly formated: {minuteStr}")
        if minute < 0 or minute > 59:
            raise ValueError(f"minutes paramerer range is [0, 59], is set to {minute}")
            
        try:
            second = float(secondStr)
        except ValueError:
            raise ValueError(f"seconds incorrectly formated: {secondStr}")
        if second < 0.0 or second >= 60.0:
            raise ValueError(f"seconds paramerer range is [0.0, 60.0), is set to {second:0.3f}")
            
        # set attributes
        self.years = year
        self.months = month
        self.days = day
        self.hours = hour
        self.minutes = minute
        self.seconds = second


@total_ordering
class zonedate(object):
    """
    Represents local time in calendar units.
    
    Public members:
      years   - Date years (integer).
      months  - Date months (integer).
      days    - Date days (integer).
      hours   - Date hours (integer).
      minutes - Date minutes (integer).
      seconds - Date seconds (float).
      gmtoff  - Seconds offset from GM (integer).
    """
    
    def __init__(self, years: int=2000, months: int=1, days: int=1, hours: int=0, minutes: int=0, seconds: float=0.0, gmtoff: int=0):
        """
        Create a zonedate object.
        
        Param: years    - Date years (integer).
        Param: months   - Date months (integer [1, 12]).
        Param: days     - Date days (integer [1, 31]).
        Param: hours    - Date hours (integer [0, 23]).
        Param: minutes  - Date minutes (integer [0, 59]).
        Param: seconds  - Date seconds (float [0.0, 60.0)).
        Param: gmtoff   - Seconds offset from GM (integer [-43200, 43200]).  
        """
        
        if years is not None:
            self.years = years
            
        if months is not None:
            if months < 1 or months > 12:
                raise ValueError(f"months paramerer range is [1, 12], is set to {months}")
            self.months = months
            
        if days is not None:
            if days < 1 or days > 31:
                raise ValueError(f"days paramerer range is [1, 31], is set to {days}")
            self.days = days
            
        if hours is not None:
            if hours < 0 or hours > 23:
                raise ValueError(f"hours paramerer range is [0, 23], is set to {hours}")
            self.hours = hours
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError(f"minutes paramerer range is [0, 59], is set to {minutes}")
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError(f"seconds paramerer range is [0.0, 60.0), is set to {seconds:0.3f}")
            self.seconds = seconds
            
        if gmtoff is None:
            gmtoff = get_gmtoff()
        self.gmtoff = gmtoff
    
    def __str__(self):
        """
        zonedate object str/print method.
        """
        
        return f"{self.years:04d}-{self.months:02d}-{self.days:02d} {self.hours:02d}:{self.minutes:02d}:{self.seconds:06.3f} [{self.gmtoff}]"
        
    def __repr__(self):
        """
        zonedate object repr string method
        """
        
        return "%s.%s(%s,%s,%s,%s,%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.years), repr(self.months), repr(self.days), repr(self.hours), repr(self.minutes), repr(self.seconds), repr(self.gmtoff))
        
    def __reduce__(self):
        """
        zonedate object pickle reduce method.
        """
        
        return (zonedate, (self.years, self.months, self.days, self.hours, self.minutes, self.seconds, self.gmtoff))
        
    def __cmp__(self, other):
        """
        zonedate comparison tests.
        """
        
        if not isinstance(other, (zonedate, date)):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        if self.to_jd() < other.to_jd():
            return -1
        elif self.to_jd() > other.to_jd():
            return 1
        else:
            return 0
            
    def __eq__(self, other):
        return True if self.__cmp__(other) == 0 else False
        
    def __lt__(self, other):
        return True if self.__cmp__(other) < 0 else False
        
    def to_date(self) -> date:
        """
        Convert local calendar time to UTC calendar time.
        Returns object of type date representing UTC time.
        """ 
        
        return zonedate_to_date(self)
        
    def to_jd(self) -> float:
        """
        Convert calendar time to Julian day.
        Returns UTC time in Julian days (float).
        """
        
        return get_julian_local_date(self)


class rst_time(object):
    """
    Represents ephemeris rist, set, and transit times.
    
    Public members:
      rise    - Rise time in UTC Julian days (float).
      set     - Set time in UTC Julian days (float).
      transit - Transit time in UTC Julian days (float).
    """
    
    def __init__(self, rise: Optional[Union[date,float]]=None, set: Optional[Union[date,float]]=None, transit: Optional[Union[date,float]]=None):
        """
        Create a rst_time object.
        
        Param: rise     - Rise time as date object or float UTC Julian days.
        Param: set      - Set time as date object or float UTC Julian days.
        Param: transit  - Transit time as date object or float UTC Julian days.
        """

        if rise is not None:
            if isinstance(rise, date):
                self.rise = rise.to_jd()
            else:
                self.rise = rise
                
        if set is not None:
            if isinstance(set, date):
                self.set = set.to_jd()
            else:
                self.set = set
                
        if transit is not None:
            if isinstance(transit, date):
                self.transit = transit.to_jd()
            else:
                self.transit = transit
                
    def __str__(self):
        """
        rst_time object str/print method.
        """
        
        return f"{self.rise:0.3f} {self.set:0.3f} {self.transit:0.3f}"
        
    def __repr__(self):
        """
        rst_time object repr string method
        """
        
        return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.rise), repr(self.set), repr(self.transit))
        
    def __reduce__(self):
        """
        rst_time object pickle reduce method.
        """
        
        return (rst_time, (self.rise, self.set, self.transit))
        
    def format(self) -> Tuple[date,date,date]:
        """
        Return a tuple (rise, set, transit) where all three are date
        objects representing the ephemeris times.
        """    
        
        return (get_date(self.rise), get_date(self.set), get_date(self.transit))


class hrz_posn(object):
    """
    Represents horizontal local position coordinates.  LWA measures azimuth angle clockwise 
    from north to east, with due north equal to 0 degrees.  Also, method for using zenith 
    angle instead of altitude angle have been added.
    
    Public members:
      az  - Position azimuth angle (float degrees).
      alt - Position altitude angle (float degrees)
    
    Members may also be accessed by subscript:
      hrz_posn[0] = az
      hrz_posn[1] = alt 
    """
    
    _astropy = None
    _distance = None
    
    def __init__(self, az: float=0.0, alt: float=0.0):
        """
        Create a hrz_posn object.
        
        Param: az   - Position azimuth angle (float degrees [0.0, 360.0), 0 = N, 90 = E).
        Param: alt  - Position altitude angle (float degrees [-90.0, 90.0]).
        """

        if az is not None:
            self.az = az
            
        if alt is not None:
            self.alt = alt
            
    @classmethod
    def from_astropy(kls, value: SkyCoord):
        if not isinstance(value, SkyCoord):
            raise TypeError("Expected an object of type SkyCoord")
            
        if value is not None and not isinstance(value.frame, AltAz):
            raise TypeError("Expected a SkyCoord in the frame of AltAz")
            
        _posn = kls()
        _posn._astropy = value
        _posn._az = value.az.deg
        _posn._alt = value.alt.deg
        try:
            _posn._distance = value.distance.to('pc').value
        except astrounits.UnitConversionError:
            pass
        return _posn
        
    @property
    def az(self) -> float:
        """
        Azimiuth in degrees.
        """
        
        return self._az
        
    @az.setter
    def az(self, value: float):
        if self.astropy is not None:
            raise RuntimeError("Cannot set az when this.astropy is not None")
            
        if value < 0.0 or value >= 360.0:
            raise ValueError(f"az paramerer range is [0.0, 360.0), is set to {value:0.3f}")
        self._az = value
        
    @property
    def alt(self) -> float:
        """
        Altitude in degrees.
        """
        
        return self._alt
        
    @alt.setter
    def alt(self, value: float):
        if self.astropy is not None:
            raise RuntimeError("Cannot set alt when this.astropy is not None")
            
        if value < -90.0 or value > 90.0:
            raise ValueError(f"alt paramerer range is [-90.0, 90.0], is set to {value:0.3f}")
        self._alt = value
        
    @property
    def distance(self) -> Optional[float]:
        """
        Distance in pc or None is the distance is unknown.
        """
        
        return self._distance
        
    @property
    def astropy(self) -> SkyCoord:
        return self._astropy
        
    def zen(self, value: Optional[float]=None) -> Union[float,None]:
        """
        If value is None, returns position zenith angle (float degrees 
        [0, 180])  Otherwise, sets the altitude according to the zenith angle
        value.
        """
        
        if value is None:
            return 90.0 - self.alt
        if value <  0.0 or value > 180.0:
            raise ValueError(f"value paramerer range is [0.0, 180.0], is set to {value:0.3f}")
        self.alt = 90.0 - value
        return None
        
    def __str__(self):
        """
        hrz_posn object print/str method.
        """
        
        return f"{self.az:0.3f} {self.alt:0.3f}"
        
    def __repr__(self):
        """
        hrz_posn object repr string method
        """
        
        return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.az), repr(self.alt))
        
    def __reduce__(self):
        """
        hrz_posn object pickle reduce method.
        """
        
        return (hrz_posn, (self.az, self.alt))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.az
        elif key == 1:
            return self.alt
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set index {key} when this.astropy is not None")
            
        if key == 0:
            self.az = value
        elif key == 1:
            self.alt = value
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __eq__(self, other):
        """
        hrz_posn equality test.
        """
        
        if not isinstance(other, hrz_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.az == other.az) and (self.alt == other.alt)
        
    def __ne__(self, other):
        """
        hrz_posn non-equality test.
        """
        
        if not isinstance(other, hrz_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.az != other.az) or (self.alt != other.alt)
        
    def to_equ(self, observer, jD):
        """
        Get equatorial/celestial coordinates from local horizontal coordinates.
        
        Param: observer - Object of type equ_posn representing the pointing's
                          position.
        Param: jD       - UTC Julian day (float).
        
        Returns object of type equ_posn representing equatorial position.
        """
        
        return get_equ_from_hrz(self, observer, jD)         
        
    def dir_cos(self) -> Tuple[float,float,float]:
        """
        Get direction cosines from horizontal coordinates.
        
        Returns: A tuple (l,m,n) of float values for the direction cosines.
                l = unit vector in E direction (azimuth = 90)
                m = unit vector in N direction (azimuth = 0)
                n = unit vector in zenith direction (zenith = 0, altitude = 90)
        """
        
        return dir_cos(self)


class equ_posn(object):
    """
    Represents equatoral/celestial position coordinates.
    
    Public members:
      ra  - Position right ascension angle (float degrees).
      dec - Position declination angle (float degrees).
    
    Members may also be accessed by subscript:
      equ_posn[0] = ra
      equ_posn[1] = dec
    """
    
    _astropy = None
    _pm_ra = None
    _pm_dec = None
    _distance = None
    
    def __init__(self, ra: float=0.0, dec: float=0.0):
        """
        Create a equ_posn object.
        
        Param: ra   - Position right ascension angle
                    Object of type hms or float degrees [0.0, 360.0).
        Param: dec  - Position declination angle
                    Object of type dms or float degrees [-90.0, 90.0].
        """
        
        if ra is not None:
            self.ra = ra
            
        if dec is not None:
            self.dec = dec
            
    @classmethod
    def from_astropy(kls, value: Union[SkyCoord,ICRS,FK4,FK5]):
        if not isinstance(value, (SkyCoord, ICRS, FK4, FK5)):
            raise TypeError("Expected an object of type SkyCoord, ICRS, FK4, or FK5")
            
        if isinstance(value, SkyCoord):
            if not isinstance(value.frame, (ICRS, FK4, FK5, PrecessedGeocentric)):
                raise TypeError("Expected a SkyCoord in the frame of ICRS, FK4, FK5, or PrecessedGeocentric")
                
        _posn = kls()
        _posn._astropy = value
        _posn._ra = value.ra.deg
        _posn._dec = value.dec.deg
        try:
            _posn._pm_ra = value.pm_ra_cosdec.to('mas/yr').value / math.cos(value.dec.rad)
            _posn._pm_dec = value.pm_dec.to('mas/yr').value
        except TypeError:
            pass
        try:
            _posn._distance = value.distance.to('pc').value
        except astrounits.UnitConversionError:
            pass
        return _posn
        
    @property
    def ra(self) -> float:
        """
        RA in degrees.
        """
        
        return self._ra
        
    @ra.setter
    def ra(self, value: Union[hms,float]):
        if self.astropy is not None:
            raise RuntimeError("Cannot set ra when this.astropy is not None")
            
        if isinstance(value, hms):
            value = value.to_deg()
        if value < 0.0 or value >= 360.0:
            raise ValueError(f"ra paramerer range is [0.0, 360.0), is set to {value:0.3f}")
        self._ra = value
        
    @property
    def dec(self) -> float:
        """
        Declination in degrees.
        """
        
        return self._dec
        
    @dec.setter
    def dec(self, value: Union[dms,float]):
        if self.astropy is not None:
            raise RuntimeError("Cannot set dec when this.astropy is not None")
            
        if isinstance(value, dms):
            value = value.to_deg()
        if value < -90.0 or value > 90.0:
            raise ValueError(f"dec paramerer range is [-90.0, 90.0], is set to {value:0.3f}")
        self._dec = value
        
    @property
    def pm_ra(self) -> Union[float,None]:
        """
        Proper motion in RA in mas/yr or None if it is unknown.
        """
        
        return self._pm_ra
        
    @property
    def pm_dec(self) -> Union[float,None]:
        """
        Proper motion in declination in mas/yr or None if it is unknown.
        """
        
        return self._pm_dec
        
    @property
    def distance(self) -> Union[float,None]:
        """
        Distance in pc or None if it is unknown.
        """
        
        return self._distance
        
    @property
    def astropy(self) -> Union[SkyCoord,ICRS,FK4,FK5]:
        return self._astropy
        
    def __str__(self):
        """
        equ_posn object str/print method.
        """
        
        return f"{self.ra:0.3f} {self.dec:0.3f}"
        
    def __repr__(self):
        """
        equ_posn object repr string method
        """
        
        return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.ra), repr(self.dec))   
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.ra
        elif key == 1:
            return self.dec
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set index {key} when this.astropy is not None")
            
        if key == 0:
            self.ra = value
        elif key == 1:
            self.dec = value
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __eq__(self, other):
        """
        equ_posn equality test.
        """
        
        if not isinstance(other, equ_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.ra == other.ra) and (self.dec == other.dec)
        
    def __ne__(self, other):
        """
        equ_posn non-equality test.
        """
        
        if not isinstance(other, equ_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.ra != other.ra) or (self.dec != other.dec)
        
    def to_hrz(self, observer: "geo_posn", jD: float) -> hrz_posn:
        """
        Get local horizontal coordinates from equatorial/celestial coordinates.
        
        Param: observer - Object of type geo_posn representing the observer's
                          location on Earth
        Param: jD       - UTC Julian day (float).
        
        Returns object of type hrz_posn representing local position.
        """
        
        return get_hrz_from_equ(self, observer, jD)
        
    def to_gal(self, jD: float) -> "gal_posn":
        """
        Get J2000 galactic coordinates from apparent equatorial coordinates.
        
        Param: jD - UTC Julian day to get position.
        
        Returns object of type gal_posn representing object's galactic 
        position.
        """
        
        equ = get_equ_prec2(self, jD, J2000_UTC_JD)
        return get_gal_from_equ(equ)
        
    def to_ecl(self, jD: float) -> "ecl_posn":
        """
        Get ecliptical coordinates from equatorial coordinates.
        
        Param: jD - UTC Julian day (float).
        
        Returns object of type ecl_posn representing eclipitcal position.
        """
        
        return get_ecl_from_equ(self, jD)
        
    def angular_separation(self, posn: "equ_posn") -> float:
        """
        Get angular separation from equatorial position.
        
        Param: posn - Object of type equ_posn representing body 2 position.
        
        Returns angular separation in degrees (float).
        """
        
        return get_angular_separation(self, posn) 
        
    def format(self) -> Tuple[hms,dms]:
        """
        Return a tuple (ra, dec) where ra is an hms object and
        dec is a dms object representing ra and dec position coordinates.
        """
            
        return (deg_to_hms(self.ra), deg_to_dms(self.dec)) 
        
    def precess(self, jD: float) -> "equ_posn":
        """
        Get position of celestial object accounting for precession.
        The equatorial coordinates should be for the J2000 epoch.
        
        Param: jD - UTC Julian day (float) to measure position.
        
        Returns: Adjusted equatorial position of object as type equ_posn. 
        """
        
        return get_equ_prec(self, jD) 
        
    def __reduce__(self):
        """
        equ_posn object pickle reduce method.
        """
            
        return (equ_posn, (self.ra,self.dec))


class gal_posn(object):
    """
    Represents galactic position coordinates.
    
    Public members:
      l - Position longitude angle (float degrees).
      b - Position latitude angle (float degrees).
    
    Members may also be accessed by subscript:
      gal_posn[0] = l
      gal_posn[1] = b
    """
    
    _astropy = None
    _pm_l = None
    _pm_b = None
    _distance = None
    
    def __init__(self, l: float=0.0, b: float=0.0):
        """
        Create a gal_posn object.
        
        Param: l - Position longitude angle. 
                Object of type dms or float degrees [-360.0, 360.0).
        Param: b - Position latitude angle. 
                Object of type dms or float degrees [-90.0, 90.0].
        """
        
        if l is not None:
            self.l = l
            
        if b is not None:
            self.b = b
            
    @classmethod
    def from_astropy(kls, value: SkyCoord):
        if not isinstance(value, SkyCoord):
            raise TypeError("Expected an object of type SkyCoord")
            
        if value is not None and not isinstance(value.frame, Galactic):
            raise TypeError("Expected a SkyCoord in the frame of Galactic")
            
        _posn = kls()
        _posn._astropy = value
        _posn._l = value.l.deg
        _posn._b = value.b.deg
        try:
            _posn._pm_l = value.pm_l_cosb.to('mas/yr').value / math.cos(value.b.rad)
            _posn._pm_b = value.pm_b.to('mas/yr').value
        except TypeError:
            pass
        try:
            _posn._distance = value.distance.to('pc').value
        except astrounits.UnitConversionError:
            pass
        return _posn
        
    @property
    def l(self) -> float:
        """
        Galactic longitude in degrees.
        """
        
        return self._l
        
    @l.setter
    def l(self, value: Union[hms,float]):
        if self.astropy is not None:
            raise RuntimeError("Cannot set l when this.astropy is not None")
            
        if isinstance(value, hms):
            value = value.to_deg()
        if value < -360.0 or value >= 360.0:
            raise ValueError(f"l paramerer range is [-360.0, 360.0), is set to {value:0.3f}")
        self._l = value
        
    @property
    def b(self) -> float:
        """
        Galactic latitude in degrees.
        """
        
        return self._b
        
    @b.setter
    def b(self, value: Union[dms,float]):
        if self.astropy is not None:
            raise RuntimeError("Cannot set b when this.astropy is not None")
            
        if isinstance(value, dms):
            value = value.to_deg()
        if value < -90.0 or value > 90.0:
            raise ValueError(f"b paramerer range is [-90.0, 90.0], is set to {value:0.3f}")
        self._b = value
        
    @property
    def pm_l(self) -> Union[float,None]:
        """
        Proper motion in Galactic longitude in mas/yr or None if it is unknown.
        """
        
        return self._pm_l
        
    @property
    def pm_b(self) -> Union[float,None]:
        """
        Proper motion in Galactic latitude in mas/yr or None if it is unknown.
        """
        
        return self._pm_b
        
    @property
    def distance(self) -> Union[float,None]:
        """
        Distance in pc or None if it is unknown.
        """
        
        return self._distance
        
    @property
    def astropy(self) -> SkyCoord:
        return self._astropy
        
    def __str__(self):
        """
        gal_posn object print/str method.
        """
        
        return f"{self.l:0.3f} {self.b:0.3f}"
        
    def __repr__(self):
        """
        gal_posn object repr string method
        """
        
        return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.l), repr(self.b))
        
    def __reduce__(self):
        """
        gal_posn object pickle reduce method.
        """
        
        return (gal_posn, (self.l, self.b))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.l
        elif key == 1:
            return self.b
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set index {key} when this.astropy is not None")
            
        if key == 0:
            self.l = value
        elif key == 1:
            self.b = value
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __eq__(self, other):
        """
        gal_posn equality test.
        """
        
        if not isinstance(other, gal_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")        
            
        return (self.l == other.l) and (self.b == other.b)
        
    def __ne__(self, other):
        """
        gal_posn non-equality test.
        """
        
        if not isinstance(other, gal_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.l != other.l) or (self.b != other.b)
        
    def to_equ(self, jD: float) -> equ_posn:
        """
        Get apparent equatorial coordinates from J2000 galactic coordinates.
        
        Param: jD - UTC Julian day to get position.
        
        Returns object of type equ_posn representing object's apparent
        equatorial position.
        """
        
        equ = get_equ_from_gal(self)
        return get_equ_prec(equ, jD)
        
    def format(self) -> Tuple[dms,dms]:
        """
        Return a tuple (lng, lat) where lng is an dms object and
        lat is a dms object representing galactic longitude and latitude
        position coordinates.
        """
        
        return (deg_to_dms(self.l), deg_to_dms(self.b))


class rect_posn(object):
    """
    Represents rectangular/Cartesian position coordinates.
    
    Public members:
      X - Position X coordinate (float).
      Y - Position Y coordinate (float).
      Z - Position Z coordinate (float).
    
    Members may also be accessed by subscript:
      rect_posn[0] = X
      rect_posn[1] = Y
      rect_posn[2] = Z
    """
    
    def __init__(self, X: float=0.0, Y: float=0.0, Z: float=0.0): 
        """
        Create a rect_posn object
        Param: X - Position X coordinate (float).
        Param: Y - Position Y coordinate (float).
        Param: Z - Position Z coordinate (float).
        """
        
        if X is not None:
            self.X = X
        if Y is not None:
            self.Y = Y
        if Z is not None:
            self.Z = Z
            
    def __str__(self):
        """
        rect_posn object str/print method.
        """
        
        return f"{self.X:0.3f} {self.Y:0.3f} {self.Z:0.3f}"
        
    def __repr__(self):
        """
        rect_posn object repr string method
        """
        
        return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.X), repr(self.Y), repr(self.Z))
        
    def __reduce__(self):
        """
        rect_posn object pickle reduce method.
        """
        
        return (rect_posn, (self.X, self.Y, self.Z))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.X
        elif key == 1:
            return self.Y
        elif key == 2:
            return self.Z
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if key == 0:
            self.X = value
        elif key == 1:
            self.Y = value
        elif key == 2:
            self.Z = value
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __eq__(self, other):
        """
        rect_posn equality test.
        """
        
        if not isinstance(other, rect_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.X == other.X) and (self.Y == other.Y) and (self.Z == other.Z)
        
        
    def __ne__(self, other):
        """
        rect_posn non-equality test.
        """
        
        if not isinstance(other, rect_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.X != other.X) or (self.Y != other.Y) or (self.Z != other.Z)


class ecl_posn(object):
    """
    Represents position as ecliptic longitude and latitude.
    
    Public members:
      lng - Position longitude coordinate (float degrees).
      lat - Position latitude coordinate (float degrees).
    
    Members may also be accessed by subscript:
      ecl_posn[0] = lng
      ecl_posn[1] = lat
    """
    
    _astropy = None
    _pm_lng = None
    _pm_lat = None
    _distance = None
    
    def __init__(self, lng: float=0.0, lat: float=0.0):
        if lng is not None:
            self.lng = lng
            
        if lat is not None:
            self.lat = lat
            
    @classmethod
    def from_astropy(kls, value: SkyCoord):
        if not isinstance(value, SkyCoord):
            raise TypeError("Expected an object of type SkyCoord")
            
        if value is not None and not isinstance(value.frame, GeocentricTrueEcliptic):
            raise TypeError("Expected a SkyCoord in the frame of GeocentricTrueEcliptic")
            
        _posn = kls()
        _posn._astropy = value
        _posn._lng = value.lon.deg
        _posn._lat = value.lat.deg
        try:
            _posn._pm_lng = value.pm_lon_coslat.to('mas/yr').value / math.cos(value.lat.rad)
            _posn._pm_lat = value.pm_lat.to('mas/yr').value
        except TypeError:
            pass
        try:
            _posn._distance = value.distance.to('pc').value
        except astrounits.UnitConversionError:
            pass
        return _posn
        
    @property
    def lng(self) -> float:
        """
        Ecliptic longitude in degees.
        """
        
        return self._lng
        
    @lng.setter
    def lng(self, value: Union[dms,float]):
        if self.astropy is not None:
            raise RuntimeError("Cannot set lng when this.astropy is not None")
            
        if isinstance(value, dms):
            value = value.to_deg()
        if value < -360.0 or value > 360.0:
            raise ValueError(f"lng parameter range is [-360.0, 360.0), is set to {value:0.3f}")
        self._lng = value
        
    @property
    def lat(self) -> float:
        """
        Ecliptic latitude in degrees.
        """
        
        return self._lat
        
    @lat.setter
    def lat(self, value: Union[dms,float]):
        if self.astropy is not None:
            raise RuntimeError("Cannot set lat when this.astropy is not None")
            
        if isinstance(value, dms):
            value = value.to_deg()
        if value < -90.0 or value > 90.0:
            raise ValueError(f"lat parameter range is [-90.0, 90.0], is set to {value:0.3f}")
        self._lat = value
        
    @property
    def pm_lng(self) -> Union[float,None]:
        """
        Proper motion in Ecliptic longitude in mas/yr or None if it is unknown.
        """
        
        return self._pm_lng
        
    @property
    def pm_lat(self) -> Union[float,None]:
        """
        Proper motion in Ecliptic latitude in mas/yr or None if it is unknown.
        """
        
        return self._pm_lat
        
    @property
    def distance(self) -> Union[float,None]:
        """
        Distance in pc or None if it is unknown.
        """
        
        return self._distance
        
    @property
    def astropy(self) -> SkyCoord:
        return self._astropy
        
    def __str__(self):
        """
        ecl_posn object print/str method.
        """
        
        return f"{self.lng:0.3f} {self.lat:0.3f}"
        
    def __repr__(self):
        """
        ecl_posn object repr string method
        """
        
        return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.lng), repr(self.lat))
        
    def __reduce__(self):
        """
        ecl_posn object pickle reduce method.
        """
        
        return (ecl_posn, (self.lng, self.lat))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.lng
        elif key == 1:
            return self.lat
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set index {key} when this.astropy is not None")
            
        if key == 0:
            self.lng = value
        elif key == 1:
            self.lat = value
        else:
            raise ValueError(f"subscript {key} out of range")
            
    def to_equ(self, jD: float) -> equ_posn:
        """
        Get equatorial coordinates from ecliptical coordinates for a given time.
        
        Param: jD - UTC Julian day (float). 
        
        Returns object of type equ_posn representing equatorial position.
        """
        
        return get_equ_from_ecl(self, jD)
        
    def format(self) -> Tuple[dms,dms]:
        """
        Return a tuple (lng, lat) where lng is an dms object and
        lat is a dms object representing longitude and latitude
        position coordinates.
        """
        
        return (deg_to_dms(self.lng), deg_to_dms(self.lat))


######################################################################
# time helper python fucntions
######################################################################

def get_gmtoff() -> int:
    """
    Get local UTC offset based on Python time module and system info.
    """ 
    
    ts = time.localtime()
    if ts[8]:
        return -time.altzone
    return -time.timezone


######################################################################
# General Conversion Functions
######################################################################

def date_to_zonedate(date: date, gmtoff: int) -> zonedate:
    """
    Convert UTC calendar time to local calendar time.
    
    Param: date   - A date object representing UTC time.
    Param: gmtoff - Seconds offset from UTC (integer -43200 to 43200).
    
    Returns object of type zonedate representing local time.
    """
    
    t0 = time.strptime(str(date), "%Y-%m-%d %H:%M:%S.%f")
    fracSec = date.seconds - int(date.seconds)
    t1 = timegm(t0) + gmtoff
    years, months, days, hours, minutes, seconds, wday, yday, dst = time.gmtime(t1)
    
    _zdate = zonedate()
    _zdate.years = years
    _zdate.months = months
    _zdate.days = days
    _zdate.hours = hours
    _zdate.minutes = minutes
    _zdate.seconds = seconds + fracSec
    _zdate.gmtoff = 1*gmtoff
    return _zdate


def zonedate_to_date(zonedate: zonedate) -> date:
    """
    Convert local calendar time to UTC calendar time.
    
    Param: zonedate - Object of type zonedate representing local time.  
    
    Returns object of type date representing UTC time.
    """
    
    t0, junk = str(zonedate).rsplit(None, 1)
    tt0 = time.strptime(t0, "%Y-%m-%d %H:%M:%S.%f")
    fracSec = zonedate.seconds - int(zonedate.seconds)
    tt1 = timegm(tt0) - zonedate.gmtoff
    years, months, days, hours, minutes, seconds, wday, yday, dst = time.gmtime(tt1)
    
    _date = date()
    _date.years = years
    _date.months = months
    _date.days = days
    _date.hours = hours
    _date.minutes = minutes
    _date.seconds = seconds + fracSec
    return _date


def rad_to_deg(radians: float) -> float:
    """
    Convert radians to degrees.
    
    Param: radians - Angle in radians (float).
    
    Returns angle in degrees (float).
    """
    
    return radians * 180.0/math.pi


def deg_to_rad(degrees: float) -> float:
    """
    Convert degres to radians.
    
    Param: degrees - Angle in degrees (float).
    
    Returns angle in radians (float).
    """
    
    return degrees * math.pi/180.0


def dms_to_rad(dms: dms) -> float:
    """
    Convert angles degrees, minutes, seconds to radians.
    
    Param: dms - Object of type dms representing angle.
    
    Returns angle in radians (float).
    """
    
    degrees = dms_to_deg(dms)
    return deg_to_rad(degrees)


def dms_to_deg(dms: dms) -> float:
    """
    Convert angles degrees, minutes, seconds to float degrees.
    
    Param: dms - Object of type dms representing angle.
    
    Returns angle in degrees (float).
    """
    
    degrees = dms.degrees + dms.minutes/60.0 + dms.seconds/3600.0
    if dms.neg:
        degrees *= -1
    return degrees


def _float_to_sexa(value):
    sgn = 1.0
    if value < 0:
        sgn = -1.0
        value *= -1
        
    d = int(value)
    m = int(value * 60) % 60
    s = (value * 3600.0) % 60.0
    
    return sgn, d, m, s


def deg_to_dms(degrees: float) -> dms:
    """
    Convert angles float degrees to degrees, minutes, seconds.
    
    Param: degrees - Angle in degrees (float). 
    
    Returns object of type dms representing angle.
    """
    
    _dms = dms()
    
    sgn,d,m,s = _float_to_sexa(degrees)
    
    _dms.degrees = d
    _dms.minutes = m
    _dms.seconds = s
    _dms.neg = True if sgn < 0 else False
    
    return _dms


def rad_to_dms(radians: float) -> dms:
    """
    Convert angles float radians to degrees, minutes, seconds.
    
    Param: radians - Angle in radians (float). 
    
    Returns object of type dms representing angle.
    """
    
    degrees = rad_to_deg(radians)
    return deg_to_dms(degrees)


def hms_to_deg(hms: hms) -> float:
    """
    Convert angles hours, minutes, seconds to float degrees.
    
    Param: hms - Object of type hms representing angle.
    
    Returns angle in degrees (float).
    """
    
    hours = hms.hours + hms.minutes/60.0 + hms.seconds/3600.0
    degrees = hours * 15.0
    return degrees


def hms_to_rad(hms: hms) -> float:
    """
    Convert angles hours, minutes, seconds to float radians.
    
    Param: hms - Object of type hms representing angle.
    
    Returns angle in radians (float).
    """
    
    degrees = hms_to_deg(hms)
    return deg_to_rad(degrees)


def deg_to_hms(degrees: float) -> hms:
    """
    Convert angles float degrees to hours, minutes, seconds.
    
    Param: degrees - Angle in degrees (float). 
    
    Returns object of type hms representing angle.
    """
    
    sgn,h,m,s = _float_to_sexa(degrees / 15.0)
    return hms(h, m, s)


def rad_to_hms(radians: float) -> hms:
    """
    Convert angles float radians to hours, minutes, seconds.
    
    Param: radians - Angle in radians (float). 
    
    Returns object of type hms representing angle.
    """
    
    degrees = rad_to_deg(radians)
    return deg_to_hms(degrees)


def add_secs_hms(hms: hms, seconds: float) -> hms:
    """
    Add seconds to time/angle hours, minutes, seconds.
    
    Param: hms      - Object of type hms representing angle.
    Param: seconds  - Seconds offset (float) to add to angle. 
    
    Returns object of type hms representing angle + offset.
    """
    
    degrees = hms_to_deg(hms)
    degrees += (seconds/3600.0) * 15.0
    return deg_to_hms(degrees)


def add_hms(source: hms, dest: hms) -> hms:
    """
    Adds time/angle hours, minutes, seconds.
    
    Param: source - Object of type hms represeting angle 1.
    Param: dest   - Object of type hms represeting angle 2.
    
    Returns object of type hms representing sum of angles.
    """
    
    degrees1 = hms_to_deg(source)
    degrees2 = hms_to_deg(dest)
    degrees1 += degrees2
    _hms = deg_to_hms(degrees1)
    
    dest.hours   = 1*_hms.hours
    dest.minutes = 1*_hms.minutes
    dest.seconds = 1*_hms.seconds
    
    return dest


def hrz_to_nswe(pos: hrz_posn) -> str:
    """
    Get cardinal/ordinal azimuth direction.
    
    Param: pos - Object of type hrz_posn giving local position.
    
    Returns string giving direction.
    """
    
    az = pos.az % 360.0
    if az < 45 or az >= 315:
        return 'N'
    elif az >= 45 and az < 135:
        return 'E'
    elif az >= 135 and az < 225:
        return 'S'
    else:
        return 'W'


def range_degrees(degrees: float) -> float:
    """
    Put angle into range [0, 360].
    
    Param: degrees - large angle (float degrees)
    
    Returns: angle in range (float degrees)
    """
    
    return degrees % 360.0


######################################################################
# General Calendar Functions
######################################################################

def get_julian_day(date: date) -> float:
    """
    Convert calendar time to Julian day.
    
    Param: date - Object of type date representing UTC time.
    
    Returns UTC time in Julian days (float).
    """
    
    s = int(date.seconds)
    us = int((date.seconds - s)*1e6)
    dt = datetime(date.years, date.months, date.days, date.hours, date.minutes, s, us)
    d = AstroTime(dt, format='datetime', scale='utc')
    return d.jd
    
    
def get_julian_local_date(zonedate: zonedate) -> float:
    """
    Convert local calendar time to Julian day.
    
    Param: zonedate - Object of type zonedate representing local time.
    
    Returns UTC time in Julian days (float).
    """
    
    _date = zonedate_to_date(zonedate)
    return get_julian_day(_date)


def get_date(jD: float) -> date:
    """
    Convert Julian day to calendar time.
    
    Param: jD - UTC time in Julian days (float).
    
    Returns object of type date representing UTC time.
    """
    
    _date = date()
    d = AstroTime(jD, format='jd', scale='utc')
    
    dt = d.datetime
    _date.years = dt.year
    _date.months = dt.month
    _date.days = dt.day
    _date.hours = dt.hour
    _date.minutes = dt.minute
    _date.seconds = dt.second + dt.microsecond/1e6
    return _date


def get_day_of_week(date: date) -> int:
    """
    Gets day of week from calendar time.
    
    Param: date - Object of type date representing UTC time.
    
    Returns day of week (0 = Sunday, 6 = Saturday).
    """
    
    jd = date.to_jd()
    return (int(round(jd)) + 1) % 7


def get_julian_from_sys() -> float:
    """
    Returns UTC Julian day (float) from system clock.
    """
    
    t0 = AstroTime.now()
    return t0.utc.jd


def get_date_from_sys() -> date:
    """
    Gets calendar time from system clock.
    
    Returns object of type date representing UTC time.
    """
    
    jD = get_julian_from_sys()
    _date = get_date(jD)
    return _date


def get_julian_from_timet(timet: float) -> float:
    """
    Gets Julian day from Unix time.
    
    Param: timet - Unix timet in seconds (float)
    
    Returns UTC Julian day (float).
    """
    
    t = AstroTime(timet, format='unix', scale='utc')
    return t.jd


def get_timet_from_julian(jD: float) -> float:
    """
    Gets Unix time from Julian day.
    
    Param: jD - UTC Julian day (float).
    
    Returns Unix timet in seconds (float).
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    return t.unix


######################################################################
# Position utility classes
######################################################################

class geo_posn(object):
    """
    Class to represent geographical position in terms of longitude, lattitude,
    and elevation.  This is a set of geodetic coordinates based on the WGS84 model.
    
    Public members:
      lng - longitude (float)
      lat - latitude (float)
      elv - elevation (float)
    
    Members may also be accessed by subscript:
      geo_posn[0] = lng
      geo_posn[1] = lat
      geo_posn[2] = elv
    """
    
    _astropy = None
    
    def __init__(self, lng: float=0.0, lat: float=0.0, elv: float=0.0):
        """
        Create a geo_posn object.
        
        Param: lng - longitude coordinate
                    Object of type dms or float degrees (-360.0, 360.0).
        Param: lat - Position latitude coordinate
                    Object of type dms or float degrees [-90.0, 90.0].
        Param: elv - elevation (float meters)
        """
        
        if lng is not None:
            self.lng = lng
            
        if lat is not None:
            self.lat = lat
            
        if elv is not None:
            self.elv = elv
            
    @classmethod
    def from_astropy(kls, value: EarthLocation):
        if not isinstance(value, EarthLocation):
            raise TypeError("Expected an object of type EarthLocation")
            
        _posn = kls()
        _posn._astropy = value
        _posn._lng = value.lon.wrap_at(360*astrounits.deg).deg
        _posn._lat = value.lat.deg
        _posn._elv = value.height.to('m').value
        return _posn
        
    @property
    def lng(self) -> float:
        return self._lng
        
    @lng.setter
    def lng(self, value: Union[dms,float]):
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set lng when this.astropy is not None")
            
        if isinstance(value, dms):
            value = value.to_deg()
        if value <= -360.0 or value >= 360.0:
            raise ValueError(f"lng parameter range is (-360.0, 360.0), is set to {value:0.3f}")
        self._lng = value
        
    @property
    def lat(self) -> float:
        return self._lat
        
    @lat.setter
    def lat(self, value: Union[dms,float]):
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set lat when this.astropy is not None")
            
        if isinstance(value, dms):
            value = value.to_deg()
        if value < -90.0 or value > 90.0:
            raise ValueError(f"lat paramerer range is [-90.0, 90.0], is set to {value:0.3f}")
        self._lat = value
    
    @property
    def elv(self) -> float:
        return self._elv
        
    @elv.setter
    def elv(self, value: float):
        if self.astropy is not None:
            raise RuntimeError(f"Cannot elv when this.astropy is not None")
            
        self._elv = value
        
    @property
    def astropy(self) -> EarthLocation:
        return self._astropy
        
    def __str__(self):
        """
        geo_posn object print/str method.
        """
        
        return f"{self.lng:0.3f} {self.lat:0.3f} {self.elv:0.3f}"
        
    def __repr__(self):
        """
        geo_posn object repr string method
        """
        
        return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.lng), repr(self.lat), repr(self.elv))
            
    def __reduce__(self):
        """
        geo_posn object pickle reduce method.
        """
        
        return (geo_posn, (self.lng, self.lat, self.elv))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.lng
        elif key == 1:
            return self.lat
        elif key == 2:
            return self.elv
        else:
            raise ValueError(f"subscript {key} out of range")
                
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if self.astropy is not None:
            raise RuntimeError(f"Cannot set index {key} when this.astropy is not None")
            
        if key == 0:
            self.lng = value
        elif key == 1:
            self.lat = value
        elif key == 2:
            self.elv = value
        else:
            raise ValueError(f"subscript {key} out of range")
        
    def __eq__(self, other):
        """
        geo_posn equality test.
        """
        
        if not isinstance(other, geo_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.lng == other.lng) and (self.lat == other.lat) and (self.elv == other.elv)
            
    def __ne__(self, other):
        """
        geo_posn non-equality test.
        """
        
        if not isinstance(other, geo_posn):
            raise TypeError(f"comparison not supported for type {type(other).__name__}")
            
        return (self.lng != other.lng) or (self.lat != other.lat) or (self.elv != other.elv)
        
    def format(self) -> Tuple[dms,dms,float]:
        """
        Return a tuple (lng, lat) where lng is an dms object and
        lat is a dms object representing longitude and latitude
        position coordinates.
        """
        
        return (deg_to_dms(self.lng), deg_to_dms(self.lat), self.elv)


######################################################################
# Transformation of Coordinates Functions
######################################################################

def get_hrz_from_equ(target: equ_posn, observer: geo_posn, jD: float) -> hrz_posn:
    """
    Get local horizontal coordinates from equatorial/celestial coordinates.
    
    Param: target   - Object of type equ_posn representing celestial position.
    Param: observer - Object of type geo_posn representing observer position.
    Param: jD       - UTC Julian day (float).
    
    Returns object of type hrz_posn representing local position.
    """
    
    try:
        elv = observer.elv
    except AttributeError:
        elv = 0.0
        
    t = AstroTime(jD, format='jd', scale='utc')
    el = EarthLocation.from_geodetic(observer.lng*astrounits.deg, observer.lat*astrounits.deg,
                                     height=elv*astrounits.m,
                                     ellipsoid='WGS84')
    sc = target.astropy
    if sc is None:
        sc = SkyCoord(target.ra*astrounits.deg, target.dec*astrounits.deg,
                      frame='fk5', equinox='J2000')
    aa = AltAz(location=el, obstime=t)
    sc = sc.transform_to(aa)
    
    return hrz_posn.from_astropy(sc)


def get_equ_from_hrz(target: hrz_posn, observer: geo_posn, jD: float) -> equ_posn:
    """
    Get equatorial/celestial coordinates from local horizontal coordinates.
    
    Param: target   - Object of type hrz_posn representing horizontal position.
    Param: observer - Object of type geo_posn representing observer position.
    Param: jD       - UTC Julian day (float).
    
    Returns object of type equ_posn representing a equatorial position in the
    astropy.coordinates.FK5 frame with equinox=J2000.
    """
    
    try:
        elv = observer.elv
    except AttributeError:
        elv = 0.0
        
    t = AstroTime(jD, format='jd', scale='utc')
    el = EarthLocation.from_geodetic(observer.lng*astrounits.deg, observer.lat*astrounits.deg,
                                     height=elv*astrounits.m,
                                     ellipsoid='WGS84')
    aa = AltAz(target.az*astrounits.deg, target.alt*astrounits.deg,
               location=el, obstime=t)
    sc = aa.transform_to(FK5(equinox='J2000'))
    
    return equ_posn.from_astropy(sc)


def get_ecl_from_rect(rect: rect_posn) -> ecl_posn:
    """
    Get ecliptical coordinates from rectangular coordinates.
    
    Param: rect - Object of type rect_posn representing position.
    
    Returns object of type ecl_posn representing ecliptical position.
    """
    
    x = rect.X
    y = rect.Y
    z = rect.Z
    
    sc = GeocentricTrueEcliptic(CartesianRepresentation(x*astrounits.au, y*astrounits.au, z*astrounits.au),
                                equinox='J2000')
    
    return ecl_posn.from_astropy(sc)


def get_equ_from_ecl(target: ecl_posn, jD: float) -> equ_posn:
    """
    Get J2000 equatorial coordinates from ecliptical coordinates for a given
    time.
    
    Param: target   - Object of type ecl_posn representing ecliptic position.
    Param: jD       - UTC Julian day (float). 
    
    Returns object of type equ_posn representing a equatorial position in the
    astropy.coordinates.FK5 frame with equinox=J2000.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    sc = target.astropy
    if sc is None:
        sc = GeocentricTrueEcliptic(target.lng*astrounits.deg, target.lat*astrounits.deg,
                                    equinox='J2000', obstime=t)
    sc = sc.transform_to(FK5(equinox='J2000'))
    
    return equ_posn.from_astropy(sc)


def get_ecl_from_equ(target: equ_posn, jD: float) -> ecl_posn:
    """
    Get ecliptical coordinates from J2000 equatorial coordinates for a given
    time.
    
    Param: target  - Object of type equ_posn representing a J2000 equatorial
                     position.
    Param: jD      - UTC Julian day (float). 
    
    Returns object of type ecl_posn representing ecliptic position in the
    astropy.coordinates.GeocentricTrueEcliptic frame with equinox=J2000 and
    obstime=jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    sc = target.astropy
    if sc is None:
        sc = SkyCoord(target.ra*astrounits.deg, target.dec*astrounits.deg,
                      frame='fk5', equinox='J2000')
    sc = sc.transform_to(GeocentricTrueEcliptic(equinox='J2000', obstime=t))
    
    return ecl_posn.from_astropy(sc)


def get_equ_from_gal(target: gal_posn) -> equ_posn:
    """
    Get J2000 equatorial coordinates from galactic coordinates.
    
    Param: target - Object of type gal_posn representing galactic position.
    
    Returns object of type equ_posn representing a equatorial position in the
    astropy.coordinates.FK5 frame with equinox=J2000.
    
    .. versionchanged:: 3.0.0
      This function now expects J2000 coordinates
    """
    
    sc = target.astropy
    if sc is None:
        sc = SkyCoord(target.l*astrounits.deg, target.b*astrounits.deg,
                      frame='galactic')
    sc = sc.transform_to(FK5(equinox='J2000'))
    
    return equ_posn.from_astropy(sc)


def get_gal_from_equ(target: equ_posn) -> gal_posn:
    """
    Get galactic coordinates from J2000 equatorial coordinates.
    
    Param: target - Object of type equ_posn representing J2000 equatorial 
                    position.
    
    Returns object of type gal_posn representing galactic position in the
    astropy.coordinates.Galactic frame.
    
    .. versionchanged:: 3.0.0
      This function now expects J2000 coordinates
    """
    
    sc = target.astropy
    if sc is None:
        sc = SkyCoord(target.ra*astrounits.deg, target.dec*astrounits.deg,
                      frame='fk5', equinox='J2000')
    sc = sc.transform_to(Galactic())
    
    return gal_posn.from_astropy(sc)


######################################################################
# Sidereal Time Functions
######################################################################

def get_apparent_sidereal_time(jD: float) -> float:
    """
    Get apparent sidereal time from Julian day.
    
    Param: jD - UTC Julian day (float).
    
    Returns GM apparent sidereal time (float hours).
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    gast = t.sidereal_time('apparent', longitude=0*astrounits.deg)
    
    return gast.hourangle


def get_mean_sidereal_time(jD: float) -> float:
    """
    Get mean sidereal time from Julian day.
    
    Param: jD - UTC Julian day (float).
    
    Returns GM mean sidereal time (float hours).
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    gmst = t.sidereal_time('mean', longitude=0*astrounits.deg)
    
    return gmst.hourangle



######################################################################
# Angular Separation Functions
######################################################################

def get_angular_separation(posn1: equ_posn, posn2: equ_posn) -> float:
    """  
    Get angular separation from equatorial positions.
    
    Param: posn1 - Object of type equ_posn representing body 1 position.
    Param: posn2 - Object of type equ_posn representing body 2 position.
    
    Returns angular separation in degrees (float).
    """
    
    sc1 = posn1.astropy
    if sc1 is None:
        sc1 = SkyCoord(posn1.ra*astrounits.deg, posn1.dec*astrounits.deg,
                       format='fk5', equinox='J2000')
    sc2 = posn2.astropy
    if sc2 is None:
        sc2 = SkyCoord(posn2.ra*astrounits.deg, posn2.dec*astrounits.deg,
                       format='fk5', equinox='J2000')
        
    sep = sc1.separation(sc2)
    
    return sep.deg


def get_rel_posn_angle(posn1: equ_posn, posn2: equ_posn):
    """
    Get relative position angle from equatorial positions.
    
    Param: posn1 - Object of type equ_posn representing body 1 position.
    Param: posn2 - Object of type equ_posn representing body 2 position.
    
    Returns position angle in degrees (float).
    """
    
    sc1 = posn1.astropy
    if sc1 is None:
        sc1 = SkyCoord(posn1.ra*astrounits.deg, posn1.dec*astrounits.deg,
                       format='fk5', equinox='J2000')
    sc2 = posn2.astropy
    if sc2 is None:
        sc2 = SkyCoord(posn2.ra*astrounits.deg, posn2.dec*astrounits.deg,
                       format='fk5', equinox='J2000')
        
    pa = sc1.position_angle(sc2)
    
    return pa.deg


######################################################################
# Apparent Position Functions
######################################################################

_DEFAULT_PROPER_MOTION = equ_posn(0.0, 0.0)

def get_apparent_posn(mean_position: equ_posn, jD: float, proper_motion: Optional[Tuple[Union[float,None],Union[float,None]]]=None) -> equ_posn:
    """
    Get apparent position of celestial object accounting for precession, nutation, 
    aberration, and optionally proper motion.
    
    Param: mean_position  - J2000 equatorial mean position of object as type 
                            equ_posn.
    Param: jD             - UTC Julian day (float) to measure position.
    Param: proper_motion  - object of type equ_posn giving object's proper motion
                            in mas/yr (optional).
    
    Returns: Apparent equatorial position in the
             astropy.coordinates.PrecessedGeocentric frame (equinox = jD; 
             obstime = jD) of object as type equ_posn.
    """
    
    if proper_motion is None:
        proper_motion = (mean_position.pm_ra, mean_position.pm_dec)
    if proper_motion[0] is None:
        proper_motion = (_DEFAULT_PROPER_MOTION[0], proper_motion[1])
    if proper_motion[1] is None:
        proper_motion = (proper_motion[0], _DEFAULT_PROPER_MOTION[1])
            
    t = AstroTime(jD, format='jd', scale='utc')
    sc = mean_position.astropy
    if sc is None:
        sc = SkyCoord(mean_position.ra*astrounits.deg, mean_position.dec*astrounits.deg,
                      pm_ra_cosdec=proper_motion[0]*math.cos(proper_motion[1]/1000/3600*math.pi/180)*astrounits.mas/astrounits.yr,  # type: ignore
                      pm_dec= proper_motion[1]*astrounits.mas/astrounits.yr,
                      frame='fk5', equinox='J2000')
    sc = sc.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
    
    return equ_posn.from_astropy(sc)


######################################################################
# Precession Functions
######################################################################

def get_equ_prec(mean_position: equ_posn, jD: float) -> equ_posn:
    """
    Get position of celestial object accounting for precession.
    Only works for converting to and from J2000 epoch.
    
    Param: mean_position - J2000 equatorial mean position of object as type equ_posn.
    Param: jD - UTC Julian day (float) to measure position.
    
    Returns: Adjusted equatorial position in the astropy.coordinates.FK5 frame
             (equinox=jD) of object as type equ_posn.
    """    
    
    sc = mean_position.astropy
    if sc is None:
        sc = SkyCoord(mean_position.ra*astrounits.deg, mean_position.dec*astrounits.deg,
                      frame='fk5', equinox='J2000')
    t = AstroTime(jD, format='jd', scale='utc')
    sc = sc.transform_to(FK5(equinox=t))
    
    return equ_posn.from_astropy(sc)


def get_equ_prec2(mean_position: equ_posn, fromJD: float, toJD: float) -> equ_posn:
    """
    Get position of celestial object accounting for precession.
    
    Param: mean_position  - equatorial first position of object as type equ_posn.
    Param: fromJD         - UTC Julian day (float) of first time.
    Param: toJD           - UTC Julian day (float) of second time.
    
    Returns: Equatorial position in the astropy.coordinates.FK5 frame of the
             object as type equ_posn converted from time 1 to time 2.
    """  
    
    t1 = AstroTime(fromJD, format='jd', scale='utc')
    sc = mean_position.astropy
    if sc is None:
        sc = SkyCoord(mean_position.ra*astrounits.deg, mean_position.dec*astrounits.deg,
                      frame='fk5', equinox=t1)
    t2 = AstroTime(toJD, format='jd', scale='utc')
    sc = sc.transform_to(FK5(equinox=t2))
    
    return equ_posn.from_astropy(sc)


######################################################################
# Rise, Set, Transit functions
######################################################################

def get_object_rst(jD: float, observer: geo_posn, target: equ_posn) -> Union[rst_time,None]:
    """
    Get rise, set, and transit times of a celstial object.
    
    Param: jD       - UTC Julian day (float) target time.
    Param: observer - object of type geo_posn giving observer position
    Param: target   - object of type equ_posn giving target equatorial position
    
    Returns: Object of type rst_time giving object's ephemeris UTC times,
            or None if the object is circumpolar.
    """

    try:
        elv = observer.elv
    except AttributeError:
        elv = 0.0
        
    _rst = rst_time()
    el = EarthLocation.from_geodetic(observer.lng*astrounits.deg, observer.lat*astrounits.deg,
                                     height=elv*astrounits.m,
                                     ellipsoid='WGS84')
    sc = target.astropy
    if sc is None:
        sc = SkyCoord(target.ra*astrounits.deg, target.dec*astrounits.deg,
                      frame='fk5', equinox='J2000')
        
    # Part 1 - Course +/1 day search for the anti-transits
    t = AstroTime(jD+np.linspace(-1, 1, 192), format='jd', scale='utc')
    aa = AltAz(location=el, obstime=t)
    aa = sc.transform_to(aa)
    antitransits = np.where(np.diff(aa.az[...].deg) < -180)[0]
    if len(antitransits) == 0:
        # No anti-transits = it never rises or never sets
        return None
        
    # Part 2 - Determine the approximate rise, set, and transit times
    #          based on when the source is up.
    above = np.where(aa.alt[antitransits[0]:antitransits[1]] >= 0)[0]
    rise = above[0] + antitransits[0]
    set = above[-1] + antitransits[0]
    transit = (rise + set) // 2
    
    # Part 3 - Refine the rough values found in the course search.
    # Note:  We are trying to emulate ephem.Observer.next_() so we want the next
    #        rising/setting/transit.
    offset = 0
    if aa.obstime[rise].jd < jD:
        offset = 1
    tr = AstroTime(aa.obstime[rise].jd+np.linspace(offset-20/1440, offset+20/1440, 41), format='jd', scale='utc')
    offset = 0
    if aa.obstime[set].jd < jD:
        offset = 1
    ts = AstroTime(aa.obstime[set].jd+np.linspace(offset-20/1440, offset+20/1440, 41), format='jd', scale='utc')
    offset = 0
    if aa.obstime[transit].jd < jD:
        offset = 1
    tt = AstroTime(aa.obstime[transit].jd+np.linspace(offset-20/1440, offset+20/1440, 41), format='jd', scale='utc')
    
    # Part 3a - Rise time via interpolation to find the zero crossing
    aa = AltAz(location=el, obstime=tr)
    aa = sc.transform_to(aa)
    rtf = interp1d(aa.alt, aa.obstime.jd)
    try:
        tr = rtf(0.0)
    except ValueError:
        return None
        
    # Part 3b - Set time via interpolation to find the zero crossing
    aa = AltAz(location=el, obstime=ts)
    aa = sc.transform_to(aa)
    stf = interp1d(aa.alt, aa.obstime.jd)
    try:
        ts = stf(0.0)
    except ValueError:
        return None
        
    # Part 3c - Transit time via fitting a cubic to the altitude as a function
    #           of time.
    aa = AltAz(location=el, obstime=tt)
    aa = sc.transform_to(aa)
    ttf = np.polyfit(aa.obstime.jd-aa.obstime.jd[np.argmax(aa.alt)], aa.alt, 3)
    tt = -ttf[1] - np.sqrt(ttf[1]**2 - 4*ttf[0]*ttf[2])
    tt /= 2*ttf[0]
    tt += aa.obstime.jd[np.argmax(aa.alt)]
    
    _rst.rise = tr
    _rst.transit = tt
    _rst.set = ts
    return _rst


SOLAR_SYSTEM_EPHEMERIS_TO_USE = 'de432s'

def _get_solar_system_rst(jD: float, observer: geo_posn, body: str) -> Union[rst_time,None]:
    """
    Get rise, set, and transit times of a solar system body.
    
    Param: jD       - UTC Julian day (float) target time.
    Param: observer - object of type geo_posn giving observer position
    Param: target   - name of the solar system body
    
    Returns: Object of type rst_time giving object's ephemeris UTC times,
            or None if the object is circumpolar.
    """

    try:
        elv = observer.elv
    except AttributeError:
        elv = 0.0
        
    _rst = rst_time()
    el = EarthLocation.from_geodetic(observer.lng*astrounits.deg, observer.lat*astrounits.deg,
                                     height=elv*astrounits.m,
                                     ellipsoid='WGS84')
        
    # Part 1 - Course +/1 day search for the anti-transits
    t = AstroTime(jD+np.linspace(-1, 1, 192), format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        sc = get_body(body, t, location=el)
        aa = AltAz(location=el, obstime=t)
        aa = sc.transform_to(aa)
    antitransits = np.where(np.diff(aa.az[...].deg) < -180)[0]
    if len(antitransits) == 0:
        # No anti-transits = it never rises or never sets
        return None
        
    # Part 2 - Determine the approximate rise, set, and transit times
    #          based on when the source is up.
    try:
        above = np.where(aa.alt[antitransits[0]:antitransits[1]] >= 0)[0]
    except IndexError:
        above = np.where(aa.alt[antitransits[0]:] >= 0)[0]
    rise = above[0] + antitransits[0]
    set = above[-1] + antitransits[0]
    transit = (rise + set) // 2
    
    # Part 3 - Refine the rough values found in the course search.
    # Note:  We are trying to emulate ephem.Observer.next_() so we want the next
    #        rising/setting/transit.
    scale = 20 if body != 'moon' else 60
    offset = 0
    if aa.obstime[rise].jd < jD:
        offset = 1
    tr = AstroTime(aa.obstime[rise].jd+np.linspace(offset-scale/1440, offset+scale/1440, 2*scale+1), format='jd', scale='utc')
    offset = 0
    if aa.obstime[set].jd < jD:
        offset = 1
    ts = AstroTime(aa.obstime[set].jd+np.linspace(offset-scale/1440, offset+scale/1440, 2*scale+1), format='jd', scale='utc')
    offset = 0
    if aa.obstime[transit].jd < jD:
        offset = 1
    tt = AstroTime(aa.obstime[transit].jd+np.linspace(offset-scale/1440, offset+scale/1440, 2*scale+1), format='jd', scale='utc')
    
    # Part 3a - Rise time via interpolation to find the zero crossing
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        sc = get_body(body, tr, location=el)
        aa = AltAz(location=el, obstime=tr)
        aa = sc.transform_to(aa)
    rtf = interp1d(aa.alt, aa.obstime.jd)
    try:
        tr = rtf(0.0)
    except ValueError:
        return None
        
    # Part 3b - Set time via interpolation to find the zero crossing
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        sc = get_body(body, ts, location=el)
        aa = AltAz(location=el, obstime=ts)
        aa = sc.transform_to(aa)
    stf = interp1d(aa.alt, aa.obstime.jd)
    try:
        ts = stf(0.0)
    except ValueError:
        return None
        
    # Part 3c - Transit time via fitting a cubic to the altitude as a function
    #           of time.
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        sc = get_body(body, tt, location=el)
        aa = AltAz(location=el, obstime=tt)
        aa = sc.transform_to(aa)
    ttf = np.polyfit(aa.obstime.jd-aa.obstime.jd[np.argmax(aa.alt)], aa.alt, 3)
    tt = -ttf[1] - np.sqrt(ttf[1]**2 - 4*ttf[0]*ttf[2])
    tt /= 2*ttf[0]
    tt += aa.obstime.jd[np.argmax(aa.alt)]
    
    _rst.rise = tr
    _rst.transit = tt
    _rst.set = ts
    return _rst


######################################################################
#  Solar Functions
######################################################################

def get_solar_equ_coords(jD: float) -> equ_posn:
    """
    Get Sun's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, and nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position that are
    in the astropy.coordinates.PrecessedGeocentric frame with equinox = jD and
    obstime = jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        b = get_body('sun', t)
        b = b.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
        
    return equ_posn.from_astropy(b)


def get_solar_rst(jD: float, observer: geo_posn) -> Union[rst_time,None]:
    """
    Get Sun's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type geo_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    return _get_solar_system_rst(jD, observer, 'sun')


######################################################################
#  Jupiter Functions
######################################################################

def get_jupiter_equ_coords(jD: float) -> equ_posn:
    """
    Get Jupiter's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position that are
    in the astropy.coordinates.PrecessedGeocentric frame with equinox = jD and
    obstime = jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        b = get_body('jupiter', t)
        b = b.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
        
    return equ_posn.from_astropy(b)


def get_jupiter_rst(jD: float, observer: geo_posn) -> Union[rst_time,None]:
    """
    Get Jupiter's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type geo_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    return _get_solar_system_rst(jD, observer, 'jupiter')


######################################################################
# Saturn Functions
######################################################################

def get_saturn_equ_coords(jD: float) -> equ_posn:
    """   
    Get Saturn's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position that are
    in the astropy.coordinates.PrecessedGeocentric frame with equinox = jD and
    obstime = jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        b = get_body('saturn', t)
        b = b.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
        
    return equ_posn.from_astropy(b)


def get_saturn_rst(jD: float, observer: geo_posn) -> Union[rst_time,None]:
    """
    Get Saturn's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type geo_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """

    return _get_solar_system_rst(jD, observer, 'saturn')


######################################################################
# Lunar Functions
######################################################################

def get_lunar_equ_coords(jD: float) -> equ_posn:
    """
    Get the Moon's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position that are
    in the astropy.coordinates.PrecessedGeocentric frame with equinox = jD and
    obstime = jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        b = get_body('moon', t)
        b = b.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
        
    return equ_posn.from_astropy(b)


def get_lunar_rst(jD: float, observer: geo_posn) -> Union[rst_time,None]:
    """
    Get the Moon's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type geo_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    return _get_solar_system_rst(jD, observer, 'moon')


######################################################################
# Venus Functions
######################################################################

def get_venus_equ_coords(jD: float) -> equ_posn:
    """
    Get Venus' apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position that are
    in the astropy.coordinates.PrecessedGeocentric frame with equinox = jD and
    obstime = jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        b = get_body('venus', t)
        b = b.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
        
    return equ_posn.from_astropy(b)


def get_venus_rst(jD: float, observer: geo_posn) -> Union[rst_time,None]:
    """
    Get Venus' rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type geo_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    return _get_solar_system_rst(jD, observer, 'venus')


######################################################################
# Mars Functions
######################################################################

def get_mars_equ_coords(jD: float) -> equ_posn:
    """
    Get Mars' apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position that are
    in the astropy.coordinates.PrecessedGeocentric frame with equinox = jD and
    obstime = jD.
    """
    
    t = AstroTime(jD, format='jd', scale='utc')
    with solar_system_ephemeris.set(SOLAR_SYSTEM_EPHEMERIS_TO_USE):
        b = get_body('mars', t)
        b = b.transform_to(PrecessedGeocentric(equinox=t, obstime=t))
        
    return equ_posn.from_astropy(b)


def get_mars_rst(jD: float, observer: geo_posn) -> Union[rst_time,None]:
    """
    Get Mars' rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type geo_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    return _get_solar_system_rst(jD, observer, 'mars')


######################################################################
# Time utility constants
######################################################################

"""
Offset in days between standary Julian day and modified Julian day.
""" 
MJD_OFFSET = 2400000.5
    
"""
Offset in days between standard Julian day and Dublin Julian day.
"""
DJD_OFFSET = 2415020.0

"""
Offset in days between UNIX time (epoch 1970/01/01) and standard Julian day.
"""
UNIX_OFFSET = 2440587.5
    
"""
The number of seconds in one day
"""
SECS_IN_DAY = 86400.0

"""
UTC Julian day of B1950.0 coordinate epoch.
"""
B1950_UTC_JD =  2433282.4235

"""
UTC Julian day of J2000.0 coordinate epoch.
"""
J2000_UTC_JD = 2451545.0

"""
Difference in seconds between TT and TAI times.
"""
TAI_TT_OFFSET = 32.184

""""
Velocity of light in meters/second.
"""
VELOCITY_OF_LIGHT = 2.99792458e8


def sec_to_jd(secs: float) -> float:
    """
    Convert seconds into julian days.
    
    Param: secs - seconds (float)
    
    Returns: Julian days as a float. 
    """
    
    return secs / SECS_IN_DAY


def jd_to_sec(jD: float) -> float:
    """
    Convert Julian days into seconds.
    
    Param: jD - Julian days (float).
    
    Returns: Seconds as a float.
    """
    
    return jD * SECS_IN_DAY    


######################################################################
# Time utility functions
######################################################################

def range_hours(hours: float) -> float:
    """
    Put an hour time value into the correct range [0.0, 24.0].
    """
    
    if (hours >= 0.0) and (hours < 24.0):
        return hours
        
    temp = int(hours / 24)
    if hours < 0.0:
        temp -= 1
    temp *= 24
    
    return hours - temp


def jd_to_mjd(jd: float) -> float:
    """
    Get modified julian day value from julian day value.
    
    Param: jd - julian day (should be >= 2400000.5) 
    
    Returns: Modified julian day.
    """
    
    if jd < MJD_OFFSET:
        raise ValueError(f"jd must be >= {MJD_OFFSET}")
        
    return (jd - MJD_OFFSET)


def mjd_to_jd(mjd: float) -> float:
    """
    Get julian day value from modified julian day value.
    
    Param: mjd - modified julian day (should be >= 0.0)
    
    Returns: Julian day.
    """
    
    if mjd < 0.0:
        raise ValueError("mjd must be >= 0.0")
        
    return (mjd + MJD_OFFSET)


FIRST_LEAP_UTC = 2441317.5

_UNIX_EPOCH_AT = AstroTime(0, format='unix', scale='utc')

def leap_secs(utcJD: float) -> float:
    """
    Get the number of leap seconds for given UTC time value.
    
    Param: utcJD - The UTC JD time.
                This should be greater than 2441317.5 (1972 JAN  1). 
    
    Returns: The number of leap seconds (float) for the UTC time.
    """
    
    if utcJD < FIRST_LEAP_UTC:
        raise ValueError(f"utcJD must be greater than {FIRST_LEAP_UTC}")
        
    t = AstroTime(utcJD, format='jd', scale='utc')
    diff_unix = t.unix
    diff_jd = (t - _UNIX_EPOCH_AT).sec
    return round(diff_jd - diff_unix, 1) + 8   # +8 for initial 10 s offset between TAI and UTC


def utc_to_tai(utcJD: float) -> float:
    """
    Get the TAI JD time value for a given UTC JD time value.
    
    Param: utcJD - The UTC JD time (float).
                This should be greater than 2441317.5 (1972 JAN  1). 
    
    Returns: The TAI JD value (float).
    """
    
    t = AstroTime(utcJD, format='jd', scale='utc')
    return t.tai.jd


def tai_to_utc(taiJD: float) -> float:
    """
    Get the UTC JD time value for a given TAI JD time value.
    
    Param: taiJD - The TAI JD time (float).
                This should be greater than 2441317.5 (1972 JAN  1).
    
    Returns: The UTC JD value (float).
    """
    
    t = AstroTime(taiJD, format='jd', scale='tai')
    return t.utc.jd


def taimjd_to_utcjd(taiMJD: float) -> float:
    """
    Get the UTC JD time value for a given TAI MJD value.
    
    Param: mjdTAI - The TAI MJD time (float).
    
    Returns: The UTC JD value (float).
    """
    
    t = AstroTime(taiMJD, format='mjd', scale='tai')
    return t.utc.jd


def utcjd_to_taimjd(utcJD: float) -> float:
    """
    Get the TAI MJD time value for a given UTC JD value.
    
    Param: utcJD - The UTC JD time (float).
    
    Returns: The TAI MJD value.
    """
    
    t = AstroTime(utcJD, format='jd', scale='utc')
    return t.tai.mjd


def unix_to_utcjd(unixTime: Union[int,float]) -> float:
    """
    Get the UTC JD time value for a given UNIX time value.
    
    Param: unixTime - the UNIX time (int/float)
    
    Returns: The UTC JD value.
    """
    
    t = AstroTime(float(unixTime), format='unix', scale='utc')
    return t.jd


def unix_to_taimjd(unixTime: Union[int,float]) -> float:
    """
    Get the TAI MJD time value for a given UNIX time value.
    
    Param: unixTime - the UNIX time (int/float)
    
    Returns: The TAI MJD value.
    """
    
    t = AstroTime(float(unixTime), format='unix', scale='utc')
    return t.tai.mjd


def utcjd_to_unix(utcJD: float) -> float:
    """
    Get UNIX time value for a given UTC JD value.
    
    Param: utcJD - The UTC JD time (float).
    
    Returns: The UNIX time
    """
    
    t = AstroTime(utcJD, format='jd', scale='utc')
    return t.unix


def taimjd_to_unix(taiMJD: float) -> float:
    """
    Get UNIX time value for a given TAI MJDvalue.
    
    Param: taiMJD - The TAI MJD time (float).
    
    Returns: The UNIX time
    """
    
    t = AstroTime(taiMJD, format='mjd', scale='tai')
    return t.utc.unix


def tai_to_tt(taiJD: float) -> float:
    """
    Get the TT JD time value for a given TAI JD time value.
    
    Param: taiJD - The TAI JD time (float).
    
    Returns: The TT JD value (float).
    """
    
    t = AstroTime(taiJD, format='jd', scale='tai')
    return t.tt.jd


def tt_to_tai(ttJD: float) -> float:
    """
    Get the TAI JD time value for a given TT JD time value.
    
    Param: ttJD - The TT JD time (float).
    
    Returns: The TAI JD value (float).
    """   
    
    t = AstroTime(ttJD, format='jd', scale='tt')
    return t.tai.jd


def utc_to_tt(utcJD: float) -> float:
    """
    Get the TT JD time value for a given UTC JD time value.
    
    Param: utcJD - The UTC JD time (float).
    
    Returns: The TT JD value (float).
    """    
    
    t = AstroTime(utcJD, format='jd', scale='utc')
    return t.tt.jd


def tt_to_utc(ttJD: float) -> float:
    """
    Get the UTC JD time value for a given TT JD time value.
    
    Param: ttJD - The TT JD time (float).
    
    Returns: The UTC JD value (float).
    """     
    
    t = AstroTime(ttJD, format='jd', scale='tt')
    return t.utc.jd


def tt_to_tdb(ttJD: float) -> float:
    """
    Get the TDB JD time value for a given TT JD time value.
    Adopted from "Astronomical Almanac Supplement", Seidelmann 1992, 2.222-1.
    
    Param: ttJD - The TT JD time (float).
    
    Returns: The TDB JD value (float).
    """    
    
    t = AstroTime(ttJD, format='jd', scale='tai')
    return t.tdb.jd  


def get_tai_from_sys() -> float:
    """
    Return the current time taken from the system clock as a
    TAI MJD (float).
    """
    
    t0 = AstroTime.now()
    return t0.tai.mjd


def hms_to_sec(hms: hms) -> float:
    """
    Convert hours, minutes, seconds to seconds.
    
    Param: hms - object of type hms representing time/angle.
    
    Returns: Seconds (float) offset of time/angle.
    """
    
    return hms.hours*3600.0 + hms.minutes*60.0 + hms.seconds


def deg_to_sec(degrees: float) -> float:
    """
    Convert longitude degrees into time seconds.
    
    Param: degrees - longitude (float degrees)
    
    Returns: Seconds (float) offset in time for longitude.
    """
    
    hours = degrees / 15.0
    return hours*3600.0


def get_local_sidereal_time(lng: float, jD: float) -> float:
    """
    Get apparent local sidereal time from Julian day.
    
    Param: lng  - longitude degrees (float), E = positive, W = negative 
    Param: jD   - UTC Julian day (float).
    
    Returns: Local mean sidereal time (float hours).
    """
    
    gast = get_apparent_sidereal_time(jD)    
    off = lng / 15.0
    return range_hours(gast + off)


######################################################################
# Pointing utility functions
######################################################################

def dir_cos(posn: hrz_posn) -> Tuple[float,float,float]:
    """
    Get direction cosines from azimuth and zenith angles.
    This function calculates the cosine values based on the LWA coordinate system:
    l = unit vector in E direction (azimuth = 90)
    m = unit vector in N direction (azimuth = 0)
    n = unit vector in zenith direction (zenith = 0, altitude = 90)
    
    Param: posn - object of type hrz_posn giving local position
    
    Returns a tuple (l,m,n) of float values for the direction cosines.
    """
    
    azRad = math.radians(posn.az)   # type: ignore
    zenRad = math.radians(posn.zen())   # type: ignore
    szen = math.sin(zenRad)
    
    l = (szen * math.sin(azRad))
    m = (szen * math.cos(azRad))
    n = math.cos(zenRad)
    
    return (l, m, n)


def get_rect_from_equ(posn: equ_posn) -> rect_posn:
    """
    Transform equatorial coordinates to rectangular coordinates.
    
    Param: posn - Object of type equ_posn giving position.
    
    Returns: Object of type rect_posn giving rectangular coordinates (normallized to 1).
    """
    
    sc = SkyCoord(posn.ra*astrounits.deg, posn.dec*astrounits.deg,
                  frame='fk5', equinox='J2000')
    sc = sc.cartesian 
    
    x = sc.x.value
    y = sc.y.value
    z = sc.z.value
    
    return rect_posn(x, y, z)


def get_equ_from_rect(posn: rect_posn) -> equ_posn:
    """
    Transform rectangular coordinates to equatorial coordinates.
    
    Param: posn - Object of type rect_posn giving position.
    
    Returns: Object of type equ_posn giving equatorial coordinates.
    """
    
    x = posn.X
    y = posn.Y
    z = posn.Z
    
    sc = SkyCoord(CartesianRepresentation(x, y, z), frame='fk5', equinox='J2000')
    
    return equ_posn.from_astropy(sc)


def get_geo_from_rect(posn: rect_posn) -> geo_posn:
    """
    Transform ECEF rectangular coordinates to geographical coordinates.
    Adapoted from "Satellite Orbits", Montenbruck and Gill 2005, 5.85 - 5.88.
    Also see gpstk ECEF::asGeodetic() method.
    
    Param: posn - object of type rect_posn giving position in m.
    
    Returns: object of type geo_posn giving geographical coordinates.
    """
    
    el = EarthLocation.from_geocentric(posn.X*astrounits.m, posn.Y*astrounits.m, posn.Z*astrounits.m)
    
    return geo_posn.from_astropy(el)


def get_rect_from_geo(posn: geo_posn) -> rect_posn:
    """
    Transform geographical coordinates to ECEF rectangular coordinates.
    Adopted from "Satellite Orbits", Montenbruck and Gill 2005, 5.83 - 5.84.
    Also see gpstk Geodetic::asECEF() method.
    
    Param: posn - object of type geo_posn giving geographical coordinates.
    
    Returns: object of type rect_posn giving ECEF position in m. 
    """
    
    el = EarthLocation.from_geodetic(posn.lng*astrounits.deg, posn.lat*astrounits.deg,
                                     height=posn.elv*astrounits.m,
                                     ellipsoid='WGS84')
    x = el.x.to('m').value
    y = el.y.to('m').value
    z = el.z.to('m').value     
    
    return rect_posn(x, y, z)


def get_precession(jD1: float, pos: equ_posn, jD2: float) -> equ_posn:
    """
    Caculate precession of equatorial coordinates from one epoch to
    another.
    
    Param: jD1 - UTC Julian day of epoch 1.
    Param: pos - object of type equ_posn giving epoch 1 position
    Param: jD2 - UTC Julian day of epoch 2.
    
    Returns: object of type equ_posn giving epoch 2 position.
    """
    
    t1 = AstroTime(jD1, format='jd', scale='utc')
    sc = SkyCoord(pos.ra*astrounits.deg, pos.dec*astrounits.deg,
                  frame='fk5', equinox=t1)
    t2 = AstroTime(jD2, format='jd', scale='utc')
    sc = sc.transform_to(FK5(equinox=t2))
    
    return equ_posn.from_astropy(sc)


def B1950_to_J2000(pos: equ_posn) -> equ_posn:
    """
    Convert B1950 epoch to J2000 epoch for equatorial coordinates.
    
    Param: pos - object of type equ_posn giving B1950 coordinates
    
    Returns: object of type equ_posn giving J2000 coordinates in the
             astropy.coordinates.FK5 frame with equinox=J2000.
    
    .. note::
        The accuracy of this function is about 0.01 degrees.
    """
    
    sc = pos.astropy
    if sc is None:
        sc = SkyCoord(pos.ra*astrounits.deg, pos.dec*astrounits.deg,
                      frame='fk4', equinox='B1950')
    sc = sc.transform_to(FK5(equinox='J2000'))
    
    return equ_posn.from_astropy(sc)


def J2000_to_B1950(pos: equ_posn) -> equ_posn:
    """
    Convert J2000 epoch to B1950 epoch for equatorial coordinates.
    
    Param: pos - object of type equ_posn giving J2000 coordinates
    
    Returns: object of type equ_posn giving B1950 coordinates in the
             astropy.coordinates.FK4 frame with equinox=B1950.
    
    .. note::
        The accuracy of this function is about 0.01 degrees.
    """   
    
    sc = pos.astropy
    if sc is None:
        sc = SkyCoord(pos.ra*astrounits.deg, pos.dec*astrounits.deg,
                      frame='fk5', equinox='J2000')
    sc = sc.transform_to(FK4(equinox='B1950'))
    
    return equ_posn.from_astropy(sc)


def resolve_name(name: str) -> equ_posn:
    """
    Given the name of an astronomical source resolve it into coordinates.
    
    Param: name - object of type str giving the name of the source
    
    Returns: object of equ_posn giving coordinates in the astropy.coordinats.ICRS
             frame with additional information about the proper motion (pm_ra and
             pm_dec; mas/yr), distance (distance; pc), and resolver service
             (resolved_by; str).
    """
    
    try:
        with urlopen('https://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame/-oxp/SNV?%s' % quote_plus(name)) as result:
            tree = ElementTree.fromstring(result.read())
            target = tree.find('Target')
            service = target.find('Resolver')   # type: ignore
            ra = service.find('jradeg') # type: ignore
            dec = service.find('jdedeg')    # type: ignore
            try:
                pm = service.find('pm') # type: ignore
            except Exception as e:
                pm = None
            try:
                plx = service.find('plx')   # type: ignore
            except Exception as e:
                plx = None
                
            service = service.attrib['name'].split('=', 1)[1]   # type: ignore
            ra = float(ra.text) # type: ignore
            dec = float(dec.text)   # type: ignore
            coordsys = 'J2000'
            if pm is not None:
                pmRA = float(pm.find('pmRA').text)  # type: ignore
                pmDec = float(pm.find('pmDE').text) # type: ignore
            else:
                pmRA = None
                pmDec = None
            if plx is not None:
                dist = float(plx.find('v').text)    # type: ignore
            else:
                dist = None
                
            if pmRA is not None:
                pmRA = pmRA*math.cos(dec*math.pi/180)*astrounits.mas/astrounits.yr
            if pmDec is not None:
                pmDec = pmDec*astrounits.mas/astrounits.yr
            if dist is not None:
                dist = dist*astrounits.pc
                
            sc = SkyCoord(ra*astrounits.deg, dec*astrounits.deg,
                          pm_ra_cosdec=pmRA, pm_dec=pmDec,
                          distance=dist,
                          frame='icrs')
            _posn = equ_posn.from_astropy(sc)
            _posn.resolved_by = service
            
    except (IOError, AttributeError, ValueError, RuntimeError) as e:
        raise RuntimeError(f"Failed to resolve source '{name}'")
        
    return _posn
