"""
Astronomical utility functions and classes based on libnova library.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import time
import math
import ephem
import numpy
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
from calendar import timegm
from datetime import datetime
from functools import total_ordering    

from lsl.common.progress import DownloadBar

from lsl.misc import telemetry
telemetry.track_module()


__version__   = '0.5'
__all__ = ['dms', 'hms', 'date', 'zonedate', 'rst_time', 'hrz_posn', 'equ_posn', 
           'gal_posn', 'rect_posn', 'lnlat_posn', 'ecl_posn', 'nutation', 
           'get_gmtoff', 'date_to_zonedate', 'zonedate_to_date', 'rad_to_deg', 
           'deg_to_rad', 'dms_to_rad', 'dms_to_deg', 'deg_to_dms', 'rad_to_dms', 
           'hms_to_deg', 'hms_to_rad', 'deg_to_hms', 'rad_to_hms', 'add_secs_hms', 
           'add_hms', 'hrz_to_nswe', 'range_degrees', 'get_julian_day', 
           'get_julian_local_date', 'get_date', 'get_day_of_week', 
           'get_julian_from_sys', 'get_date_from_sys', 'get_julian_from_timet', 
           'get_timet_from_julian', 'get_hrz_from_equ', 'get_equ_from_hrz', 
           'get_ecl_from_rect', 'get_equ_from_ecl', 'get_ecl_from_equ', 
           'get_equ_from_gal', 'get_gal_from_equ', 'get_equ2000_from_gal', 
           'get_gal_from_equ2000', 'get_apparent_sidereal_time', 
           'get_mean_sidereal_time', 'get_angular_separation', 'get_rel_posn_angle', 
           'get_apparent_posn', 'get_equ_prec', 'get_equ_prec2', 'get_nutation', 
           'get_equ_nut', 'get_equ_aber', 'get_equ_pm', 'get_object_rst', 
           'get_solar_equ_coords', 'get_solar_rst', 'get_jupiter_equ_coords', 
           'get_jupiter_rst', 'get_saturn_equ_coords', 'get_saturn_rst', 
           'get_lunar_equ_coords', 'get_lunar_rst', 'get_venus_equ_coords', 
           'get_venus_rst', 'get_mars_equ_coords', 'get_mars_rst', 'sec_to_jd', 
           'jd_to_sec',  'range_hours', 'jd_to_mjd', 'mjd_to_jd', 'leap_secs', 
           'utc_to_tai', 'tai_to_utc', 'taimjd_to_utcjd', 'utcjd_to_taimjd', 
           'unix_to_utcjd', 'unix_to_taimjd', 'utcjd_to_unix', 'taimjd_to_unix', 
           'tai_to_tt', 'tt_to_tai', 'utc_to_tt', 'tt_to_utc', 'tt_to_tdb', 
           'get_tai_from_sys', 'hms_to_sec', 'deg_to_sec', 'get_local_sidereal_time', 
           'geo_posn', 'dir_cos', 'get_rect_from_equ', 'get_equ_from_rect', 
           'get_geo_from_rect', 'get_rect_from_geo', 'get_precession', 'B1950_to_J2000', 
           'J2000_to_B1950']
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
    
    def __init__(self, neg = False, degrees = 0, minutes = 0, seconds = 0.0):
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
                raise ValueError("degrees parameter range is [0, 359], is set to %d" % degrees)
            self.degrees = degrees
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError("minutes parameter range is [0, 59], is set to %d" % minutes)
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError("seconds parameter range is [0.0, 60.0), is set to %0.3f" % seconds)
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
            
        sec = "%02.2f" % self.seconds   	
        
        return "%s%02d %02d %s" % (sign, self.degrees, self.minutes, sec.zfill(5))
        
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
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
                
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
        
    def to_deg(self):
        """
        Convert angles degrees, minutes, seconds to float degrees.
        Returns angle in degrees (float).
        """
        
        return dms_to_deg(self)
        
    def to_hms(self):
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
    
    def __init__(self, hours = 0, minutes = 0, seconds = 0.0):
        """
        Create a hms object.
        
        Param: hours    - Angle/time hours (integer [0, 23]).
        Param: minutes  - Angle/time minutes (integer [0, 59]).
        Param: seconds  - Angle/time seconds (float [0.0, 60.0)).
        """
        
        if hours is not None:
            if hours < 0 or hours > 23:
                raise ValueError("hours paramerer range is [0, 23], is set to %d" % hours)
            self.hours = hours
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError("minutes paramerer range is [0, 59], is set to %d" % minutes)
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError("seconds paramerer range is [0.0, 60.0), is set to %0.3f" % seconds)
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
        
        sec = "%02.2f" % self.seconds
        
        return "%02d %02d %s" % (self.hours, self.minutes, sec.zfill(5))
        
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
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
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
        
    def to_deg(self):
        """
        Convert angles hours, minutes, seconds to float degrees.
        Returns angle in degrees (float).
        """
        
        return hms_to_deg(self)
        
    def to_dms(self):
        """
        Convert angle hours, minutes, seconds to degrees, minutes, seconds.
        Returns: object of type dms representing angle.
        """
        
        return deg_to_dms(self.to_deg())
        
    def to_sec(self):
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
    
    def __init__(self, years = 2000, months = 1, days = 1, hours = 0, minutes = 0, seconds = 0.0):
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
                raise ValueError("months paramerer range is [1, 12], is set to %d" % months)
            self.months = months
            
        if days is not None:
            if days < 1 or days > 31:
                raise ValueError("days paramerer range is [1, 31], is set to %d" % days)
            self.days = days
            
        if hours is not None:
            if hours < 0 or hours > 23:
                raise ValueError("hours paramerer range is [0, 23], is set to %d" % hours)
            self.hours = hours
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError("minutes paramerer range is [0, 59], is set to %d" % minutes)
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError("degrees paramerer range is [0.0, 60.0), is set to %0.3f" % seconds)
            self.seconds = seconds
            
    def __str__(self):
        """
        date object print/str method.
        """
        
        sec = "%02.3f" % self.seconds
        
        return "%04d-%02d-%02d %02d:%02d:%s" % (self.years, self.months, self.days, self.hours, self.minutes, sec.zfill(6))
        
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
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
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
        
    def to_zone(self, gmtoff = None):
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
    
    def to_jd(self):
        """
        Convert calendar time to Julian day.
        Returns UTC time in Julian days (float).
        """
        
        return get_julian_day(self)
        
    def load(self, dateStr, timeStr):
        """
        Load date object from formatted strings.
        
        Param: dateStr - A string of format YYYY-MM-DD giving date.
        Param: timeStr - A string of format HH:MM:SS.S giving time.
        """
        
        # parse date string
        try:
            (yearStr, monthStr, dayStr) = dateStr.split('-')
        except ValueError:
            raise ValueError("date incorrectly formated: %s" % dateStr) 
            
        try:
            year = int(yearStr)
        except ValueError:
            raise ValueError("year incorrectly formated: %s" % yearStr)
            
        try:
            month = int(monthStr)
        except ValueError:
            raise ValueError("month incorrectly formated: %s" % monthStr)
        if month < 1 or month > 12:
            raise ValueError("months paramerer range is [1, 12], is set to %d" % month)
            
        try:
            day = int(dayStr)
        except ValueError:
            raise ValueError("day incorrectly formated: %s" % dayStr)
        if day < 1 or day > 31:
            raise ValueError("days paramerer range is [1, 31], is set to %d" % day)
            
        # parse time sting 
        try:    
            (hourStr, minuteStr, secondStr) = timeStr.split(':') 
        except ValueError:
            raise ValueError("time incorrectly formated: %s" % timeStr)         
            
        try:
            hour = int(hourStr)
        except ValueError:
            raise ValueError("hour incorrectly formated: %s" % hourStr)
        if hour < 0 or hour > 23:
            raise ValueError("hours paramerer range is [0, 23], is set to %d" % hour)
            
        try:
            minute = int(minuteStr)
        except ValueError:
            raise ValueError("minutes incorrectly formated: %s" % minuteStr)
        if minute < 0 or minute > 59:
            raise ValueError("minutes paramerer range is [0, 59], is set to %d" % minute)
            
        try:
            second = float(secondStr)
        except ValueError:
            raise ValueError("seconds incorrectly formated: %s" % secondStr)
        if second < 0.0 or second >= 60.0:
            raise ValueError("degrees paramerer range is [0.0, 60.0), is set to %0.3f" % second)
            
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
    
    def __init__(self, years = 2000, months = 1, days = 1, hours = 0, minutes = 0, seconds = 0.0, gmtoff = 0):
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
                raise ValueError("months paramerer range is [1, 12], is set to %d" % months)
            self.months = months
            
        if days is not None:
            if days < 1 or days > 31:
                raise ValueError("days paramerer range is [1, 31], is set to %d" % days)
            self.days = days
            
        if hours is not None:
            if hours < 0 or hours > 23:
                raise ValueError("hours paramerer range is [0, 23], is set to %d" % hours)
            self.hours = hours
            
        if minutes is not None:
            if minutes < 0 or minutes > 59:
                raise ValueError("minutes paramerer range is [0, 59], is set to %d" % minutes)
            self.minutes = minutes
            
        if seconds is not None:
            if seconds < 0.0 or seconds >= 60.0:
                raise ValueError("degrees paramerer range is [0.0, 60.0), is set to %0.3f" % seconds)
            self.seconds = seconds
            
        if gmtoff is None:
            gmtoff = get_gmtoff()
        self.gmtoff = gmtoff
    
    def __str__(self):
        """
        zonedate object str/print method.
        """
        
        sec = "%02.3f" % self.seconds 
        
        return "%04d-%02d-%02d %02d:%02d:%s [%d]" % (self.years, self.months, self.days, self.hours, self.minutes, sec.zfill(6), self.gmtoff)
        
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
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
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
        
    def to_date(self):
        """
        Convert local calendar time to UTC calendar time.
        Returns object of type date representing UTC time.
        """ 
        
        return zonedate_to_date(self)
        
    def to_jd(self):
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
    
    def __init__(self, rise = None, set = None, transit = None):
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
        
        return "%0.3f %0.3f %0.3f" % (self.rise, self.set, self.transit)
        
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
        
    def format(self):
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
    
    def __init__(self, az = 0.0, alt = 0.0):
        """
        Create a hrz_posn object.
        
        Param: az   - Position azimuth angle (float degrees [0.0, 360.0), 0 = N, 90 = E).
        Param: alt  - Position altitude angle (float degrees [-90.0, 90.0]).
        """

        if az is not None:
            if az < 0.0 or az >= 360.0:
                raise ValueError("az paramerer range is [0.0, 360.0), is set to %0.3f" % az)
            self.az = az
            
        if alt is not None:
            if alt < -90.0 or alt > 90.0:
                raise ValueError("alt paramerer range is [-90.0, 90.0], is set to %0.3f" % alt)
            self.alt = alt
            
    def zen(self, value = None):
        """
        If value is None, returns position zenith angle (float degrees 
        [0, 180])  Otherwise, sets the altitude according to the zenith angle
        value.
        """
        
        if value is None:
            return 90.0 - self.alt
        if value <  0.0 or value > 180.0:
            raise ValueError("value paramerer range is [0.0, 180.0], is set to %0.3f" % value)
        self.alt = 90.0 - value 
        
    def __str__(self):
        """
        hrz_posn object print/str method.
        """
        
        return "%0.3f %0.3f" % (self.az, self.alt)
        
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
            raise ValueError("subscript %d out of range" % key)
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if key == 0:
            self.az = value
        elif key == 1:
            self.alt = value
        else:
            raise ValueError("subscript %d out of range" % key)
            
    def __eq__(self, other):
        """
        hrz_posn equality test.
        """
        
        if not isinstance(other, hrz_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.az == other.az) and (self.alt == other.alt)
        
    def __ne__(self, other):
        """
        hrz_posn non-equality test.
        """
        
        if not isinstance(other, hrz_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.az != other.az) or (self.alt != other.alt)
        
    def to_equ(self, observer, jD):
        """
        Get equatorial/celestial coordinates from local horizontal coordinates.
        
        Param: observer - Object of type lnlat_posn representing observer position.
        Param: jD       - UTC Julian day (float).
        
        Returns object of type equ_posn representing equatorial position.
        """
        
        return get_equ_from_hrz(self, observer, jD)         
        
    def dir_cos(self):
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
    
    def __init__(self, ra = 0.0, dec = 0.0):
        """
        Create a equ_posn object.
        
        Param: ra   - Position right ascension angle
                    Object of type hms or float degrees [0.0, 360.0).
        Param: dec  - Position declination angle
                    Object of type dms or float degrees [-90.0, 90.0].
        """
        
        if ra is not None:
            if isinstance(ra, hms):
                ra = ra.to_deg()
            if ra < 0.0 or ra >= 360.0:
                raise ValueError("ra paramerer range is [0.0, 360.0), is set to %0.3f" % ra)
            self.ra = ra
            
        if dec is not None:
            if isinstance(dec, dms):
                dec = dec.to_deg()
            if dec < -90.0 or dec > 90.0:
                raise ValueError("dec paramerer range is [-90.0, 90.0], is set to %0.3f" % dec)
            self.dec = dec
            
    def __str__(self):
        """
        equ_posn object str/print method.
        """
        
        return "%0.3f %0.3f" % (self.ra, self.dec)
        
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
            raise ValueError("subscript %s out of range" % key)
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if key == 0:
            self.ra = value
        elif key == 1:
            self.dec = value
        else:
            raise ValueError("subscript %s out of range" % key)
            
    def __eq__(self, other):
        """
        equ_posn equality test.
        """
        
        if not isinstance(other, equ_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.ra == other.ra) and (self.dec == other.dec)
        
    def __ne__(self, other):
        """
        equ_posn non-equality test.
        """
        
        if not isinstance(other, equ_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.ra != other.ra) or (self.dec != other.dec)
        
    def to_hrz(self, observer, jD):
        """
        Get local horizontal coordinates from equatorial/celestial coordinates.
        
        Param: observer - Object of type lnlat_posn representing observer position.
        Param: jD       - UTC Julian day (float).
        
        Returns object of type hrz_posn representing local position.
        """
        
        return get_hrz_from_equ(self, observer, jD)
        
    def to_gal(self, jD):
        """
        Get J2000 galactic coordinates from apparent equatorial coordinates.
        
        Param: jD - UTC Julian day to get position.
        
        Returns object of type gal_posn representing object's galactic 
        position.
        """
        
        equ = get_equ_prec2(self, jD, J2000_UTC_JD)
        return get_gal_from_equ2000(equ)
        
    def to_ecl(self, jD):
        """
        Get ecliptical coordinates from equatorial coordinates.
        
        Param: jD - UTC Julian day (float).
        
        Returns object of type lnlat_posn representing eclipitcal position.
        """
        
        return get_ecl_from_equ(self, jD)
        
    def angular_separation(self, posn):
        """
        Get angular separation from equatorial position.
        
        Param: posn - Object of type equ_posn representing body 2 position.
        
        Returns angular separation in degrees (float).
        """
        
        return get_angular_separation(self, posn) 
        
    def format(self):
        """
        Return a tuple (ra, dec) where ra is an hms object and
        dec is a dms object representing ra and dec position coordinates.
        """
            
        return (deg_to_hms(self.ra), deg_to_dms(self.dec)) 
        
    def precess(self, jD):
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
    
    def __init__(self, l = 0.0, b = 0.0):
        """
        Create a gal_posn object.
        
        Param: l - Position longitude angle. 
                Object of type dms or float degrees [0.0, 360.0).
        Param: b - Position latitude angle. 
                Object of type dms or float degrees [-90.0, 90.0].
        """
        
        if l is not None:
            if isinstance(l, dms):
                l = l.to_deg()
            if l < -360.0 or l >= 360.0:
                raise ValueError("l parameter range is [-360.0, 360.0), is set to %0.3f" % l)
            self.l = l
            
        if b is not None:
            if isinstance(b, dms):
                b = b.to_deg()
            if b < -90.0 or b > 90.0:
                raise ValueError("b paramerer range is [-90.0, 90.0], is set to %0.3f" % b)
            self.b = b
            
    def __str__(self):
        """
        gal_posn object print/str method.
        """
        
        return "%0.3f %0.3f" % (self.l, self.b)
        
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
            raise ValueError("subscript %s out of range" % key)
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if key == 0:
            self.l = value
        elif key == 1:
            self.b = value
        else:
            raise ValueError("subscript %s out of range" % key)
            
    def __eq__(self, other):
        """
        gal_posn equality test.
        """
        
        if not isinstance(other, gal_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)        
            
        return (self.l == other.l) and (self.b == other.b)
        
    def __ne__(self, other):
        """
        gal_posn non-equality test.
        """
        
        if not isinstance(other, gal_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.l != other.l) or (self.b != other.b)
        
    def to_equ(self, jD):
        """
        Get apparent equatorial coordinates from J2000 galactic coordinates.
        
        Param: jD - UTC Julian day to get position.
        
        Returns object of type equ_posn representing object's apparent
        equatorial position.
        """
        
        equ = get_equ2000_from_gal(self)
        return get_equ_prec(equ, jD)
        
    def format(self):
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
    
    def __init__(self, X = 0.0, Y = 0.0, Z = 0.0): 
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
        
        return "%0.3f %0.3f %0.3f" % (self.X, self.Y, self.Z)
        
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
            raise ValueError("subscript %s out of range" % key)
            
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
            raise ValueError("subscript %s out of range" % key)
            
    def __eq__(self, other):
        """
        rect_posn equality test.
        """
        
        if not isinstance(other, rect_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.X == other.X) and (self.Y == other.Y) and (self.Z == other.Z)
        
        
    def __ne__(self, other):
        """
        rect_posn non-equality test.
        """
        
        if not isinstance(other, rect_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.X != other.X) or (self.Y != other.Y) or (self.Z != other.Z)


class lnlat_posn(object):
    """
    Represents position coordinates in latitude and longitude.
    When representing a geographical location, the longitude is negative
    when measured west of GM and positive is measured east of GM.
    
    Public members:
      lng - Position longitude coordinate (float degrees).
      lat - Position latitude coordinate (float degrees).
    
    Members may also be accessed by subscript:
      lnlat_posn[0] = lng
      lnlat_posn[1] = lat
    """
    
    def __init__(self, lng = 0.0, lat = 0.0):
        """
        Create a lnlat_posn object.
        
        Param: lng - Position longitude coordinate
                    Object of type dms or float degrees (-360.0, 360.0).
        Param: lat - Position latitude coordinate
                    Object of type dms or float degrees [-90.0, 90.0].
        """
                
        if lng is not None:
            if isinstance(lng, dms):
                lng = lng.to_deg()
            if lng <= -360.0 or lng >= 360.0:
                raise ValueError("lng parameter range is (-360.0, 360.0), is set to %0.3f" % lng)
            self.lng = lng
            
        if lat is not None:
            if isinstance(lat, dms):
                lat = lat.to_deg()
            if lat < -90.0 or lat > 90.0:
                raise ValueError("lat paramerer range is [-90.0, 90.0], is set to %0.3f" % lat)
            self.lat = lat
            
    def __str__(self):
        """
        lnlat_posn object print/str method.
        """
        
        return "%0.3f %0.3f" % (self.lng, self.lat)
        
    def __repr__(self):
        """
        lnlat_posn object repr string method
        """
        
        return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.lng), repr(self.lat))
        
    def __reduce__(self):
        """
        lnlat_posn object pickle reduce method.
        """
        
        return (lnlat_posn, (self.lng, self.lat))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 0:
            return self.lng
        elif key == 1:
            return self.lat
        else:
            raise ValueError("subscript %s out of range" % key)
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if key == 0:
            self.lng = value
        elif key == 1:
            self.lat = value
        else:
            raise ValueError("subscript %s out of range" % key)
            
    def __eq__(self, other):
        """
        lnlat_posn equality test.
        """
        
        if not isinstance(other, lnlat_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.lng == other.lng) and (self.lat == other.lat)
        
    def __ne__(self, other):
        """
        lnlat_posn non-equality test.
        """
        
        if not isinstance(other, lnlat_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return (self.lng != other.lng) or (self.lat != other.lat)
        
    def format(self):
        """
        Return a tuple (lng, lat) where lng is an dms object and
        lat is a dms object representing longitude and latitude
        position coordinates.
    """
        
        return (deg_to_dms(self.lng), deg_to_dms(self.lat))


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
    
    def __init__(self, lng = 0.0, lat = 0.0):
        self.lng = lng
        self.lat = lat
        
    def to_equ(self, jD):
        """
        Get equatorial coordinates from ecliptical coordinates for a given time.
        
        Param: jD - UTC Julian day (float). 
        
        Returns object of type equ_posn representing equatorial position.
        """
        
        return get_equ_from_ecl(self, jD)
        
    def __reduce__(self):
        """
        ecl_posn object pickle reduce method.
        """
        
        return (ecl_posn, (self.lng, self.lat))
        
    def format(self):
        """
        Return a tuple (lng, lat) where lng is an dms object and
        lat is a dms object representing longitude and latitude
        position coordinates.
        """
        
        return (deg_to_dms(self.lng), deg_to_dms(self.lat))


class nutation(object):
    """
    Provides nutation information in longitude and obliquity.
    
    Public members:
      longitude - Nutation in longitude (float degrees).
      obliquity - Nutation in ecliptic obliquity (float degrees).
      ecliptic - Obliquity of the ecliptic (float degrees).
    """
    
    def __init__(self, longitude = 0.0, obliquity = 0.0, ecliptic = 0.0):
        """
        Create a nutation object.
        
        Param: longitude  - Nutation in longitude.
                            Object of type dms or float degrees (-360.0, 360.0).
        Param: obliquity  - Nutation in obliquity.
                            Object of type dms or float degrees [-90.0, 90.0].
        Param: ecliptic   - Obliquity of the ecliptic.
                            Object of type dms or float degrees [-90.0, 90.0].
        """
        
        if longitude is not None:
            if isinstance(longitude, dms):
                longitude = longitude.to_deg()
            if longitude <= -360.0 or longitude >= 360.0:
                raise ValueError("longitude parameter range is (-360.0, 360.0), is set to %0.3f" % longitude)
            self.longitude = longitude
            
        if obliquity is not None:
            if isinstance(obliquity, dms):
                obliquity = obliquity.to_deg()
            if obliquity < -90.0 or obliquity > 90.0:
                raise ValueError("obliquity paramerer range is [-90.0, 90.0], is set to %0.3f" % obliquity)
            self.obliquity = obliquity
            
        if ecliptic is not None:
            if isinstance(ecliptic, dms):
                ecliptic = ecliptic.to_deg()
            if ecliptic < -90.0 or ecliptic > 90.0:
                raise ValueError("ecliptic paramerer range is [-90.0, 90.0], is set to %0.3f" % ecliptic)
            self.ecliptic = ecliptic  
            
    def __str__(self):
        """
        nutation object print/str method.
        """
        
        return "%0.3f %0.3f %0.3f" % (self.longitude, self.obliquity, self.ecliptic) 
        
    def __repr__(self):
        """
        nutation object repr string method
        """
        
        return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.longitude), repr(self.obliquity), repr(self.ecliptic))
        
    def __reduce__(self):
        """
        nutation object pickle reduce method.
        """
        
        return (nutation, (self.longitude, self.obliquity, self.ecliptic))
        
    def format(self):
        """
        Return a tuple (lng, obl, ecl) where lng is an dms object,
        obl is a dms object, and ecl is a dms object representing nutation
        in longitude and obliquity, and obliquity of the ecliptic.
        """
        
        return (deg_to_dms(self.longitude), deg_to_dms(self.obliquity), deg_to_dms(self.ecliptic))


######################################################################
# time helper python fucntions
######################################################################

def get_gmtoff():
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

def date_to_zonedate(date, gmtoff):
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


def zonedate_to_date(zonedate):
    """
    Convert local calendar time to UTC calendar time.
    
    Param: zonedate - Object of type zonedate representing local time.  
    
    Returns object of type date representing UTC time.
    """
    
    t0, junk = str(zonedate).rsplit(None, 1)
    t0 = time.strptime(t0, "%Y-%m-%d %H:%M:%S.%f")
    fracSec = zonedate.seconds - int(zonedate.seconds)
    t1 = timegm(t0) - zonedate.gmtoff
    years, months, days, hours, minutes, seconds, wday, yday, dst = time.gmtime(t1)
    
    _date = date()
    _date.years = years
    _date.months = months
    _date.days = days
    _date.hours = hours
    _date.minutes = minutes
    _date.seconds = seconds + fracSec
    return _date


def rad_to_deg(radians):
    """
    Convert radians to degrees.
    
    Param: radians - Angle in radians (float).
    
    Returns angle in degrees (float).
    """
    
    return radians * 180.0/math.pi


def deg_to_rad(degrees):
    """
    Convert degres to radians.
    
    Param: degrees - Angle in degrees (float).
    
    Returns angle in radians (float).
    """
    
    return degrees * math.pi/180.0


def dms_to_rad(dms):
    """
    Convert angles degrees, minutes, seconds to radians.
    
    Param: dms - Object of type dms representing angle.
    
    Returns angle in radians (float).
    """
    
    degrees = dms_to_rad(dms)
    return deg_to_rad(degrees)


def dms_to_deg(dms):
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


def deg_to_dms(degrees):
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


def rad_to_dms(radians):
    """
    Convert angles float radians to degrees, minutes, seconds.
    
    Param: radians - Angle in radians (float). 
    
    Returns object of type dms representing angle.
    """
    
    degrees = rad_to_deg(radians)
    return deg_to_dms(degrees)


def hms_to_deg(hms):
    """
    Convert angles hours, minutes, seconds to float degrees.
    
    Param: hms - Object of type hms representing angle.
    
    Returns angle in degrees (float).
    """
    
    hours = hms.hours + hms.minutes/60.0 + hms.seconds/3600.0
    degrees = hours * 15.0
    return degrees


def hms_to_rad(hms):
    """
    Convert angles hours, minutes, seconds to float radians.
    
    Param: hms - Object of type hms representing angle.
    
    Returns angle in radians (float).
    """
    
    degrees = hms_to_deg(hms)
    return deg_to_rad(degrees)


def deg_to_hms(degrees):
    """
    Convert angles float degrees to hours, minutes, seconds.
    
    Param: degrees - Angle in degrees (float). 
    
    Returns object of type hms representing angle.
    """
    
    sgn,h,m,s = _float_to_sexa(degrees / 15.0)
    return hms(h, m, s)


def rad_to_hms(radians):
    """
    Convert angles float radians to hours, minutes, seconds.
    
    Param: radians - Angle in radians (float). 
    
    Returns object of type hms representing angle.
    """
    
    degrees = rad_to_deg(radians)
    return deg_to_hms(degrees)


def add_secs_hms(hms, seconds):
    """
    Add seconds to time/angle hours, minutes, seconds.
    
    Param: hms      - Object of type hms representing angle.
    Param: seconds  - Seconds offset (float) to add to angle. 
    
    Returns object of type hms representing angle + offset.
    """
    
    degrees = hms_to_deg(hms)
    degrees += (seconds/3600.0) * 15.0
    return deg_to_hms(degrees)


def add_hms(source, dest):
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
    
    dest.degrees = 1*_hms.hours
    dest.minutes = 1*_hms.minutes
    dest.seconds = 1*_hms.seconds
    
    return dest


def hrz_to_nswe(pos):
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


def range_degrees(degrees):
    """
    Put angle into range [0, 360].
    
    Param: degrees - large angle (float degrees)
    
    Returns: angle in range (float degrees)
    """
    
    return degrees % 360.0


######################################################################
# General Calendar Functions
######################################################################

def get_julian_day(date):
    """
    Convert calendar time to Julian day.
    
    Param: date - Object of type date representing UTC time.
    
    Returns UTC time in Julian days (float).
    """
    
    _date = ephem.Date("%i/%02i/%02i %02i:%02i:%09.6f" % (date.years, date.months, date.days, date.hours, date.minutes, date.seconds))
    jd = float(_date)+DJD_OFFSET
    return jd
    
    
def get_julian_local_date(zonedate):
    """
    Convert local calendar time to Julian day.
    
    Param: zonedate - Object of type zonedate representing local time.
    
    Returns UTC time in Julian days (float).
    """
    
    _date = zonedate_to_date(zonedate)
    return get_julian_day(_date)


def get_date(jD):
    """
    Convert Julian day to calendar time.
    
    Param: jD - UTC time in Julian days (float).
    
    Returns object of type date representing UTC time.
    """
    
    _date = date()
    d = ephem.Date(jD-DJD_OFFSET)
    
    years,months,days,hours,minutes,seconds = d.tuple()
    _date.years = years
    _date.months = months
    _date.days = days
    _date.hours = hours
    _date.minutes = minutes
    _date.seconds = seconds
    return _date


def get_day_of_week(date):
    """
    Gets day of week from calendar time.
    
    Param: date - Object of type date representing UTC time.
    
    Returns day of week (0 = Sunday, 6 = Saturday).
    """
    
    jd = date.to_jd()
    return (int(round(jd)) + 1) % 7


def get_julian_from_sys():
    """
    Returns UTC Julian day (float) from system clock.
    """
    
    t0 = time.time()
    return unix_to_utcjd(t0)


def get_date_from_sys():
    """
    Gets calendar time from system clock.
    
    Returns object of type date representing UTC time.
    """
    
    jD = get_julian_from_sys()
    _date = get_date(jD)
    return _date


def get_julian_from_timet(timet):
    """
    Gets Julian day from Unix time.
    
    Param: timet - Unix timet in seconds (integer)
    
    Returns UTC Julian day (float).
    """
    
    jD = float(timet) / SECS_IN_DAY + UNIX_OFFSET
    return jD


def get_timet_from_julian(jD):
    """
    Gets Unix time from Julian day.
    
    Param: jD - UTC Julian day (float).
    
    Returns Unix timet in seconds (integer).
    """
    
    timet = int((jD - UNIX_OFFSET) * SECS_IN_DAY)
    return timet


######################################################################
# Transformation of Coordinates Functions
######################################################################

def get_hrz_from_equ(target, observer, jD):
    """
    Get local horizontal coordinates from equatorial/celestial coordinates.
    
    Param: target   - Object of type equ_posn representing celestial position.
    Param: observer - Object of type lnlat_posn representing observer position.
    Param: jD       - UTC Julian day (float).
    
    Returns object of type hrz_posn representing local position.
    """
    
    _posn = hrz_posn()
    b = ephem.FixedBody()
    b._ra = deg_to_rad(target.ra)
    b._dec = deg_to_rad(target.dec)
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
    b.compute(o)
    az = rad_to_deg(b.az)
    alt = rad_to_deg(b.alt)
    
    _posn.az = az
    _posn.alt = alt
    return _posn


def get_equ_from_hrz(target, observer, jD):
    """
    Get equatorial/celestial coordinates from local horizontal coordinates.
    
    Param: target   - Object of type hrz_posn representing horizontal position.
    Param: observer - Object of type lnlat_posn representing observer position.
    Param: jD       - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    ra,dec = o.radec_of(deg_to_rad(target.az), deg_to_rad(target.alt))
    ra = rad_to_deg(ra)
    dec = rad_to_deg(dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


def get_ecl_from_rect(rect):
    """
    Get ecliptical coordinates from rectangular coordinates.
    
    Param: rect - Object of type rect_posn representing position.
    
    Returns object of type lnlat_posn representing ecliptical position.
    """
    
    _posn = ecl_posn()
    _posn.lng = 1*rect.lng
    _posn.lat = 1*rect.lat
    return _posn


def get_equ_from_ecl(target, jD):
    """
    Get equatorial coordinates from ecliptical coordinates for a given time.
    
    Param: target   - Object of type lnlat_posn representing ecliptic position.
    Param: jD       - UTC Julian day (float). 
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    ecl = ephem.Ecliptic(deg_to_rad(target.lng), deg_to_rad(target.lat), epoch=ephem.B1950)
    equ = ephem.Equatorial(ecl)
    ra = rad_to_deg(equ.ra)
    dec = rad_to_deg(equ.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


def get_ecl_from_equ(target, jD):
    """
    Get ecliptical coordinates from equatorial coordinates for a given time.
    
    Param: target  - Object of type equ_posn representing equatorial position.
    Param: jD       - UTC Julian day (float). 
    
    Returns object of type ecl_posn representing ecliptic position.
    """
    
    _posn = ecl_posn()
    equ = ephem.Equatorial(deg_to_rad(target.ra), deg_to_rad(target.dec), epoch=ephem.B1950)
    ecl = ephem.Ecliptic(equ)
    l = rad_to_deg(ecl.lon)
    b = rad_to_deg(ecl.lat)
    
    _posn.lng = l
    _posn.lat = b
    return _posn    


def get_equ_from_gal(target):
    """
    Get B1950 equatorial coordinates from galactic coordinates.
    
    Param: target - Object of type gal_posn representing galactic position.
    
    Returns object of type equ_posn representing B1950 equatorial position.
    """
    
    _posn = equ_posn()
    gal = ephem.Galactic(deg_to_rad(target.l), deg_to_rad(target.b))
    equ = ephem.Equatorial(gal)
    equ = ephem.Equatorial(equ, epoch=ephem.B1950)
    ra = rad_to_deg(equ.ra)
    dec = rad_to_deg(equ.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


def get_gal_from_equ(target):
    """
    Get galactic coordinates from B1950 equatorial coordinates.
    
    Param: target - Object of type equ_posn representing B1950 equatorial 
                    position.
    
    Returns object of type gal_posn representing galactic position.
    """
    
    _posn = gal_posn()
    equ = ephem.Equatorial(deg_to_rad(target.ra), deg_to_rad(target.dec), epoch=ephem.B1950)
    ecl = ephem.Galactic(equ)
    l = rad_to_deg(ecl.lon)
    b = rad_to_deg(ecl.lat)
    
    _posn.l = l
    _posn.b = b
    return _posn


def get_equ2000_from_gal(target):
    """
    Get J2000 equatorial coordinates from galactic coordinates.
    
    Param: target - Object of type gal_posn representing galactic position.
    
    Returns object of type equ_posn representing J2000 equatorial position.
    """
    
    _posn = equ_posn()
    gal = ephem.Galactic(deg_to_rad(target.l), deg_to_rad(target.b))
    equ = ephem.Equatorial(gal, epoch=ephem.J2000)
    ra = rad_to_deg(equ.ra)
    dec = rad_to_deg(equ.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


def get_gal_from_equ2000(target):
    """
    Get galactic coordinates from J2000 equatorial coordinates.
    
    Param: target - Object of type equ_posn representing J2000 equatorial 
                    position.
    
    Returns object of type gal_posn representing galactic position.
    """
    
    _posn = gal_posn()
    equ = ephem.Equatorial(deg_to_rad(target.ra), deg_to_rad(target.dec), epoch=ephem.J2000)
    ecl = ephem.Galactic(equ)
    l = rad_to_deg(ecl.lon)
    b = rad_to_deg(ecl.lat)
    
    _posn.l = l
    _posn.b = b
    return _posn


######################################################################
# Sidereal Time Functions
######################################################################

def get_apparent_sidereal_time(jD):
    """
    Get apparent sidereal time from Julian day.
    
    Param: jD - UTC Julian day (float).
    
    Returns GM apparent sidereal time (float hours).
    
    From: http://aa.usno.navy.mil/faq/docs/GAST.php
    """
    
    gmst = get_mean_sidereal_time(jD)
    
    D = jD - 2451545.0
    Omega = 125.04 - 0.052954*D
    L = 280.47 + 0.98565*D
    epsilon = 23.4393 - 0.0000004*D
    deltaPhi = -0.000319*math.sin(deg_to_rad(Omega)) - 0.000024*math.sin(2*deg_to_rad(L))
    
    eqeq = deltaPhi*math.cos(deg_to_rad(epsilon))
    
    gast = gmst + eqeq
    return gast


def get_mean_sidereal_time(jD):
    """
    Get mean sidereal time from Julian day.
    
    Param: jD - UTC Julian day (float).
    
    Returns GM mean sidereal time (float hours).
    
    From: http://aa.usno.navy.mil/faq/docs/GAST.php
    """
    
    D = jD - 2451545.0
    D0 = (int(jD) - 0.5) - 2451545.0
    H = (D - D0) * 24.0
    T = D / 36525.0
    
    gmst =  6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T*T
    gmst %= 24
    return gmst



######################################################################
# Angular Separation Functions
######################################################################

def get_angular_separation(posn1, posn2):
    """  
    Get angular separation from equatorial positions.
    
    Param: posn1 - Object of type equ_posn representing body 1 position.
    Param: posn2 - Object of type equ_posn representing body 2 position.
    
    Returns angular separation in degrees (float).
    """
    
    sep = ephem.separation((deg_to_rad(posn1.ra), deg_to_rad(posn1.dec)), (deg_to_rad(posn2.ra), deg_to_rad(posn2.dec)))
    sep = rad_to_deg(sep)
    
    return sep


def get_rel_posn_angle(posn1, posn2):
    """
    Get relative position angle from equatorial positions.
    
    Param: posn1 - Object of type equ_posn representing body 1 position.
    Param: posn2 - Object of type equ_posn representing body 2 position.
    
    Returns position angle in degrees (float).
    
    Based on dpav.f from SLALIB.
    """
    
    d1 = dir_cos(posn1)
    d2 = dir_cos(posn2)
    
    w1 = math.sqrt( d1[0]**2 + d1[1]**2 + d1[2]**2 )
    if w1 != 0:
        d1[0] /= w1
        d1[1] /= w1
        d1[2] /= w1
    w2 = math.sqrt( d2[0]**2 + d2[1]**2 + d2[2]**2 )
    if w2 != 0:
        d2[0] /= w2
        d2[1] /= w2
        d2[2] /= w2
        
    sq = d2[1]*d1[0] - d2[0]*d1[1]
    cq = d2[2]*(d1[0]**2+d1[1]**2) - d1[2]*(d2[0]*d1[0]+d2[1]*d1[1])
    if sq == 0 and cq == 0:
        cq = 1.0
    ang = math.atan2(sq, cq)
    ang = rad_to_deg(ang)
    return ang


######################################################################
# Apparent Position Functions
######################################################################

_DEFAULT_PROPER_MOTION = equ_posn(0.0, 0.0)

def get_apparent_posn(mean_position, jD, proper_motion = None):
    """
    Get apparent position of celestial object accounting for precession, nutation, 
    aberration, and optionally proper motion.
    
    Param: mean_position  - J2000 equatorial mean position of object as type 
                            equ_posn.
    Param: jD             - UTC Julian day (float) to measure position.
    Param: proper_motion  - object of type equ_posn giving object's proper motion
                            (optional).
    
    Returns: Apparent equatorial position of object as type equ_posn.
    """
    
    if proper_motion is None:
        proper_motion = _DEFAULT_PROPER_MOTION  
        
    _posn = get_equ_pm(mean_position, proper_motion, jD)
    _posn = get_equ_aber(_posn, jD)
    return get_equ_prec(_posn, jD)


######################################################################
# Precession Functions
######################################################################

def get_equ_prec(mean_position, jD):
    """
    Get position of celestial object accounting for precession.
    Only works for converting to and from J2000 epoch.
    
    Param: mean_position - J2000 equatorial mean position of object as type equ_posn.
    Param: jD - UTC Julian day (float) to measure position.
    
    Returns: Adjusted equatorial position of object as type equ_posn.
    """    
    
    _posn = equ_posn() 
    equ = ephem.Equatorial(deg_to_rad(mean_position.ra), deg_to_rad(mean_position.dec), epoch=ephem.J2000)
    equ = ephem.Equatorial(equ, epoch=jD-DJD_OFFSET)
    ra = rad_to_deg(equ.ra)
    dec = rad_to_deg(equ.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


def get_equ_prec2(mean_position, fromJD, toJD):
    """
    Get position of celestial object accounting for precession.
    
    Param: mean_position  - equatorial first position of object as type equ_posn.
    Param: fromJD         - UTC Julian day (float) of first time.
    Param: toJD           - UTC Julian day (float) of second time.
    
    Returns: Equatorial position of object as type equ_posn converted from
            time 1 to time 2.
    """  
    
    _posn = equ_posn() 
    equ = ephem.Equatorial(deg_to_rad(mean_position.ra), deg_to_rad(mean_position.dec), epoch=fromJD-DJD_OFFSET)
    equ = ephem.Equatorial(equ, epoch=toJD-DJD_OFFSET)
    ra = rad_to_deg(equ.ra)
    dec = rad_to_deg(equ.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


######################################################################
# Nutation Functions
######################################################################

def get_nutation(jD):
    """
    Get nutation corrections for a given time.
    
    Param: jD - UTC Julian day (float) to measure nutation.
    
    Returns: Nutation corrections as object of type nutation.
    
    Based on the nutate.pro and co_nutate.pro from the AstroIDL
    library.
    """
    
    # form time in Julian centuries from 1900.0
    t = (jD - 2451545.0) / 36525.
    
    # Mean elongation of the Moon
    coeff1 = numpy.array([297.85036, 445267.111480, -0.0019142, 1.0/189474])
    d = deg_to_rad(numpy.polyval(coeff1[::-1], t))
    d %= (2*math.pi)
    
    # Sun's mean anomaly
    coeff2 = numpy.array([357.52772, 35999.050340, -0.0001603, -1.0/3e5])
    M = deg_to_rad(numpy.polyval(coeff2[::-1], t))
    M %= (2*math.pi)
    
    # Moon's mean anomaly
    coeff3 = numpy.array([134.96298, 477198.867398, 0.0086972, 1.0/5.625e4])
    Mprime = deg_to_rad(numpy.polyval(coeff3[::-1], t))
    Mprime %= (2*math.pi)
    
    # Moon's argument of latitude
    coeff4 = numpy.array([93.27191, 483202.017538, -0.0036825, -1.0/3.27270e5])
    F = deg_to_rad(numpy.polyval(coeff4[::-1], t)) 
    F %= (2*math.pi)
    
    # Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
    # measured from the mean equinox of the date
    coeff5 = numpy.array([125.04452, -1934.136261, 0.0020708, 1.0/4.5e5])
    omega = deg_to_rad(numpy.polyval(coeff5[::-1], t))
    omega %= (2*math.pi)
    
    d_lng = numpy.array([0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,2,0,2,0,0,-2,0,2,0,0,-2,
                    0,-2,0,0,2,-2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,
                    0,-2,-2,0,-1,-2,1,0,0,-1,0,0,2,0,2])
                    
    m_lng = numpy.array([0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    2,0,2,1,0,-1,0,0,0,1,1,-1,0,0,0,0,0,0,-1,-1,0,0,0,1,0,0,
                    1,0,0,0,-1,1,-1,-1,0,-1])
                    
    mp_lng = numpy.array([0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,2,2,1,0,0,
                    -1,0,-1,0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,
                    1,0,0,1,-2,1,1,1,-1,3,0])
                    
    f_lng = numpy.array([0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2,2,2,0,0,
                    2,0,0,0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,
                    0,0,2,2,2,2])
                    
    om_lng = numpy.array([1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,2,0,0,1,0,1,
                    2,1,1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,
                    0,2,2,2,2])
                    
    sin_lng = numpy.array([-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, 
                        -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 
                        26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, 
                        -7, 6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3])
                        
    sdelt = numpy.array([-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1, 
                    0,0,0.1, 0,-0.1,0,0,0,0,0,0,0,0,0,0, -0.1, 0, 0.1,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) 
                    
    cos_lng = numpy.array([92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,0,-70,-53,0, 
                        -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,7,6,0,5,3,
                        -3,0,3,3,0,-3,-3,3,3,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
                        
    cdelt = numpy.array([8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
                    
    # Sum the periodic terms 
    arg = d_lng*d + m_lng*M + mp_lng*Mprime + f_lng*F + om_lng*omega
    sarg = numpy.sin(arg)
    carg = numpy.cos(arg)
    lng = 0.0001*( (sdelt*t + sin_lng)*sarg ).sum()
    obl = 0.0001*( (cdelt*t + cos_lng)*carg ).sum()

    T = (jD -2451545.0) / 36525.0
    ecl = 23.4392911*3600. - 46.8150*T - 0.00059*T*T + 0.001813*T*T*T
    ecl = (ecl + obl)/3600
    
    _nut = nutation()
    _nut.longitude = lng / 3600.0
    _nut.obliquity = obl / 3600.0
    _nut.ecliptic = ecl
    return _nut


def get_equ_nut(position, jD):
    """
    Get the position of a celesital object accounting for nutation.
    
    Param: mean_position  -  Equatorial position of object as type 
                             equ_posn.
    Param: jD             - UTC Julian day (float) to measure nutation.
    
    Returns: Adjusted equatorial position of object as type equ_posn.
    
    Based on the AstroIDL co_nutate.pro procedure
    """    
    
    # Get the nutation
    nut = get_nutation(jD)
    
    # Convert RA/dec into cartesian
    ra  = deg_to_rad(position.ra)
    dec = deg_to_rad(position.dec)
    x = numpy.cos(dec) * numpy.cos(ra)
    y = numpy.cos(dec) * numpy.sin(ra)
    z = numpy.sin(dec)
    
    # Apply the nutation
    ecl = deg_to_rad(nut.ecliptic)
    obl = deg_to_rad(nut.obliquity)
    lng = deg_to_rad(nut.longitude)
    
    x2 = x - (y*numpy.cos(ecl)*lng + z*numpy.sin(ecl)*lng)
    y2 = y + (x*numpy.cos(ecl)*lng - z*obl)
    z2 = z + (x*numpy.sin(ecl)*lng + y*obl)
    
    # Back to RA/dec
    r = numpy.sqrt(x2**2 + y2**2 + z2**2)
    xyproj = numpy.sqrt(x2**2 + y2**2)
    
    ra, dec = 0.0, 0.0
    if xyproj == 0.0 and z != 0.0:
        ra = 0.0
        dec = numpy.arcsin(y2/x2)
    if xyproj != 0.0:
        ra = numpy.arctan2(y2, x2)
        dec = numpy.arcsin(z2/r)
        
    # Create the output object and update it
    _posn = equ_posn()
    _posn.ra = rad_to_deg(ra)
    _posn.dec = rad_to_deg(dec)
    
    return _posn


######################################################################
# Aberration Functions
######################################################################


def get_equ_aber(mean_position, jD): 
    """
    Get position of celestial object accounting for aberration.
    
    Param: mean_position  - J2000 equatorial mean position of object as type 
                            equ_posn.
    Param: jD             - UTC Julian day (float) to measure aberration.
    
    Returns: Adjusted equatorial position of object as type equ_posn.
    
    Based on the libnova ln_get_equ_aber() function.
    """    
    
    _posn = equ_posn()
    # speed of light in 10-8 au per day
    c = 17314463350.0
    
    # calc T
    T = (jD - 2451545.0) / 36525.0
    
    # calc planetary perturbutions
    L2 = 3.1761467 + 1021.3285546 * T
    L3 = 1.7534703 + 628.3075849 * T
    L4 = 6.2034809 + 334.0612431 * T
    L5 = 0.5995464 + 52.9690965 * T
    L6 = 0.8740168 + 21.329909095 * T
    L7 = 5.4812939 + 7.4781599 * T
    L8 = 5.3118863 + 3.8133036 * T
    LL = 3.8103444 + 8399.6847337 * T
    D = 5.1984667 + 7771.3771486 * T
    MM = 2.3555559 + 8328.6914289 * T
    F = 1.6279052 + 8433.4661601 * T
    
    X = 0.0
    Y = 0.0
    Z = 0.0
    
    # terms
    TERMS = 36
    
    arguments = [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
                [0, 2, 0, -1, 0, 0, 0, 0, 0, 0, 0],
                [0, 3, -8, 3, 0, 0, 0, 0, 0, 0, 0],
                [0, 5, -8, 3, 0, 0, 0, 0, 0, 0, 0],
                [2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0],
                [0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 3, 0, -2, 0, 0, 0, 0, 0, 0, 0],
                [1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
                [2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 3, -2, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 2, -1, 0],
                [8, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [8, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
                [3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0],
                [3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, -2, 0, 0]]
                
    x_coefficients = [[-1719914, -2, -25, 0],
                    [6434, 141, 28007, -107],
                    [715, 0, 0, 0],
                    [715, 0, 0, 0],
                    [486, -5, -236, -4],
                    [159, 0, 0, 0],
                    [0, 0, 0, 0],
                    [39, 0, 0, 0],
                    [33, 0, -10, 0],
                    [31, 0, 1, 0],
                    [8, 0, -28, 0],
                    [8, 0, -28, 0],
                    [21, 0, 0, 0],
                    [-19, 0, 0, 0],
                    [17, 0, 0, 0],
                    [16, 0, 0, 0],
                    [16, 0, 0, 0],
                    [11, 0, -1, 0],
                    [0, 0, -11, 0],
                    [-11, 0, -2, 0],
                    [-7, 0, -8, 0],
                    [-10, 0, 0, 0],
                    [-9, 0, 0, 0], 
                    [-9, 0, 0, 0],
                    [0, 0, -9, 0],
                    [0, 0, -9, 0],
                    [8, 0, 0, 0],
                    [8, 0, 0, 0], 
                    [-4, 0, -7, 0],
                    [-4, 0, -7, 0],
                    [-6, 0, -5, 0],
                    [-1, 0, -1, 0],
                    [4, 0, -6, 0],
                    [0, 0, -7, 0],
                    [5, 0, -5, 0],
                    [5, 0, 0, 0]]
                    
    y_coefficients = [[25, -13, 1578089, 156],
                    [25697, -95, -5904, -130],
                    [6, 0, -657, 0], 
                    [0, 0, -656, 0],
                    [-216, -4, -446, 5],
                    [2, 0, -147, 0],
                    [0, 0, 26, 0],
                    [0, 0, -36, 0],
                    [-9, 0, -30, 0],
                    [1, 0, -28, 0],
                    [25, 0, 8, 0],
                    [-25, 0, -8, 0],
                    [0, 0, -19, 0],
                    [0, 0, 17, 0],
                    [0, 0, -16, 0],
                    [0, 0, 15, 0],
                    [1, 0, -15, 0],
                    [-1, 0, -10, 0],
                    [-10, 0, 0, 0],
                    [-2, 0, 9, 0],
                    [-8, 0, 6, 0], 
                    [0, 0, 9, 0], 
                    [0, 0, -9, 0], 
                    [0, 0, -8, 0],
                    [-8, 0, 0, 0],
                    [8, 0, 0, 0],
                    [0, 0, -8, 0],
                    [0, 0, -7, 0], 
                    [-6, 0, -4, 0],
                    [6, 0, -4, 0],
                    [-4, 0, 5, 0],
                    [-2, 0, -7, 0],
                    [-5, 0, -4, 0],
                    [-6, 0, 0, 0], 
                    [-4, 0, -5, 0],
                    [0, 0, -5, 0]]
                    
    z_coefficients = [[10, 32, 684185, -358],
                    [11141, -48, -2559, -55],
                    [-15, 0, -282, 0],
                    [0, 0, -285, 0],
                    [-94, 0, -193, 0],
                    [-6, 0, -61, 0],
                    [0, 0, 59, 0],
                    [0, 0, 16, 0],
                    [-5, 0, -13, 0],
                    [0, 0, -12, 0],
                    [11, 0, 3, 0],
                    [-11, 0, -3, 0],
                    [0, 0, -8, 0],
                    [0, 0, 8, 0],
                    [0, 0, -7, 0],
                    [1, 0, 7, 0],
                    [-3, 0, -6, 0],
                    [-1, 0, 5, 0],
                    [-4, 0, 0, 0],
                    [-1, 0, 4, 0],
                    [-3, 0, 3, 0],
                    [0, 0, 4, 0],
                    [0, 0, -4, 0],
                    [0, 0, -4, 0],
                    [-3, 0, 0, 0],
                    [3, 0, 0, 0],
                    [0, 0, -3, 0],
                    [0, 0, -3, 0],
                    [-3, 0, 2, 0],
                    [3, 0, -2, 0],
                    [-2, 0, 2, 0],
                    [1, 0, -4, 0],
                    [-2, 0, -2, 0],
                    [-3, 0, 0, 0],
                    [-2, 0, -2, 0],
                    [0, 0, -2, 0]]
                    
    # Sum the terms
    for i in range(TERMS):
        A = arguments[i][0]*L2 + arguments[i][1]*L3 + arguments[i][2]*L4 + arguments[i][3]*L5 + arguments[i][4]*L6 + \
                arguments[i][5]*L7 + arguments[i][6]*L8 + arguments[i][7]*LL + arguments[i][8]*D + arguments[i][9]*MM + \
                arguments[i][10]*F
        X += (x_coefficients[i][0] + x_coefficients[i][1]*T) * math.sin(A) + (x_coefficients[i][2] + x_coefficients[i][3]*T) * math.cos(A)
        Y += (y_coefficients[i][0] + y_coefficients[i][1]*T) * math.sin(A) + (y_coefficients[i][2] + y_coefficients[i][3]*T) * math.cos(A)
        Z += (z_coefficients[i][0] + z_coefficients[i][1]*T) * math.sin(A) + (z_coefficients[i][2] + z_coefficients[i][3]*T) * math.cos(A)
        
    # Equ 22.4
    mean_ra = deg_to_rad(mean_position.ra)
    mean_dec = deg_to_rad(mean_position.dec)
    
    delta_ra = (Y * math.cos(mean_ra) - X * math.sin(mean_ra)) / (c * math.cos(mean_dec))
    delta_dec = (X * math.cos(mean_ra) + Y * math.sin(mean_ra)) * math.sin(mean_dec) - Z * math.cos(mean_dec)
    delta_dec /= -c
    
    _posn.ra = rad_to_deg(mean_ra + delta_ra)
    _posn.dec = rad_to_deg(mean_dec + delta_dec)
    return _posn


######################################################################
# Proper Motion Functions
######################################################################   

def get_equ_pm(mean_position, proper_motion, jD):
    """
    Adjusts equatorial position of a stellar object accouting for proper motion.
    
    Param: mean_position - J2000 equatorial mean position of object as type equ_posn.
    Param: proper_motion - Object of type equ_posn giving object's proper motion
                           (units are deg/year).
    Param: jD - UTC Julian day (float) to measure position.
    
    Returns: Adjusted equatorial position of object as type equ_posn.
    
    Based on pm.f, dcs2c.f, and dcc2s.f from SLALIB.
    """
    
    _posn = equ_posn()
    ra = ephem.degrees(str(mean_position.ra))
    dec = ephem.degrees(str(mean_position.dec))
    pmRA = ephem.degrees(str(proper_motion.ra))
    pmDec = ephem.degrees(str(proper_motion.dec))
    
    p = [math.cos(ra)*math.cos(dec), 
        math.sin(ra)*math.cos(dec), 
        math.sin(dec)]
    em = [-pmRA*p[1] - pmDec*math.cos(ra)*math.sin(dec), 
        pmRA*p[0]  - pmDec*math.sin(ra)*math.sin(dec), 
        pmDec*math.cos(dec)]
        
    t = jD - (ephem.J2000+DJD_OFFSET)
    p[0] += em[0]*t
    p[1] += em[1]*t
    p[2] += em[2]*t
    
    r = math.sqrt(p[0]**2 + p[1]**2)
    if r == 0:
        ra = 0.0 % (2*math.pi)
    else:
        ra = math.atan2(p[1], p[0])
    if p[2] == 0:
        dec = 0.0
    else:
        dec = math.atan2(p[2], r)
        
    ra = rad_to_deg(ra)
    dec = rad_to_deg(dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


######################################################################
# Rise, Set, Transit functions
######################################################################

def get_object_rst(jD, observer, target):
    """
    Get rise, set, and transit times of a celstial object.
    
    Param: jD       - UTC Julian day (float) target time.
    Param: observer - object of type lnlat_posn giving observer position
    Param: target   - object of type equ_posn giving target equatorial position
    
    Returns: Object of type rst_time giving object's ephemeris UTC times,
            or None if the object is circumpolar.
    """

    _rst = rst_time()
    b = ephem.FixedBody()
    b._ra = deg_to_rad(target.ra)
    b._dec = deg_to_rad(target.dec)
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
    
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


######################################################################
#  Solar Functions
######################################################################

def get_solar_equ_coords(jD):
    """
    Get Sun's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, and nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    b = ephem.Sun()
    b.compute(jD-DJD_OFFSET)
    
    _posn.ra = rad_to_deg(b.g_ra)
    _posn.dec = rad_to_deg(b.g_dec)
    return _posn


def get_solar_rst(jD, observer):
    """
    Get Sun's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type lnlat_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar..
    """
    
    _rst = rst_time()
    b = ephem.Sun()
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


######################################################################
#  Jupiter Functions
######################################################################

def get_jupiter_equ_coords(jD):
    """
    Get Jupiter's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    b = ephem.Jupiter()
    b.compute(jD-DJD_OFFSET)
    
    _posn.ra = rad_to_deg(b.g_ra)
    _posn.dec = rad_to_deg(b.g_dec)
    return _posn


def get_jupiter_rst(jD, observer):
    """
    Get Jupiter's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type lnlat_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    _rst = rst_time()
    b = ephem.Jupiter()
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


######################################################################
# Saturn Functions
######################################################################

def get_saturn_equ_coords(jD):
    """   
    Get Saturn's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    b = ephem.Saturn()
    b.compute(jD-DJD_OFFSET)
    
    _posn.ra = rad_to_deg(b.g_ra)
    _posn.dec = rad_to_deg(b.g_dec)
    return _posn


def get_saturn_rst(jD, observer):
    """
    Get Saturn's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type lnlat_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """

    _rst = rst_time()
    b = ephem.Saturn()
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


######################################################################
# Lunar Functions
######################################################################

def get_lunar_equ_coords(jD):
    """
    Get the Moon's apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    b = ephem.Moon()
    b.compute(jD-DJD_OFFSET)
    
    _posn.ra = rad_to_deg(b.g_ra)
    _posn.dec = rad_to_deg(b.g_dec)
    return _posn


def get_lunar_rst(jD, observer):
    """
    Get the Moon's rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type lnlat_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    _rst = rst_time()
    b = ephem.Moon()
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


######################################################################
# Venus Functions
######################################################################

def get_venus_equ_coords(jD):
    """
    Get Venus' apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    b = ephem.Venus()
    b.compute(jD-DJD_OFFSET)
    
    _posn.ra = rad_to_deg(b.g_ra)
    _posn.dec = rad_to_deg(b.g_dec)
    return _posn


def get_venus_rst(jD, observer):
    """
    Get Venus' rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type lnlat_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    _rst = rst_time()
    b = ephem.Venus()
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


######################################################################
# Mars Functions
######################################################################

def get_mars_equ_coords(jD):
    """
    Get Mars' apparent equatorial coordinates from Julian day.
    Accounts for aberration and precession, but not nutation.
    
    Param: jD - UTC Julian day (float).
    
    Returns object of type equ_posn representing equatorial position.
    """
    
    _posn = equ_posn()
    b = ephem.Mars()
    b.compute(jD-DJD_OFFSET)
    
    _posn.ra = rad_to_deg(b.g_ra)
    _posn.dec = rad_to_deg(b.g_dec)
    return _posn


def get_mars_rst(jD, observer):
    """
    Get Mars' rise, transit, set times from Julian day.
    
    Param: jD       - UTC Julian day (float).
    Param: observer - Object of type lnlat_posn representing observer position.
    
    Returns Object of type rst_time represeting UTC ephemeris times,
            or None if the object is circumpolar.
    """
    
    _rst = rst_time()
    b = ephem.Mars()
    
    o = ephem.Observer()
    o.date = jD-DJD_OFFSET
    o.lon = deg_to_rad(observer.lng)
    o.lat = deg_to_rad(observer.lat)
    try:
        o.elev = observer.elv
    except:
        pass
        
    b.compute(o)
    
    try:
        _rst.rise = o.next_rising(b)+DJD_OFFSET
        _rst.transit = o.next_transit(b)+DJD_OFFSET
        _rst.set = o.next_setting(b)+DJD_OFFSET
        return _rst
    except (ephem.NeverUpError, ephem.AlwaysUpError):
        return None


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


def sec_to_jd(secs):
    """
    Convert seconds into julian days.
    
    Param: secs - seconds (float)
    
    Returns: Julian days as a float. 
    """
    
    return secs / SECS_IN_DAY


def jd_to_sec(jD):
    """
    Convert Julian days into seconds.
    
    Param: jD - Julian days (float).
    
    Returns: Seconds as a float.
    """
    
    return jD * SECS_IN_DAY    


######################################################################
# This code is run at import time to load the UTC/TAI almanac data 
# The information is provided from the USNO data file available from
# site location http://maia.usno.navy.mil/ser7/tai-utc.dat.
######################################################################

# Create the cache directory
if not os.path.exists(os.path.join(os.path.expanduser('~'), '.lsl')):
    os.mkdir(os.path.join(os.path.expanduser('~'), '.lsl'))
_CACHE_DIR = os.path.join(os.path.expanduser('~'), '.lsl', 'astro_cache')
if not os.path.exists(_CACHE_DIR):
    os.mkdir(_CACHE_DIR)

######################################################################
#
# Create empty UTC->leap second list.  The function _parse_tai_file()
# will populate this list with tuples of three elements.  The three
# tuple elements (utc_jd, leap_sec, leap_jd) are:
#
#   utc_jd    - The UTC Julian day time at which a leap second 
#               adjustmentwas made
#   leap_sec  - The cumulative number of leap seconds to add to the
#               UTC value to get TAI
#   leap_jd   - The 'leap_sec' value converted to Julian days
#
######################################################################

FIRST_LEAP_UTC = 2441317.5

_LEAP_SEC_LIST = []

def _parse_tai_file():
    # get path to almanac data file
    download = True
    datName = os.path.join(_CACHE_DIR, 'Leap_Second.dat')
    if not os.path.exists(datName):
        from lsl.common.paths import DATA
        oldName = os.path.join(DATA, 'astro', 'Leap_Second.dat')
        with open(oldName, 'rb') as oh:
            with open(datName, 'wb') as dh:
                dh.write(oh.read())
                
    # check for expiration
    with open(datName, 'r') as datFile:
        for l in datFile:
            if l.find('File expires on') != -1:
                expDate = l.strip().rstrip().split('on ', 1)[1]
                expDate = datetime.strptime(expDate, "%d %B %Y")
                if datetime.utcnow() < expDate:
                    download = False
                    
    # download as needed
    if download:
        url = "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat"
        
        print("Downloading %s" % url)
        lsFH = urlopen(url, timeout=120)
        meta = lsFH.info()
        try:
            remote_size = int(meta.getheaders("Content-Length")[0])
        except AttributeError:
            remote_size = 1
        pbar = DownloadBar(max=remote_size)
        
         -Length")[0]))
        while True:
            new_data = lsFH.read(32768)
            if len(new_data) == 0:
                break
            pbar.inc(len(new_data))
            try:
                data += new_data
            except NameError:
                data = new_data
            sys.stdout.write(pbar.show()+'\r')
            sys.stdout.flush()
        lsFH.close()
        sys.stdout.write(pbar.show()+'\n')
        sys.stdout.flush()
        
        with open(datName, 'wb') as fh:
            fh.write(data)
            
    # read tai-utc.dat file to get conversion info
    with open(datName, 'r') as datFile:
        lineNum = 0
        for l in datFile:
            if l.startswith('#'):
                continue
            elif len(l) < 3:
                continue
                
            # split
            utcMJD, day, month, year, leapSec = l.split(None, 4)
            
            # get UTC JD of leap second boundaries
            try:
                utcJD = float(utcMJD) + MJD_OFFSET
            except ValueError:
                raise RuntimeError("line %d of %s file not correctly formatted" % (lineNum, datName))
                
            # only get values prior to UTC JD 2441317.5 (1972 JAN  1)
            if utcJD < FIRST_LEAP_UTC:
                lineNum += 1
                continue
                
            # get leap second asjustment value
            try:
                leapSec = float(leapSec)
            except ValueError:
                raise RuntimeError("line %d of %s file not correctly formatted" % (lineNum, datName))
                
            # add entry to list
            _LEAP_SEC_LIST.append((utcJD, leapSec, sec_to_jd(leapSec)))
            lineNum += 1


_parse_tai_file()


######################################################################
# Time utility functions
######################################################################

def range_hours(hours):
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


def jd_to_mjd(jd):
    """
    Get modified julian day value from julian day value.
    
    Param: jd - julian day (should be >= 2400000.5) 
    
    Returns: Modified julian day.
    """
    
    if jd < MJD_OFFSET:
        raise ValueError("jd must be >= %f" % MJD_OFFSET)
        
    return (jd - MJD_OFFSET)


def mjd_to_jd(mjd):
    """
    Get julian day value from modified julian day value.
    
    Param: mjd - modified julian day (should be >= 0.0)
    
    Returns: Julian day.
    """
    
    if mjd < 0.0:
        raise ValueError("mjd must be >= 0.0")
        
    return (mjd + MJD_OFFSET)


def leap_secs(utcJD):
    """
    Get the number of leap seconds for given UTC time value.
    
    Param: utcJD - The UTC JD time.
                This should be greater than 2441317.5 (1972 JAN  1). 
    
    Returns: The number of leap seconds (float) for the UTC time.
    """
    
    # check the last entry first since it is likely the time is
    # current or future
    p = _LEAP_SEC_LIST[-1]
    if utcJD >= p[0]:
        return p[1]
        
    if utcJD < FIRST_LEAP_UTC:
        raise ValueError("utcJD must be greater than %f" % FIRST_LEAP_UTC)
        
    # search the conversion list for the UTC JD range
    p = _LEAP_SEC_LIST[0]
    for e in _LEAP_SEC_LIST[1:]:
        if utcJD < e[0]:
            return p[1]
        p = e
        
    # the time is after the last conversion range entry
    p = _LEAP_SEC_LIST[-1]
    return p[1]    


def utc_to_tai(utcJD):
    """
    Get the TAI JD time value for a given UTC JD time value.
    
    Param: utcJD - The UTC JD time (float).
                This should be greater than 2441317.5 (1972 JAN  1). 
    
    Returns: The TAI JD value (float).
    """
    
    # check the last entry first since it is likely the time is
    # current or future
    p = _LEAP_SEC_LIST[-1]
    if utcJD >= p[0]:
        return utcJD + p[2]
        
    if utcJD < FIRST_LEAP_UTC:
        raise ValueError("utcJD must be greater than %f" % FIRST_LEAP_UTC)
        
    # search the conversion list for the UTC JD range
    p = _LEAP_SEC_LIST[0]
    for e in _LEAP_SEC_LIST[1:]:
        if utcJD < e[0]:
            return utcJD + p[2]
        p = e
        
    # the time is after the last conversion range entry
    p = _LEAP_SEC_LIST[-1]
    return utcJD + p[2]


def tai_to_utc(taiJD):
    """
    Get the UTC JD time value for a given TAI JD time value.
    
    Param: taiJD - The TAI JD time (float).
                This should be greater than 2441317.5 (1972 JAN  1).
    
    Returns: The UTC JD value (float).
    """
    
    # check the last entry first since it is likely the time is
    # current or future
    p = _LEAP_SEC_LIST[-1]
    tai = p[0] + p[2]
    if taiJD >= tai:
        return taiJD - p[2]
        
    p = _LEAP_SEC_LIST[0]
    firstTai = p[0] + p[2]
    if taiJD < firstTai:
        raise ValueError("taiJD must be greater than %f" % firstTai)
        
    # search the conversion list for the TAI JD range
    p = _LEAP_SEC_LIST[0]
    for e in _LEAP_SEC_LIST[1:]:
        tai = e[0] + e[2]
        if taiJD < tai:
            return taiJD - p[2] 
        p = e
        
    # the time is after the last conversion range entry
    p = _LEAP_SEC_LIST[-1]
    return taiJD - p[2]


def taimjd_to_utcjd(taiMJD):
    """
    Get the UTC JD time value for a given TAI MJD value.
    
    Param: mjdTAI - The TAI MJD time (float).
    
    Returns: The UTC JD value (float).
    """
    
    jd = mjd_to_jd(taiMJD)
    return tai_to_utc(jd)


def utcjd_to_taimjd(utcJD):
    """
    Get the TAI MJD time value for a given UTC JD value.
    
    Param: utcJD - The UTC JD time (float).
    
    Returns: The TAI MJD value.
    """
    
    tai = utc_to_tai(utcJD)
    return jd_to_mjd(tai)


def unix_to_utcjd(unixTime):
    """
    Get the UTC JD time value for a given UNIX time value.
    
    Param: unixTime - the UNIX time (int/float)
    
    Returns: The UTC JD value.
    """
    
    utcJD = float(unixTime) / SECS_IN_DAY + UNIX_OFFSET
    return utcJD


def unix_to_taimjd(unixTime):
    """
    Get the TAI MJD time value for a given UNIX time value.
    
    Param: unixTime - the UNIX time (int/float)
    
    Returns: The TAI MJD value.
    """
    
    utcJD = unix_to_utcjd(unixTime)
    taiMJD = utcjd_to_taimjd(utcJD)
    return taiMJD


def utcjd_to_unix(utcJD):
    """
    Get UNIX time value for a given UTC JD value.
    
    Param: utcJD - The UTC JD time (float).
    
    Returns: The UNIX time
    """
    
    unixTime = (utcJD - UNIX_OFFSET) * SECS_IN_DAY
    return unixTime


def taimjd_to_unix(taiMJD):
    """
    Get UNIX time value for a given TAI MJDvalue.
    
    Param: taiMJD - The TAI MJD time (float).
    
    Returns: The UNIX time
    """
    
    utcJD = taimjd_to_utcjd(taiMJD)
    unixTime = utcjd_to_unix(utcJD)
    return unixTime


def tai_to_tt(taiJD):
    """
    Get the TT JD time value for a given TAI JD time value.
    
    Param: taiJD - The TAI JD time (float).
    
    Returns: The TT JD value (float).
    """
    
    return taiJD + sec_to_jd(TAI_TT_OFFSET)


def tt_to_tai(ttJD):
    """
    Get the TAI JD time value for a given TT JD time value.
    
    Param: ttJD - The TT JD time (float).
    
    Returns: The TAI JD value (float).
    """   
    
    return ttJD - sec_to_jd(TAI_TT_OFFSET)


def utc_to_tt(utcJD):
    """
    Get the TT JD time value for a given UTC JD time value.
    
    Param: utcJD - The UTC JD time (float).
    
    Returns: The TT JD value (float).
    """    
    
    return tai_to_tt(utc_to_tai(utcJD))


def tt_to_utc(ttJD):
    """
    Get the UTC JD time value for a given TT JD time value.
    
    Param: ttJD - The TT JD time (float).
    
    Returns: The UTC JD value (float).
    """     
    
    return tai_to_utc(tt_to_tai(ttJD))


def tt_to_tdb(ttJD):
    """
    Get the TDB JD time value for a given TT JD time value.
    Adopted from "Astronomical Almanac Supplement", Seidelmann 1992, 2.222-1.
    
    Param: ttJD - The TT JD time (float).
    
    Returns: The TDB JD value (float).
    """    
    
    g = math.radians(357.53 + (0.9856003 * (ttJD - J2000_UTC_JD)))
    return (ttJD + (0.001658 * math.sin(g)) + (0.000014 * math.sin(2.0 * g)))  


def get_tai_from_sys():
    """
    Return the current time taken from the system clock as a
    TAI MJD (float).
    """
    
    return utcjd_to_taimjd(get_julian_from_sys())


def hms_to_sec(hms):
    """
    Convert hours, minutes, seconds to seconds.
    
    Param: hms - object of type hms representing time/angle.
    
    Returns: Seconds (float) offset of time/angle.
    """
    
    return hms.hours*3600.0 + hms.minutes*60.0 + hms.seconds


def deg_to_sec(degrees):
    """
    Convert longitude degrees into time seconds.
    
    Param: degrees - longitude (float degrees)
    
    Returns: Seconds (float) offset in time for longitude.
    """
    
    hours = degrees / 15.0
    return hours*3600.0


def get_local_sidereal_time(lng, jD):
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
# Position utility classes
######################################################################

class geo_posn(lnlat_posn):
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
    
    def __init__(self, lng = 0.0, lat = 0.0, elv = 0.0):
        """
        Create a geo_posn object.
        
        Param: lng - longitude (float degrees)
        Param: lat - latitude (float degrees)
        Param: elv - elevation (float meters)
        """
        
        lnlat_posn.__init__(self, lng, lat)
        
        self.elv = elv
        
    def __str__(self):
        """
        Get string representation of geo_posn object.
        """
        
        return lnlat_posn.__str__(self) + " %0.3f" % self.elv
        
    def __reduce__(self):
        """
        geo_posn object pickle reduce method.
        """
        
        return (geo_posn, (self.lng, self.lat, self.elv))
        
    def __repr__(self):
        """
        geo_posn object repr string method.
        """
        
        return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.lng), repr(self.lat), repr(self.elv))
        
    def __getitem__(self, key):
        """
        Get members by subscript index.
        """
        
        if key == 2:
            return self.elv
        else:
            return lnlat_posn.__getitem__(self, key)
            
    def __setitem__(self, key, value):
        """
        Set members by subscript index.
        """
        
        if key == 2:
            self.elv = value
        else:
            lnlat_posn.__setitem__(self, key, value)
            
    def __eq__(self, other):
        """
        geo_posn equality operation.
        """
        
        if not isinstance(other, geo_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return lnlat_posn.__eq__(self, other) and (self.elv == other.elv)
        
    def __ne__(self, other):
        """
        geo_posn non-equality operation.
        """
        
        if not isinstance(other, geo_posn):
            raise TypeError("comparison not supported for type %s" % type(other).__name__)
            
        return lnlat_posn.__ne__(self, other) or (self.elv != other.elv)


######################################################################
# Pointing utility functions
######################################################################

def dir_cos(posn):
    """
    Get direction cosines from azimuth and zenith angles.
    This function calculates the cosine values based on the LWA coordinate system:
    l = unit vector in E direction (azimuth = 90)
    m = unit vector in N direction (azimuth = 0)
    n = unit vector in zenith direction (zenith = 0, altitude = 90)
    
    Param: posn - object of type hrz_posn giving local position
    
    Returns a tuple (l,m,n) of float values for the direction cosines.
    """
    
    azRad = math.radians(posn.az)
    zenRad = math.radians(posn.zen())
    szen = math.sin(zenRad)
    
    l = (szen * math.sin(azRad))
    m = (szen * math.cos(azRad))
    n = math.cos(zenRad)
    
    return (l, m, n)


def get_rect_from_equ(posn):
    """
    Transform equatorial coordinates to rectangular coordinates.
    
    Param: posn - Object of type equ_posn giving position.
    
    Returns: Object of type rect_posn giving rectangular coordinates (normallized to 1).
    """
    
    raRad = math.radians(posn.ra)
    decRad = math.radians(posn.dec)
    cdec = math.cos(decRad) 
    
    x = (math.cos(raRad) * cdec)
    y = (math.sin(raRad) * cdec)
    z = math.sin(decRad)
    
    return rect_posn(x, y, z)


def get_equ_from_rect(posn):
    """
    Transform rectangular coordinates to equatorial coordinates.
    
    Param: posn - Object of type rect_posn giving position.
    
    Returns: Object of type equ_posn giving equatorial coordinates.
    """
    
    x = posn.X
    y = posn.Y
    z = posn.Z
    
    t = math.sqrt((x * x) + (y * y))
    ra = math.degrees(math.atan2(y, x))
    dec = math.degrees(math.atan2(z, t))
    
    return equ_posn(range_degrees(ra), dec)


def get_geo_from_rect(posn):
    """
    Transform ECEF rectangular coordinates to geographical coordinates.
    Adapoted from "Satellite Orbits", Montenbruck and Gill 2005, 5.85 - 5.88.
    Also see gpstk ECEF::asGeodetic() method.
    
    Param: posn - object of type rect_posn giving position.
    
    Returns: object of type geo_posn giving geographical coordinates.
    """
    
    x = posn.X
    y = posn.Y
    z = posn.Z
    
    e2 = 6.69437999014e-3
    a = 6378137.0
    
    p = math.sqrt((x * x) + (y * y))
    dz = (e2 * z)
    
    while(True):
        sz = (z + dz)
        slat =  sz / math.sqrt((x * x) + (y * y) + (sz * sz))
            
        N = a / math.sqrt(1.0 - (e2 * slat * slat))  
        
        olddz = dz
        dz = (N * e2 * slat) 
            
        dzerr = math.fabs(dz - olddz)  
        if dzerr < 1.0e-9:
            break 
            
    sz = (z + dz)  
    lon = math.atan2(y, x)
    lat = math.atan2(sz, p)
    h = math.sqrt((x * x) + (y * y) + (sz * sz)) - N 
    
    lat = math.degrees(lat)
    lon = math.degrees(lon)  
    
    return geo_posn(range_degrees(lon), lat, h)


def get_rect_from_geo(posn):
    """
    Transform geographical coordinates to ECEF rectangular coordinates.
    Adopted from "Satellite Orbits", Montenbruck and Gill 2005, 5.83 - 5.84.
    Also see gpstk Geodetic::asECEF() method.
    
    Param: posn - object of type geo_posn giving geographical coordinates.
    
    Returns: object of type rect_posn giving ECEF position. 
    """
    
    lon = math.radians(posn.lng)
    lat = math.radians(posn.lat)
    h = posn.elv     
    
    a = 6378137.0
    e2 = 6.69437999014e-3
    
    clat = math.cos(lat)
    slat = math.sin(lat)
    
    rad_cur  = a / math.sqrt(1.0 - e2*slat*slat)
    
    x = (rad_cur + h) * clat * math.cos(lon)
    y = (rad_cur + h) * clat * math.sin(lon)
    z = ((1.0 - e2) * rad_cur + h) * slat
    
    return rect_posn(x, y, z)


def get_precession(jD1, pos, jD2):
    """
    Caculate precession of equatorial coordinates from one epoch to
    another.
    
    Param: jD1 - UTC Julian day of epoch 1.
    Param: pos - object of type equ_posn giving epoch 1 position
    Param: jD2 - UTC Julian day of epoch 2.
    
    Returns: object of type equ_posn giving epoch 2 position.
    """
    
    rect = get_rect_from_equ(pos)
    
    # the precession function time paramters should be in TDB
    # the UTC->TDB conversion, however, currently does not support
    # times before 1972
    # this might result in a small error in position
    (xp, yp, zp) = _precession(jD1, (rect.X, rect.Y, rect.Z), jD2) 
    
    rect.X = xp
    rect.Y = yp
    rect.Z = zp
    
    return get_equ_from_rect(rect)


def _precession(tjd1, pos, tjd2):
    """
    Precesses equatorial rectangular coordinates from one epoch to
    another.  The coordinates are referred to the mean equator and
    equinox of the two respective epochs.
    
    Adapoted from NOVAS-C library Version 2.0.1, function precession(). 
    
    CREDITS:
    Astronomical Applications Dept, U.S. Naval Observatory
    
    REFERENCES:
    Explanatory Supplement to AE and AENA (1961); pp. 30-34.
    Lieske, J., et al. (1977). Astron. & Astrophys. 58, 1-16.
    Lieske, J. (1979). Astron. & Astrophys. 73, 282-284.
    Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, pp. 1197-1210.
    Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
    Subroutines"; USNO internal document dated 20 Oct 1988; revised 15 Mar 1990.
    
    INPUT
    ARGUMENTS:
    tjd1 (double)
    TDB Julian date of first epoch.
    pos[3] (double)
    Position vector, geocentric equatorial rectangular coordinates,
    referred to mean equator and equinox of first epoch.
    tjd2 (double)
    TDB Julian date of second epoch.
    
    RETURNED
    VALUE:
    pos2[3] (double)
    Position vector, geocentric equatorial rectangular coordinates,
    referred to mean equator and equinox of second epoch.
    """
    
    T0 = 2451545.00000000
    RAD2SEC = 206264.806247096355
    
    #
    #  't' and 't1' below correspond to Lieske's "big T" and "little t".
    #
    t = (tjd1 - T0) / 36525.0
    t1 = (tjd2 - tjd1) / 36525.0
    t02 = t * t
    t2 = t1 * t1
    t3 = t2 * t1
    
    #
    #  'zeta0', 'zee', 'theta' below correspond to Lieske's "zeta-sub-a",
    #  "z-sub-a", and "theta-sub-a".
    #
    zeta0 = (2306.2181 + 1.39656 * t - 0.000139 * t02) * t1 \
        + (0.30188 - 0.000344 * t) * t2 + 0.017998 * t3
        
    zee = (2306.2181 + 1.39656 * t - 0.000139 * t02) * t1 \
        + (1.09468 + 0.000066 * t) * t2 + 0.018203 * t3
        
    theta = (2004.3109 - 0.85330 * t - 0.000217 * t02) * t1 \
        + (-0.42665 - 0.000217 * t) * t2 - 0.041833 * t3
        
    zeta0 /= RAD2SEC
    zee /= RAD2SEC
    theta /= RAD2SEC
    
    #
    #  Precalculate trig terms.
    #
    cz0 = math.cos (zeta0)
    sz0 = math.sin (zeta0)
    ct = math.cos (theta)
    st = math.sin (theta)
    cz = math.cos (zee)
    sz = math.sin (zee)
    
    #
    #  Precession rotation matrix follows.
    #
    xx =  cz0 * ct * cz - sz0 * sz
    yx = -sz0 * ct * cz - cz0 * sz
    zx = -st * cz
    xy = cz0 * ct * sz + sz0 * cz
    yy = -sz0 * ct * sz + cz0 * cz
    zy = -st * sz
    xz = cz0 * st
    yz = -sz0 * st
    zz = ct
    
    #
    #  Perform rotation.
    #
    xr = xx * pos[0] + yx * pos[1] + zx * pos[2]
    yr = xy * pos[0] + yy * pos[1] + zy * pos[2]
    zr = xz * pos[0] + yz * pos[1] + zz * pos[2]
    
    return (xr, yr, zr)


def B1950_to_J2000(pos):
    """
    Convert B1950 epoch to J2000 epoch for equatorial coordinates.
    
    Param: pos - object of type equ_posn giving B1950 coordinates
    
    Returns: object of type equ_posn giving J2000 coordinates.
    
    .. note::
        The accuracy of this function is about 0.01 degrees.
    """
    
    _posn = equ_posn()
    coord = ephem.Equatorial(deg_to_rad(pos.ra), deg_to_rad(pos.dec), epoch=ephem.B1950)
    coord = ephem.Equatorial(coord, epoch=ephem.J2000)
    ra = rad_to_deg(coord.ra)
    dec = rad_to_deg(coord.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn


def J2000_to_B1950(pos):
    """
    Convert J2000 epoch to B1950 epoch for equatorial coordinates.
    
    Param: pos - object of type equ_posn giving J2000 coordinates
    
    Returns: object of type equ_posn giving B1950 coordinates.
    
    .. note::
        The accuracy of this function is about 0.01 degrees.
    """   
    
    _posn = equ_posn()
    coord = ephem.Equatorial(deg_to_rad(pos.ra), deg_to_rad(pos.dec), epoch=ephem.J2000)
    coord = ephem.Equatorial(coord, epoch=ephem.B1950)
    ra = rad_to_deg(coord.ra)
    dec = rad_to_deg(coord.dec)
    
    _posn.ra = ra
    _posn.dec = dec
    return _posn
