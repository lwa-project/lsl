"""
Module that contains commmon routines for the SDF and IDF modules.
"""

import os
import re
import math
import warnings
from datetime import datetime, timedelta, timezone

try:
    import zoneinfo
except ImportError:
    from backports import zoneinfo

from astropy.time import Time as AstroTime

from lsl.astro import date as astroDate, get_date as astroGetDate
from lsl.common.stations import lwa1
from lsl.common.color import colorfy

__version__ = '0.1'
__all__ = ['pid_print', 'render_file_size', 'render_bandwidth', 'parse_time']


def pid_print(*args, **kwds):
    print(f"[{os.getpid()}]", *args, **kwds)


def render_file_size(size):
    """Convert a file size in bytes to a easy-to-use string."""
    
    units = 'B'
    if size >= 1024**4:
        size /= 1024.0**4
        units = 'TB'
    elif size >= 1024**3:
        size /= 1024.0**3
        units = 'GB'
    elif size >= 1024**2:
        size /= 1024.0**2
        units = 'MB'
    elif size >= 1024**1:
        size /= 1024.0**1
        units = 'kB'
        
    return f"{size:.2f} {units}"


def render_bandwidth(filter, filter_codes):
    """Convert a filter number to an easy-to-use string."""
    
    bw = filter_codes[filter]
    units = 'Hz'
    if bw > 1e6:
        bw /= 1e6
        units = 'MHz'
    elif bw > 1e3:
        bw /= 1e3
        units = 'kHz'
        
    return f"{bw:.3f} {units}"


def _get_equinox_equation(jd):
    """
    Compute the equation of the equinoxes (nutation in right ascension) in 
    hours for the specified Julian Date.
    
    From:
    http://aa.usno.navy.mil/faq/docs/GAST.php
    """
    
    # Get the number of days since January 1, 2000 @ 12:00 UT
    D = jd - 2451545.0
    
    # Compute the obliquity
    epsilon = 23.4393 - 0.0000004*D
    
    # Compute the mean longitude of the Sun
    L = 280.47 + 0.98565*D
    
    # Compute the longitude of the Moon's ascending node
    Omega = 125.04 - 0.052954*D
    
    # The nutation in the longitude (hours)
    deltaPsi = -0.000319*math.sin(Omega*math.pi/180.0) - 0.000024*math.sin(2*L*math.pi/180.0)
    
    # Return the equation of the equinoxes
    return deltaPsi * math.cos(epsilon*math.pi/180.0)


_dtRE = re.compile(r'^((?P<tz>[A-Z]{2,3}) )?(?P<year>\d{4})[ -/]((?P<month>\d{1,2})|(?P<mname>[A-Za-z]{3}))[ -/](?P<day>\d{1,2})[ T](?P<hour>\d{1,2}):(?P<minute>\d{1,2}):(?P<second>\d{1,2}(\.\d{1,6})?) ?(?P<tzOffset>[-+]\d{1,2}:?\d{1,2})?$')
_EST = zoneinfo.ZoneInfo('US/Eastern')
_CST = zoneinfo.ZoneInfo('US/Central')
_MST = zoneinfo.ZoneInfo('US/Mountain')
_PST = zoneinfo.ZoneInfo('US/Pacific')


def parse_time(s, station=lwa1):
    """
    Given a time zone-aware datetime instance or a string in the format of 
    (UTC) YYYY MM DD HH:MM:SS.SSS, return the corresponding UTC datetime object.
    This function goes a little beyond what datetime.strptime does in the 
    since that it handle both integer and float seconds as well as does the 
    appropriate rounding to get millisecond precision.
    
    .. versionchanged:: 2.0.0
        Added support for astropy.time.Time instances
    
    .. versionchanged:: 1.2.0
        Renamed the 'site' keyword to 'station'
    
    .. versionchanged:: 1.0.0
        Renamed to parse_time()
        Added support for timezone-aware datetime instances
    """
    
    if isinstance(s, AstroTime):
        s = s.utc.datetime
        
    if isinstance(s, datetime):
        if s.tzinfo is None:
            raise ValueError("Only time zone-aware datetime instances are supported.")
            
        # Round the microsecond value to milliseconds
        us = s.microsecond
        us = int(round(us/1000.0))*1000
        dtObject = s.replace(microsecond=us)
        
    else:
        mtch = _dtRE.match(s)
        if mtch is None:
            raise ValueError(f"Unparsable time string: '{s}'")
        else:
            year = int(mtch.group('year'))
            day = int(mtch.group('day'))
            
            hour = int(mtch.group('hour'))
            minute = int(mtch.group('minute'))
            second = math.floor(float(mtch.group('second')))
            microsecond = int(round((float(mtch.group('second')) - second)*1000.0))*1000
            second = int(second)
            
            if mtch.group('mname') is None:
                month = int(mtch.group('month'))
            else:
                monthName = mtch.group('mname').lower()
                if monthName == 'jan':
                    month = 1
                elif monthName == 'feb':
                    month = 2
                elif monthName == 'mar':
                    month = 3
                elif monthName == 'apr':
                    month = 4
                elif monthName == 'may':
                    month = 5
                elif monthName == 'jun':
                    month = 6
                elif monthName == 'jul':
                    month = 7
                elif monthName == 'aug':
                    month = 8
                elif monthName == 'sep':
                    month = 9
                elif monthName == 'oct':
                    month = 10
                elif monthName == 'nov':
                    month = 11
                elif monthName == 'dec':
                    month = 12
                else:
                    raise ValueError(f"Unknown month abbreviation: '{monthName}'")
                    
            if mtch.group('tz') is None and mtch.group('tzOffset') is None:
                tz = timezone.utc
            elif mtch.group('tzOffset') is not None:
                tzOffsetSign = 1
                if mtch.group('tzOffset')[0] == '-':
                    tzOffsetSign = -1
                tzOffsetHours = int( mtch.group('tzOffset').replace(':', '')[1:3] )
                tzOffsetMinutes = int( mtch.group('tzOffset').replace(':', '')[3:5] )
                
                tzOffsetSeconds = tzOffsetHours*3600 + tzOffsetMinutes*60
                tzOffsetSeconds *= tzOffsetSign
                tz = timezone(timedelta(seconds=tzOffsetSeconds))
            else:
                tzName = mtch.group('tz')
                if tzName in ['UT', 'UTC']:
                    tz = _UTC
                elif tzName in ['EST', 'EDT']:
                    tz = _EST
                elif tzName in ['CST', 'CDT']:
                    tz = _CST
                elif tzName in ['MST', 'MDT']:
                    tz = _MST
                elif tzName in ['PST', 'PDT']:
                    tz = _PST
                elif tzName in ['LST',]:
                    tz = 'LST'
                else:
                    ## Exhaustive search through pytz.  This may yield strange matches...
                    warnings.warn(colorfy("{{%%yellow Entering zoneinfo search mode for '%s'}}" % tzName), RuntimeWarning)
                    
                    tzFound = False
                    tzNormal = datetime(year, month, day)
                    for tzi in list(zoneinfo.available_timezones()):
                        tz = zoneinfo.ZoneInfo(tzi)
                        try:
                            cTZName = tzNormal.replace(tzinfo=tz).tzname()
                        except TypeError:
                            continue
                        if cTZName == tzName:
                            tzFound = True
                            break
                            
                    if not tzFound:
                        raise ValueError(f"Unknown time zone: '{tzName}'")
                        
            if tz == 'LST':
                # Deal with sidereal times...
                #
                # NOTE:
                # The RMS on this method is ~0.4 seconds over the years 
                # 2000 to 2100.  This should be "good enough" for scheduling
                # purposes.
                
                # Get the position of the observer on the Earth and the Julian 
                # Date of midnight UT for the day we want to map LST to
                dt = astroDate(year, month, day, 0, 0, 0)
                jd = dt.to_jd()
                
                # Get the LST in hours
                LST = hour + minute/60.0 + (second + microsecond/1e6)/3600.0
                
                # Get the Greenwich apparent ST for LST using the longitude of 
                # the site.  The site longitude is stored as radians, so convert
                # to hours first.
                GAST = LST - station.long*12/math.pi
                
                # Get the Greenwich mean ST by removing the equation of the 
                # equinoxes (or some approximation thereof)
                GMST = GAST - _get_equinox_equation(jd)
                
                # Get the value of D0, days since January 1, 2000 @ 12:00 UT, 
                # and T, the number of centuries since the year 2000.  The value
                # of T isn't terribly important but it is nice to include
                D0 = jd - 2451545.0
                T = D0 / 36525.0
                
                # Solve for the UT hour for this LST and map onto 0 -> 24 hours
                # From: http://aa.usno.navy.mil/faq/docs/GAST.php
                H  = GMST - 6.697374558 - 0.06570982441908*D0 - 0.000026*T**2
                H /= 1.002737909350795
                while H < 0:
                    H += 24/1.002737909350795
                while H > 24:
                    H -= 24/1.002737909350795
                    
                # Get the full Julian Day that this corresponds to
                jd += H/24.0
                
                # Convert the JD back to a time and extract the relevant 
                # quantities needed to build a datetime instance
                dt = astroGetDate(jd)
                
                tz = _UTC
                year = dt.years
                month = dt.months
                day = dt.days
                hour = dt.hours
                minute = dt.minutes
                second = int(dt.seconds)
                microsecond = int(round((dt.seconds - second)*1e6))
                ## Trim the microsecond down to the millisecond level
                microsecond = int(round(microsecond/1000.0))*1000
                
            # Localize as the appropriate time zone
            dtObject = datetime(year, month, day, hour, minute, second, microsecond), tzinfo=tz)
            
    # Return as UTC
    return dtObject.astimezone(timezone.utc)
