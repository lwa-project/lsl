# -*- coding: utf-8 -*-

"""
Astronomical utility functions and classes based on libnova library.
"""

######################################################################
# $Id: astro.py 84 2010-05-19 22:24:00Z dwood $
######################################################################


import time
import math

from lsl import libnova


__revision__  = "$ Revision: 86 $"
__version__   = "dev"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


######################################################################
# version helper functions
######################################################################


def get_libnova_version():
  """
  Get version of libnova C library in use.
  
  Returns: A tuple of version numbers for libnova C library.
  """
    
  vs = libnova.ln_get_version()
  vs = vs.split('.')
  return (int(vs[0]), int(vs[1]), int(vs[2]))
    

def _hexversion(version):
  """
  Convert tuple of version numbers to flat integer for comparisons.
  """
  
  return (version[0] << 16) + (version[1] << 8) + version[2]
    

_MIN_LIBNOVA_VER = _hexversion((0, 12, 0))
_MAX_LIBNOVA_VER = _hexversion((0, 13, 0))

_CUR_LIBNOVA_VER = _hexversion(get_libnova_version())

if (_CUR_LIBNOVA_VER < _MIN_LIBNOVA_VER) or (_CUR_LIBNOVA_VER > _MAX_LIBNOVA_VER):
  raise ImportError("libnova version %s not supported" % libnova.ln_get_version())
    
    

######################################################################
# python class wrappers for libnova C structures
######################################################################


class dms(libnova.ln_dms):
  """
  Wrapper for libnova ln_dms structure.
  Represents angles in degrees, minutes, seconds.
  
  Public members:
    neg     - True if measured west of GM, south of EQ
              False if measured east of GM, north of EQ.
    degrees - Angle degrees (integer).
    minutes - Angle minutes (integer).
    seconds - Angle seconds (float).
  """
  
	
  def __init__(self, neg = None, degrees = None, minutes = None, seconds = None):
    """
    Create a dms object.
    
    Param: neg      - True if measured west of GM, south of EQ
                      False if measured east of GM, north of EQ.
    Param: degrees  - Angle degrees (integer [0, 359)).
    Param: minutes  - Angle minutes (integer [0, 59]).
    Param: seconds  - Angle seconds (float [0.0, 60.0)).
    """
          
    libnova.ln_dms.__init__(self)

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
        
    return "%s.%s(%s,%s,%s,%s)" % (type(self).__module__,type(self).__name__,
      repr(self.neg), repr(self.degrees), repr(self.minutes), repr(self.seconds))


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
          
    return cmp(self.to_deg(), o)


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
		


class hms(libnova.ln_hms):
  """
  Wrapper for libnova ln_hms structure.
  Represents times/angles in hours, minutes, seconds.
  
  Public members:
    hours - Angle/time hours (integer).
    minutes - Angle/time minutes (integer).
    seconds - Angle/time seconds (float).
  """
  
	
  def __init__(self, hours = None, minutes = None, seconds = None):
    """
    Create a hms object.
    
    Param: hours    - Angle/time hours (integer [0, 23]).
    Param: minutes  - Angle/time minutes (integer [0, 59]).
    Param: seconds  - Angle/time seconds (float [0.0, 60.0)).
    """
        
    libnova.ln_hms.__init__(self)

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
        
    return "%s.%s(%s,%s,%s)" % (type(self).__module__,type(self).__name__,
      repr(self.hours), repr(self.minutes), repr(self.seconds))           
		
		
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
        
    return cmp(self.to_deg(), o)


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
        
            

class date(libnova.ln_date):
  """
  Wrapper for libnova ln_date structure.
  Represents UT time in calendar units.
  
  Public members:
    years   - Date years (integer).
    months  - Date months (integer).
    days    - Date days (integer).
    hours   - Date hours (integer).
    minutes - Date minutes (integer).
    seconds - Date seconds (float).
  """


  def __init__(self, years = None, months = None, days = None, hours = None,
                     minutes = None, seconds = None):
    """
    Create a date object.
    
    Param: years    - Date years (integer).
    Param: months   - Date months (integer [1, 12]).
    Param: days     - Date days (integer [1, 31]).
    Param: hours    - Date hours (integer [0, 23]).
    Param: minutes  - Date minutes (integer [0, 59]).
    Param: seconds  - Date seconds (float [0.0, 60.0)).
    """

    libnova.ln_date.__init__(self)

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

    return "%04d-%02d-%02d %02d:%02d:%s" % (self.years, self.months, self.days, 
      self.hours, self.minutes, sec.zfill(6))
    
    
  def __repr__(self):
    """
    date object repr string method
    """
        
    return "%s.%s(%s,%s,%s,%s,%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.years), repr(self.months), repr(self.days), repr(self.hours),
      repr(self.minutes), repr(self.seconds))
    
    
  def __reduce__(self):
    """
    date object pickle reduce method.
    """
        
    return (date, (self.years, self.months, self.days, self.hours, 
      self.minutes, self.seconds))
            
            
  def __cmp__(self, other):
    """
    date comparison tests.
    """
        
    if not isinstance(other, (date, zonedate)):
      raise TypeError("comparison not supported for type %s" % type(other).__name__)
        
    return cmp(self.to_jd(), other.to_jd())   
        
    
  def to_zone(self, gmtoff = None):
    """
    Convert UTC calendar time to local calendar time.
    
    Param: gmtoff - local seconds offset from UTC (integer -43200 to 43200).
                    if set to None, the time module offset value is used for current location 
    
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



class zonedate(libnova.ln_zonedate):
  """
  Wrapper for libnova ln_zonedate structure.
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
	
  
  def __init__(self, years = None, months = None, days = None, hours = None,
                     minutes = None, seconds = None, gmtoff = None):
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

    libnova.ln_zonedate.__init__(self)
            
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
    
    return "%04d-%02d-%02d %02d:%02d:%s [%d]" % (self.years, self.months, self.days, 
      self.hours, self.minutes, sec.zfill(6), self.gmtoff)
        
     
  def __repr__(self):
    """
    zonedate object repr string method
    """
    
    return "%s.%s(%s,%s,%s,%s,%s,%s,%s)" % (type(self).__module__,type(self).__name__,
      repr(self.years), repr(self.months), repr(self.days), repr(self.hours), 
      repr(self.minutes), repr(self.seconds), repr(self.gmtoff))
     
     
  def __reduce__(self):
    """
    zonedate object pickle reduce method.
    """
        
    return (zonedate, (self.years, self.months, self.days, self.hours,
      self.minutes, self.seconds, self.gmtoff))
            
    
  def __cmp__(self, other):
    """
    zonedate comparison tests.
    """
        
    if not isinstance(other, (zonedate, date)):
      raise TypeError("comparison not supported for type %s" % type(other).__name__)
    
    return cmp(self.to_jd(), other.to_jd())
    
       
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
            


class rst_time(libnova.ln_rst_time):
  """
  Wrapper for libnova ln_rst_time structure.
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

    libnova.ln_rst_time.__init__(self)

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
        
    return "%s.%s(%s,%s,%s)" % (type(self).__module__,type(self).__name__,
      repr(self.rise), repr(self.set), repr(self.transit))
    
    
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



        
class hrz_posn(libnova.ln_hrz_posn):
  """
  Wrapper for libnova ln_hrz_posn structure.
  Represents horizontal local position coordinates.  The original libnova convention has been
  modified for the LWA wrapper.  libnova measures azimuth angle clockwise from south to
  west, with due south equal to 0 degrees.  LWA measures azimuth angle clockwise from north
  to east, with due north equal to 0 degrees.  Also, method for using zenith angle instead
  of altitude angle have been added.
    
  Public members:
    az  - Position azimuth angle (float degrees).
    alt - Position altitude angle (float degrees)
        
  Members may also be accessed by subscript:
    hrz_posn[0] = az
    hrz_posn[1] = alt    
  """

  def __init__(self, az = None, alt = None):
    """
    Create a hrz_posn object.
    
    Param: az   - Position azimuth angle (float degrees [0.0, 360.0), 0 = N, 90 = E).
    Param: alt  - Position altitude angle (float degress [-90.0, 90.0]).
    """

    libnova.ln_hrz_posn.__init__(self)

    if az is not None:
      if az < 0.0 or az >= 360.0:
        raise ValueError("az paramerer range is [0.0, 360.0), is set to %0.3f" % az)
      self.az = az
            
    if alt is not None:
      if alt < -90.0 or alt > 90.0:
        raise ValueError("alt paramerer range is [-90.0, 90.0], is set to %0.3f" % alt)
      self.alt = alt


  def __setattr__(self, name, value):
    """
    Returns position azimuth angle (float degrees [0.0, 360.0), 0 = N, 90 = E).
    """
    
    if name == 'az':
      value = range_degrees(value + 180.0)
    libnova.ln_hrz_posn.__setattr__(self, name, value)

        
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
        
    return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.az), repr(self.alt))


  def __reduce__(self):
    """
    hrz_posn object pickle reduce method.
    """

    return (hrz_posn, (self.az, self.alt))


  def __getattribute__(self, name):
    """
    Set position azimuth angle (float degrees [0.0, 360.0), 0 = N, 90 = E).
    """

    value = libnova.ln_hrz_posn.__getattribute__(self, name)
    if name == 'az':
      value = range_degrees(value - 180.0)
    return value
        
        
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
               


class equ_posn(libnova.ln_equ_posn):
  """
  Wrapper for libnova ln_equ_posn structure.
  Represents equatoral/celestial position coordinates.
    
  Public members:
    ra  - Position right ascension angle (float degrees).
    dec - Position declination angle (float degrees).
        
  Members may also be accessed by subscript:
    equ_posn[0] = ra
    equ_posn[1] = dec
  """

  def __init__(self, ra = None, dec = None):
    """
    Create a equ_posn object.
    
    Param: ra   - Position right ascension angle
                  Object of type hms or float degrees [0.0, 360.0).
    Param: dec  - Position declination angle
                  Object of type dms or float degrees [-90.0, 90.0].
    """

    libnova.ln_equ_posn.__init__(self)

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
        
    return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.ra), repr(self.dec))   
            
            
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



class gal_posn(libnova.ln_gal_posn):
  """
  Wrapper for libnova ln_gal_posn structure.
  Represents galactic position coordinates.
    
  Public members:
    l - Position longitude angle (float degrees).
    b - Position latitude angle (float degrees).
        
  Members may also be accessed by subscript:
    gal_posn[0] = l
    gal_posn[1] = b
  """

  def __init__(self, l = None, b = None):
    """
    Create a gal_posn object.
    
    Param: l - Position longitude angle. 
               Object of type dms or float degrees [0.0, 360.0).
    Param: b - Position latitude angle. 
               Object of type dms or float degrees [-90.0, 90.0].
    """

    libnova.ln_gal_posn.__init__(self)

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
        
    return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.l), repr(self.b))
        
    
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

	

class rect_posn(libnova.ln_rect_posn):
  """
  Wrapper for libnova ln_rect_posn structure.
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

  def __init__(self, X = None, Y = None, Z = None): 
    """
    Create a rect_posn object
    Param: X - Position X coordinate (float).
    Param: Y - Position Y coordinate (float).
    Param: Z - Position Z coordinate (float).
    """

    libnova.ln_rect_posn.__init__(self)
            
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
        
    return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.X), repr(self.Y), repr(self.Z))

    
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
        
        

class lnlat_posn(libnova.ln_lnlat_posn):
  """
  Wrapper for libnova ln_lnlat_posn structure.
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
    
  def __init__(self, lng = None, lat = None):
    """
    Create a lnlat_posn object.
    
    Param: lng - Position longitude coordinate
                 Object of type dms or float degress (-360.0, 360.0).
    Param: lat - Position latitude coordinate
                 Object of type dms or float degrees [-90.0, 90.0].
    """
        
    libnova.ln_lnlat_posn.__init__(self)
            
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
        
    return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.lng), repr(self.lat))
        
        
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
    


class ecl_posn(lnlat_posn):
  """
  Represents position as ecliptic longitude and latitude.
    
  Public members:
    lng - Position longitude coordinate (float degrees).
    lat - Position latitude coordinate (float degrees).
        
  Members may also be accessed by subscript:
    ecl_posn[0] = lng
    ecl_posn[1] = lat
  """
    
    
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
        
        
        
class nutation(libnova.ln_nutation):
  """
  Wrapper for libnova ln_nutation structure.
  Provides nutation information in longitude and obliquity.
  
  Public members:
    longitude - Nutation in longitude (float degrees).
    obliquity - Nutation in ecliptic obliquity (float degrees).
    ecliptic - Obliquity of the ecliptic (float degrees).
  """
    
    
  def __init__(self, longitude = None, obliquity = None, ecliptic = None):
    """
    Create a nutation object.
    
    Param: longitude  - Nutation in longitude.
                        Object of type dms or float degress (-360.0, 360.0).
    Param: obliquity  - Nutation in obliquity.
                        Object of type dms or float degrees [-90.0, 90.0].
    Param: ecliptic   - Obliquity of the ecliptic.
                        Object of type dms or float degrees [-90.0, 90.0].
    """
        
    libnova.ln_nutation.__init__(self)
        
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
        
    return "%s.%s(%s,%s,%s)" % (type(self).__module__, type(self).__name__,
      repr(self.longitude), repr(self.obliquity), repr(ecliptic))
    
     
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
        
    return (deg_to_dms(self.longitude), deg_to_dms(self.obliquity),
      deg_to_dms(self.ecliptic))  
        
        
        

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
# python fucntion wrappers for libnova C General Conversion Functions
######################################################################



def date_to_zonedate(date, gmtoff):
  """
  Wrapper for for libnova ln_date_to_zonedate() function.
  Convert UTC calendar time to local calendar time.
  
  Param: date   - A date object representing UTC time.
  Param: gmtoff - Seconds offset from UTC (integer -43200 to 43200).
  
  Returns object of type zonedate representing local time.
  """
	
  _zdate = zonedate()
  libnova.ln_date_to_zonedate(date, _zdate, gmtoff)
  return _zdate
	
	
	
def zonedate_to_date(zonedate):
  """
  Wrapper for for libnova ln_zonedate_to_date() function.
  Convert local calendar time to UTC calendar time.
  
  Param: zonedate - Object of type zonedate representing local time.  
  
  Returns object of type date representing UTC time.
  """
	
  _date = date()
  libnova.ln_zonedate_to_date(zonedate, _date)
  return _date	


def rad_to_deg(radians):
  """
  Wrapper for libnova ln_rad_to_deg() function.
  Convert radians to degress.
  
  Param: radians - Angle in radians (float).
  
  Returns angle in degress (float).
  """
	
  return libnova.ln_rad_to_deg(radians)
	
	
def deg_to_rad(degrees):
  """
  Wrapper for libnova ln_deg_to_rad() function.
  Convert degress to radians.
  
  Param: degrees - Angle in degrees (float).
  
  Returns angle in radians (float).
  """
	
  return libnova.ln_deg_to_rad(degrees)
	

def dms_to_rad(dms):
  """
  Wrapper for libnova ln_dms_to_rad() function.
  Convert angles degrees, minutes, seconds to radians.
  
  Param: dms - Object of type dms representing angle.
  
  Returns angle in radians (float).
  """
	
  return libnova.ln_dms_to_rad(dms)
	

def dms_to_deg(dms):
  """
  Wrapper for libnova ln_dms_to_deg() function.
  Convert angles degrees, minutes, seconds to float degrees.
  
  Param: dms - Object of type dms representing angle.
  
  Returns angle in degrees (float).
  """
	
  return libnova.ln_dms_to_deg(dms)


def deg_to_dms(degrees):
  """
  Wrapper for libnova ln_deg_to_dms() function.
  Convert angles float degrees to degrees, minutes, seconds.
  
  Param: degrees - Angle in degrees (float). 
  
  Returns object of type dms representing angle.
  """
	
  _dms = dms()
  libnova.ln_deg_to_dms(degrees, _dms)
  return _dms


def rad_to_dms(radians):
  """
  Wrapper for libnova ln_rad_to_dms() function.
  Convert angles float radians to degrees, minutes, seconds.
  
  Param: radians - Angle in radians (float). 
  
  Returns object of type dms representing angle.
  """
	
  _dms = dms()
  libnova.ln_rad_to_dms(radians, _dms)
  return _dms
	

def hms_to_deg(hms):
  """
  Wrapper for libnova ln_hms_to_deg() function.
  Convert angles hours, minutes, seconds to float degrees.
  
  Param: hms - Object of type hms representing angle.
  
  Returns angle in degrees (float).
  """
	
  return libnova.ln_hms_to_deg(hms)
	
	
def hms_to_rad(hms):
  """
  Wrapper for libnova ln_hms_to_rad() function.
  Convert angles hours, minutes, seconds to float radians.
  
  Param: hms - Object of type hms representing angle.
  
  Returns angle in radians (float).
  """
	
  return libnova.ln_hms_to_rad(hms)	
	

def deg_to_hms(degrees):
  """
  Wrapper for libnova ln_deg_to_hms() function.
  Convert angles float degrees to hours, minutes, seconds.
  
  Param: degrees - Angle in degrees (float). 
  
  Returns object of type hms representing angle.
  """
	
  _hms = hms()
  libnova.ln_deg_to_hms(degrees, _hms)
  return _hms


def rad_to_hms(radians):
  """
  Wrapper for libnova ln_rad_to_hms() function.
  Convert angles float radians to hours, minutes, seconds.
  
  Param: radians - Angle in radians (float). 
  
  Returns object of type hms representing angle.
  """
	
  _hms = hms()
  libnova.ln_rad_to_hms(radians, _hms)
  return _hms


def add_secs_hms(hms, seconds):
  """
  Wrapper for libnova ln_add_secs_hms() function.
  Add seconds to time/angle hours, minutes, seconds.
  
  Param: hms      - Object of type hms representing angle.
  Param: seconds  - Seconds offset (float) to add to angle. 
  
  Returns object of type hms representing angle + offset.
  """
	
  libnova.ln_add_secs_hms(hms, seconds)
  return hms
	
	
def add_hms(source, dest):
  """
  Wrapper for libnova ln_add_hms() function.
  Adds time/angle hours, minutes, seconds.
  
  Param: source - Object of type hms represeting angle 1.
  Param: dest   - Object of type hms represeting angle 2.
  
  Returns object of type hms representing sum of angles.
  """
	
  libnova.ln_add_hms(source, dest)
  return dest


def hrz_to_nswe(pos):
  """
  Wrapper for libnova ln_hrz_to_nswe() function.
  Get cardinal/ordinal azimuth direction.
  
  Param: pos - Object of type hrz_posn giving local position.
  
  Returns string giving direction.
  """
    
  return libnova.ln_hrz_to_nswe(pos)
    	

def range_degrees(degrees):
  """
  Put angle into range [0, 360].
  
  Param: degrees - large angle (float degrees)
  
  Returns: angle in range (float degrees)
  """

  return libnova.ln_range_degrees(degrees)		



######################################################################
# python fucntion wrappers for libnova C General Calendar Functions
######################################################################



def get_julian_day(date):
  """
  Wrapper for libnova ln_get_julian_day() function.
  Convert calendar time to Julian day.
  
  Param: date - Object of type date representing UTC time.
  
  Returns UTC time in Julian days (float).
  """
	
  return libnova.ln_get_julian_day(date)
    
    
def get_julian_local_date(zonedate):
  """
  Wrapper for libnova ln_get_julian_local_date() function.
  Convert local calendar time to Julian day.
  
  Param: zonedate - Object of type zonedate representing local time.
  
  Returns UTC time in Julian days (float).
  """
    
  return libnova.ln_get_julian_local_date(zonedate)   


def get_date(jD):
  """
  Wrapper for libnova ln_get_date() function.
  Convert Julian day to calendar time.
  
  Param: jD - UTC time in Julian days (float).
  
  Returns object of type date representing UTC time.
  """
	
  _date = date()
  libnova.ln_get_date(jD, _date)
  return _date
	
	
	
def get_day_of_week(date):
  """
  Wrapper for libnova ln_get_day_of_week() function.
  Gets day of week from calendar time.
  
  Param: date - Object of type date representing UTC time.
  
  Returns day of week (0 = Sunday, 6 = Saturday).
  """
	
  return libnova.ln_get_day_of_week(date)
	
		
	
def get_julian_from_sys():
  """
  Wrapper for libnova ln_get_julian_from_sys() function.
  
  Returns UTC Julian day (float) from system clock.
  """
	
  return libnova.ln_get_julian_from_sys()
	
	

def get_date_from_sys():
  """
  Wrapper for libnova ln_get_date_from_sys() function.
  Gets calendar time from system clock.
  
  Returns object of type date representing UTC time.
  """
	
  _date = date()
  libnova.ln_get_date_from_sys(_date)
  return _date



def get_julian_from_timet(time_):
  """
  Wrapper for libnova ln_get_julian_from_timet() function.
  Gets Julian day from Unix time.
  
  Param: time_ - Unix timet in seconds (integer)
  
  Returns UTC Julian day (float).
  """
	
  tt = libnova.new_time_t()
  libnova.time_t_assign(tt, time_)
  jD = libnova.ln_get_julian_from_timet(tt)
  libnova.delete_time_t(tt)
  return jD
	
	
	
def get_timet_from_julian(jD):
  """
  Wrapper for libnova ln_get_timet_from_julian() function.
  Gets Unix time from Julian day.
  
  Param: jD - UTC Julian day (float).
  
  Returns Unix timet in seconds (integer).
  """
	
  tt = libnova.new_time_t()
  libnova.ln_get_timet_from_julian(jD, tt)
  time_ = libnova.time_t_value(tt)
  libnova.delete_time_t(tt)
  return time_



######################################################################
# python fucntion wrappers for libnova C Transformation of Coordinates 
# Functions
######################################################################


def get_hrz_from_equ(object_, observer, jD):
  """
  Wrapper for libnova ln_get_hrz_from_equ() function.
  Get local horizontal coordinates from equatorial/celestial coordinates.
  
  Param: object_  - Object of type equ_posn representing celestial position.
  Param: observer - Object of type lnlat_posn representing observer position.
  Param: jD       - UTC Julian day (float).
  
  Returns object of type hrz_posn representing local position.
  """

  _posn = hrz_posn()
  libnova.ln_get_hrz_from_equ(object_, observer, jD, _posn)
  return _posn
    
    
def get_equ_from_hrz(object_, observer, jD):
  """
  Wrapper for libnova ln_get_equ_from_hrz() function.
  Get equatorial/celestial coordinates from local horizontal coordinates.
  
  Param: object_  - Object of type hrz_posn representing horizontal position.
  Param: observer - Object of type lnlat_posn representing observer position.
  Param: jD       - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """
    
  _posn = equ_posn()
  libnova.ln_get_equ_from_hrz(object_, observer, jD, _posn)
  return _posn     


def get_ecl_from_rect(rect):
  """
  Wrapper for libnova ln_get_ecl_from_rect() function.
  Get ecliptical coordinates from rectangular coordinates.
  
  Param: rect - Object of type rect_posn representing position.
  
  Returns object of type lnlat_posn representing ecliptical position.
  """
    
  _posn = lnlat_posn()
  libnova.ln_get_ecl_from_rect(rect, _posn)
  return _posn
    
    
def get_equ_from_ecl(object_, jD):
  """
  Wrapper for libnova ln_get_equ_from_ecl() function.
  Get equatorial coordinates from ecliptical coordinates for a given time.
  
  Param: object_  - Object of type lnlat_posn representing ecliptic position.
  Param: jD       - UTC Julian day (float). 
  
  Returns object of type equ_posn representing equatorial position.
  """
    
  _posn = equ_posn()
  libnova.ln_get_equ_from_ecl(object_, jD, _posn)
  return _posn
    
    
def get_ecl_from_equ(object_, jD):
  """
  Wrapper for libnova ln_ecl_from_equ() function.
  Get ecliptical coordinates from equatorial coordinates for a given time.
  
  Param: object_  - Object of type equ_posn representing equatorial position.
  Param: jD       - UTC Julian day (float). 
  
  Returns object of type ecl_posn representing ecliptic position.
  """
    
  _posn = ecl_posn()
  libnova.ln_get_ecl_from_equ(object_, jD, _posn)
  return _posn    
         

def get_equ_from_gal(object_):
  """
  Wrapper for libnova ln_get_equ_from_gal() function.
  Get B1950 equatorial coordinates from galactic coordinates.
  
  Param: object_ - Object of type gal_posn representing galactic position.
  
  Returns object of type equ_posn representing B1950 equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_equ_from_gal(object_, _posn)
  return _posn


def get_gal_from_equ(object_):
  """
  Wrapper for libnova ln_gal_from_equ() function.
  Get galactic coordinates from B1950 equatorial coordinates.
  
  Param: object_ - Object of type equ_posn representing B1950 equatorial 
                   position.
  
  Returns object of type gal_posn representing galactic position.
  """

  _posn = gal_posn()
  libnova.ln_get_gal_from_equ(object_, _posn)
  return _posn
        
        
_LN_GET_EQU2000_FROM_GAL_MIN = _hexversion((0, 12, 0))        
    
def get_equ2000_from_gal(object_):
  """
  Wrapper for libnova ln_get_equ2000_from_gal() function.
  Get J2000 equatorial coordinates from galactic coordinates.
  
  Param: object_ - Object of type gal_posn representing galactic position.
  
  Returns object of type equ_posn representing J2000 equatorial position.
  """
    
  if _CUR_LIBNOVA_VER < _LN_GET_EQU2000_FROM_GAL_MIN:
        
    _posn = get_equ_from_gal(object_)
    return B1950_to_J2000(_posn)
        
  else:

    _posn = equ_posn()
    libnova.ln_get_equ2000_from_gal(object_, _posn)
    return _posn
    
    
_LN_GET_GAL_FROM_EQU2000_MIN = _hexversion((0, 12, 0))    
    
def get_gal_from_equ2000(object_):
  """
  Wrapper for libnova ln_gal_from_equ2000() function.
  Get galactic coordinates from J2000 equatorial coordinates.
  
  Param: object_ - Object of type equ_posn representing J2000 equatorial 
                   position.
  
  Returns object of type gal_posn representing galactic position.
  """
    
  if _CUR_LIBNOVA_VER < _LN_GET_GAL_FROM_EQU2000_MIN:
   
    _posn = J2000_to_B1950(object_)
    return get_gal_from_equ(_posn)
        
  else:

    _posn = gal_posn()
    libnova.ln_get_gal_from_equ2000(object_, _posn)
    return _posn



######################################################################
# python fucntion wrappers for libnova C Sidereal Time Functions
######################################################################


def get_apparent_sidereal_time(jD):
  """
  Wrapper for libnova ln_get_apparent_sidereal_time() function.
  Get apparent sidereal time from Julian day.
  
  Param: jD - UTC Julian day (float).
  
  Returns GM apparent sidereal time (float hours).
  """
    
  return libnova.ln_get_apparent_sidereal_time(jD)


def get_mean_sidereal_time(jD):
  """
  Wrapper for libnova ln_get_mean_sidereal_time() function.
  Get mean sidereal time from Julian day.
  
  Param: jD - UTC Julian day (float).
  
  Returns GM mean sidereal time (float hours).
  """
    
  return libnova.ln_get_mean_sidereal_time(jD)



######################################################################
# python fucntion wrappers for libnova C Angular Separation Functions
######################################################################


def get_angular_separation(posn1, posn2):
  """  
  Wrapper for libnova ln_get_angular_separation() function.
  Get angular separation from equatorial positions.
  
  Param: posn1 - Object of type equ_posn representing body 1 position.
  Param: posn2 - Object of type equ_posn representing body 2 position.
  
  Returns angular separation in degrees (float).
  """
    
  return libnova.ln_get_angular_separation(posn1, posn2)
     


def get_rel_posn_angle(posn1, posn2):
  """
  Wrapper for libnova ln_get_rel_posn_angle() function.
  Get relative position angle from equatorial positions.
  
  Param: posn1 - Object of type equ_posn representing body 1 position.
  Param: posn2 - Object of type equ_posn representing body 2 position.
  
  Returns position angle in degrees (float).
  """
    
  return libnova.ln_get_rel_posn_angle(posn1, posn2)
    
    
######################################################################
# python fucntion wrappers for libnova C Apparent Position Functions
######################################################################


_DEFAULT_PROPER_MOTION = equ_posn(0.0, 0.0)

_LN_GET_APPARENT_POSN_FIX_MIN = _hexversion((0, 12, 0))
_LN_GET_APPARENT_POSN_FIX_MAX = _hexversion((0, 12, 1))


def get_apparent_posn(mean_position, jD, proper_motion = None):
  """
  Wrapper for libnova ln_get_apparent_posn() function.
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
    
  # libnova ln_get_equ_pm is broken for version 0.12.0
  # internally, ln_get_apparent_posn() calls this function, breaking this also
  # use copy of version 0.11.0 code
    
  if (_CUR_LIBNOVA_VER >= _LN_GET_APPARENT_POSN_FIX_MIN) and \
    (_CUR_LIBNOVA_VER <= _LN_GET_APPARENT_POSN_FIX_MAX):
    
    _posn = get_equ_pm(mean_position, proper_motion, jD)
    _posn = get_equ_aber(_posn, jD)
    return get_equ_prec(_posn, jD)
    
  # let ln_get_apparent_posn() do the work
    
  else:
    
    _posn = equ_posn() 
    libnova.ln_get_apparent_posn(mean_position, proper_motion, jD, _posn)
    _posn.ra = range_degrees(_posn.ra)
    return _posn
               
    
######################################################################
# python fucntion wrappers for libnova C Precession Functions
######################################################################


def get_equ_prec(mean_position, jD):
  """
  Wrapper for libnova ln_get_equ_prec() function.
  Get position of celestial object accounting for precession.
  Only works for converting to and from J2000 epoch.
  
  Param: mean_position - J2000 equatorial mean position of object as type equ_posn.
  Param: jD - UTC Julian day (float) to measure position.
  
  Returns: Adjusted equatorial position of object as type equ_posn.
  """    
    
  _posn = equ_posn()    
  libnova.ln_get_equ_prec(mean_position, jD, _posn)
  _posn.ra = range_degrees(_posn.ra)
  return _posn
    
    
_LN_GET_EQU_PREC_2_MIN = _hexversion((0, 12, 0))

def get_equ_prec2(mean_position, fromJD, toJD):
  """
  Wrapper for libnova ln_get_equ_prec2() function.
  Get position of celestial object accounting for precession.
  If libnova has a version lower than 0.12.0, the ln_get_equ_prec_2()
  function is not available and the NOVAS C library adoptation 
  of precession() is called instead.
  
  Param: mean_position  - equatorial first position of object as type equ_posn.
  Param: fromJD         - UTC Julian day (float) of first time.
  Param: toJD           - UTC Julian day (float) of second time.
  
  Returns: Equatorial position of object as type equ_posn converted from
           time 1 to time 2.
  """  
    
  if _CUR_LIBNOVA_VER < _LN_GET_EQU_PREC_2_MIN:
    return get_precession(fromJD, mean_position, toJD)
    
  _posn = equ_posn()
  libnova.ln_get_equ_prec2(mean_position, fromJD, toJD, _posn)
  return _posn
    
    
    
######################################################################
# python fucntion wrappers for libnova C Nutation Functions
######################################################################


def get_nutation(jD):
  """
  Wrapper for libnova ln_get_nutation() function.
  Get nutation corrections for a given time.
  
  Param: jD - UTC Julian day (float) to measure nutation.
  
  Returns: Nutation corrections as object of type nutation.
  """
    
  _nut = nutation()
  libnova.ln_get_nutation(jD, _nut)
  return _nut
    

  
######################################################################
# python fucntion wrappers for libnova C Aberration Functions
######################################################################


def get_equ_aber(mean_position, jD): 
  """
  Wrapper for libnova ln_get_equ_aber() function.
  Get position of celestial object accounting for aberration.
  
  Param: mean_position  - J2000 equatorial mean position of object as type 
                          equ_posn.
  Param: jD             - UTC Julian day (float) to measure aberration.
  
  Returns: Adjusted equatorial position of object as type equ_posn.
  """    
    
  _posn = equ_posn()
  libnova.ln_get_equ_aber(mean_position, jD, _posn)
  return _posn

    
######################################################################
# python fucntion wrappers for libnova C Proper Motion Functions
######################################################################   
    

_LN_GET_EQU_PM_FIX_MIN = _hexversion((0, 12, 0))
_LN_GET_EQU_PM_FIX_MAX = _hexversion((0, 12, 1))

def get_equ_pm(mean_position, proper_motion, jD):
  """
  Wrapper for libnova ln_get_equ_pm() function.
  Adjusts equatorial position of a stellar object accouting for proper motion.
  
  Param: mean_position - J2000 equatorial mean position of object as type equ_posn.
  Param: proper_motion - Object of type equ_posn giving object's proper motion
                         (units are deg/year).
  Param: jD - UTC Julian day (float) to measure position.
  
  Returns: Adjusted equatorial position of object as type equ_posn.
  """
    
  # libnova ln_get_equ_pm is broken for version 0.12.0
  # use copy of version 0.11.0 code
   
  if (_CUR_LIBNOVA_VER >= _LN_GET_EQU_PM_FIX_MIN) and \
    (_CUR_LIBNOVA_VER <= _LN_GET_EQU_PM_FIX_MAX):
        
     T = (jD - J2000_UTC_JD) / 365.25
	
     # change original ra and dec to radians
            
     mean_ra = math.radians(mean_position.ra)
     mean_dec = math.radians(mean_position.dec)

     # calc proper motion
	
     mean_ra += T * math.radians(proper_motion.ra)
     mean_dec += T * math.radians(proper_motion.dec)
	
     # change to degrees 
	    
     mean_ra = range_degrees(math.degrees(mean_ra))
     mean_dec = math.degrees(mean_dec)
            
     return equ_posn(mean_ra, mean_dec)
   
    
  # otherwise, let ln_get_equ_pm do the work
    
  else:
    
    _posn = equ_posn()
    libnova.ln_get_equ_pm(mean_position, proper_motion, jD, _posn)
    return _posn
    
    

######################################################################
# python fucntion wrappers for libnova C Rise, Set, Transit functions
######################################################################


def get_object_rst(jD, observer, object_):
  """
  Wrapper for libnova ln_get_object_rst() function.
  Get rise, set, and transit times of a celstial object.
  
  Param: jD       - UTC Julian day (float) target time.
  Param: observer - object of type lnlat_posn giving observer position
  Param: object_  - object of type equ_posn giving target equatorial position
  
  Returns: Object of type rst_time giving object's ephemeris UTC times,
           or None if the object is circumpolar.
  """

  _rst = rst_time()
  if libnova.ln_get_object_rst(jD, observer, object_, _rst) != 0:
    return None
  return _rst

    


######################################################################
# python fucntion wrappers for libnova C Solar Functions
######################################################################


def get_solar_equ_coords(jD):
  """
  Wrapper for libnova ln_get_solar_equ_coords() function.
  Get Sun's apparent equatorial coordinates from Julian day.
  Accounts for aberration and precession, and nutation.
  
  Param: jD - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_solar_equ_coords(jD, _posn)
  return get_equ_prec(_posn, jD)


def get_solar_rst(jD, observer):
  """
  Wrapper for libnova ln_get_solar_rst() function.
  Get Sun's rise, transit, set times from Julian day.
  
  Param: jD       - UTC Julian day (float).
  Param: observer - Object of type lnlat_posn representing observer position.
  
  Returns Object of type rst_time represeting UTC ephemeris times,
          or None if the object is circumpolar..
  """

  _rst = rst_time()
  if libnova.ln_get_solar_rst(jD, observer, _rst) != 0:
    return None
  return _rst
    



######################################################################
# python fucntion wrappers for libnova C Jupiter Functions
######################################################################


def get_jupiter_equ_coords(jD):
  """
  Wrapper for libnova ln_get_jupiter_equ_coords() function.
  Get Jupiter's apparent equatorial coordinates from Julian day.
  Accounts for aberration and precession, but not nutation.
  
  Param: jD - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_jupiter_equ_coords(jD, _posn)
  return get_equ_prec(_posn, jD)


def get_jupiter_rst(jD, observer):
  """
  Get Jupiter's rise, transit, set times from Julian day.
  
  Param: jD       - UTC Julian day (float).
  Param: observer - Object of type lnlat_posn representing observer position.
  
  Returns Object of type rst_time represeting UTC ephemeris times,
          or None if the object is circumpolar.
  """

  _rst = rst_time()
  if libnova.ln_get_jupiter_rst(jD, observer, _rst) != 0:
    return None
  return _rst



######################################################################
# python fucntion wrappers for libnova C Saturn Functions
######################################################################


def get_saturn_equ_coords(jD):
  """
  Wrapper for libnova ln_get_saturn_equ_coords() function.    
  Get Saturn's apparent equatorial coordinates from Julian day.
  Accounts for aberration and precession, but not nutation.
  
  Param: jD - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_saturn_equ_coords(jD, _posn)
  return get_equ_prec(_posn, jD)


def get_saturn_rst(jD, observer):
  """
  Get Saturn's rise, transit, set times from Julian day.
  
  Param: jD       - UTC Julian day (float).
  Param: observer - Object of type lnlat_posn representing observer position.
  
  Returns Object of type rst_time represeting UTC ephemeris times,
          or None if the object is circumpolar.
  """

  _rst = rst_time()
  if libnova.ln_get_saturn_rst(jD, observer, _rst) != 0:
    return None 
  return _rst



######################################################################
# python fucntion wrappers for libnova C Lunar Functions
######################################################################


def get_lunar_equ_coords(jD):
  """
  Wrapper for libnova ln_get_lunar_equ_coords() function.
  Get the Moon's apparent equatorial coordinates from Julian day.
  Accounts for aberration and precession, but not nutation.
  
  Param: jD - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_lunar_equ_coords(jD, _posn)
  return get_equ_prec(_posn, jD)


def get_lunar_rst(jD, observer):
  """
  Get the Moon's rise, transit, set times from Julian day.
  
  Param: jD       - UTC Julian day (float).
  Param: observer - Object of type lnlat_posn representing observer position.
  
  Returns Object of type rst_time represeting UTC ephemeris times,
          or None if the object is circumpolar.
  """

  _rst = rst_time()
  if libnova.ln_get_lunar_rst(jD, observer, _rst) != 0:
    return None
  return _rst



######################################################################
# python fucntion wrappers for libnova C Venus Functions
######################################################################


def get_venus_equ_coords(jD):
  """
  Wrapper for libnova ln_get_venus_equ_coords() function.
  Get Venus' apparent equatorial coordinates from Julian day.
  Accounts for aberration and precession, but not nutation.
  
  Param: jD - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_venus_equ_coords(jD, _posn)
  return get_equ_prec(_posn, jD)


def get_venus_rst(jD, observer):
  """
  Get Venus' rise, transit, set times from Julian day.
  
  Param: jD       - UTC Julian day (float).
  Param: observer - Object of type lnlat_posn representing observer position.
  
  Returns Object of type rst_time represeting UTC ephemeris times,
          or None if the object is circumpolar.
  """

  _rst = rst_time()
  if libnova.ln_get_venus_rst(jD, observer, _rst) != 0:
    return None
  return _rst
    
    
    
    
######################################################################
# python fucntion wrappers for libnova C Mars Functions
######################################################################


def get_mars_equ_coords(jD):
  """
  Wrapper for libnova ln_get_mars_equ_coords() function.
  Get Mars' apparent equatorial coordinates from Julian day.
  Accounts for aberration and precession, but not nutation.
  
  Param: jD - UTC Julian day (float).
  
  Returns object of type equ_posn representing equatorial position.
  """

  _posn = equ_posn()
  libnova.ln_get_mars_equ_coords(jD, _posn)
  return get_equ_prec(_posn, jD)


def get_mars_rst(jD, observer):
  """
  Get Mars' rise, transit, set times from Julian day.
  
  Param: jD       - UTC Julian day (float).
  Param: observer - Object of type lnlat_posn representing observer position.
  
  Returns Object of type rst_time represeting UTC ephemeris times,
          or None if the object is circumpolar.
  """

  _rst = rst_time()
  if libnova.ln_get_mars_rst(jD, observer, _rst) != 0:
    return None
  return _rst



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
B1950_UTC_JD =  libnova.B1950

"""
UTC Julian day of J2000.0 coordinate epoch.
"""
J2000_UTC_JD = libnova.JD2000

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


# lookup location of TAI/UTC almanac data file in location 
# <package_dir>/refdata/astro/tai-utc.dat

def _parse_tai_file():

  import os
  from lsl.common.paths import data as dataPath

  # get path to almanac data file
     
  datName = os.path.join(dataPath, 'astro', 'tai-utc.dat') 
  if not os.path.exists(datName):
    raise RuntimeError("file %s not found" % datName)        

  # read tai-utc.dat file to get conversion info

  datFile = open(datName, 'r')
  datLines = datFile.readlines()
  datFile.close()

  lineNum = 0
  for l in datLines:

    # get UTC JD of leap second boundaries
    
    try:
      utcJD = float(l[16:26])
    except ValueError:
      raise RuntimeError("line %d of %s file not correctly formatted" % (lineNum, datName))
        
    # only get values prior to UTC JD 2441317.5 (1972 JAN  1)
    
    if utcJD < FIRST_LEAP_UTC:
      lineNum += 1
      continue    
        
    # get leap second asjustment value
    
    try:
      leapSec = float(l[36:48])
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
  
  #utcJD = float(unixTime) / SECS_IN_DAY + UNIX_OFFSET
  utcJD = get_julian_from_timet(unixTime)
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
  
  #unixTime = (utcJD - UNIX_OFFSET) * SECS_IN_DAY
  unixTime = get_timet_from_julian(utcJD)
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
    
  return ((hms.hours * 60.0 * 60.0) + (hms.minutes * 60.0) + hms.seconds)  



def deg_to_sec(degrees):
  """
  Convert longitude degrees into time seconds.
  
  Param: degrees - longitude (float degrees)
  
  Returns: Seconds (float) offset in time for longitude.
  """
    
  return (degrees * 3600.0) / 15.0
 
    
    
def get_local_sidereal_time(lng, jD):
  """
  Get apparent local sidereal time from Julian day.
  
  Param: lng  - longitude degrees (float), E = positive, W = negative 
  Param: jD   - UTC Julian day (float).
  
  Returns: Local apparent sidereal time (float hours).
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

  def __init__(self, lng = None, lat = None, elv = 0.0):
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
        
    return (geo_posn,(self.lng,self.lat,self.elv))


  def __repr__(self):
    """
    geo_posn object repr string method.
    """
        
    return "%s.%s(%s,%s,%s)" % (type(self).__module__,type(self).__name__,repr(self.lng),repr(self.lat),repr(self.elv))


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
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97,
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

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

    xx =  cz0 * ct * cz - sz0 * sz;
    yx = -sz0 * ct * cz - cz0 * sz;
    zx = -st * cz;
    xy = cz0 * ct * sz + sz0 * cz;
    yy = -sz0 * ct * sz + cz0 * cz;
    zy = -st * sz;
    xz = cz0 * st;
    yz = -sz0 * st;
    zz = ct;

    #
    #  Perform rotation.
    #

    xr = xx * pos[0] + yx * pos[1] + zx * pos[2];
    yr = xy * pos[0] + yy * pos[1] + zy * pos[2];
    zr = xz * pos[0] + yz * pos[1] + zz * pos[2];

    return (xr, yr, zr)
                    

def B1950_to_J2000(pos):
  """
  Convert B1950 epoch to J2000 epoch for equatorial coordinates.
  
  Param: pos - object of type equ_posn giving B1950 coordinates
  
  Returns: object of type equ_posn giving J2000 coordinates.
  """

  return get_equ_prec2(pos, B1950_UTC_JD, J2000_UTC_JD)


def J2000_to_B1950(pos):
  """
  Convert J2000 epoch to B1950 epoch for equatorial coordinates.
  
  Param: pos - object of type equ_posn giving J2000 coordinates
  
  Returns: object of type equ_posn giving B1950 coordinates.
  """   

  return get_equ_prec2(pos, J2000_UTC_JD, B1950_UTC_JD)
    
    

######################################################################
# LWA_Common.astro module informal unit test / demo
# give NRL observation data at current system clock time
######################################################################


if __name__ == '__main__':
    
    # get libnova version
    
    ver = get_libnova_version()
    print("Using libnova version %d.%d.%d" % (ver[0], ver[1], ver[2]))

    # input NRL longitude and latitude
    
    lng = dms(True, 77, 1, 40)
    lat = dms(False, 38, 49, 30)  
    nrl_lnlat = lnlat_posn(lng, lat)
    
    print('---------------------------------------------------------------')
    print('NRL location')
    print('---------------------------------------------------------------')
    print('Longitude:         %s (%0.3f)' % (lng, nrl_lnlat.lng))
    print('Latitude:          %s (%0.3f)' % (lat, nrl_lnlat.lat)) 
        
    # calculate offset from GM
    
    nrl_gm_hms = deg_to_hms(-nrl_lnlat.lng)
    nrl_gm_off = nrl_gm_hms.to_sec()
    print('GM Offset:         %s (%0.3f)' % (nrl_gm_hms, nrl_gm_off))
    
    # get current UTC time from system clock
    
    utc = get_julian_from_sys()
    utcd = get_date(utc)
    lcld = utcd.to_zone()
    unixt = get_timet_from_julian(utc)
    
    # caculate sidereal times
    
    gm_sid = get_apparent_sidereal_time(utc)
    nrl_sid = get_local_sidereal_time(nrl_lnlat.lng, utc)
        
    print('---------------------------------------------------------------')
    print('Current time')
    print('---------------------------------------------------------------')
    print('UTC time:          %s (%0.3f)' % (utcd, utc))
    print('NRL local time:    %s' % str(lcld))
    print('GM sidereal time:  %0.3f' % gm_sid)
    print('NRL sidereal time: %0.3f' % nrl_sid)
    print('UNIX time:         %d' % unixt)   
    
    # calculate nutation
    
    nut = get_nutation(utc)
    (nut_lng, nut_obl, nut_ecl) = nut.format()
    
    print('---------------------------------------------------------------')
    print('Nutation')
    print('---------------------------------------------------------------')
    print('Longitude nutation: %s (%0.4f)' % (nut_lng, nut.longitude))
    print('Obliquity nutation: %s (%0.4f)' % (nut_obl, nut.obliquity))
    print('Ecliptic obliquity: %s (%0.4f)' % (nut_ecl, nut.ecliptic))
    
    # calculate Solar phenomena
    
    sun_rst = get_solar_rst(utc, nrl_lnlat)
    (sun_utc_rise, sun_utc_set, sun_utc_trans) = sun_rst.format()
    sun_lcl_rise = sun_utc_rise.to_zone()
    sun_lcl_trans = sun_utc_trans.to_zone()
    sun_lcl_set = sun_utc_set.to_zone()
    
    sun_equ = get_solar_equ_coords(utc) 
    (sun_ra, sun_dec) = sun_equ.format()
    sun_hrz = sun_equ.to_hrz(nrl_lnlat, utc)
    sun_ecl = sun_equ.to_ecl(utc)
    (sun_lng, sun_lat) = sun_ecl.format()
    
    print('---------------------------------------------------------------')
    print('Sun')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (sun_ra, sun_equ.ra))
    print('DEC:               %s (%0.3f)' % (sun_dec, sun_equ.dec)) 
    print('Ecl longitude:     %s (%0.3f)' % (sun_lng, sun_ecl.lng))
    print('Ecl latitude:      %s (%0.3f)' % (sun_lat, sun_ecl.lat))             
    print('Rise:              %s (%0.3f) [%s]' % (sun_utc_rise, sun_rst.rise, sun_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (sun_utc_trans, sun_rst.transit, sun_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (sun_utc_set, sun_rst.set, sun_lcl_set))
    print('Azimuth:           %0.3f %s' % (sun_hrz.az, hrz_to_nswe(sun_hrz)))
    print('Altitude:          %0.3f' % sun_hrz.alt)
    print('Zenith:            %0.3f' % sun_hrz.zen()) 
     
    
    # calculate Lunar phenomena
    
    moon_rst = get_lunar_rst(utc, nrl_lnlat)
    (moon_utc_rise, moon_utc_set, moon_utc_trans) = moon_rst.format()
    moon_lcl_rise = moon_utc_rise.to_zone()
    moon_lcl_trans = moon_utc_trans.to_zone()
    moon_lcl_set = moon_utc_set.to_zone()
    
    moon_equ = get_lunar_equ_coords(utc) 
    (moon_ra, moon_dec) = moon_equ.format()
    moon_hrz = moon_equ.to_hrz(nrl_lnlat, utc)
    moon_ecl = moon_equ.to_ecl(utc)
    (moon_lng, moon_lat) = moon_ecl.format()
    
    moon_sun_ang = sun_equ.angular_separation(moon_equ)
    
    print('---------------------------------------------------------------')
    print('Moon')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (moon_ra, moon_equ.ra))
    print('DEC:               %s (%0.3f)' % (moon_dec, moon_equ.dec)) 
    print('Ecl longitude:     %s (%0.3f)' % (moon_lng, moon_ecl.lng))
    print('Ecl latitude:      %s (%0.3f)' % (moon_lat, moon_ecl.lat))             
    print('Rise:              %s (%0.3f) [%s]' % (moon_utc_rise, moon_rst.rise, moon_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (moon_utc_trans, moon_rst.transit, moon_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (moon_utc_set, moon_rst.set, moon_lcl_set))
    print('Azimuth:           %0.3f %s' % (moon_hrz.az, hrz_to_nswe(moon_hrz)))
    print('Altitude:          %0.3f' % moon_hrz.alt)
    print('Zenith:            %0.3f' % moon_hrz.zen())
    print('Sun angle:         %0.3f' % moon_sun_ang) 
    
    
    # calculate Venus phenomena
    
    venus_rst = get_venus_rst(utc, nrl_lnlat)
    (venus_utc_rise, venus_utc_set, venus_utc_trans) = venus_rst.format()
    venus_lcl_rise = venus_utc_rise.to_zone()
    venus_lcl_trans = venus_utc_trans.to_zone()
    venus_lcl_set = venus_utc_set.to_zone()
    
    venus_equ = get_venus_equ_coords(utc) 
    (venus_ra, venus_dec) = venus_equ.format()
    venus_hrz = venus_equ.to_hrz(nrl_lnlat, utc)
    venus_ecl = venus_equ.to_ecl(utc)
    (venus_lng, venus_lat) = venus_ecl.format()
    
    venus_sun_ang = sun_equ.angular_separation(venus_equ)
    venus_elong = sun_ecl.lng - venus_ecl.lng
    if venus_elong < 0:
        venus_dir = 'E'
    else:
        venus_dir = 'W'
    
    print('---------------------------------------------------------------')
    print('Venus')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (venus_ra, venus_equ.ra))
    print('DEC:               %s (%0.3f)' % (venus_dec, venus_equ.dec)) 
    print('Ecl longitude:     %s (%0.3f)' % (venus_lng, venus_ecl.lng))
    print('Ecl latitude:      %s (%0.3f)' % (venus_lat, venus_ecl.lat))             
    print('Rise:              %s (%0.3f) [%s]' % (venus_utc_rise, venus_rst.rise, venus_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (venus_utc_trans, venus_rst.transit, venus_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (venus_utc_set, venus_rst.set, venus_lcl_set))
    print('Azimuth:           %0.3f %s' % (venus_hrz.az, hrz_to_nswe(venus_hrz)))
    print('Altitude:          %0.3f' % venus_hrz.alt)
    print('Zenith:            %0.3f' % venus_hrz.zen())
    print('Sun angle:         %0.3f' % venus_sun_ang)
    print('Elongation:        %s %s (%0.3f)' % (deg_to_dms(venus_elong), venus_dir, venus_elong))
    
    
    # calculate Mars phenomena
    
    mars_rst = get_mars_rst(utc, nrl_lnlat)
    (mars_utc_rise, mars_utc_set, mars_utc_trans) = mars_rst.format()
    mars_lcl_rise = mars_utc_rise.to_zone()
    mars_lcl_trans = mars_utc_trans.to_zone()
    mars_lcl_set = mars_utc_set.to_zone()
    
    mars_equ = get_mars_equ_coords(utc) 
    (mars_ra, mars_dec) = mars_equ.format()
    mars_hrz = mars_equ.to_hrz(nrl_lnlat, utc)
    mars_ecl = mars_equ.to_ecl(utc)
    (mars_lng, mars_lat) = mars_ecl.format()
    
    mars_sun_ang = sun_equ.angular_separation(mars_equ)
    
    print('---------------------------------------------------------------')
    print('Mars')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (mars_ra, mars_equ.ra))
    print('DEC:               %s (%0.3f)' % (mars_dec, mars_equ.dec)) 
    print('Ecl longitude:     %s (%0.3f)' % (mars_lng, mars_ecl.lng))
    print('Ecl latitude:      %s (%0.3f)' % (mars_lat, mars_ecl.lat))             
    print('Rise:              %s (%0.3f) [%s]' % (mars_utc_rise, mars_rst.rise, mars_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (mars_utc_trans, mars_rst.transit, mars_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (mars_utc_set, mars_rst.set, mars_lcl_set))
    print('Azimuth:           %0.3f %s' % (mars_hrz.az, hrz_to_nswe(mars_hrz)))
    print('Altitude:          %0.3f' % mars_hrz.alt)
    print('Zenith:            %0.3f' % mars_hrz.zen())
    print('Sun angle:         %0.3f' % mars_sun_ang)
    
    
    # calculate Jupiter phenomena
    
    jupiter_rst = get_jupiter_rst(utc, nrl_lnlat)
    (jupiter_utc_rise, jupiter_utc_set, jupiter_utc_trans) = jupiter_rst.format()
    jupiter_lcl_rise = jupiter_utc_rise.to_zone()
    jupiter_lcl_trans = jupiter_utc_trans.to_zone()
    jupiter_lcl_set = jupiter_utc_set.to_zone()
    
    jupiter_equ = get_jupiter_equ_coords(utc) 
    (jupiter_ra, jupiter_dec) = jupiter_equ.format()
    jupiter_hrz = jupiter_equ.to_hrz(nrl_lnlat, utc)
    jupiter_ecl = jupiter_equ.to_ecl(utc)
    (jupiter_lng, jupiter_lat) = jupiter_ecl.format()
    
    jupiter_sun_ang = sun_equ.angular_separation(jupiter_equ)
    
    print('---------------------------------------------------------------')
    print('Jupiter')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (jupiter_ra, jupiter_equ.ra))
    print('DEC:               %s (%0.3f)' % (jupiter_dec, jupiter_equ.dec)) 
    print('Ecl longitude:     %s (%0.3f)' % (jupiter_lng, jupiter_ecl.lng))
    print('Ecl latitude:      %s (%0.3f)' % (jupiter_lat, jupiter_ecl.lat))             
    print('Rise:              %s (%0.3f) [%s]' % (jupiter_utc_rise, jupiter_rst.rise, jupiter_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (jupiter_utc_trans, jupiter_rst.transit, jupiter_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (jupiter_utc_set, jupiter_rst.set, jupiter_lcl_set))
    print('Azimuth:           %0.3f %s' % (jupiter_hrz.az, hrz_to_nswe(jupiter_hrz)))
    print('Altitude:          %0.3f' % jupiter_hrz.alt)
    print('Zenith:            %0.3f' % jupiter_hrz.zen())
    print('Sun angle:         %0.3f' % jupiter_sun_ang)
    
    
    # calculate Saturn phenomena
    
    saturn_rst = get_saturn_rst(utc, nrl_lnlat)
    (saturn_utc_rise, saturn_utc_set, saturn_utc_trans) = saturn_rst.format()
    saturn_lcl_rise = saturn_utc_rise.to_zone()
    saturn_lcl_trans = saturn_utc_trans.to_zone()
    saturn_lcl_set = saturn_utc_set.to_zone()
    
    saturn_equ = get_saturn_equ_coords(utc) 
    (saturn_ra, saturn_dec) = saturn_equ.format()
    saturn_hrz = saturn_equ.to_hrz(nrl_lnlat, utc)
    saturn_ecl = saturn_equ.to_ecl(utc)
    (saturn_lng, saturn_lat) = saturn_ecl.format()
    
    saturn_sun_ang = sun_equ.angular_separation(saturn_equ)
    
    print('---------------------------------------------------------------')
    print('Saturn')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (saturn_ra, saturn_equ.ra))
    print('DEC:               %s (%0.3f)' % (saturn_dec, saturn_equ.dec)) 
    print('Ecl longitude:     %s (%0.3f)' % (saturn_lng, saturn_ecl.lng))
    print('Ecl latitude:      %s (%0.3f)' % (saturn_lat, saturn_ecl.lat))             
    print('Rise:              %s (%0.3f) [%s]' % (saturn_utc_rise, saturn_rst.rise, saturn_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (saturn_utc_trans, saturn_rst.transit, saturn_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (saturn_utc_set, saturn_rst.set, saturn_lcl_set))
    print('Azimuth:           %0.3f %s' % (saturn_hrz.az, hrz_to_nswe(saturn_hrz)))
    print('Altitude:          %0.3f' % saturn_hrz.alt)
    print('Zenith:            %0.3f' % saturn_hrz.zen())
    print('Sun angle:         %0.3f' % saturn_sun_ang)
    
    # calculate SgrA phenomena
    
    sgra_j2000_equ = equ_posn(hms(17, 42, 48.1), dms(True, 28, 55, 8))
    sgra_equ = get_apparent_posn(sgra_j2000_equ, utc)
    sgra_rst = get_object_rst(utc, nrl_lnlat, sgra_equ)
    (sgra_utc_rise, sgra_utc_set, sgra_utc_trans) = sgra_rst.format()
    sgra_lcl_rise = sgra_utc_rise.to_zone()
    sgra_lcl_trans = sgra_utc_trans.to_zone()
    sgra_lcl_set = sgra_utc_set.to_zone()
    
    (sgra_ra, sgra_dec) = sgra_equ.format()
    sgra_hrz = sgra_equ.to_hrz(nrl_lnlat, utc)
    sgra_gal = sgra_equ.to_gal(utc)
    (sgra_l, sgra_b) = sgra_gal.format()
    
    sgra_sun_ang = sun_equ.angular_separation(sgra_equ)
    
    print('---------------------------------------------------------------')
    print('SgrA')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (sgra_ra, sgra_equ.ra))
    print('DEC:               %s (%0.3f)' % (sgra_dec, sgra_equ.dec)) 
    print('Gal longitude:     %s (%0.3f)' % (sgra_l, sgra_gal.l))
    print('Gal latitude:      %s (%0.3f)' % (sgra_b, sgra_gal.b))           
    print('Rise:              %s (%0.3f) [%s]' % (sgra_utc_rise, sgra_rst.rise, sgra_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (sgra_utc_trans, sgra_rst.transit, sgra_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (sgra_utc_set, sgra_rst.set, sgra_lcl_set))
    print('Azimuth:           %0.3f %s' % (sgra_hrz.az, hrz_to_nswe(sgra_hrz)))
    print('Altitude:          %0.3f' % sgra_hrz.alt)
    print('Zenith:            %0.3f' % sgra_hrz.zen())
    print('Sun angle:         %0.3f' % sgra_sun_ang)
    
    # calculate CasA phenomena
    
    casa_j2000_equ = equ_posn(hms(23, 23, 22.7), dms(False, 58, 49, 16))
    casa_equ = get_apparent_posn(casa_j2000_equ, utc)
    
    (casa_ra, casa_dec) = casa_equ.format()
    casa_hrz = casa_equ.to_hrz(nrl_lnlat, utc)
    casa_gal = casa_equ.to_gal(utc)
    (casa_l, casa_b) = casa_gal.format()
    
    casa_sun_ang = sun_equ.angular_separation(casa_equ)
    
    print('---------------------------------------------------------------')
    print('CasA')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (casa_ra, casa_equ.ra))
    print('DEC:               %s (%0.3f)' % (casa_dec, casa_equ.dec)) 
    print('Gal longitude:     %s (%0.3f)' % (casa_l, casa_gal.l))
    print('Gal latitude:      %s (%0.3f)' % (casa_b, casa_gal.b)) 
    print('Azimuth:           %0.3f %s' % (casa_hrz.az, hrz_to_nswe(casa_hrz)))
    print('Altitude:          %0.3f' % casa_hrz.alt)
    print('Zenith:            %0.3f' % casa_hrz.zen())
    print('Sun angle:         %0.3f' % casa_sun_ang)
         

    # calculate CygA phenomena
    
    cyga_j2000_equ = equ_posn(hms(19, 59, 27.8), dms(False, 40, 44, 2))
    cyga_equ = get_apparent_posn(cyga_j2000_equ, utc)
    cyga_rst = get_object_rst(utc, nrl_lnlat, cyga_equ)
    (cyga_utc_rise, cyga_utc_set, cyga_utc_trans) = cyga_rst.format()
    cyga_lcl_rise = cyga_utc_rise.to_zone()
    cyga_lcl_trans = cyga_utc_trans.to_zone()
    cyga_lcl_set = cyga_utc_set.to_zone()
    
    (cyga_ra, cyga_dec) = cyga_equ.format()
    cyga_hrz = cyga_equ.to_hrz(nrl_lnlat, utc)
    cyga_gal = cyga_equ.to_gal(utc)
    (cyga_l, cyga_b) = cyga_gal.format()
    
    cyga_sun_ang = sun_equ.angular_separation(cyga_equ)
    
    print('---------------------------------------------------------------')
    print('CygA')
    print('---------------------------------------------------------------')
    print('RA:                %s (%0.3f)' % (cyga_ra, cyga_equ.ra))
    print('DEC:               %s (%0.3f)' % (cyga_dec, cyga_equ.dec)) 
    print('Gal longitude:     %s (%0.3f)' % (cyga_l, cyga_gal.l))
    print('Gal latitude:      %s (%0.3f)' % (cyga_b, cyga_gal.b))           
    print('Rise:              %s (%0.3f) [%s]' % (cyga_utc_rise, cyga_rst.rise, cyga_lcl_rise))
    print('Transit:           %s (%0.3f) [%s]' % (cyga_utc_trans, cyga_rst.transit, cyga_lcl_trans))
    print('Set:               %s (%0.3f) [%s]' % (cyga_utc_set, cyga_rst.set, cyga_lcl_set))
    print('Azimuth:           %0.3f %s' % (cyga_hrz.az, hrz_to_nswe(cyga_hrz)))
    print('Altitude:          %0.3f' % cyga_hrz.alt)
    print('Zenith:            %0.3f' % cyga_hrz.zen())
    print('Sun angle:         %0.3f' % cyga_sun_ang)
    
    
