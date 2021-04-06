"""
Module that contains all of the relevant class to build up a representation 
of a session definition file as defined in MCS0030v5.  The hierarchy of classes
is:
  * Project - class that holds all of the information about the project (including
              the observer) and one or more sessions.  Technically, a SD file has 
              only one session but this approach allows for the generation of 
              multiple SD files from a single Project object.
  * Observer - class that hold the observer's name and numeric ID
  * Session - class that holds all of the observations associated with a particular 
              DP output.  
  * Observations - class that hold information about a particular observation.  It
                   includes a variety of attributes that are used to convert human-
                   readable inputs to SDF data values.  The observation class is 
                   further subclasses into:
                     - TBW - class for TBW observations
                     - TBN - class for TBN observations
                     - DRX - class for general DRX observation, with sub-classes:
                       * Solar - class for solar tracking
                       * Jovian - class for Jovian tracking
                     - Stepped - class for stepped observations
  * BeamStep - class that holds the information about a particular step in a Stepped
               Observation
    
All of the classes, except for Stepped and BeamStep, are complete and functional.  In 
addition, most class contain 'validate' attribute functions that can be used to 
determine if the project/session/observation are valid or not given the constraints of
the DP system.

In addition to providing the means for creating session definition files from scratch, 
this module also includes a simple parser for SD files.

.. versionchanged:: 2.0.0
    Added support for astropy.time.Time and astropy.coordinates.Angle instances

.. versionchanged:: 1.0.0
    Added the get_observation_start_stop() function.
    Renamed parse_timeString() to parse_time()
    parse_time() can now accept dates/times as timezone-aware datetime instances
    Observations can now be initialized with durations as timedelta instances
    Observations can now be initialized with RA/dec/az/alt as ephem.hours and 
    ephem.degrees instances
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import re
import copy
import math
import pytz
import ephem
import weakref
import warnings
from datetime import datetime, timedelta

from astropy import units as astrounits
from astropy.time import Time as AstroTime
from astropy.coordinates import Angle as AstroAngle

from lsl.transform import Time
from lsl.astro import utcjd_to_unix, MJD_OFFSET, DJD_OFFSET
from lsl.astro import date as astroDate, get_date as astroGetDate
from lsl.common.color import colorfy

from lsl.common.mcs import LWA_MAX_NSTD, datetime_to_mjdmpm, mjdmpm_to_datetime
from lsl.common.dp import freq_to_word, word_to_freq
from lsl.common.stations import lwa1
from lsl.reader.tbn import FILTER_CODES as TBNFilters
from lsl.reader.drx import FILTER_CODES as DRXFilters
from lsl.reader.tbw import FRAME_SIZE as TBWSize
from lsl.reader.tbn import FRAME_SIZE as TBNSize
from lsl.reader.drx import FRAME_SIZE as DRXSize

from lsl.config import LSL_CONFIG
OBSV_CONFIG = LSL_CONFIG.view('observing')

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '1.2'
__all__ = ['UCF_USERNAME_RE', 'Observer', 'ProjectOffice', 'Project', 'Session', 'Observation', 'TBW', 'TBN', 'DRX', 'Solar', 'Jovian', 'Stepped', 'BeamStep', 'parse_sdf',  'get_observation_start_stop', 'is_valid']


_dtRE = re.compile(r'^((?P<tz>[A-Z]{2,3}) )?(?P<year>\d{4})[ -/]((?P<month>\d{1,2})|(?P<mname>[A-Za-z]{3}))[ -/](?P<day>\d{1,2})[ T](?P<hour>\d{1,2}):(?P<minute>\d{1,2}):(?P<second>\d{1,2}(\.\d{1,6})?) ?(?P<tzOffset>[-+]\d{1,2}:?\d{1,2})?$')
_UTC = pytz.utc
_EST = pytz.timezone('US/Eastern')
_CST = pytz.timezone('US/Central')
_MST = pytz.timezone('US/Mountain')
_PST = pytz.timezone('US/Pacific')
_DRSUCapacityTB = 10
# Factors for computing the time it takes to read out a TBW from the number 
# of samples
_TBW_TIME_SCALE = 196000
_TBW_TIME_GAIN = 2500


# UCF Username RE
UCF_USERNAME_RE = re.compile(r'ucfuser:[ \t]*(?P<username>[a-zA-Z0-9_]+)(\/(?P<subdir>[a-zA-Z0-9\/\+\-_]+))?')


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
            raise ValueError("Unparsable time string: '%s'" % s)
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
                    raise ValueError("Unknown month abbreviation: '%s'" % monthName)
                    
            if mtch.group('tz') is None and mtch.group('tzOffset') is None:
                tz = _UTC
            elif mtch.group('tzOffset') is not None:
                tzOffsetSign = 1
                if mtch.group('tzOffset')[0] == '-':
                    tzOffsetSign = -1
                tzOffsetHours = int( mtch.group('tzOffset').replace(':', '')[1:3] )
                tzOffsetMinutes = int( mtch.group('tzOffset').replace(':', '')[3:5] )
                
                tz = pytz.FixedOffset(tzOffsetSign*(tzOffsetHours*60+tzOffsetMinutes), {})
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
                    import warnings
                    warnings.warn("Entering pytz search mode for '%s'" % tzName, RuntimeWarning)
                    
                    tzFound = False
                    tzNormal = datetime(year, month, day)
                    for tzi in pytz.common_timezones[::-1]:
                        tz = pytz.timezone(tzi)
                        try:
                            cTZName = tz.tzname(tzNormal, is_dst=False)
                        except TypeError:
                            cTZName = tz.tzname(tzNormal)
                        if cTZName == tzName:
                            tzFound = True
                            break
                            
                    if not tzFound:
                        raise ValueError("Unknown time zone: '%s'" % tzName)
                        
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
            dtObject = tz.localize(datetime(year, month, day, hour, minute, second, microsecond))
            
    # Return as UTC
    return dtObject.astimezone(_UTC)


class _TypedParentList(list):
    """
    Sub-class of list that restricts the list's contents to certain object 
    types.  Plus, it allows the entries to have a _parent reference that points
    back to who owns the list.
    """
    
    def __init__(self, allowed_types, parent=None, iterable=None):
        if not isinstance(allowed_types, (list, tuple)):
            allowed_types = [allowed_types,]
        self.allowed_types = tuple(allowed_types)
        if parent is not None:
            self.parent = weakref.proxy(parent)
        else:
            self.parent = None
            
        if iterable is not None:
            if not all([isinstance(value, self.allowed_types) for value in iterable]):
                raise TypeError("Expected one of: %s" % (', '.join([t.__name__ for t in self.allowed_types])))
            list.extend(self, iterable)
            
    def append(self, value):
        if isinstance(value, self.allowed_types):
            if self.parent is not None:
                value._parent = self.parent
            list.append(self, value)
        else:
            raise TypeError("Expected one of: %s" % (', '.join([t.__name__ for t in self.allowed_types])))
            
    def extend(self, values):
        if all([isinstance(value, self.allowed_types) for value in values]):
            if self.parent is not None:
                for value in values:
                    value._parent = self.parent
            list.extend(self, values)
            
    def insert(self, index, value):
        if isinstance(value, self.allowed_types):
            if self.parent is not None:
                value._parent = self.parent
            list.insert(self, index, value)
        else:
            raise TypeError("Expected one of: %s" % (', '.join([t.__name__ for t in self.allowed_types])))
            
    def __setitem__(self, index, value):
        if isinstance(value, self.allowed_types):
            if self.parent is not None:
                value._parent = self.parent
            list.__setitem__(self, index, value)
        else:
            raise TypeError("Expected one of: %s" % (', '.join([t.__name__ for t in self.allowed_types])))
            

class Observer(object):
    """Class to hold information about an observer."""
    
    def __init__(self, name, id, first=None, last=None):
        self.name = name
        self.first = first
        self.last = last
        self.id = int(id)
        
    @classmethod
    def autofilled(cls):
        name = OBSV_CONFIG.get('observer_name')
        id = OBSV_CONFIG.get('observer_id')
        if name is None or id is None:
            raise RuntimeError("Auto-fill values for the observer cannot be loaded from the configuration file")
        return cls(name, id)
        
    def join_name(self):
        if self.first != '':
            self.name = ', '.join([self.last, self.first])
        else:
            self.name = self.last
        
    def split_name(self):
        try:
            self.last, self.first = self.name.split(', ', 1)
        except ValueError:
            self.last = self.name
            self.first = ''


class ProjectOffice(object):
    """Class to hold comments from the LWA object office.  This class isn't really 
    needed to create SD files, but it is helpful for parsing SD files."""
    
    def __init__(self, project=None, sessions=None, observations=None):
        self.project = project
        if sessions is None:
            self.sessions = []
        else:
            self.sessions = sessions
        if observations is None:
            self.observations = []
        else:
            self.observations = observations


class Project(object):
    """
    Class to hold all the information about a specific session for a 
    project/proposal.
    
    .. versionchanged:: 1.2.1
        Added a new writeto() method to directly write the SDF to a file.
    """
    
    def __init__(self, observer, name, id, sessions=None, comments=None, project_office=None):
        if not isinstance(observer, Observer):
            raise TypeError("Expected 'observer' to be an Observer")
        self.observer = observer
        self.name = name
        self.id = id
        self.comments = comments
        self.sessions = _TypedParentList(Session, self)
        if sessions is not None:
            if isinstance(sessions, (list, tuple)):
                self.sessions.extend(sessions)
            else:
                self.sessions.append(sessions)
        if project_office is None:
            self.project_office = ProjectOffice()
        else:
            if not isinstance(project_office, ProjectOffice):
                raise TypeError("Expected 'project_office' to be a ProjectOffice")
            self.project_office = project_office
            
    @classmethod
    def autofilled(cls, sessions=None, comments=None, project_office=None):
        observer = Observer.autofilled()
        name = OBSV_CONFIG.get('project_name')
        id = OBSV_CONFIG.get('project_id')
        if name is None or id is None:
            raise RuntimeError("Auto-fill values for the project cannot be loaded from the configuration file")
        return cls(observer, name, id, sessions=session, comments=comments, project_office=project_office)
        
    def update(self):
        """Update the various sessions that are part of this project."""
        
        for ses in self.sessions:
            ses.update()
            
    def validate(self, verbose=False):
        """Examine all of the sessions and all of their observations to check
        for validity.  If everything is valid, return True.  Otherwise, return
        False."""
        
        self.update()
        
        failures = 0
        sessionCount = 1
        if len(self.id) > 8:
            if verbose:
                print("[%i] Project ID is too long" % (os.getpid(),))
            failures += 1
            
        for session in self.sessions:
            if verbose:
                print("[%i] Validating session %i" % (os.getpid(), sessionCount))
            if not session.validate(verbose=verbose):
                failures += 1
                
            if session.station != self.sessions[0].station:
                print("[%i] Session station mis-match" % os.getpid())
                failures += 1
                
            sessionCount += 1
            
        if failures == 0:
            return True
        else:
            return False
            
    @staticmethod
    def _render_file_size(size):
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
            
        return "%.2f %s" % (size, units)
        
    @staticmethod
    def _render_bandwidth(filter, filter_codes):
        """Convert a filter number to an easy-to-use string."""
        
        if filter_codes[filter] > 1e6:
            return "%.3f MHz" % (filter_codes[filter]/1e6,)
        elif filter_codes[filter] > 1e3:
            return "%.3f kHz" % (filter_codes[filter]/1e3,)
        else:
            return "%.3f Hz" % (filter_codes[filter],)
            
    def append(self, newSession):
        """Add a new Session to the list of sessions."""
        
        self.sessions.append(newSession)
        
    def render(self, session=0, verbose=False):
        """Create a session definition file that corresponds to the specified 
        session.  Returns the SD file's contents as a string."""
        
        if not self.validate(verbose=verbose) :
            raise RuntimeError("Invalid session/observation parameters.  Aborting.")
        if session >= len(self.sessions):
            raise IndexError("Invalid session index")
        
        self.sessions[session].update()
        self.sessions[session].observations.sort()
        
        ses = self.sessions[session]
        try:
            # Try to pull out the project office comments about the session
            pos = self.project_office.sessions[session]
        except:
            pos = None
        try:
            # Try to pull out the project office comments about the observations
            poo = self.project_office.observations[session]
        except:
            poo = []
        # Enforce that the number of project office observation comments match the
        # actual number of observations
        while (len(ses.observations) - len(poo)) > 0:
            poo.append(None)
            
        # Combine the session comments together in an intelligent fashion
        ## Observer comments
        if ses.ucf_username is not None:
            clean = ''
            if ses.comments:
                clean = UCF_USERNAME_RE.sub('', ses.comments)
            ses.comments = 'ucfuser:%s' % ses.ucf_username
            if len(clean) > 0:
                ses.comments += ';;%s' % clean
        ## Project office comments, including the data return method
        if pos != 'None' and pos is not None:
            pos = 'Requested data return method is %s;;%s' % (ses.dataReturnMethod, pos)
            
        ## PI Information
        output = ""
        output += "PI_ID            %s\n" % (self.observer.id,)
        output += "PI_NAME          %s\n" % (self.observer.name,)
        output += "\n"
        
        ## Project Information
        output += "PROJECT_ID       %s\n" % (self.id,)
        output += "PROJECT_TITLE    %s\n" % (self.name,)
        output += "PROJECT_REMPI    %s\n" % (self.comments[:4090] if self.comments else 'None provided',)
        output += "PROJECT_REMPO    %s\n" % (self.project_office.project,)
        output += "\n"
        
        ## Session Information
        output += "SESSION_ID       %s\n" % (ses.id,)
        output += "SESSION_TITLE    %s\n" % ('None provided' if ses.name is None else ses.name,)
        output += "SESSION_REMPI    %s\n" % (ses.comments[:4090] if ses.comments else 'None provided',)
        output += "SESSION_REMPO    %s\n" % ("Requested data return method is %s" % ses.dataReturnMethod if pos == 'None' or pos is None else pos,)
        if ses.configuration_authority != 0:
            output += "SESSION_CRA      %i\n" % (ses.configuration_authority,)
        if ses.drx_beam != -1:
            output += "SESSION_DRX_BEAM %i\n" % (ses.drx_beam,)
        if ses.spcSetup[0] != 0 and ses.spcSetup[1] != 0:
            output += "SESSION_SPC      %i %i%s\n" % (ses.spcSetup[0], ses.spcSetup[1], '' if ses.spcMetatag == None else ses.spcMetatag)
        for component in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
            if ses.recordMIB[component] != -1:
                output += "SESSION_MRP_%s  %i\n" % (component, ses.recordMIB[component])
        for component in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
            if ses.updateMIB[component] != -1:
                output += "SESSION_MUP_%s  %i\n" % (component, ses.updateMIB[component])
        if ses.include_mcssch_log:
            output += "SESSION_LOG_SCH  %i\n" % (ses.include_mcssch_log,)
        if ses.include_mcsexe_log:
            output += "SESSION_LOG_EXE  %i\n" % (ses.include_mcsexe_log,)
        if ses.include_station_smib:
            output += "SESSION_INC_SMIB %i\n" % (ses.include_station_smib,)
        if ses.include_station_design:
            output += "SESSION_INC_DES  %i\n" % (ses.include_station_design,)
        output += "\n"
        
        ## Observations
        for i,obs in enumerate(ses.observations):
            obsID = i + 1
            
            output += "OBS_ID           %i\n" % (obsID,)
            output += "OBS_TITLE        %s\n" % (obs.name if obs.name else 'None provided',)
            output += "OBS_TARGET       %s\n" % (obs.target if obs.target else 'None provided',)
            output += "OBS_REMPI        %s\n" % (obs.comments[:4090] if obs.comments else 'None provided',)
            output += "OBS_REMPO        %s\n" % ("Estimated data volume for this observation is %s" % self._render_file_size(obs.dataVolume) if poo[i] == 'None' or poo[i] == None else poo[i],)
            output += "OBS_START_MJD    %i\n" % (obs.mjd,)
            output += "OBS_START_MPM    %i\n" % (obs.mpm,)
            output += "OBS_START        %s\n" % (obs.start.strftime("%Z %Y/%m/%d %H:%M:%S") if isinstance(obs.start, datetime) else obs.start,)
            output += "OBS_DUR          %i\n" % (obs.dur,)
            output += "OBS_DUR+         %s\n" % (obs.duration,)
            output += "OBS_MODE         %s\n" % (obs.mode,)
            if obs.beamDipole is not None:
                output += "OBS_BDM          %i %6.4f %6.4f %s\n" % (tuple(obs.beamDipole))
            if obs.mode == 'TBN':
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (self._render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_RADEC':
                output += "OBS_RA           %.9f\n" % (obs.ra,)
                output += "OBS_DEC          %+.9f\n" % (obs.dec,)
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (self._render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_SOL':
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (self._render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_JOV':
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (self._render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'STEPPED':
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (self._render_bandwidth(obs.filter, obs.filter_codes),)
                output += "OBS_STP_N        %i\n" % (len(obs.steps),)
                output += "OBS_STP_RADEC    %i\n" % (obs.steps[0].is_radec,)
                for j,step in enumerate(obs.steps):
                    stpID = j + 1
                    
                    output += "OBS_STP_C1[%i]      %.9f\n" % (stpID, step.c1)
                    output += "OBS_STP_C2[%i]      %+.9f\n" % (stpID, step.c2)
                    output += "OBS_STP_T[%i]       %i\n" % (stpID, step.dur)
                    output += "OBS_STP_FREQ1[%i]   %i\n" % (stpID, step.freq1)
                    output += "OBS_STP_FREQ1+[%i]  %.9f MHz\n" % (stpID, step.frequency1/1e6)
                    output += "OBS_STP_FREQ2[%i]   %i\n" % (stpID, step.freq2)
                    output += "OBS_STP_FREQ2+[%i]  %.9f MHz\n" % (stpID, step.frequency2/1e6)
                    output += "OBS_STP_B[%i]       %s\n" % (stpID, step.beam)
                    if step.beam == 'SPEC_DELAYS_GAINS':
                        for k,delay in enumerate(step.delays):
                            dlyID = k + 1
                            
                            output += "OBS_BEAM_DELAY[%i][%i] %i\n" % (stpID, dlyID, delay)
                        for k,gain in enumerate(step.gains):
                            gaiID = k + 1
                            
                            output += "OBS_BEAM_GAIN[%i][%i][1][1] %i\n" % (stpID, gaiID, gain[0][0])
                            output += "OBS_BEAM_GAIN[%i][%i][1][2] %i\n" % (stpID, gaiID, gain[0][1])
                            output += "OBS_BEAM_GAIN[%i][%i][2][1] %i\n" % (stpID, gaiID, gain[1][0])
                            output += "OBS_BEAM_GAIN[%i][%i][2][2] %i\n" % (stpID, gaiID, gain[1][1])
            ## FEE power settings
            if all(j == obs.fee_power[0] for j in obs.fee_power):
                ### All the same
                if obs.fee_power[0][0] != -1 and obs.fee_power[0][1] != -1:
                    output += "OBS_FEE[%i][1]  %i\n" % (0, obs.fee_power[0][0])
                    output += "OBS_FEE[%i][2]  %i\n" % (0, obs.fee_power[0][1])
            else:
                ### Some different
                for j,fee in enumerate(obs.fee_power):
                    feeID = j + 1
                    
                    if fee[0] != -1:
                        output += "OBS_FEE[%i][1]  %i\n" % (feeID, fee[0])
                    if fee[1] != -1:
                        output += "OBS_FEE[%i][2]  %i\n" % (feeID, fee[1])
            ## ASP filter setting
            if all(j == obs.asp_filter[0] for j in obs.asp_filter):
                ### All the same
                if obs.asp_filter[0] != -1:
                    output += "OBS_ASP_FLT[%i]  %i\n" % (0, obs.asp_filter[0])
            else:
                ### Some different
                for j,flt in enumerate(obs.asp_filter):
                    fltID = j + 1
                    
                    if flt != -1:
                        output += "OBS_ASP_FLT[%i]  %i\n" % (fltID, flt)
            ## First attenuator setting
            if all(j == obs.asp_atten_1[0] for j in obs.asp_atten_1):
                ### All the same
                if obs.asp_atten_1[0] != -1:
                    output += "OBS_ASP_AT1[%i]  %i\n" % (0, obs.asp_atten_1[0])
            else:
                ### Some different
                for j,at1 in enumerate(obs.asp_atten_1):
                    at1ID = j + 1
                    
                    if at1 != -1:
                        output += "OBS_ASP_AT1[%i]  %i\n" % (at1ID, at1)
            ## Second attenuator setting
            if all(j == obs.asp_atten_2[0] for j in obs.asp_atten_2):
                ### All the same
                if obs.asp_atten_2[0] != -1:
                    output += "OBS_ASP_AT2[%i]  %i\n" % (0, obs.asp_atten_2[0])
            else:
                ### Some different
                for j,at2 in enumerate(obs.asp_atten_2):
                    at2ID = j + 1
                    
                    if at2 != -1:
                        output += "OBS_ASP_AT2[%i]  %i\n" % (at2ID, at2)
            ## Second attenuator setting
            if all(j == obs.asp_atten_split[0] for j in obs.asp_atten_split):
                ### All the same
                if obs.asp_atten_split[0] != -1:
                    output += "OBS_ASP_ATS[%i]  %i\n" % (0, obs.asp_atten_split[0])
            else:
                ### Some different
                for j,ats in enumerate(obs.asp_atten_split):
                    atsID = j + 1
                    
                    if ats != -1:
                        output += "OBS_ASP_ATS[%i]  %i\n" % (atsID, ats)
            ## TBW settings
            if obs.mode == 'TBW':
                output += "OBS_TBW_BITS     %i\n" % (obs.bits,)
                output += "OBS_TBW_SAMPLES  %i\n" % (obs.samples,)
            ## TBN gain
            elif obs.mode == 'TBN':
                if obs.gain != -1:
                    output += "OBS_TBN_GAIN     %i\n" % (obs.gain,)
            ## DRX gain
            else:
                if obs.gain != -1:
                    output += "OBS_DRX_GAIN     %i\n" % (obs.gain,)
            output += "\n"
            
        return output
        
    def writeto(self, filename, session=0, verbose=False, overwrite=False):
        """Create a session definition file that corresponds to the specified 
        session and write it to the provided filename."""
        
        if os.path.exists(filename) and not overwrite:
            raise RuntimeError("'%s' already exists" % filename)
            
        output = self.render(session=session, verbose=verbose)
        with open(filename, 'w') as fh:
            fh.write(output)


class Observation(object):
    """
    Class to hold the specifics of an observations.  It currently
    handles TBW, TBN, TRK_RADEC, TRK_SOL, TRK_JOV, and Stepped
    
    .. versionchanged:: 1.0.0
        Added support for RA/dec values as ephem.hours/ephem.degrees instances
    """
    
    _parent = None
    
    id = 1
    dur = 0
    
    def __init__(self, name, target, start, duration, mode, ra, dec, frequency1, frequency2, filter, gain=-1, max_snr=False, comments=None):
        self.name = name
        self.target = target
        self.ra = ra
        self.dec = dec
        self.start = start
        self.duration = duration
        self.mode = mode
        self.beamDipole = None
        self.frequency1 = float(frequency1)
        self.frequency2 = float(frequency2)
        self.filter = int(filter)
        self.max_snr = bool(max_snr)
        self.comments = comments
        
        self.beam = None
        self.dataVolume = None
        
        self.fee_power = [[-1,-1] for n in range(LWA_MAX_NSTD)]
        self.asp_filter = [-1 for n in range(LWA_MAX_NSTD)]
        self.asp_atten_1 = [-1 for n in range(LWA_MAX_NSTD)]
        self.asp_atten_2 = [-1 for n in range(LWA_MAX_NSTD)]
        self.asp_atten_split = [-1 for n in range(LWA_MAX_NSTD)]

        self.gain = int(gain)
        
        self.update()
        
    def __str__(self):
        """Return a nice string to describe the observation."""
        
        return "%s Obs. of '%s':\n Start %s\n Duration %s\n Filter: %i\n Frequency: %.3f; %.3f Hz\n RA: %.6f hr\n Dec. %.6f d\n" % (self.mode, self.name, self.start, self.duration, self.filter, self.frequency1, self.frequency2, self.ra, self.dec)
        
    def update(self):
        """Update the computed parameters from the string values."""
        
        self.beam = self.get_beam_type()
        self.dataVolume = self.estimate_bytes()
        
    @property
    def start(self):
        """Start time."""
        
        utc = mjdmpm_to_datetime(self.mjd, self.mpm)
        return utc.strftime("UTC %Y/%m/%d %H:%M:%S.%f")
        
    @start.setter
    def start(self, value):
        utc = parse_time(value)
        self.mjd, self.mpm = datetime_to_mjdmpm(utc)
        
    @property
    def duration(self):
        """Duration in seconds."""
        
        s, ms = self.dur//1000, (self.dur%1000)/1000.0
        h = s // 3600
        m = (s // 60) % 60
        s = s % 60
    
        return "%i:%02i:%06.3f" % (h, m, s+ms)
        
    @duration.setter
    def duration(self, value):
        if isinstance(value, str):
            fields = value.split(':')
            s = float(fields.pop())
            try:
                m = int(fields.pop(), 10)
            except IndexError:
                m = 0
            try:
                h = int(fields.pop(), 10)
            except IndexError:
                h = 0
            seconds = h*3600 + m*60 + int(s)
            ms = int(round((s - int(s))*1000))
            if ms >= 1000:
                seconds += 1
                ms -= 1000
            seconds = seconds + ms/1000.0
            
        elif isinstance(value, timedelta):
            seconds = value.days*86400 + value.seconds
            ms = int(round(value.microseconds/1000.0))/1000.0
            seconds = seconds + ms
            
        elif isinstance(value, astrounits.quantity.Quantity):
            seconds = seconds.to('s').value
            
        else:
            seconds = value
            
        self.dur = int(round(seconds*1000))
        
    @property
    def ra(self):
        """Target RA (J2000)."""
        
        return self._ra
        
    @ra.setter
    def ra(self, value):
        if isinstance(value, ephem.Angle):
            value = value * 12.0/math.pi
        elif isinstance(value, AstroAngle):
            value = value.to('hourangle').value
        elif isinstance(value, str):
            value = AstroAngle(value).to('hourangle').value
        if value < 0.0 or value >= 24.0:
            raise ValueError("Invalid value for RA '%.6f' hr" % value)
        self._ra = value
        
    @property
    def dec(self):
        """Target dec. (J2000)."""
        
        return self._dec
        
    @dec.setter
    def dec(self, value):
        if isinstance(value, ephem.Angle):
            value = value * 180.0/math.pi
        elif isinstance(value, AstroAngle):
            value = value.to('deg').value
        elif isinstance(value, str):
            value = AstroAngle(value).to('deg').value
        if value < -90.0 or value > 90.0:
            raise ValueError("Invalid value for dec. '%+.6f' deg" % value)
        self._dec = value
        
    @property
    def frequency1(self):
        """Tuning 1 frequency in Hz."""
        
        return word_to_freq(self.freq1)
        
    @frequency1.setter
    def frequency1(self, value):
        if isinstance(value, astrounits.quantity.Quantity):
            value = value.to('Hz').value
        self.freq1 = freq_to_word(float(value))
        
    @property
    def frequency2(self):
        """Tuning 2 frequency in Hz."""
        
        return word_to_freq(self.freq2)
        
    @frequency2.setter
    def frequency2(self, value):
        if isinstance(value, astrounits.quantity.Quantity):
            value = value.to('Hz').value
        self.freq2 = freq_to_word(float(value))
        
    def get_beam_type(self):
        """Return a valid value for beam type based on whether maximum S/N beam 
        forming has been requested."""
        
        if self.max_snr:
            return 'MAX_SNR'
        else:
            return 'SIMPLE'
    
    def estimate_bytes(self):
        """Place holder for functions that return the estimate size of the data
        set being defined by the observation."""
        
        raise NotImplementedError
        
    @property
    def fixed_body(self):
        """Place holder for functions that return ephem.Body objects (or None)
        that define the pointing center of the observation."""
        
        return None
        
    @property
    def target_visibility(self):
        """Place holder for functions that return the fractional visibility of the 
        target during the observation period."""
        
        return 1.0
    
    def validate(self, verbose=False):
        """Place holder for functions that evaluate the observation and return True 
        if it is valid, False otherwise."""
        
        raise NotImplementedError
        
    def _validate_asp(self, verbose=False):
        """Evaulate the FEE and ASP options associated with an observation and
        return True if valid, False otherwise."""
        
        station = lwa1
        if self._parent is not None:
            station = self._parent.station
        is_dp = station.interface.backend == 'dp'
        nstand = station.interface.get_module('mcs').LWA_MAX_NSTD
                
        failures = 0
        # FEE
        if len(self.fee_power) < nstand:
            failures += 1
            if verbose:
                print("[%i] Error: Invalid number of FEE power settings (%i != %i)" % (os.getpid(), len(self.fee_power), nstand))
        for f,fee in enumerate(self.fee_power):
            if not isinstance(fee, (tuple, list)):
                failures += 1
                if verbose:
                    print("[%i] Error: Expected a tuple or list for the FEE %i power setting" % (os.getpid(), f))
                continue
            if len(fee) != 2:
                failures += 1
                if verbose:
                    print("[%i] Error: Invalid number of polarizations on FEE %i (%i != 2)" % (os.getpid(), f, len(fee)))
                continue
            for p in (0, 1):
                if fee[p] not in (-1, 0, 1):
                    failures += 1
                    if verbose:
                        print("[%i] Error: Invalid power setting on FEE %i, polarization %i '%i'" % (os.getpid(), f, p, fee[p]))
                        
        # ASP
        ## Filter
        if len(self.asp_filter) < nstand:
            failures += 1
            if verbose:
                print("[%i] Error: Invalid number of ASP filter settings (%i < %i)" % (os.getpid(), len(self.asp_filter), nstand))
        for f,filt in enumerate(self.asp_filter):
            if is_dp and filt > 3:
                warnings.warn("ASP filter %i is degenerate with %i for DP-based stations" % (filt, filt-4), RuntimeWarning)
                
            if filt not in (-1, 0, 1, 2, 3, 4, 5):
                failures += 1
                if verbose:
                    print("[%i] Error: Invalid ASP filter setting on stand %i '%i'" % (os.getpid(), f, filt))
        ## AT1/AT2/ATS
        if len(self.asp_atten_1) < nstand:
            failures += 1
            if verbose:
                print("[%i] Error: Invalid number of ASP attenuator 1 settings (%i < %i)" % (os.getpid(), len(self.asp_atten_1), nstand))
        for f,atten in enumerate(self.asp_atten_1):
            if atten < -1 or atten > 15:
                failures += 1
                if verbose:
                    print("[%i] Error: Invalid ASP attenuator 1 setting on stand %i '%i'" % (os.getpid(), f, atten))
        if len(self.asp_atten_2) < nstand:
            failures += 1
            if verbose:
                print("[%i] Error: Invalid number of ASP attenuator 2 settings (%i < %i)" % (os.getpid(), len(self.asp_atten_2), nstand))
        for f,atten in enumerate(self.asp_atten_2):
            if atten < -1 or atten > 15:
                failures += 1
                if verbose:
                    print("[%i] Error: Invalid ASP attenuator 2 setting on stand %i '%i'" % (os.getpid(), f, atten))
        if len(self.asp_atten_split) < nstand:
            failures += 1
            if verbose:
                print("[%i] Error: Invalid number of ASP attenuator split settings (%i < %i)" % (os.getpid(), len(self.asp_atten_split), nstand))
        for f,atten in enumerate(self.asp_atten_split):
            if atten < -1 or atten > 15:
                failures += 1
                if verbose:
                    print("[%i] Error: Invalid ASP attenuator split setting on stand %i '%i'" % (os.getpid(), f, atten))
                    
        # Any failures indicates a bad FEE/ASP configuration
        if failures == 0:
            return True
        else:
            return False
            
    def __eq__(self, other):
        if isinstance(other, Observation):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf == startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __ne__(self, other):
        if isinstance(other, Observation):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf != startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __gt__(self, other):
        if isinstance(other, Observation):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf > startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __ge__(self, other):
        if isinstance(other, Observation):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf >= startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __lt__(self, other):
        if isinstance(other, Observation):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf < startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __le__(self, other):
        if isinstance(other, Observation):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf <= startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)


class TBW(Observation):
    """Sub-class of Observation specifically for TBW observations.  It features a
    reduced number of parameters needed to setup the observation and provides extra
    information about the number of data bits and the number of samples.
    
    .. note::
        TBW read-out times in ms are calculated using (samples/196000+1)*5000 per
        MCS
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string)
     * integer number of samples
    
    Optional Keywords:
     * bits - number of data bits (4 or 12)
     * comments - comments about the observation
    """
    
    def __init__(self, name, target, start, samples, bits=12, comments=None):
        self.samples = int(samples)
        self.bits = int(bits)
        
        duration = (self.samples / _TBW_TIME_SCALE + 1)*_TBW_TIME_GAIN
        durStr = '%02i:%02i:%06.3f' % (int(duration/1000.0)/3600, int(duration/1000.0)%3600/60, duration/1000.0%60)
        Observation.__init__(self, name, target, start, durStr, 'TBW', 0.0, 0.0, 0.0, 0.0, 1, comments=comments)
        
    def estimate_bytes(self):
        """Estimate the data volume for the specified type and duration of 
        observations.  For TBW:
        
            bytes = samples / samplesPerFrame * 1224 bytes * 260 stands
        """
        
        SamplesPerFrame = 400
        if self.bits == 4:
            SamplesPerFrame = 1200
        nFrames = self.samples / SamplesPerFrame
        nBytes = nFrames * TBWSize * LWA_MAX_NSTD
        return nBytes
        
    def validate(self, verbose=False):
        """Evaluate the observation and return True if it is valid, False
        otherwise."""
        
        self.update()
        
        failures = 0
        # Basic - Sample size and data bits agreement
        if self.bits not in [4, 12]:
            if verbose:
                print("[%i] Error: Invalid number of data bits '%i'" % (os.getpid(), self.bits))
            failures += 1
        if self.bits == 12 and self.samples > 12000000:
            if verbose:
                print("[%i] Error: Invalid number of samples for 12-bit data (%i > 12000000)" % (os.getpid(), self.samples))
            failures += 1
        if self.bits == 4 and self.samples > 36000000:
            if verbose:
                print("[%i] Error: Invalid number of samples for 4-bit data (%i > 36000000)" % (os.getpid(), self.samples))
            failures += 1
            
        # Advanced - Data Volume
        if self.dataVolume >= (_DRSUCapacityTB*1024**4):
            if verbose:
                print("[%i] Error: Data volume exceeds %i TB DRSU limit" % (os.getpid(), _DRSUCapacityTB))
            failures += 1
            
        # Advanced - ASP
        failures += not self._validate_asp(verbose=verbose)
        
        # Any failures indicates a bad observation
        if failures == 0:
            return True
        else:
            return False


class TBN(Observation):
    """Sub-class of Observation specifically for TBN observations.   It features a
    reduced number of parameters needed to setup the observation.
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * observation duration (HH:MM:SS.SSS string or timedelta instance)
     * observation frequency (Hz)
     * integer filter code
    
    Optional Keywords:
     * comments - comments about the observation
    """
    
    def __init__(self, name, target, start, duration, frequency, filter, gain=-1, comments=None):
        self.filter_codes = TBNFilters
        Observation.__init__(self, name, target, start, duration, 'TBN', 0.0, 0.0, frequency, 0.0, filter, gain=gain, comments=comments)
        
    def estimate_bytes(self):
        """Estimate the data volume for the specified type and duration of 
        observations.  For TBN:
        
            bytes = duration * sample_rate / 512 * 1048 bytes * 260 stands * 2 pols.
        """
        
        try:
            nFrames = self.dur/1000.0 * self.filter_codes[self.filter] / 512
        except KeyError:
            nFrames = 0
        nBytes = nFrames * TBNSize * LWA_MAX_NSTD * 2
        return nBytes
        
    def validate(self, verbose=False):
        """Evaluate the observation and return True if it is valid, False
        otherwise.
        
        ..note::
            This version of sdf allows for TBN tuning between 5 and 93 MHz.
        """
        
        self.update()
        
        station = lwa1
        if self._parent is not None:
            station = self._parent.station
        backend = station.interface.get_module('backend')
        
        failures = 0
        # Basic - Duration, frequency, and filter code values
        if self.dur < 1:
            if verbose:
                print("[%i] Error: Specified a duration of length zero" % os.getpid())
            failures += 1
        if self.freq1 < backend.TBN_TUNING_WORD_MIN or self.freq1 > backend.TBN_TUNING_WORD_MAX:
            if verbose:
                print("[%i] Error: Specified frequency is outside of DP tuning range" % os.getpid())
            failures += 1
        if self.filter not in [1, 2, 3, 4, 5, 6, 7]:
            if verbose:
                print("[%i] Error: Invalid filter code '%i'" % (os.getpid(), self.filter))
            failures += 1
            
        # Advanced - Data Volume
        if self.dataVolume >= (_DRSUCapacityTB*1024**4):
            if verbose:
                print("[%i] Error: Data volume exceeds %i TB DRSU limit" % (os.getpid(), _DRSUCapacityTB))
            failures += 1
            
        # Advanced - ASP
        failures += not self._validate_asp(verbose=verbose)
        
        # Any failures indicates a bad observation
        if failures == 0:
            return True
        else:
            return False


class DRX(Observation):
    """Sub-class of Observation specifically for DRX-style observations.
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * observation duration (HH:MM:SS.SSS string or timedelta instance)
     * observation RA in hours, J2000.0 or ephem.hours instance
     * observation Dec in degrees, J2000.0 or ephem.hours instance
     * observation tuning frequency 1 (Hz)
     * observation tuning frequency 1 (Hz)
     * integer filter code
    
    Optional Keywords:
     * max_snr - specifies if maximum signal-to-noise beam forming is to be used
                 (default = False)
     * comments - comments about the observation
    """
    
    filter_codes = DRXFilters
    
    def __init__(self, name, target, start, duration, ra, dec, frequency1, frequency2, filter, gain=-1, max_snr=False, comments=None):
        Observation.__init__(self, name, target, start, duration, 'TRK_RADEC', ra, dec, frequency1, frequency2, filter, gain=gain, max_snr=max_snr, comments=comments)
        
    def set_beamdipole_mode(self, stand, beam_gain=0.04, dipole_gain=1.0, pol='X'):
        """Convert the current observation to a 'beam-dipole mode' 
        observation with the specified stand.  Setting the stand to zero
        will disable the 'beam-dipole mode' for this observation'.
        
        Keywords:
         * beam_gain - BAM gain to use for each dipole in the beam
                     default: 0.04; range: 0.0 to 1.0
         * dipole_gain - BAM gain to use for the single dipole
                        default: 1.0; range: 0.0 to 1.0
         * pol - Polarization to record  default: "X"
        """
        
        # Validate
        if stand < 0 or stand > LWA_MAX_NSTD:
            raise ValueError("Stand number %i is out of range: 0 <= stand <= %i" % (stand, LWA_MAX_NSTD))
        if beam_gain < 0.0 or beam_gain > 1.0:
            raise ValueError("Beam BAM gain is out of range: 0.0 <= beam_gain <= 1.0")
        if dipole_gain < 0.0 or dipole_gain > 1.0:
            raise ValueError("Dipole BAM gain is out of range: 0.0 <= dipole_gain <= 1.0")
        if pol.upper() not in ('X', 'Y'):
            raise ValueError("Unknown polarization.  Valid values are 'X' and 'Y'")
        
        # Go
        if stand == 0:
            ## Disable beam-dipole mode
            self.beamDipole = None
        else:
            ## Stand -> DP Stand
            station = lwa1
            if self._parent is not None:
                station = self._parent.station
                
            for ant in station.antennas:
                if ant.stand.id == stand:
                    dpStand = (ant.digitizer+1)/2
                    
            self.beamDipole = [dpStand, beam_gain, dipole_gain, pol.upper()]
            
    def estimate_bytes(self):
        """Estimate the data volume for the specified type and duration of 
        observations.  For DRX:
        
            bytes = duration * sample_rate / 4096 * 4128 bytes * 2 tunings * 2 pols.
        """
        
        try:
            sample_rate = self.filter_codes[self.filter]
        except KeyError:
            sample_rate = 0.0
            
        data_rate = DRXSize * 4 * sample_rate / 4096
        if self._parent is not None:
            nchan, nwin = self._parent.spcSetup
            nprod = 4 if self._parent.spcMetatag in ('{Stokes=IQUV}',) else 2
            SPCSize = 76 + nchan*2*nprod*4
            try:
                data_rate = SPCSize * sample_rate / (nchan*nwin)
            except ZeroDivisionError:
                pass
                
        return self.dur/1000.0 * data_rate
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the observation is 
        pointed.  None if the observation mode is an all-sky mode."""
        
        pnt = ephem.FixedBody()
        pnt._ra = self.ra / 12.0 * math.pi
        pnt._dec = self.dec / 180.0 * math.pi
        pnt._epoch = ephem.J2000
        return pnt
        
    @property
    def target_visibility(self):
        """Return the fractional visibility of the target during the observation 
        period."""
        
        station = lwa1
        if self._parent is not None:
            station = self._parent.station
            
        pnt = self.fixed_body
        
        vis = 0
        cnt = 0
        dt = 0.0
        max_alt = 0.0
        while dt <= self.dur/1000.0:
            station.date = self.mjd + (self.mpm/1000.0 + dt)/3600/24.0 + MJD_OFFSET - DJD_OFFSET
            pnt.compute(station)
            max_alt = max([max_alt, pnt.alt])
            
            cnt += 1
            if pnt.alt > 0:
                vis += 1
                
            dt += 300.0
            
        if max_alt < 20*math.pi/180:
            #warnings.warn("Maximum altitude for this observation is %.1f degrees" % (max_alt*180/math.pi))
            pass
            
        return float(vis)/float(cnt)
        
    def validate(self, verbose=False):
        """Evaluate the observation and return True if it is valid, False
        otherwise."""
        
        self.update()
        
        station = lwa1
        if self._parent is not None:
            station = self._parent.station
        backend = station.interface.get_module('backend')
        
        failures = 0
        # Basic - Duration, frequency, and filter code values
        if self.dur < 1:
            if verbose:
                print("[%i] Error: Specified a duration of length zero" % os.getpid())
            failures += 1
        if self.freq1 < backend.DRX_TUNING_WORD_MIN or self.freq1 > backend.DRX_TUNING_WORD_MAX:
            if verbose:
                print("[%i] Error: Specified frequency for tuning 1 is outside of DP tuning range" % os.getpid())
            failures += 1
        if (self.freq2 < backend.DRX_TUNING_WORD_MIN or self.freq2 > backend.DRX_TUNING_WORD_MAX) and self.freq2 != 0:
            if verbose:
                print("[%i] Error: Specified frequency for tuning 2 is outside of DP tuning range" % os.getpid())
            failures += 1
        if self.filter not in [1, 2, 3, 4, 5, 6, 7]:
            if verbose:
                print("[%i] Error: Invalid filter code '%i'" % (os.getpid(), self.filter))
            failures += 1
            
        # Advanced - Target Visibility
        if self.target_visibility < 1.0:
            if verbose:
                print("[%i] Error: Target is only above the horizon for %.1f%% of the observation" % (os.getpid(), self.target_visibility*100.0))
            failures += 1
            
        # Advanced - Data Volume
        if self.dataVolume >= (_DRSUCapacityTB*1024**4):
            if verbose:
                print("[%i] Error: Data volume exceeds %i TB DRSU limit" % (os.getpid(), _DRSUCapacityTB))
            failures += 1
            
        # Advanced - ASP
        failures += not self._validate_asp(verbose=verbose)
        
        # Any failures indicates a bad observation
        if failures == 0:
            return True
        else:
            return False


class Solar(DRX):
    """Sub-class of DRX specifically for Solar DRX observations.   It features a
    reduced number of parameters needed to setup the observation.
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * observation duration (HH:MM:SS.SSS string or timedelta instance)
     * observation tuning frequency 1 (Hz)
     * observation tuning frequency 1 (Hz)
     * integer filter code
    
    Optional Keywords:
     * max_snr - specifies if maximum signal-to-noise beam forming is to be used
                 (default = False)
     * comments - comments about the observation
    """
    
    def __init__(self, name, target, start, duration, frequency1, frequency2, filter, gain=-1, max_snr=False, comments=None):
        Observation.__init__(self, name, target, start, duration, 'TRK_SOL', 0.0, 0.0, frequency1, frequency2, filter, gain=gain, max_snr=max_snr, comments=comments)
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the observation is 
        pointed.  None if the observation mode is either TBN or TBW."""
        
        return ephem.Sun()


class Jovian(DRX):
    """Sub-class of DRX specifically for Jovian DRX observations.   It features a
    reduced number of parameters needed to setup the observation.
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * observation duration (HH:MM:SS.SSS string or timedelta instance)
     * observation tuning frequency 1 (Hz)
     * observation tuning frequency 1 (Hz)
     * integer filter code
    
    Optional Keywords:
     * max_snr - specifies if maximum signal-to-noise beam forming is to be used
                 (default = False)
     * comments - comments about the observation
    """
    
    def __init__(self, name, target, start, duration, frequency1, frequency2, filter, gain=-1, max_snr=False, comments=None):
        Observation.__init__(self, name, target, start, duration, 'TRK_JOV', 0.0, 0.0, frequency1, frequency2, filter, gain=gain, max_snr=max_snr, comments=comments)
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the observation is 
        pointed.  None if the observation mode is either TBN or TBW."""
        
        return ephem.Jupiter()


class Stepped(Observation):
    """Sub-class of Observation for dealing with STEPPED-mode observations.  It 
    features a reduced number of parameters needed to setup the observation and added
    support for the individual steps.
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * integer filter code
    
    Optional Keywords:
     * steps - array of BeamStep objects that specify the different steps
     * comments - comments about the observation
    """
    
    def __init__(self, name, target, start, filter, steps=None, is_radec=True, gain=-1, comments=None):
        self.is_radec = bool(is_radec)
        self.steps = _TypedParentList(BeamStep, self)
        if steps is not None:
            if isinstance(steps, (tuple, list)):
                self.steps.extend(steps)
            else:
                self.steps.append(steps)
        self.filter_codes = DRXFilters
        Observation.__init__(self, name, target, start, 'please_dont_warn_me', 'STEPPED', 0.0, 0.0, 0.0, 0.0, filter, gain=gain, max_snr=False, comments=comments)
        
    def update(self):
        """Update the computed parameters from the string values."""
        
        self.beam = self.get_beam_type()
        
        disabledBeamDipole = False
        duration = 0
        for step in self.steps:
            step.update()
            duration += step.dur
            
            ## Disable beam-dipole mode for STEPPED-mode observations that 
            ## use custom delays and gains
            if step.delays is not None and step.gains is not None:
                if not disabledBeamDipole:
                    self.set_beamdipole_mode(0)
                    disabledBeamDipole = True
                    
        self.dur = duration
        self.dataVolume = self.estimate_bytes()
        
    @property
    def duration(self):
        """Parse the list of BeamStep objects to get the total observation 
        duration as the number of seconds in that period."""
        
        duration = 0
        for step in self.steps:
            duration += step.dur
            
        s, ms = duration//1000, (duration%1000)/1000.0
        h = s // 3600
        m = (s // 60) % 60
        s = s % 60
    
        return "%i:%02i:%06.3f" % (h, m, s+ms)
        
    @duration.setter
    def duration(self, value):
        if value != 'please_dont_warn_me':
            warnings.warn("The duration of a STEPPED observation can only be changed by adjusting the step durations", RuntimeWarning)
            
    def append(self, newStep):
        """Add a new BeamStep step to the list of steps."""
        
        self.steps.append(newStep)
        self.update()
        
    def set_beamdipole_mode(self, stand, beam_gain=0.04, dipole_gain=1.0, pol='X'):
        """Convert the current observation to a 'beam-dipole mode' 
        observation with the specified stand.  Setting the stand to zero
        will disable the 'beam-dipole mode' for this observation'.
        
        Keywords:
         * beam_gain - BAM gain to use for each dipole in the beam
                      default: 0.04; range: 0.0 to 1.0
         * dipole_gain - BAM gain to use for the single dipole
                        default: 1.0; range: 0.0 to 1.0
         * pol - Polarization to record  default: "X"
        """
        
        # Validate
        if stand < 0 or stand > LWA_MAX_NSTD:
            raise ValueError("Stand number %i is out of range: 0 <= stand <= %i" % (stand, LWA_MAX_NSTD))
        if beam_gain < 0.0 or beam_gain > 1.0:
            raise ValueError("Beam BAM gain is out of range: 0.0 <= beam_gain <= 1.0")
        if dipole_gain < 0.0 or dipole_gain > 1.0:
            raise ValueError("Dipole BAM gain is out of range: 0.0 <= dipole_gain <= 1.0")
        if pol.upper() not in ('X', 'Y'):
            raise ValueError("Unknown polarization.  Valid values are 'X' and 'Y'")
        
        # Go
        if stand == 0:
            ## Disable beam-dipole mode
            self.beamDipole = None
        else:
            ## Stand -> DP Stand
            station = lwa1
            if self._parent is not None:
                station = self._parent.station
                
            for ant in station.antennas:
                if ant.stand.id == stand:
                    dpStand = (ant.digitizer+1)/2
                    
            self.beamDipole = [dpStand, beam_gain, dipole_gain, pol.upper()]
            
    def estimate_bytes(self):
        """Estimate the data volume for the specified type and duration of 
        observations.  For DRX:
        
            bytes = duration * sample_rate / 4096 * 4128 bytes * 2 tunings * 2 pols.
        """
        
        dur = 0
        for step in self.steps:
            dur += step.dur
            
        try:
            sample_rate = self.filter_codes[self.filter]
        except KeyError:
            sample_rate = 0.0
            
        data_rate = DRXSize * 4 * sample_rate / 4096
        if self._parent is not None:
            nchan, nwin = self._parent.spcSetup
            nprod = 4 if self._parent.spcMetatag in ('{Stokes=IQUV}',) else 2
            SPCSize = 76 + nchan*2*nprod*4
            try:
                data_rate = SPCSize * sample_rate / (nchan*nwin)
            except ZeroDivisionError:
                pass
                
        return dur/1000.0 * data_rate
        
    @property
    def target_visibility(self):
        """Return the fractional visibility of the target during the observation 
        period."""
        
        station = lwa1
        if self._parent is not None:
            station = self._parent.station
            
        pnt = self.fixed_body
        
        vis = 0
        cnt = 0
        relStart = 0
        max_alt = 0.0
        for step in self.steps:
            if step.is_radec:
                pnt = step.fixed_body
                
                dt = 0.0
                while dt <= self.dur/1000.0:
                    station.date = self.mjd + (relStart/1000.0 + self.mpm/1000.0 + dt)/3600/24.0 + MJD_OFFSET - DJD_OFFSET
                    pnt.compute(station)
                    max_alt = max([max_alt, pnt.alt])
                    
                    cnt += 1
                    if pnt.alt > 0:
                        vis += 1
                        
                    dt += 300.0
            else:
                max_alt = max([max_alt, step.c2])
                
                cnt += 1
                if step.c2 > 0:
                    vis += 1
            
            relStart += step.dur
            
        if max_alt < 20*math.pi/180:
            #warnings.warn("Maximum altitude for this observation is %.1f degrees" % (max_alt*180/math.pi))
            pass
            
        return float(vis)/float(cnt)
        
    def validate(self, verbose=False):
        """Evaluate the observation and return True if it is valid, False
        otherwise."""
        
        self.update()
        
        failures = 0
        # Basic - filter setup
        if self.filter not in [1, 2, 3, 4, 5, 6, 7]:
            if verbose:
                print("[%i] Error: Invalid filter code '%i'" % (os.getpid(), self.filter))
            failures += 1
            
        # Basic - steps
        stepCount = 1
        for step in self.steps:
            if verbose:
                print("[%i] Validating step %i" % (os.getpid(), stepCount))
            if not step.validate(verbose=verbose):
                failures += 1
            if step.is_radec != self.is_radec:
                if verbose:
                    print("[%i] Error: Step is not of the same coordinate type as observation" % os.getpid())
                failures += 1
                
            stepCount += 1
            
        # Advanced - Target Visibility
        if self.target_visibility < 1.0:
            if verbose:
                print("[%i] Error: Target steps only above the horizon for %.1f%% of the observation" % (os.getpid(), self.target_visibility*100.0))
            failures += 1
            
        # Advanced - Data Volume
        if self.dataVolume >= (_DRSUCapacityTB*1024**4):
            if verbose:
                print("[%i] Error: Data volume exceeds %i TB DRSU limit" % (os.getpid(), _DRSUCapacityTB))
            failures += 1
            
        # Advanced - ASP
        failures += not self._validate_asp(verbose=verbose)
        
        # Any failures indicates a bad observation
        if failures == 0:
            return True
        else:
            return False


class BeamStep(object):
    """Class for holding all of the information (pointing center, tuning frequencies, 
    etc.)associated with a particular step.  
    
    Required Keywords:
     * pointing coordinate 1 (RA [hours] or azimuth [degrees] or ephem.hours/ephem.degrees 
       instance)
     * pointing coordinate 2 (dec or elevation/altitude [degrees] or ephem.degrees instance)
     * observation duration (HH:MM:SS.SSS string or timedelta instance)
     * observation tuning frequency 1 (Hz)
     * observation tuning frequency 1 (Hz)
    
    Optional Keywords:
     * is_radec - whether the coordinates are in RA/Dec or Az/El pairs (default=RA/Dec)
     * max_snr - specifies if maximum signal-to-noise beam forming is to be used
                 (default = False)
     * spec_delays - 520 list of delays to apply for each antenna
     * spec_gains - 260 by 2 by 2 list of gains ([[XY, XY], [YX, YY]]) to apply for each antenna
    
    .. note::
        If `spec_delays` is specified, `spec_gains` must also be specified.
        Specifying both `spec_delays` and `spec_gains` overrides the `max_snr` keyword.
    
    .. versionchanged:: 1.0.0
        Added support for azimuth/altitude and RA/dec values as ephem.hours/ephem.degrees 
        instances
    """
    
    def __init__(self, c1, c2, duration, frequency1, frequency2, is_radec=True, max_snr=False, spec_delays=None, spec_gains=None):
        self.is_radec = bool(is_radec)
        if self.is_radec:
            convFactor = 12.0/math.pi
        else:
            convFactor = 180.0/math.pi
        self.c1 = float(c1) * (convFactor if type(c1).__name__ == 'Angle' else 1.0)
        self.c2 = float(c2) * (180.0/math.pi if type(c2).__name__ == 'Angle' else 1.0)
        if isinstance(duration, timedelta):
            # Make sure the number of microseconds agree with milliseconds
            us = int(round(duration.microseconds/1000.0))*1000
            duration = timedelta(days=duration.days, seconds=duration.seconds, microseconds=us)
        self.duration = duration
        self.frequency1 = float(frequency1)
        self.frequency2 = float(frequency2)
        self.max_snr = bool(max_snr)
        self.delays = spec_delays
        self.gains = spec_gains
        
        self.beam = None
        
        self.update()
        
    def __str__(self):
        c1s = "RA" if self.is_radec else "Az"
        c2s = "Dec" if self.is_radec else "Alt"
        return "Step of %s %.3f, %s %.3f for %s at %.3f and %.3f MHz" % (c1s, self.c1, c2s, self.c2, self.duration, self.frequency1/1e6, self.frequency2/1e6)
        
    @property
    def c1(self):
        """Coordinate 1 - hours (J2000) if RA, degrees if azimuth."""
        
        return self._c1
        
    @c1.setter
    def c1(self, value):
        if isinstance(value, ephem.Angle):
            if self.is_radec:
                value = value * 12.0/math.pi
            else:
                value = value * 180.0/math.pi
        elif isinstance(value, AstroAngle):
            if self.is_radec:
                value = value.to('hourangle').value
            else:
                value = value.to('deg').value
        elif isinstance(value, str):
            value = AstroAngle(value)
            if self.is_radec:
                value = value.to('hourangle').value
            else:
                value = value.to('deg').value
        if self.is_radec:
            if value < 0.0 or value >=24.0:
                raise ValueError("Invalid value for RA '%.6f' hr" % value)
        else:
            if value < 0.0 or value >= 360.0:
                raise ValueError("Invalid value for azimuth '%.6f' deg" % value)
        self._c1 = value
        
    @property
    def c2(self):
        """Coordinate 2 - degrees (J2000) if dec., degrees if elevation."""
        
        return self._c2
        
    @c2.setter
    def c2(self, value):
        if isinstance(value, ephem.Angle):
            value = value * 180.0/math.pi
        elif isinstance(value, AstroAngle):
            value = value.to('deg').value
        elif isinstance(value, str):
            value = AstroAngle(value).to('deg').value
        if self.is_radec:
            if value < -90.0 or value > 90.0:
                raise ValueError("Invalid value for dec. '%.6f' deg" % value)
        else:
            if value < 0.0 or value > 90.0:
                raise ValueError("Invalid value for elevation '%.6f' deg" % value)
        self._c2 = value
        
    @property
    def duration(self):
        """Duration in seconds."""
        
        s, ms = self.dur//1000, (self.dur%1000)/1000.0
        h = s // 3600
        m = (s // 60) % 60
        s = s % 60
    
        return "%i:%02i:%06.3f" % (h, m, s+ms)
        
    @duration.setter
    def duration(self, value):
        if isinstance(value, str):
            fields = value.split(':')
            s = float(fields.pop())
            try:
                m = int(fields.pop(), 10)
            except IndexError:
                m = 0
            try:
                h = int(fields.pop(), 10)
            except IndexError:
                h = 0
            seconds = h*3600 + m*60 + int(s)
            ms = int(round((s - int(s))*1000))
            if ms >= 1000:
                seconds += 1
                ms -= 1000
            seconds = seconds + ms/1000.0
            
        elif isinstance(value, timedelta):
            seconds = value.days*86400 + value.seconds
            ms = int(round(value.microseconds/1000.0))/1000.0
            seconds = seconds + ms
            
        elif isinstance(value, astrounits.quantity.Quantity):
            seconds = seconds.to('s').value
            
        else:
            seconds = value
        self.dur = int(round(seconds*1000))
        
    @property
    def frequency1(self):
        """Tuning 1 frequency in Hz."""
        
        return word_to_freq(self.freq1)
        
    @frequency1.setter
    def frequency1(self, value):
        if isinstance(value, astrounits.quantity.Quantity):
            value = value.to('Hz').value
        self.freq1 = freq_to_word(float(value))
        
    @property
    def frequency2(self):
        """Tuning 2 frequency in Hz."""
        
        return word_to_freq(self.freq2)
        
    @frequency2.setter
    def frequency2(self, value):
        if isinstance(value, astrounits.quantity.Quantity):
            value = value.to('Hz').value
        self.freq2 = freq_to_word(float(value))
        
    def update(self):
        """
        Update the settings.
        """
        
        self.beam = self.get_beam_type()
        
    def get_beam_type(self):
        """Return a valid value for beam type based on whether maximum S/N beam 
        forming has been requested."""
        
        if self.delays is not None and self.gains is not None:
            return 'SPEC_DELAYS_GAINS'
        else:
            if self.max_snr:
                return 'MAX_SNR'
            else:
                return 'SIMPLE'
                
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the observation is 
        pointed.  None if the observation mode is either TBN or TBW."""
        
        if self.is_radec:
            pnt = ephem.FixedBody()
            pnt._ra = self.c1 / 12.0 * math.pi
            pnt._dec = self.c2 / 180.0 * math.pi
            pnt._epoch = ephem.J2000
            
        else:
            pnt = None
            
        return pnt
            
    def validate(self, verbose=False):
        """Evaluate the step and return True if it is valid, False otherwise."""
        
        self.update()
        
        station = lwa1
        if self._parent is not None:
            if self._parent._parent is not None:
                station = self._parent._parent.station
        mandc = station.interface.get_module('mcs')
        backend = station.interface.get_module('backend')
        
        failures = 0
        # Basic - Delay and gain settings are correctly configured
        if self.delays is not None:
            if len(self.delays) != 2*mandc.LWA_MAX_NSTD:
                failures += 1
                if verbose:
                    print("[%i] Error: Specified delay list had the wrong number of antennas" % os.getpid())
            if self.gains is None:
                failures += 1
                if verbose:
                    print("[%i] Error: Delays specified but gains were not" % os.getpid())
        if self.gains is not None:
            if len(self.gains) != mandc.LWA_MAX_NSTD:
                failures += 1
                if verbose:
                    print("[%i] Error: Specified gain list had the wrong number of stands" % os.getpid())
            for g,gain in enumerate(self.gains):
                if len(gain) != 2:
                    failures += 1
                    if verbose:
                        print("[%i] Error: Expected a 2x2 matrix of gain values for stand %i" % (os.getpid(), g))
                else:
                    if len(gain[0]) != 2 or len(gain[1]) != 2:
                        failures += 1
                        if verbose:
                            print("[%i] Error: Expected a 2x2 matrix of gain values for stand %i" % (os.getpid(), g))
            if self.delays is None:
                failures += 1
                if verbose:
                    print("[%i] Error: Gains specified but delays were not" % os.getpid())
        # Basic - Observation time
        if self.dur < 5:
            if verbose:
                print("[%i] Error: step dwell time (%i ms) is too short" % (os.getpid(), self.dur))
            failures += 1
        # Basic - Frequency and filter code values
        if self.freq1 < backend.DRX_TUNING_WORD_MIN or self.freq1 > backend.DRX_TUNING_WORD_MAX:
            if verbose:
                print("[%i] Error: Specified frequency for tuning 1 is outside of DP tuning range" % os.getpid())
            failures += 1
        if (self.freq2 < backend.DRX_TUNING_WORD_MIN or self.freq2 > backend.DRX_TUNING_WORD_MAX) and self.freq2 != 0:
            if verbose:
                print("[%i] Error: Specified frequency for tuning 2 is outside of DP tuning range" % os.getpid())
            failures += 1
        # Any failures indicates a bad observation
        if failures == 0:
            return True
        else:
            return False


class Session(object):
    """Class to hold all of the observations in a session."""
    
    _allowed_modes = (TBW, TBN, DRX, Stepped)
    
    _parent = None
    
    def __init__(self, name, id, observations=None, data_return_method='DRSU', comments=None, station=lwa1):
        self.name = name
        self.id = int(id)
        self.observations = _TypedParentList(self._allowed_modes, self)
        if observations is not None:
            if isinstance(observations, (tuple, list)):
                self.observations.extend(observations)
            else:
                self.observations.append(observations)
        self.dataReturnMethod = data_return_method
        self.ucf_username = None
        self.comments = comments
        
        self.configuration_authority = 0
        self.drx_beam = -1
        self.spcSetup = [0, 0]
        self.spcMetatag = None
        
        self.recordMIB = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
        self.updateMIB = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
        
        self.include_mcssch_log = False
        self.include_mcsexe_log = False
        
        self.include_station_smib = False
        self.include_station_design = False
        
        self.station = station
        
    @property
    def station(self):
        return self._station
        
    @station.setter
    def station(self, value):
        """
        Update the station used by the project for source computations.
        
        .. versionadded:: 1.2.0
        """
        
        if value.interface.sdf != 'lsl.common.sdf':
            raise ValueError("Incompatible station: expected %s, got %s" % \
                             (value.interface.sdf, 'lsl.common.sdf'))
            
        self._station = value
        self.update()
        
    def append(self, newObservation):
        """Add a new Observation to the list of observations."""
        
        self.observations.append(newObservation)
        
    @property
    def spectrometer_channels(self):
        """Number of spectrometer channesl to output, 0 is disables."""
        
        return self.spcSetup[0]
        
    @spectrometer_channels.setter
    def spectrometer_channels(self, value):
        """Set the number of spectrometer channels to generate, 0 to disable."""
        
        value = int(value)
        if value not in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192):
            raise ValueError("Invalid DR spectrometer channel count '%i'" % value)
        self.spcSetup[0] = value
        
    @property
    def spectrometer_integration(self):
        """Number of FFT windows per spectrometer integration, 0 to disable."""
        
        return self.spcSetup[1]
        
    @spectrometer_integration.setter
    def spectrometer_integration(self, value):
        """Set the number of FFT window per spectrometer integration to use, 
        0 to disable."""
        
        value = int(value)
        if value not in (384, 768, 1536, 3072, 6144, 12288, 24576, 49152, 98304, 196608):
            raise ValueError("Invalid DR spectrometer integration count '%i'" % value)
        self.spcSetup[1] = value
        
    @property
    def spectrometer_metatag(self):
        """Spectrometer polarization selection."""
        
        return self.spcMetatag
        
    @spectrometer_metatag.setter
    def spectrometer_metatag(self, value):
        """Set the spectrometer metatag, '' to disable."""
        
        if value == '':
            value = None
        elif value is not None:
            if value[0] != '{':
                value = '{'+value
            if value[-1] != '}':
                value = value+'}'
        if value not in (None, '', 
                         '{Stokes=XXYY}', '{Stokes=CRCI}', '{Stokes=XXCRCIYY}', 
                         '{Stokes=I}', '{Stokes=IV}', '{Stokes=IQUV}'):
            raise ValueError("Invalid DR spectrometer mode '%s'" % value)
        self.spcMetatag = value
        
    def set_mib_record_interval(self, component, interval):
        """Set the record interval for one of the level-1 subsystems (ASP, DP_, etc.) to
        a particular value in minutes.  A KeyError is raised if an invalid sub-system is
        specified.
        
        Special Values are:
          * -1 = use the MCS default interval
          * 0 = never record the MIB entries (the entries are still updated, however)
        """
        
        if component not in self.recordMIB.keys():
               raise KeyError("Unknown subsystem '%s'" % component)
        self.recordMIB[component] = int(interval)
        
    def set_mib_update_interval(self, component, interval):
        """Set the update interval for one of the level-1 subsystems (ASP, DP_, etc.) to 
        a particular value in minutes.  A KeyError is raised if an invalid sub-system is
        specified.
        
        Special Values are:
         * -1 = use the MCS default interval
         * 0 = request no updates to the MIB entries
        """
        
        if component not in self.updateMIB.keys():
               raise KeyError("Unknown subsystem '%s'" % component)
        self.updateMIB[component] = int(interval)
        
    @property
    def data_return_method(self):
        return self.dataReturnMethod
        
    @data_return_method.setter
    def data_return_method(self, method):
        """Set the data return method for the session.  Valid values are: UCF, DRSU, and 
        'USB Harddrives'."""
        
        if method not in ('UCF', 'DRSU', 'USB Harddrives'):
            raise ValueError("Unknown data return method: %s" % method)
            
        self.dataReturnMethod = method
        
    def update(self):
        """Update the various observations in the session."""
        
        for obs in self.observations:
            obs.update()
            
    def validate(self, verbose=False):
        """Examine all of the observations associated with the session to check
        for validity.  If everything is valid, return True.  Otherwise, return
        False."""
        
        self.update()
        
        backend = self.station.interface.get_module('backend')
        
        failures = 0
        totalData = 0.0
        if self.id < 1 or self.id > 9999:
            if verbose:
                print("[%i] Error: Invalid session ID number '%i'" % (os.getpid(), self.id))
            failures += 1
            
        if self.configuration_authority < 0 or self.configuration_authority > 65535:
            if verbose:
                print("[%i] Error: Invalid configuraton request authority '%i'" % (os.getpid(), self.configuration_authority))
            failures += 1
        if self.drx_beam != -1 and self.drx_beam not in list(range(1, backend.DRX_BEAMS_MAX+1)):
            if verbose:
                print("[%i] Error: Invalid beam number '%i'" % (os.getpid(), self.drx_beam))
            failures += 1
        for key in list(self.recordMIB.keys()):
            if self.recordMIB[key] < -1:
                if verbose:
                    print("[%i] Error: Invalid recording interval for '%s' MIB entry '%i'" % (os.getpid(), key, self.recordMIB[key]))
                failures += 1
            if self.updateMIB[key] < -1:
                if verbose:
                    print("[%i] Error: Invalid update interval for '%s' MIB entry '%i'" % (os.getpid(), key, self.updateMIB[key]))
                failures += 1
                
        if self.spcSetup[0] > 0 or self.spcSetup[1] > 0 or self.spcMetatag not in (None, ''):
            if self.spcSetup[0] not in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192):
                if verbose:
                    print("[%i] Error: Invalid DR spectrometer channel count '%i'" % (os.getpid(), self.spcSetup[0]))
                failures += 1
            if self.spcSetup[1] not in (384, 768, 1536, 3072, 6144, 12288, 24576, 49152, 98304, 196608):
                if verbose:
                    print("[%i] Error: Invalid DR spectrometer integration count '%i'" % (os.getpid(), self.spcSetup[1]))
                failures += 1
            if self.spcMetatag not in (None, '', '{Stokes=XXYY}', '{Stokes=IQUV}', '{Stokes=IV}'):
                if verbose:
                    print("[%i] Error: Invalid DR spectrometer mode '%s'" % (os.getpid(), self.spcMetatag))
                failures += 1
            if len(self.observations) > 0:
                if self.observations[0].mode in ('TBW', 'TBN'):
                    if verbose:
                        print("[%i] Error: DR spectrometer incompatible with '%s'" % (os.getpid(), self.observations[0].mode))
                    failures += 1
                    
        observationCount = 1
        for obs in self.observations:
            if verbose:
                print("[%i] Validating observation %i" % (os.getpid(), observationCount))
            
            if not obs.validate(verbose=verbose):
                failures += 1
            totalData += obs.dataVolume
            
            observationCount += 1

        # Make sure that the observations don't overlap
        sObs = self.observations
        
        for i in range(len(sObs)):
            maxOverlaps = 1
            overlaps = []
            nOverlaps = 0

            for j in range(len(sObs)):
                if verbose and i != j:
                    print("[%i] Checking for overlap between observations %i and %i" % (os.getpid(), i+1, j+1))

                cStart = int(sObs[j].mjd)*24*3600*1000 + int(sObs[j].mpm)
                cStop = cStart + int(sObs[j].dur)
                pStart = int(sObs[i].mjd)*24*3600*1000 + int(sObs[i].mpm)
                pStop = pStart + int(sObs[i].dur)

                if pStart >= cStart and pStart < cStop:
                    nOverlaps += 1
                    
                    if i != j:
                        overlaps.append(j)
            
            if nOverlaps > maxOverlaps:
                if verbose:
                    print("[%i] Error: Observation %i overlaps with %s" % (os.getpid(), i+1, ','.join(["%i" % (j+1) for j in overlaps])))
                failures += 1
            
        if totalData >= (_DRSUCapacityTB*1024**4):
            if verbose:
                print("[%i] Error: Total data volume for session exceeds %i TB DRSU limit" % (os.getpid(), _DRSUCapacityTB,))
            failures += 1
        
        if failures == 0:
            return True
        else:
            return False
            
    def __eq__(self, other):
        if isinstance(other, Session):
            self.observations.sort()
            other.observations.sort()
            
            startSelf = self.observations[0].mjd + self.observations[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.observations[0].mjd + other.observations[0].mpm / (1000.0*3600.0*24.0)
            return startSelf == startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __ne__(self, other):
        if isinstance(other, Session):
            self.observations.sort()
            other.observations.sort()
            
            startSelf = self.observations[0].mjd + self.observations[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.observations[0].mjd + other.observations[0].mpm / (1000.0*3600.0*24.0)
            return startSelf != startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __gt__(self, other):
        if isinstance(other, Session):
            self.observations.sort()
            other.observations.sort()
            
            startSelf = self.observations[0].mjd + self.observations[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.observations[0].mjd + other.observations[0].mpm / (1000.0*3600.0*24.0)
            return startSelf > startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __ge__(self, other):
        if isinstance(other, Session):
            self.observations.sort()
            other.observations.sort()
            
            startSelf = self.observations[0].mjd + self.observations[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.observations[0].mjd + other.observations[0].mpm / (1000.0*3600.0*24.0)
            return startSelf >= startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __lt__(self, other):
        if isinstance(other, Session):
            self.observations.sort()
            other.observations.sort()
            
            startSelf = self.observations[0].mjd + self.observations[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.observations[0].mjd + other.observations[0].mpm / (1000.0*3600.0*24.0)
            return startSelf < startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)
            
    def __le__(self, other):
        if isinstance(other, Session):
            self.observations.sort()
            other.observations.sort()
            
            startSelf = self.observations[0].mjd + self.observations[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.observations[0].mjd + other.observations[0].mpm / (1000.0*3600.0*24.0)
            return startSelf <= startOther
        else:
            raise TypeError("Unsupported type: '%s'" % type(other).__name__)


def _parse_create_obs_object(obs_temp, beam_temps=None, verbose=False):
    """Given a obs_temp dictionary of observation parameters and, optionally, a list of
    beam_temp step parameters, return a complete Observation object corresponding to 
    those values."""
    
    # If the observation ID is 0, do nothing.
    if obs_temp['id'] == 0:
        return None
        
    if beam_temps is None:
        beam_temps = []
        
    # Create a time string for the start time in UTC.  This is a little tricky 
    # because of the rounding to the nearest millisecond which has to be done
    # to the datetime object.
    start = Time(obs_temp['mjd'] + obs_temp['mpm'] / 1000.0 / 3600.0 / 24.0, format='MJD').utc_py_date
    start += timedelta(microseconds=(int(round(start.microsecond/1000.0)*1000.0)-start.microsecond))
    utcString = start.strftime("UTC %Y %m %d %H:%M:%S.%f")
    utcString = utcString[:-3]
    
    # Build up a string representing the observation duration.  For TBW observations 
    # this needs to be wrapped in a try...expect statement to catch errors.
    try:
        dur = obs_temp['duration']
        dur = float(dur) / 1000.0
        durString = '%02i:%02i:%06.3f' % (dur/3600.0, (dur%3600.0)/60.0, dur%60.0)
    except:
        pass
        
    # Convert the frequencies from "tuning words" to Hz
    f1 = word_to_freq(obs_temp['freq1'])
    f2 = word_to_freq(obs_temp['freq2'])
    
    # Get the mode and run through the various cases
    mode = obs_temp['mode']
    if verbose:
        print("[%i] Obs %i is mode %s" % (os.getpid(), obs_temp['id'], mode))
        
    if mode == 'TBW':
        obsOut = TBW(obs_temp['name'], obs_temp['target'], utcString, 12000000, comments=obs_temp['comments'])
        obsOut.bits = 1*obs_temp['tbwBits']
        if obs_temp['tbwSamples'] > 0:
            obsOut.samples = 1*obs_temp['tbwSamples']
        else:
            obsOut.samples = 12000000 if obsOut.bits == 12 else 36000000
    elif mode == 'TBN':
        obsOut = TBN(obs_temp['name'], obs_temp['target'], utcString, durString, f1, obs_temp['filter'], gain=obs_temp['gain'], comments=obs_temp['comments'])
    elif mode == 'TRK_RADEC':
        obsOut = DRX(obs_temp['name'], obs_temp['target'], utcString, durString, obs_temp['ra'], obs_temp['dec'], f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'TRK_SOL':
        obsOut = Solar(obs_temp['name'], obs_temp['target'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'TRK_JOV':
        obsOut = Jovian(obs_temp['name'], obs_temp['target'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'STEPPED':
        if verbose:
            print("[%i] -> found %i steps" % (os.getpid(), len(beam_temps)))
            
        obsOut = Stepped(obs_temp['name'], obs_temp['target'], utcString, obs_temp['filter'], is_radec=obs_temp['stpRADec'], steps=[], gain=obs_temp['gain'], comments=obs_temp['comments'])
        for beam_temp in beam_temps:
            try:
                dur = beam_temp['duration']
                dur = float(dur) / 1000.0
                durString = '%02i:%02i:%06.3f' % (dur/3600.0, (dur%3600.0)/60.0, dur%60.0)
            except:
                pass
            
            f1 = word_to_freq(beam_temp['freq1'])
            f2 = word_to_freq(beam_temp['freq2'])
            
            if beam_temp['delays'] is not None:
                if len(beam_temp['delays']) != 2*LWA_MAX_NSTD:
                    raise RuntimeError("Invalid number of delays for custom beamforming")
            if beam_temp['gains'] is not None:
                if len(beam_temp['gains']) != LWA_MAX_NSTD:
                    raise RuntimeError("Invalid number of gains for custom beamforming")
                    
            obsOut.append( BeamStep(beam_temp['c1'], beam_temp['c2'], durString, f1, f2, obs_temp['stpRADec'], beam_temp['MaxSNR'], beam_temp['delays'], beam_temp['gains']) )
    else:
        raise RuntimeError("Invalid mode encountered: %s" % mode)
        
    # Set the beam-dipole mode information (if applicable)
    if obs_temp['beamDipole'] is not None:
        obsOut.beamDipole = obs_temp['beamDipole']
        
    # Set the ASP/FEE values
    obsOut.fee_power = copy.deepcopy(obs_temp['obsFEE'])
    obsOut.asp_filter = copy.deepcopy(obs_temp['aspFlt'])
    obsOut.asp_atten_1 = copy.deepcopy(obs_temp['aspAT1'])
    obsOut.asp_atten_2 = copy.deepcopy(obs_temp['aspAT2'])
    obsOut.asp_atten_split = copy.deepcopy(obs_temp['aspATS'])
    
    # Force the observation to be updated
    obsOut.update()
    
    # Return the newly created Observation object
    return obsOut


def parse_sdf(filename, verbose=False):
    """
    Given a filename, read the file's contents into the SDM instance and return
    that instance.
    """
    
    # Create the keyword regular expression to deal with various indicies included 
    # in the keywords
    kwdRE = re.compile(r'(?P<keyword>[A-Z_0-9\+]+)(\[(?P<id1>[0-9]+?)\])?(\[(?P<id2>[0-9]+?)\])?(\[(?P<id3>[0-9]+?)\])?(\[(?P<id4>[0-9]+?)\])?')
    
    # Create the metatag regular expression to deal with spectrometer mode settings
    metaRE = re.compile(r'\{.*\}')
    
    # Create empty objects to get things started.  Values will get filled in as they
    # are found in the file
    po = ProjectOffice()
    observer = Observer('observer_name', 0)
    project = Project(observer, 'project_name', 'project_id', project_office=po)
    session = Session('session_name', 0, observations=[])
    project.sessions = [session,]
    project.project_office.sessions = []
    project.project_office.observations = [[],]
    
    obs_temp = {'id': 0, 'name': '', 'target': '', 'ra': 0.0, 'dec': 0.0, 'start': '', 'duration': '', 'mode': '', 
                'beamDipole': None, 'freq1': 0, 'freq2': 0, 'filter': 0, 'MaxSNR': False, 'comments': None, 
                'stpRADec': True, 'tbwBits': 12, 'tbwSamples': 0, 'gain': -1, 
                'obsFEE': [[-1,-1] for n in range(LWA_MAX_NSTD)], 
                'aspFlt': [-1 for n in range(LWA_MAX_NSTD)], 'aspAT1': [-1 for n in range(LWA_MAX_NSTD)], 
                'aspAT2': [-1 for n in range(LWA_MAX_NSTD)], 'aspATS': [-1 for n in range(LWA_MAX_NSTD)]}
    beam_temp = {'id': 0, 'c1': 0.0, 'c2': 0.0, 'duration': 0, 'freq1': 0, 'freq2': 0, 'MaxSNR': False, 'delays': None, 'gains': None}
    beam_temps = []
    
    # Loop over the file
    with open(filename, 'r') as fh:
        for line in fh:
            # Trim off the newline character and skip blank lines
            line = line.replace('\n', '')
            if len(line) == 0 or line.isspace():
                continue
            
            # Split into a keyword, value pair and run it through the regular expression
            # to deal with any indicies present
            try:
                keywordSection, value = line.split(None, 1)
            except:
                continue
            
            mtch = kwdRE.match(keywordSection)
            keyword = mtch.group('keyword')
        
            ids = [-1, -1, -1, -1]
            for i in range(4):
                try:
                    ids[i] = int(mtch.group('id%i' % (i+1)))
                except TypeError:
                    pass
                
            # Skip over the observer comment lines (denoted by a plus sign at the end) 
            # of the keyword
            if keyword[-1] == '+':
                continue
            
            # Observer Info
            if keyword == 'PI_ID':
                project.observer.id = int(value)
                continue
            if keyword == 'PI_NAME':
                project.observer.name = value
                project.observer.split_name()
                continue
            
            # Project/Proposal Info
            if keyword == 'PROJECT_ID':
                project.id = value
                continue
            if keyword == 'PROJECT_TITLE':
                project.name = value
                continue
            if keyword == 'PROJECT_REMPI':
                project.comments = value
                continue
            if keyword == 'PROJECT_REMPO':
                project.project_office.project = value
                continue
            
            # Session Info
            if keyword == 'SESSION_ID':
                project.sessions[0].id = int(value)
                continue
            if keyword == 'SESSION_TITLE':
                project.sessions[0].name = value
                continue
            if keyword == 'SESSION_REMPI':
                mtch = UCF_USERNAME_RE.search(value)
                if mtch is not None:
                    project.sessions[0].ucf_username = mtch.group('username')
                    if mtch.group('subdir') is not None:
                        project.sessions[0].ucf_username = os.path.join(project.sessions[0].ucf_username, mtch.group('subdir'))
                project.sessions[0].comments = value
                continue
            if keyword == 'SESSION_REMPO':
                project.project_office.sessions.append(None)
                parts = value.split(';;', 1)
                first = parts[0]
                try:
                    second = parts[1]
                except IndexError:
                    second = ''
                
                if first[:31] == 'Requested data return method is':
                    # Catch for project office comments that are data return related
                    project.sessions[0].dataReturnMethod = first[32:]
                    project.project_office.sessions[0] = second
                else:
                    # Catch for standard (not data related) project office comments
                    project.project_office.sessions[0] = value
                continue
            if keyword == 'SESSION_CRA':
                project.sessions[0].configuration_authority = int(value)
                continue
            if keyword[0:12] == 'SESSION_MRP_':
                component = keyword[12:]
                project.sessions[0].recordMIB[component] = int(value)
                continue
            if keyword[0:12] == 'SESSION_MUP_':
                component = keyword[12:]
                project.sessions[0].updateMIB[component] = int(value)
                continue
            if keyword == 'SESSION_LOG_SCH':
                project.sessions[0].include_mcssch_log = bool(value)
                continue
            if keyword == 'SESSION_LOG_EXE':
                project.sessions[0].include_mcsexe_log = bool(value)
                continue
            if keyword == 'SESSION_INC_SMIB':
                project.sessions[0].include_station_smib = bool(value)
                continue
            if keyword == 'SESSION_INC_DES':
                project.sessions[0].include_station_design = bool(value)
                continue
            if keyword == 'SESSION_DRX_BEAM':
                project.sessions[0].drx_beam = int(value)
                continue
            if keyword == 'SESSION_SPC':
                # Remove the ' marks
                value = value.replace("'", "")
                # Excise the metatags
                mtch = metaRE.search(value)
                if mtch is not None:
                    metatag = mtch.group(0)
                    value = metaRE.sub('', value)
                else:
                    metatag = None
                    
                project.sessions[0].spcSetup = [int(i) for i in value.lstrip().rstrip().split(None, 1)]
                project.sessions[0].spcMetatag = metatag
                # If the input field is '' the value of spcSetup is [].  This
                # isn't good for the SDF render so reset [] to [0, 0]
                if project.sessions[0].spcSetup == []:
                    project.sessions[0].spcSetup = [0, 0]
                    project.sessions[0].spcMetatag = None
                continue
                
            # Observation Info
            if keyword == 'OBS_ID':
                if obs_temp['id'] != 0:
                    project.sessions[0].observations.append( _parse_create_obs_object(obs_temp, beam_temps=beam_temps, verbose=verbose) )
                    beam_temp = {'id': 0, 'c1': 0.0, 'c2': 0.0, 'duration': 0, 'freq1': 0, 'freq2': 0, 'MaxSNR': False, 'delays': None, 'gains': None}
                    beam_temps = []
                obs_temp['id'] = int(value)
                project.project_office.observations[0].append( None )
            
                if verbose:
                    print("[%i] Started obs %i" % (os.getpid(), int(value)))
                
                continue
            if keyword == 'OBS_TITLE':
                obs_temp['name'] = value
                continue
            if keyword == 'OBS_TARGET':
                obs_temp['target'] = value
                continue
            if keyword == 'OBS_REMPI':
                obs_temp['comments'] = value
                continue
            if keyword == 'OBS_REMPO':
                project.project_office.observations[0][-1] = value
                continue
            if keyword == 'OBS_START_MJD':
                obs_temp['mjd'] = int(value)
                continue
            if keyword == 'OBS_START_MPM':
                obs_temp['mpm'] = int(value)
                continue
            if keyword == 'OBS_DUR':
                obs_temp['duration'] = int(value)
                continue
            if keyword == 'OBS_MODE':
                obs_temp['mode'] = value
                continue
            if keyword == 'OBS_BDM':
                # Remove the ' marks
                value = value.replace("'", "")
                try:
                    stand, beam_gain, dipole_gain, pol = value.lstrip().rstrip().split(None, 3)
                    obs_temp['beamDipole'] = [int(stand), float(beam_gain), float(dipole_gain), pol]
                except ValueError:
                    pass
            if keyword == 'OBS_RA':
                obs_temp['ra'] = float(value)
                continue
            if keyword == 'OBS_DEC':
                obs_temp['dec'] = float(value)
                continue
            if keyword == 'OBS_B':
                if value != 'SIMPLE':
                    obs_temp['MaxSNR'] = True
                continue
            if keyword == 'OBS_FREQ1':
                obs_temp['freq1'] = int(value)
                continue
            if keyword == 'OBS_FREQ2':
                obs_temp['freq2'] = int(value)
                continue
            if keyword == 'OBS_BW':
                obs_temp['filter'] = int(value)
                continue
            if keyword == 'OBS_STP_RADEC':
                obs_temp['stpRADec'] = bool(int(value))
                continue
            
            # Individual Stepped Beam Observations - This is a bit messy because of
            # trying to keep up when a new step is encountered.  This adds in some 
            # overhead to all of the steps.
            if keyword == 'OBS_STP_C1':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['c1'] = float(value)
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['c1'] = float(value)
                continue
            
            if keyword == 'OBS_STP_C2':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['c2'] = float(value)
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['c2'] = float(value)
                continue
            
            if keyword == 'OBS_STP_T':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['duration'] = int(value)
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['duration'] = int(value)
                continue
            
            if keyword == 'OBS_STP_FREQ1':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['freq1'] = int(value)
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['freq1'] = int(value)
                continue
            
            if keyword == 'OBS_STP_FREQ2':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['freq2'] = int(value)
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    beam_temps[-1]['freq2'] = int(value)
                continue
            
            if keyword == 'OBS_STP_B':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                
                    if value in ('MAX_SNR', '2'):
                        beam_temps[-1]['MaxSNR'] = True
                    
                    elif value in ('SPEC_DELAYS_GAINS', '3'):
                        beam_temps[-1]['delays'] = []
                        beam_temps[-1]['gains'] = []
                        for bdi in range(2*LWA_MAX_NSTD):
                            beam_temps[-1]['delays'].append( 0 )
                            if bdi < LWA_MAX_NSTD:
                                beam_temps[-1]['gains'].append( [[0, 0], [0, 0]] )
                            
                    else:
                        beam_temps[-1]['MaxSNR'] = False
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    
                    if value in ('MAX_SNR', '2'):
                        beam_temps[-1]['MaxSNR'] = True
                    
                    elif value in ('SPEC_DELAYS_GAINS', '3'):
                        beam_temps[-1]['delays'] = []
                        beam_temps[-1]['gains'] = []
                        for bdi in range(2*LWA_MAX_NSTD):
                            beam_temps[-1]['delays'].append( 0 )
                            if bdi < LWA_MAX_NSTD:
                                beam_temps[-1]['gains'].append( [[0, 0], [0, 0]] )
                            
                    else:
                        beam_temps[-1]['MaxSNR'] = False
                continue
            
            if keyword == 'OBS_BEAM_DELAY':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    try:
                        beam_temps[-1]['delays'][ids[1]-1] = int(value)
                    except IndexError:
                        raise RuntimeError("Invalid index encountered when parsing OBS_BEAM_DELAY")
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    try:
                        beam_temps[-1]['delays'][ids[1]-1] = int(value)
                    except IndexError:
                        raise RuntimeError("Invalid index encountered when parsing OBS_BEAM_DELAY")
                continue
            
            if keyword == 'OBS_BEAM_GAIN':
                if len(beam_temps) == 0:
                    beam_temps.append( copy.deepcopy(beam_temp) )
                    beam_temps[-1]['id'] = ids[0]
                    try:
                        beam_temps[-1]['gains'][ids[1]-1][ids[2]-1][ids[3]-1] = int(value)
                    except IndexError:
                        pass
                else:
                    if beam_temps[-1]['id'] != ids[0]:
                        beam_temps.append( copy.deepcopy(beam_temps[-1]) )
                        beam_temps[-1]['id'] = ids[0]
                    try:
                        beam_temps[-1]['gains'][ids[1]-1][ids[2]-1][ids[3]-1] = int(value)
                    except IndexError:
                        pass
                continue
            
            # Session wide settings at the end of the observations
            if keyword == 'OBS_FEE':
                if ids[0] == 0:
                    for n in range(len(obs_temp['obsFEE'])):
                        obs_temp['obsFEE'][n][ids[1]-1] = int(value)
                else:
                    obs_temp['obsFEE'][ids[0]-1][ids[1]-1] = int(value)
                continue
            if keyword == 'OBS_ASP_FLT':
                if ids[0] == 0:
                    for n in range(len(obs_temp['aspFlt'])):
                        obs_temp['aspFlt'][n] = int(value)
                else:
                    obs_temp['aspFlt'][ids[0]-1] = int(value)
                continue
            if keyword == 'OBS_ASP_AT1':
                if ids[0] == 0:
                    for n in range(len(obs_temp['aspAT1'])):
                        obs_temp['aspAT1'][n] = int(value)
                else:
                    obs_temp['aspAT1'][ids[0]-1] = int(value)
                continue
            if keyword == 'OBS_ASP_AT2':
                if ids[0] == 0:
                    for n in range(len(obs_temp['aspAT2'])):
                        obs_temp['aspAT2'][n] = int(value)
                else:
                    obs_temp['aspAT2'][ids[0]-1] = int(value)
                continue
            if keyword == 'OBS_ASP_ATS':
                if ids[0] == 0:
                    for n in range(len(obs_temp['aspATS'])):
                        obs_temp['aspATS'][n] = int(value)
                else:
                    obs_temp['aspATS'][ids[0]-1] = int(value)
                continue
            if keyword == 'OBS_TBW_BITS':
                obs_temp['tbwBits'] = int(value)
                continue
            if keyword == 'OBS_TBW_SAMPLES':
                obs_temp['tbwSamples'] = int(value)
                continue
            if keyword == 'OBS_TBN_GAIN':
                obs_temp['gain'] = int(value)
                continue
            if keyword == 'OBS_DRX_GAIN':
                obs_temp['gain'] = int(value)
                continue
            
            # Keywords that might indicate this is for ADP-based stations/actually an IDF
            if keyword in ('OBS_TBF_SAMPLES', 'RUN_ID'):
                raise RuntimeError("Invalid keyword encountered: %s" % keyword)
            
        # Create the final observation
        if obs_temp['id'] != 0:
            project.sessions[0].observations.append( _parse_create_obs_object(obs_temp, beam_temps=beam_temps, verbose=verbose) )
            beam_temps = []
            
    # Return the project
    return project


def get_observation_start_stop(obs):
    """
    Given an observation, get the start and stop times (returned as a two-
    element tuple of UTC datetime instances).
    
    .. versionadded:: 1.0.0
    """
    
    # UNIX timestamp for the start
    tStart = utcjd_to_unix(obs.mjd + MJD_OFFSET)
    tStart += obs.mpm / 1000.0
    
    # UNIX timestamp for the stop
    tStop = tStart +  obs.dur / 1000.0
    
    # Conversion to a timezone-aware datetime instance
    tStart = _UTC.localize( datetime.utcfromtimestamp(tStart) )
    tStop  = _UTC.localize( datetime.utcfromtimestamp(tStop ) )
    
    # Make sure we have an integer number of milliseconds
    ## Start
    us = tStart.microsecond
    us = int(round(us/1000.0))*1000
    tStart = tStart.replace(microsecond=us)
    ## Stop
    us = tStop.microsecond
    us = int(round(us/1000.0))*1000
    tStop = tStop.replace(microsecond=us)
    
    # Return
    return tStart, tStop


def is_valid(filename, verbose=False):
    """
    Given a filename, see if it is valid SDF file or not.
    
    .. versionadded:: 1.2.0
    """
    
    passes = 0
    failures = 0
    try:
        proj = parse_sdf(filename)
        passes += 1
        if verbose:
            print(colorfy("Parser - {{%green OK}}"))
            
        valid = proj.validate()
        if valid:
            passes += 1
            if verbose:
                print(colorfy("Validator - {{%green OK}}"))
        else:
            failures += 1
            if verbose:
                print(colorfy("Validator -{{%red {{%bold FAILED}}}}"))
                
    except IOError as e:
        raise e
    except:
        failures += 1
        if verbose:
            print(colorfy("Parser - {{%red {{%bold FAILED}}}}"))
            
    if verbose:
        print("---")
        print("%i passed / %i failed" % (passes, failures))
        
    return False if failures else True
