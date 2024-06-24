"""
Module that contains all of the relevant class to build up a representation 
of a interferometer definition file.  The hierarchy of classes is:
  * Project - class that holds all of the information about the project (including
              the observer) and one or more runs.  Technically, a ID file has 
              only one run but this approach allows for the generation of 
              multiple SD files from a single Project object.
  * Observer - class that hold the observer's name and numeric ID
  * Run - class that holds all of the scans associated and the associated
          correlator setup for a with a particular interferometer run
  * Scan - class that hold information about a particular scan.  It
           includes a variety of attributes that are used to convert human-
           readable inputs to SDF data values.  The scan class is 
           further subclasses into:
             - DRX - class for general DRX scan, with sub-classes:
               * Solar - class for solar tracking
               * Jovian - class for Jovian tracking
    
Most class contain 'validate' attribute functions that can be used to determine if the 
project/run/scan are valid or not given the constraints of
the ADP system.

In addition to providing the means for creating interferometer definition files from 
scratch, this module also includes a simple parser for ID files.

.. versionchanged:: 2.0.0
    Added support for astropy.time.Time and astropy.coordinates.Angle instances

.. versionadded:: 1.2.4
"""

import os
import re
import copy
import math
import pytz
import ephem
from functools import total_ordering
from datetime import datetime, timedelta

from astropy import units as astrounits
from astropy.coordinates import Angle as AstroAngle

from lsl.transform import Time
from lsl.astro import utcjd_to_unix, MJD_OFFSET, DJD_OFFSET
from lsl.common._sdf_utils import *
from lsl.common.color import colorfy

from lsl.common.mcsADP import datetime_to_mjdmpm, mjdmpm_to_datetime
from lsl.common.adp import freq_to_word, word_to_freq
from lsl.common.stations import LWAStation, get_full_stations, lwa1
from lsl.reader.drx import FILTER_CODES as DRXFilters
from lsl.reader.drx import FRAME_SIZE as DRXSize
from lsl.common.sdf import UCF_USERNAME_RE, Observer
from lsl.common import sdf, sdfADP

from lsl.config import LSL_CONFIG
OBSV_CONFIG = LSL_CONFIG.view('observing')

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['UCF_USERNAME_RE', 'Observer', 'ProjectOffice', 'Project', 'Run', 'Scan', 'DRX', 'Solar', 'Jovian', 'parse_idf',  'get_scan_start_stop', 'is_valid', '__version__']


_UTC = pytz.utc
_DRSUCapacityTB = 10

_MAX_ALT_PHASE_CENTERS = 10


class ProjectOffice(object):
    """Class to hold comments from the LWA object office.  This class isn't really 
    needed to create ID files, but it is helpful for parsing ID files."""
    
    def __init__(self, project=None, runs=None, scans=None):
        self.project = project
        if runs is None:
            self.runs = []
        else:
            self.runs = runs
        if scans is None:
            self.scans = []
        else:
            self.scans = scans


class Project(object):
    """
    Class to hold all the information about a specific interferometer run for a 
    project/proposal.
    """
    
    def __init__(self, observer, name, id, runs=None, comments=None, project_office=None):
        if not isinstance(observer, Observer):
            raise TypeError("Expected 'observer' to be an Observer")
        self.observer = observer
        self.name = name
        self.id = id
        self.comments = comments
        self.runs = sdf._TypedParentList(Run, self)
        if runs is not None:
            if isinstance(runs, (list, tuple)):
                self.runs.extend(runs)
            else:
                self.runs.append(runs)
        if project_office is None:
            self.project_office = ProjectOffice()
        else:
            if not isinstance(project_office, ProjectOffice):
                raise TypeError("Expected 'project_office' to be a ProjectOffice")
            self.project_office = project_office
            
    def __str__(self):
        return f"{self.id}: {self.name} with {len(self.runs)} run(s) for {str(self.observer)}"
        
    @classmethod
    def autofilled(cls, runs=None, comments=None, project_office=None):
        observer = Observer.autofilled()
        name = OBSV_CONFIG.get('project_name')
        id = OBSV_CONFIG.get('project_id')
        if name is None or id is None:
            raise RuntimeError("Auto-fill values for the project cannot be loaded from the configuration file")
        return cls(observer, name, id, runs=runs, comments=comments, project_office=project_office)
        
    def update(self):
        """Update the various runs that are part of this project."""
        
        for ses in self.runs:
            ses.update()
            
    def validate(self, verbose=False):
        """Examine all of the runs and all of their scans to check
        for validity.  If everything is valid, return True.  Otherwise, return
        False."""
        
        self.update()
        
        failures = 0
        runCount = 1
        if len(self.id) > 8:
            if verbose:
                pid_print("Project ID is too long")
            failures += 1
            
        for run in self.runs:
            if verbose:
                pid_print(f"Validating run {runCount}")
            if not run.validate(verbose=verbose):
                failures += 1
                
            runCount += 1
            
        if failures == 0:
            return True
        else:
            return False
            
    def append(self, newRun):
        """Add a new run to the list of runs."""
        
        self.runs.append(newRun)
        
    def render(self, run=0, verbose=False):
        """Create a run definition file that corresponds to the specified 
        run.  Returns the ID file's contents as a string."""
        
        if not self.validate(verbose=verbose) :
            raise RuntimeError("Invalid run/scan parameters.  Aborting.")
        if run >= len(self.runs):
            raise IndexError("Invalid run index")
        
        self.runs[run].update()
        self.runs[run].scans.sort()
            
        ses = self.runs[run]
        try:
            # Try to pull out the project office comments about the run
            pos = self.project_office.runs[run]
        except (TypeError, IndexError):
            pos = None
        try:
            # Try to pull out the project office comments about the scans
            poo = self.project_office.scans[run]
        except (TypeError, IndexError):
            poo = []
        # Enforce that the number of project office scan comments match the
        # actual number of scans
        while (len(ses.scans) - len(poo)) > 0:
            poo.append(None)
            
        # Combine the run comments together in an intelligent fashion
        ## Observer comments
        if ses.ucf_username is not None:
            clean = ''
            if ses.comments:
                clean = UCF_USERNAME_RE.sub('', ses.comments)
            ses.comments = f"ucfuser:{ses.ucf_username}"
            if len(clean) > 0:
                ses.comments += ';;%s' % clean
        ## Project office comments, including the data return method
        if pos != 'None' and pos is not None:
            pos = f"Requested data return method is {ses.dataReturnMethod};;{pos}"
            
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
        
        ## Run Information
        output += "RUN_ID           %s\n" % (ses.id,)
        output += "RUN_TITLE        %s\n" % ('None provided' if ses.name is None else ses.name,)
        output += "RUN_STATIONS     %s\n" % (','.join([station.id for station in ses.stations]),)
        output += "RUN_CHANNELS     %i\n" % (ses.correlator_channels,)
        output += "RUN_INTTIME      %.3f\n" % (ses.correlator_inttime,)
        output += "RUN_BASIS        %s\n" % (ses.correlator_basis,)
        output += "RUN_REMPI        %s\n" % (ses.comments[:4090] if ses.comments else 'None provided',)
        output += "RUN_REMPO        %s\n" % (f"Requested data return method is {ses.dataReturnMethod}" if pos == 'None' or pos is None else pos[:4090],)
        output += "\n"
                    
        ## Scans
        for i,obs in enumerate(ses.scans):
            obsID = i + 1
            
            output += "SCAN_ID          %i\n" % (obsID,)
            output += "SCAN_TARGET      %s\n" % (obs.target,)
            output += "SCAN_INTENT      %s\n" % (obs.intent,)
            output += "SCAN_REMPI       %s\n" % (obs.comments[:4090] if obs.comments else 'None provided',)
            output += "SCAN_REMPO       %s\n" % ("Estimated raw data volume for this scan is %s per station; %s total" % (render_file_size(obs.dataVolumeStation), render_file_size(obs.dataVolume)) if poo[i] == 'None' or poo[i] is None else poo[i],)
            output += "SCAN_START_MJD   %i\n" % (obs.mjd,)
            output += "SCAN_START_MPM   %i\n" % (obs.mpm,)
            output += "SCAN_START       %s\n" % (obs.start.strftime("%Z %Y/%m/%d %H:%M:%S") if isinstance(obs.start, datetime) else obs.start,)
            output += "SCAN_DUR         %i\n" % (obs.dur,)
            output += "SCAN_DUR+        %s\n" % (obs.duration,)
            output += "SCAN_MODE        %s\n" % (obs.mode,)
            if obs.mode == 'TRK_RADEC':
                output += "SCAN_RA          %.9f\n" % (obs.ra,)
                output += "SCAN_DEC         %+.9f\n" % (obs.dec,)
                if obs.pm[0] != 0.0 or obs.pm[1] != 0.0:
                    output += "SCAN_PM_RA       %+.1f\n" % (obs.pm[0],)
                    output += "SCAN_PM_DEC      %+.1f\n" % (obs.pm[1],)
            output += "SCAN_FREQ1       %i\n" % (obs.freq1,)
            output += "SCAN_FREQ1+      %.9f MHz\n" % (obs.frequency1/1e6,)
            output += "SCAN_FREQ2       %i\n" % (obs.freq2,)
            output += "SCAN_FREQ2+      %.9f MHz\n" % (obs.frequency2/1e6,)
            output += "SCAN_BW          %i\n" % (obs.filter,)
            output += "SCAN_BW+         %s\n" % (render_bandwidth(obs.filter, obs.FILTER_CODES),)
            ## Alternate phase centers
            if len(obs.alt_phase_centers) > 0:
                output += "SCAN_ALT_N             %i\n" % (len(obs.alt_phase_centers),)
                for j,phase_center in enumerate(obs.alt_phase_centers):
                    output += "SCAN_ALT_TARGET[%i]    %s\n" % (j+1, phase_center.target)  
                    output += "SCAN_ALT_INTENT[%i]    %s\n" % (j+1, phase_center.intent) 
                    output += "SCAN_ALT_RA[%i]        %.9f\n" % (j+1, phase_center.ra)  
                    output += "SCAN_ALT_DEC[%i]       %+.9f\n" % (j+1, phase_center.dec)
                    if phase_center.pm[0] != 0.0 or phase_center.pm[1] != 0.0:
                        output += "SCAN_ALT_PM_RA[%i]       %+.1f\n" % (j+1, phase_center.pm[0])
                        output += "SCAN_ALT_PM_DEC[%i]      %+.1f\n" % (j+1, phase_center.pm[1])
                        
            ## ASP filter setting
            if obs.asp_filter != -1:
                output += "SCAN_ASP_FLT     %i\n" % (obs.asp_filter,)
            ## DRX gain
            if obs.gain != -1:
                output += "SCAN_DRX_GAIN    %i\n" % (obs.gain,)
            output += "\n"
            
        return output
        
    def writeto(self, filename, run=0, verbose=False, overwrite=False):
        """Create a run definition file that corresponds to the specified 
        run and write it to the provided filename."""
        
        if os.path.exists(filename) and not overwrite:
            raise RuntimeError(f"'{filename}' already exists")
            
        output = self.render(run=run, verbose=verbose)
        with open(filename, 'w') as fh:
            fh.write(output)
            
    def generate_sdfs(self, starting_session_id=1, run=0, verbose=False):
        """Convert the ID file into a collection of `lsl.common.sdfADP.Project` instances
        that can be used to write SD files."""
        
        if not self.validate(verbose=verbose) :
            raise RuntimeError("Invalid run/scan parameters.  Aborting.")
        if run >= len(self.runs):
            raise IndexError("Invalid run index")
        
        self.runs[run].update()
        self.runs[run].scans.sort()   
        ses = self.runs[run]
        
        # Build the SDFs
        ## Setup the common information
        start = mjdmpm_to_datetime(ses.scans[0].mjd, ses.scans[0].mpm)
        new_observer = sdf.Observer(copy.deepcopy(self.observer.name), copy.deepcopy(self.observer.id))
        ## Go!
        sdfs = []
        for i,station in enumerate(ses.stations):
            ### Session
            sdfmod = sdf if station.interface.sdf == 'lsl.common.sdf' else sdfADP
            
            session = sdfmod.Session("%s - %s (%i of %i)" % (ses.name, station.id, i+1, len(ses.stations)), 
                                     starting_session_id, observations=[], station=station)
            session.drx_beam = 1
            session.ucf_username = 'eLWA/%s_%s_%s_%04i' % (self.id, start.strftime('%y%m%d'), start.strftime('%H%M'), ses.id)
            session.data_return_method = 'UCF'
            session.include_station_smib = True
            
            ## Project Office
            new_projoff = sdfmod.ProjectOffice(project=copy.deepcopy(self.project_office.project), 
                                               sessions=copy.deepcopy([self.project_office.runs[run],]), 
                                               observations=copy.deepcopy([self.project_office.scans[run],]))
            
            ## Observations
            for o,obs in enumerate(ses.scans):
                obs_start = mjdmpm_to_datetime(obs.mjd, obs.mpm)
                if isinstance(obs, DRX):
                    ### Apply the proper motion to generate the SDFs
                    delta_epoch = (obs.mjd + obs.mpm/1000.0/86400.0 - 51544.5) / 365.25
                    ra = obs.ra + delta_epoch * obs.pm[0]/math.cos(obs.dec*math.pi/180)/1000.0/3600.0/15.0
                    dec = obs.dec + delta_epoch * obs.pm[1]/1000.0/3600.0
                    comments = ''
                    if comments is not None:
                        comments += obs.comments
                    if obs.pm[0] != 0.0 or obs.pm[1] != 0.0:
                        comments = comments+";;Applied proper motion of %+.1f mas/yr in RA and %+.1f mas/yr in dec" % (obs.pm[0], obs.pm[1])
                        
                    new_obs = sdfmod.DRX(obs.intent, obs.target, _UTC.localize(obs_start), obs.duration, 
                                         ra, dec, 
                                         obs.frequency1, obs.frequency2, obs.filter, 
                                         gain=obs.gain, max_snr=False, comments=comments)
                elif isinstance(obs, Solar):
                    new_obs = sdfmod.Solar(obs.intent, obs.target, _UTC.localize(obs_start), obs.duration, 
                                           obs.frequency1, obs.frequency2, obs.filter, 
                                           gain=obs.gain, max_snr=False, comments=obs.comments)
                elif isinstance(obs, Jovian):
                    new_obs = sdfmod.Jovian(obs.intent, obs.target, _UTC.localize(obs_start), obs.duration, 
                                            obs.frequency1, obs.frequency2, obs.filter, 
                                            gain=obs.gain, max_snr=False, comments=obs.comments)
                else:
                    raise RuntimeError("This should never happen")
                session.append(new_obs)
                
                npc = len(obs.alt_phase_centers)
                for p,phase_center in enumerate(reversed(obs.alt_phase_centers)):
                    cid = npc - p
                    alt_t, alt_i, alt_r, alt_d = phase_center.target, phase_center.intent, phase_center.ra, phase_center.dec
                    alt_pr, alt_pd = phase_center.pm
                    
                    ### Apply the proper motion to generate the SDFs
                    delta_epoch = (obs.mjd + obs.mpm/1000.0/86400.0 - 51544.5) / 365.25
                    alt_r = alt_r + delta_epoch * alt_pr/math.cos(alt_d*math.pi/180)/1000.0/3600.0/15.0
                    alt_d = alt_d + delta_epoch * alt_pd/1000.0/3600.0
                    
                    try:
                        new_projoff.observations[0][o]
                    except IndexError:
                        new_projoff.observations[0][o] = None
                    if new_projoff.observations[0][o] is None:
                        new_projoff.observations[0][o] = ''
                    new_projoff.observations[0][o] = "alttarget%i:%s;;altintent%i:%s;;altra%i:%.9f;;altdec%i:%+.9f;;%s" % (cid, alt_t, cid, alt_i, cid, alt_r, cid, alt_d, new_projoff.observations[0][o])
                    
            ## Project
            project = sdfmod.Project(new_observer, "%s - %s (%i of %i)" % (self.name, station.id, i+1, len(ses.stations)), 
                                     copy.deepcopy(self.id), sessions=[session,], comments=copy.deepcopy(self.comments), 
                                     project_office=new_projoff)
            if project.project_office.sessions[0] is None:
               project.project_office.sessions[0] = ''
            project.project_office.sessions[0] = "corrchannels:%i;;corrinttime:%.3f;;corrbasis:%s;;origuser:%s;;origreturn:%s;;%s" % (ses.correlator_channels, ses.correlator_inttime, ses.correlator_basis, ses.ucf_username, ses.data_return_method.lower(), project.project_office.sessions[0])
            
            ## Save an increment the session ID
            sdfs.append(project)
            starting_session_id += 1
            
        # Done
        return sdfs


@total_ordering
class Run(object):
    """Class to hold all of the scans in an interferometer run."""
    
    def __init__(self, name, id, scans=None, data_return_method='DRSU', comments=None, correlator_channels=512, correlator_inttime=1.0, correlator_basis='linear', stations=get_full_stations()):
        self.name = name
        self.id = int(id)
        self.scans = sdf._TypedParentList(Scan, self)
        if scans is not None:
            if isinstance(scans, (tuple, list)):
                self.scans.extend(scans)
            else:
                self.scans.append(scans)
        self.data_return_method = data_return_method
        self.ucf_username = None
        self.comments = comments
        
        self.correlator_inttime = correlator_inttime
        self.correlator_channels = correlator_channels
        self.correlator_basis = correlator_basis
        self.stations = sdf._TypedParentList(LWAStation, None, stations)
        
    def __str__(self):
        return f"{self.id}: {self.name} with {len(self.scans)} scans and correlator setup:\n  channels: {self.correlator_channels}\n  int. time: {self.correlator_inttime}\n  basis: {self.correlator_basis}\n  stations: {' '.join([s.id for s in self.stations])}"
        
    @property
    def stations(self):
        """List of LWA stations to use in the interferometer."""
        
        return self._stations
        
    @stations.setter
    def stations(self, value):
        """
        Update the stations used by the project for source computations.
        """
        
        value = sdf._TypedParentList(LWAStation, None, value)
        if len(value) < 2:
            raise ValueError("Need at least two stations to form an interferometer")
        self._stations = value
        self.update()
        
    def append(self, newScan):
        """Add a new Scan to the list of scans."""
        
        self.scans.append(newScan)
        
    @property
    def correlator_channels(self):
        """Number of correlator channels to use."""
        
        return self.corr_channels
        
    @correlator_channels.setter
    def correlator_channels(self, value):
        """Set the number of correlator channels to generate."""
        
        value = int(value)
        if value < 16 or value > 32768 or value % 2:
            raise ValueError(f"Invalid correlator channel count '{value}'")
        self.corr_channels = value
        
    @property
    def correlator_inttime(self):
        """Correlator integration time in seconds."""
        
        return self.corr_inttime
        
    @correlator_inttime.setter
    def correlator_inttime(self, value):
        """Set the number of correlatro integration time to use."""
        
        value = float(value)
        if value < 0.1 or value > 10.0:
            raise ValueError(f"Invalid integration time '{value:.3f}'")
        self.corr_inttime = value
        
    @property
    def correlator_basis(self):
        """Correlator output polarization basis."""
        return self.corr_basis
        
    @correlator_basis.setter
    def correlator_basis(self, value):
        """Set the correlator output polarization basis."""
        
        if value.lower() not in ('linear', 'circular', 'stokes'):
            raise ValueError(f"Unknown polarization basis: {value}")
        self.corr_basis = value
        
    @property
    def data_return_method(self):
        return self.dataReturnMethod
        
    @data_return_method.setter
    def data_return_method(self, method):
        """Set the data return method for the run.  Valid values are: UCF, DRSU, and 
        'USB Harddrives'."""
        
        if method not in ('UCF', 'DRSU', 'USB Harddrives'):
            raise ValueError(f"Unknown data return method: {method}")
            
        self.dataReturnMethod = method
        
    def update(self):
        """Update the various scans in the run."""
        
        for obs in self.scans:
            obs.update()
            
    def validate(self, verbose=False):
        """Examine all of the scans associated with the run to check
        for validity.  If everything is valid, return True.  Otherwise, return
        False."""
        
        self.update()
        
        failures = 0
        totalData = 0.0
        if self.id < 1 or self.id > 9999:
            if verbose:
                pid_print(f"Error: Invalid run ID number '{self.id}'")
            failures += 1
            
        station_count = {}
        for station in self.stations:
            try:
                station_count[station.id] += 1
            except KeyError:
                station_count[station.id] = 1
        for station in station_count:
            if station_count[station] != 1:
                if verbose:
                    pid_print(f"Error: Station '{station}' is included {station_count[station]} times")
                failures += 1
                
        scanCount = 1
        for obs in self.scans:
            if verbose:
                pid_print(f"Validating scan {scanCount}")
                
            if not obs.validate(verbose=verbose):
                failures += 1
            totalData += obs.dataVolume
                
            if scanCount > 1:
                if obs.filter != self.scans[scanCount-2].filter:
                    if verbose:
                        pid_print(f"Error: Filter code changes at scan {scanCount}")
                    failures += 1
                    
            scanCount += 1
            
        # Make sure that the scans don't overlap
        sObs = self.scans
        
        for i in range(len(sObs)):
            maxOverlaps = 1
            overlaps = []
            nOverlaps = 0

            for j in range(len(sObs)):
                if verbose and i != j:
                    pid_print(f"Checking for overlap between scans {i+1} and {j+1}")

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
                    pid_print(f"Error: Scan {i+1} overlaps with "+(','.join(["%i" % (j+1) for j in overlaps])))
                failures += 1
            
        if totalData >= (len(self.stations)*_DRSUCapacityTB*1024**4):
            if verbose:
                pid_print(f"Error: Total data volume for run exceeds per-station {_DRSUCapacityTB} TB DRSU limit")
            failures += 1
        
        if failures == 0:
            return True
        else:
            return False
            
    def __eq__(self, other):
        if isinstance(other, Run):
            self.scans.sort()
            other.scans.sort()
            
            startSelf = self.scans[0].mjd + self.scans[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.scans[0].mjd + other.scans[0].mpm / (1000.0*3600.0*24.0)
            return startSelf == startOther
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __lt__(self, other):
        if isinstance(other, Run):
            self.scans.sort()
            other.scans.sort()
            
            startSelf = self.scans[0].mjd + self.scans[0].mpm / (1000.0*3600.0*24.0)
            startOther = other.scans[0].mjd + other.scans[0].mpm / (1000.0*3600.0*24.0)
            return startSelf < startOther
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")


@total_ordering
class Scan(object):
    """
    Class to hold the specifics of a scans.  It currently
    handles TRK_RADEC, TRK_SOL, and TRK_JOV.
    """
    
    _parent = None
    
    id = 1
    FILTER_CODES = DRXFilters

    def __init__(self, target, intent, start, duration, mode, ra, dec, frequency1, frequency2, filter, gain=-1, pm=[0.0, 0.0], comments=None):
        self.target = target
        self.intent = intent
        self.ra = ra
        self.dec = dec
        self.pm = pm
        self.start = start
        self.duration = duration
        self.mode = mode
        self.beamDipole = None
        self.frequency1 = float(frequency1)
        self.frequency2 = float(frequency2)
        self.filter = int(filter)
        self.comments = comments
        
        self.dataVolumeStation = None
        self.dataVolume = None
        
        self.asp_filter = -1
        
        self.gain = int(gain)
        
        self.alt_phase_centers = sdf._TypedParentList(AlternatePhaseCenter, self)
        
        self.update()
        
    def __str__(self):
        """Return a nice string to describe the scan."""
        
        return "%s Scan of '%s':\n Intent %s\n Start %s\n Duration %s\n Filter: %i\n Frequency: %.3f; %.3f\n RA: %.6f hr @ %+.3f mas/yr\n Dec. %.6f d @ +%.3f mas/yr\n" % (self.mode, self.target, self.intent, self.start, self.duration, self.filter, self.frequency1, self.frequency2, self.ra, self.pm[0], self.dec, self.pm[1])
        
    def update(self):
        """Update the computed parameters from the string values."""
        
        stations = [lwa1,]
        if self._parent is not None:
            stations = self._parent.stations
            
        self.dataVolumeStation = self.estimate_bytes()
        self.dataVolume = self.dataVolumeStation*len(stations)
        
        # Update the associated alternate phase centers
        for phase_center in self.alt_phase_centers:
            phase_center.update()
            
    @property
    def intent(self):
        """Obsevational intent."""
        
        return self._intent
        
    @intent.setter
    def intent(self, value):
        value = value.lower()
        if value not in ('fluxcal', 'phasecal', 'target', 'dummy'):
            raise ValueError(f"Invalid scan intent '{value}'")
        self._intent = value
        
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
            raise ValueError(f"Invalid value for RA '{value:.6f}' hr")
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
            raise ValueError(f"Invalid value for dec. '{value:+.6f}' deg")
        self._dec = value
        
    @property
    def pm(self):
        """Target proper motion in mas/yr."""
        
        return self._pm
        
    @pm.setter
    def pm(self, value):
        if value is None:
            self._pm = [0.0, 0.0]
        elif isinstance(value, (tuple, list)):
            if len(value) == 2:
                self._pm = [value[0], value[1]]
            else:
                raise ValueError("Expected a two-element tuple of list of proper motion values in mas/yr")
        else:
            raise ValueError("Expected a two-element tuple of list of proper motion values in mas/yr")
            
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
        
    def add_alt_phase_center(self, target_or_apc, intent=None, ra=None, dec=None, pm=None):
        """Add an alternate phase center to the scan."""
        
        if isinstance(target_or_apc, AlternatePhaseCenter):
            apc = target_or_apc
        else:
            apc = AlternatePhaseCenter(target_or_apc, intent, ra, dec, pm=pm)
        self.alt_phase_centers.append(apc)
        
    def estimate_bytes(self):
        """Estimate the data volume for the specified type and duration of 
        scans.  For DRX:
        
            bytes = duration * sample_rate / 4096 * 4128 bytes * 2 tunings * 2 pols.
        """
        
        try:
            nFrames = self.dur/1000.0 * self.FILTER_CODES[self.filter] / 4096
        except KeyError:
            nFrames = 0
        nBytes = nFrames * DRXSize * 4
        return nBytes
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the scan is 
        pointed.  None if the scan mode is TBN."""
        
        pnt = ephem.FixedBody()
        pnt.name = self.target
        pnt._ra = self.ra / 12.0 * math.pi
        pnt._dec = self.dec / 180.0 * math.pi
        pnt._pmra = self.pm[0]
        pnt._pmdec = self.pm[1]
        pnt._epoch = ephem.J2000
        return pnt
        
    @property
    def target_visibility(self):
        """Return the fractional visibility of the target during the scan 
        period."""
        
        stations = [lwa1,]
        if self._parent is not None:
            stations = self._parent.stations
            
        vis_list = []
        max_alt = 0.0
        for station in stations:
            lwa = station.get_observer()
            pnt = self.fixed_body
            
            vis = 0
            cnt = 0
            dt = 0.0
            while dt <= self.dur/1000.0:
                lwa.date = self.mjd + (self.mpm/1000.0 + dt)/3600/24.0 + MJD_OFFSET - DJD_OFFSET
                pnt.compute(lwa)
                max_alt = max([max_alt, pnt.alt])
                
                cnt += 1
                if pnt.alt > 0:
                    vis += 1
                    
                dt += 300.0
                
            vis_list.append(float(vis)/float(cnt))
            
        if max_alt < 20*math.pi/180:
            #warnings.warn("Maximum altitude for this scan is %.1f degrees" % (max_alt*180/math.pi))
            pass
            
        return min(vis_list)
        
    def validate(self, verbose=False):
        """Evaluate the scan and return True if it is valid, False otherwise."""
        
        self.update()
        
        stations = [lwa1,]
        if self._parent is not None:
            stations = self._parent.stations
        tuning_min = max([station.interface.get_module('backend').DRX_TUNING_WORD_MIN for station in stations])
        tuning_max = min([station.interface.get_module('backend').DRX_TUNING_WORD_MAX for station in stations])
        
        failures = 0
        # Basic - Intent, duration, frequency, and filter code values
        if self.dur < 1:
            if verbose:
                pid_print("Error: Specified a duration of length zero")
            failures += 1
        if self.freq1 < tuning_min or self.freq1 > tuning_max:
            if verbose:
                pid_print("Error: Specified frequency for tuning 1 is outside of LWA tuning range")
            failures += 1
        if (self.freq2 < tuning_min or self.freq2 > tuning_max) and self.freq2 != 0:
            if verbose:
                pid_print("Error: Specified frequency for tuning 2 is outside of LWA tuning range")
            failures += 1
        if self.filter not in [1, 2, 3, 4, 5, 6, 7]:
            if verbose:
                pid_print(f"Error: Invalid filter code '{self.filter}'")
            failures += 1
            
        # Advanced - Target Visibility
        if self.target_visibility < 1.0:
            if verbose:
                pid_print(f"Error: Target is only above the horizon for {self.target_visibility*100.0:.1f}% of the scan")
            failures += 1
            
        # Advanced - alternate phase centers
        if len(self.alt_phase_centers) > _MAX_ALT_PHASE_CENTERS:
            if verbose:
                pid_print("Error: too many alternate phase centers defined")
            failures += 1
        for j,phase_center in enumerate(self.alt_phase_centers):
            if not phase_center.validate(verbose=verbose):
                if verbose:
                    pid_print(f"Error: invalid alternate phase center {j+1}")
                failures += 1
                
            ## Closeness to pointing center
            pnt = self.fixed_body
            alt_pnt = phase_center.fixed_body
            
            lwa = stations[0].get_observer()
            lwa.date = self.mjd + (self.mpm/1000.0 + self.dur/1000.0/2.0)/3600/24.0 + MJD_OFFSET - DJD_OFFSET
            pnt.compute(lwa)
            alt_pnt.compute(lwa)
            
            alt_sep = ephem.separation(pnt, alt_pnt) * 180/math.pi
            
            beam = 2.0*74e6/max([self.frequency1, self.frequency2])
            if alt_sep > beam/2.0:
                if verbose:
                    pid_print(f"Error: alternate phase center {j+1} is {alt_sep:.1f} degrees from pointing center")
                failures += 1
                
        # Advanced - Data Volume
        if self.dataVolumeStation >= (_DRSUCapacityTB*1024**4):
            if verbose:
                pid_print(f"Error: Data volume exceeds {_DRSUCapacityTB} TB DRSU limit")
            failures += 1
            
        # Any failures indicates a bad scan
        if failures == 0:
            return True
        else:
            return False
            
    def __eq__(self, other):
        if isinstance(other, Scan):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf == startOther
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")
            
    def __lt__(self, other):
        if isinstance(other, Scan):
            startSelf = self.mjd + self.mpm / (1000.0*3600.0*24.0)
            startOther = other.mjd + other.mpm / (1000.0*3600.0*24.0)
            return startSelf < startOther
        else:
            raise TypeError(f"Unsupported type: '{type(other).__name__}'")


class DRX(Scan):
    """
    Required Arguments:
     * scan target
     * scan intent
     * scan start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * scan duration (HH:MM:SS.SSS string or timedelta instance)
     * scan RA in hours, J2000.0 or ephem.hours instance
     * scan Dec in degrees, J2000.0 or ephem.hours instance
     * scan tuning frequency 1 (Hz)
     * scan tuning frequency 1 (Hz)
     * integer filter code
    
    Optional Keywords:
     * comments - comments about the scan
    """
    
    alt_phase_centers = []
    
    def __init__(self, target, intent, start, duration, ra, dec, frequency1, frequency2, filter, gain=-1, pm=[0.0, 0.0], comments=None):
        Scan.__init__(self, target, intent, start, duration, 'TRK_RADEC', ra, dec, frequency1, frequency2, filter, gain=gain, pm=pm, comments=comments)


class Solar(Scan):
    """Sub-class of DRX specifically for Solar DRX scans.   It features a
    reduced number of parameters needed to setup the scan.
    
    Required Arguments:
     * scan target
     * scan intent
     * scan start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * scan duration (HH:MM:SS.SSS string or timedelta instance)
     * scan tuning frequency 1 (Hz)
     * scan tuning frequency 1 (Hz)
     * integer filter code
    
    Optional Keywords:
     * comments - comments about the scan
    """
    
    def __init__(self, target, intent, start, duration, frequency1, frequency2, filter, gain=-1, comments=None):
        Scan.__init__(self, target, intent, start, duration, 'TRK_SOL', 0.0, 0.0, frequency1, frequency2, filter, gain=gain, comments=comments)
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the scan is 
        pointed.  None if the scan mode is TBN."""
        
        return ephem.Sun()


class Jovian(Scan):
    """Sub-class of DRX specifically for Jovian DRX scans.   It features a
    reduced number of parameters needed to setup the scan.
    
    Required Arguments:
     * scan target
     * scan intent
     * scan start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string or timezone-
       aware datetime instance)
     * scan duration (HH:MM:SS.SSS string or timedelta instance)
     * scan tuning frequency 1 (Hz)
     * scan tuning frequency 1 (Hz)
     * integer filter code
    
    Optional Keywords:
     * comments - comments about the scan
    """
    
    def __init__(self, target, intent, start, duration, frequency1, frequency2, filter, gain=-1, comments=None):
        Scan.__init__(self, target, intent, start, duration, 'TRK_JOV', 0.0, 0.0, frequency1, frequency2, filter, gain=gain, comments=comments)
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the scan is 
        pointed.  None if the scan mode is TBN."""
        
        return ephem.Jupiter()


class AlternatePhaseCenter(object):
    """Class to hold an alternate phase center for a scan."""
    
    def __init__(self, target, intent, ra, dec, pm=[0.0, 0.0]):
        self.target = target
        self.intent = intent
        self.ra = ra
        self.dec = dec
        self.pm = pm
        
    @property
    def intent(self):
        """Observational intent."""
        
        return self._intent
        
    @intent.setter
    def intent(self, value):
        value = value.lower()
        if value not in ('fluxcal', 'phasecal', 'target', 'dummy'):
            raise ValueError(f"Invalid scan intent '{value}'")
        self._intent = value
        
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
            raise ValueError(f"Invalid value for RA '{value:.6f}' hr")
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
            raise ValueError(f"Invalid value for dec. '{value:+.6f}' deg")
        self._dec = value
        
    @property
    def pm(self):
        """Target proper motion in mas/yr."""
        
        return self._pm
        
    @pm.setter
    def pm(self, value):
        if value is None:
            self._pm = [0.0, 0.0]
        elif isinstance(value, (tuple, list)):
            if len(value) == 2:
                self._pm = [value[0], value[1]]
            else:
                raise ValueError("Expected a two-element tuple of list of proper motion values in mas/yr")
        else:
            raise ValueError("Expected a two-element tuple of list of proper motion values in mas/yr")
            
    def update(self):
        """Update the computed parameters from the string values."""
        
        pass
        
    @property
    def fixed_body(self):
        """Return an ephem.Body object corresponding to where the scan is 
        pointed.  None if the scan mode is TBN."""
        
        pnt = ephem.FixedBody()
        pnt.name = self.target
        pnt._ra = self.ra / 12.0 * math.pi
        pnt._dec = self.dec / 180.0 * math.pi
        pnt._pmra = self.pm[0]
        pnt._pmdec = self.pm[1]
        pnt._epoch = ephem.J2000
        return pnt
        
    @property
    def target_visibility(self):
        """Return the fractional visibility of the target during the scan 
        period."""
        
        mjd, mpm, dur = 0, 0, 0
        stations = [lwa1,]
        if self._parent is not None:
            mjd = self._parent.mjd
            mpm = self._parent.mpm
            dur = self._parent.dur
            if self._parent._parent is not None:
                stations = self._parent._parent.stations
                
        vis_list = []
        max_alt = 0.0
        for station in stations:
            lwa = station.get_observer()
            pnt = self.fixed_body
            
            vis = 0
            cnt = 0
            dt = 0.0
            while dt <= dur/1000.0:
                lwa.date = mjd + (mpm/1000.0 + dt)/3600/24.0 + MJD_OFFSET - DJD_OFFSET
                pnt.compute(lwa)
                max_alt = max([max_alt, pnt.alt])
                
                cnt += 1
                if pnt.alt > 0:
                    vis += 1
                    
                dt += 300.0
                
            vis_list.append(float(vis)/float(cnt))
            
        if max_alt < 20*math.pi/180:
            #warnings.warn("Maximum altitude for this scan is %.1f degrees" % (max_alt*180/math.pi))
            pass
            
        return min(vis_list)
        
    def validate(self, verbose=False):
        """Basic validation of the pointing, that's it."""
        
        failures = 0
        
        ## Advanced - Target Visibility
        if self.target_visibility < 1.0:
            if verbose:
                pid_print(f"Error: Target is only above the horizon for {self.target_visibility*100.0:.1f}% of the scan")
            failures += 1
            
        # Any failures indicates a bad alternate phase center
        if failures == 0:
            return True
        else:
            return False


def _parse_create_scan_object(obs_temp, alt_temps=[], verbose=False):
    """Given a obs_temp dictionary of scan parameters, return a complete Scan object 
    corresponding to those values."""
    
    # If the scan ID is 0, do nothing.
    if obs_temp['id'] == 0:
        return None
        
    # Create a time string for the start time in UTC.  This is a little tricky 
    # because of the rounding to the nearest millisecond which has to be done
    # to the datetime object.
    start = Time(obs_temp['mjd'] + obs_temp['mpm'] / 1000.0 / 3600.0 / 24.0, format='MJD').utc_py_date
    start += timedelta(microseconds=(int(round(start.microsecond/1000.0)*1000.0)-start.microsecond))
    utcString = start.strftime("UTC %Y %m %d %H:%M:%S.%f")
    utcString = utcString[:-3]
    
    # Build up a string representing the scan duration.
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
        pid_print(f"Scan {obs_temp['id']} is mode {mode}")
        
    if mode == 'TRK_RADEC':
        obsOut = DRX(obs_temp['target'], obs_temp['intent'], utcString, durString, obs_temp['ra'], obs_temp['dec'], f1, f2, obs_temp['filter'], gain=obs_temp['gain'], pm=obs_temp['pm'], comments=obs_temp['comments'])
    elif mode == 'TRK_SOL':
        obsOut = Solar(obs_temp['target'], obs_temp['intent'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], comments=obs_temp['comments'])
    elif mode == 'TRK_JOV':
        obsOut = Jovian(obs_temp['target'], obs_temp['intent'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], comments=obs_temp['comments'])
    else:
        raise RuntimeError(f"Invalid mode encountered: {mode}")
        
    # Add in the alternate phase centers
    if obs_temp['nAlt'] != len(alt_temps):
        raise RuntimeError("Mis-match between SCAN_ALT_N and the number of alternate phase centers")
    for alt_temp in alt_temps:
        obsOut.add_alt_phase_center(alt_temp['target'], alt_temp['intent'], alt_temp['ra'], alt_temp['dec'], pm=alt_temp['pm'])
        
    # Set the ASP/FEE values
    obsOut.asp_filter = copy.deepcopy(obs_temp['aspFlt'])
    
    # Force the scan to be updated
    obsOut.update()
    
    # Return the newly created Scan object
    return obsOut


def parse_idf(filename, verbose=False):
    """
    Given a filename, read the file's contents into the IDF instance and return
    that instance.
    """
    
    # Create the keyword regular expression to deal with various indicies included 
    # in the keywords
    kwdRE = re.compile(r'(?P<keyword>[A-Z_0-9\+]+)(\[(?P<id1>[0-9]+?)\])?(\[(?P<id2>[0-9]+?)\])?(\[(?P<id3>[0-9]+?)\])?(\[(?P<id4>[0-9]+?)\])?')
    
    # Create empty objects to get things started.  Values will get filled in as they
    # are found in the file
    po = ProjectOffice()
    observer = Observer('observer_name', 0)
    project = Project(observer, 'project_name', 'project_id', project_office=po)
    run = Run('run_name', 0, scans=[])
    project.runs = [run,]
    project.project_office.runs = []
    project.project_office.scans = [[],]
    
    obs_temp = {'id': 0, 'target': '', 'intent': '', 'ra': 0.0, 'dec': 0.0, 'pm':[0.0, 0.0], 'start': '', 'duration': '', 'mode': '', 
                'freq1': 0, 'freq2': 0, 'filter': 0, 'comments': None, 'nAlt': 0, 'gain': -1, 
                'aspFlt': -1}
    alt_temp = {'id': 0, 'target': '', 'intent': '', 'ra': 0.0, 'dec': 0.0, 'pm':[0.0, 0.0]}
    alt_temps = []
    
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
            except ValueError:
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
            
            # Run Info
            if keyword == 'RUN_ID':
                project.runs[0].id = int(value)
                continue
            if keyword == 'RUN_TITLE':
                project.runs[0].name = value
                continue
            if keyword == 'RUN_STATIONS':
                use_stations = []
                possible = get_full_stations()
                for field in value.split(','):
                    field = field.strip().rstrip()
                    if field.lower() == 'all':
                        use_stations = copy.deepcopy(possible)
                        break
                    else:
                        for station in possible:
                            if station.id == field:
                                use_stations.append(station)
                                break
                project.runs[0].stations  = use_stations
            if keyword == 'RUN_CHANNELS':
                project.runs[0].correlator_channels = int(value)
                continue
            if keyword == 'RUN_INTTIME':
                project.runs[0].correlator_inttime = float(value)
                continue
            if keyword == 'RUN_BASIS':
                project.runs[0].correlator_basis = value
                continue
            if keyword == 'RUN_REMPI':
                mtch = UCF_USERNAME_RE.search(value)
                if mtch is not None:
                    project.runs[0].ucf_username = mtch.group('username')
                    if mtch.group('subdir') is not None:
                        project.runs[0].ucf_username = os.path.join(project.runs[0].ucf_username, mtch.group('subdir'))
                project.runs[0].comments = UCF_USERNAME_RE.sub('', value)
                continue
            if keyword == 'RUN_REMPO':
                project.project_office.runs.append(None)
                parts = value.split(';;', 1)
                first = parts[0]
                try:
                    second = parts[1]
                except IndexError:
                    second = ''
                
                if first[:31] == 'Requested data return method is':
                    # Catch for project office comments that are data return related
                    project.runs[0].dataReturnMethod = first[32:]
                    project.project_office.runs[0] = second
                else:
                    # Catch for standard (not data related) project office comments
                    project.project_office.runs[0] = value
                continue
            
            # Scan Info
            if keyword == 'SCAN_ID':
                if obs_temp['id'] != 0:
                    project.runs[0].scans.append( _parse_create_scan_object(obs_temp, alt_temps=alt_temps, verbose=verbose) )
                    alt_temp = {'id': 0, 'target': '', 'intent': '', 'ra': 0.0, 'dec': 0.0, 'pm':[0.0, 0.0]}
                    alt_temps = []
                obs_temp['id'] = int(value)
                project.project_office.scans[0].append( None )
            
                if verbose:
                    pid_print(f"Started scan {value}")
                
                continue
            if keyword == 'SCAN_TARGET':
                obs_temp['target'] = value
                continue
            if keyword == 'SCAN_INTENT':
                obs_temp['intent'] = value
                continue
            if keyword == 'SCAN_REMPI':
                obs_temp['comments'] = value
                continue
            if keyword == 'SCAN_REMPO':
                project.project_office.scans[0][-1] = value
                continue
            if keyword == 'SCAN_START_MJD':
                obs_temp['mjd'] = int(value)
                continue
            if keyword == 'SCAN_START_MPM':
                obs_temp['mpm'] = int(value)
                continue
            if keyword == 'SCAN_DUR':
                obs_temp['duration'] = int(value)
                continue
            if keyword == 'SCAN_MODE':
                obs_temp['mode'] = value
                continue
            if keyword == 'SCAN_RA':
                obs_temp['ra'] = float(value)
                continue
            if keyword == 'SCAN_DEC':
                obs_temp['dec'] = float(value)
                continue
            if keyword == 'SCAN_PM_RA':
                obs_temp['pm'][0] = float(value)
                continue
            if keyword == 'SCAN_PM_DEC':
                obs_temp['pm'][1] = float(value)
                continue
            if keyword == 'SCAN_FREQ1':
                obs_temp['freq1'] = int(value)
                continue
            if keyword == 'SCAN_FREQ2':
                obs_temp['freq2'] = int(value)
                continue
            if keyword == 'SCAN_BW':
                obs_temp['filter'] = int(value)
                continue
            
            # Alternate phase centers
            if keyword == 'SCAN_ALT_N':
                obs_temp['nAlt'] = int(value)
                continue
            if keyword == 'SCAN_ALT_TARGET':
                if len(alt_temps) == 0:
                    alt_temps.append( copy.deepcopy(alt_temp) )
                    alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['target'] = value
                else:
                    if alt_temps[-1]['id'] != ids[0]:
                        alt_temps.append( copy.deepcopy(alt_temps[-1]) )
                        alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['target'] = value
                continue
            if keyword == 'SCAN_ALT_INTENT':
                if len(alt_temps) == 0:
                    alt_temps.append( copy.deepcopy(alt_temp) )
                    alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['intent'] = value
                else:
                    if alt_temps[-1]['id'] != ids[0]:
                        alt_temps.append( copy.deepcopy(alt_temp) )
                        alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['intent'] = value
                continue
            if keyword == 'SCAN_ALT_RA':
                if len(alt_temps) == 0:
                    alt_temps.append( copy.deepcopy(alt_temp) )
                    alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['ra'] = float(value)
                else:
                    if alt_temps[-1]['id'] != ids[0]:
                        alt_temps.append( copy.deepcopy(alt_temp) )
                        alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['ra'] = float(value)
                continue
            if keyword == 'SCAN_ALT_DEC':
                if len(alt_temps) == 0:
                    alt_temps.append( copy.deepcopy(alt_temp) )
                    alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['dec'] = float(value)
                else:
                    if alt_temps[-1]['id'] != ids[0]:
                        alt_temps.append( copy.deepcopy(alt_temp) )
                        alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['dec'] = float(value)
                continue
            if keyword == 'SCAN_ALT_PM_RA':
                if len(alt_temps) == 0:
                    alt_temps.append( copy.deepcopy(alt_temp) )
                    alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['pm'][0] = float(value)
                else:
                    if alt_temps[-1]['id'] != ids[0]:
                        alt_temps.append( copy.deepcopy(alt_temp) )
                        alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['pm'][0] = float(value)
                continue
            if keyword == 'SCAN_ALT_PM_DEC':
                if len(alt_temps) == 0:
                    alt_temps.append( copy.deepcopy(alt_temp) )
                    alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['pm'][1] = float(value)
                else:
                    if alt_temps[-1]['id'] != ids[0]:
                        alt_temps.append( copy.deepcopy(alt_temp) )
                        alt_temps[-1]['id'] = ids[0]
                    alt_temps[-1]['pm'][1] = float(value)
                continue
            # Run wide settings at the end of the scans
            if keyword == 'SCAN_ASP_FLT':
                obs_temp['aspFlt'] = int(value)
                continue
            if keyword == 'SCAN_DRX_GAIN':
                obs_temp['gain'] = int(value)
                continue
            
            # Keywords that might indicate this is a SDF
            if keyword in ('SESSION_ID', 'SESSION_DRX_BEAM'):
                raise RuntimeError(f"Invalid keyword encountered: {keyword}")
            
        # Create the final scan
        if obs_temp['id'] != 0:
            project.runs[0].scans.append( _parse_create_scan_object(obs_temp, alt_temps=alt_temps, verbose=verbose) )
            
    # Return the project
    return project


def get_scan_start_stop(obs):
    """
    Given a scan, get the start and stop times (returned as a two-
    element tuple of UTC datetime instances).
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
    Given a filename, see if it is valid IDF file or not.
    """
    
    passes = 0
    failures = 0
    try:
        proj = parse_idf(filename)
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
                print(colorfy("Validator - {{%red {{%bold FAILED}}}}"))
                
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
