"""
Module that contains all of the relevant class to build up a representation 
of a session definition file as defined in MCS0030v5 and updated for LWA-SV.  The 
hierarchy of classes is:
  * Project - class that holds all of the information about the project (including
              the observer) and one or more sessions.  Technically, a SD file has 
              only one session but this approach allows for the generation of 
              multiple SD files from a single Project object.
  * Observer - class that hold the observer's name and numeric ID
  * Session - class that holds all of the observations associated with a particular 
              ADP output.  
  * Observations - class that hold information about a particular observation.  It
                   includes a variety of attributes that are used to convert human-
                   readable inputs to SDF data values.  The observation class is 
                   further subclasses into:
                     - TBF - class for TBF observations
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
the ADP system.

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

import os
import re
import copy
import math
import pytz
from datetime import datetime, timedelta

from lsl.transform import Time
from lsl.astro import utcjd_to_unix, MJD_OFFSET
from lsl.common.color import colorfy

from lsl.common.mcsADP import LWA_MAX_NSTD
from lsl.common.adp import word_to_freq, fC
from lsl.common.stations import lwasv
from lsl.reader.drx import FILTER_CODES as DRXFilters
from lsl.reader.tbf import FRAME_SIZE as TBFSize

from lsl.common._sdf_utils import *
from lsl.common.sdf import Observer, ProjectOffice
from lsl.common.sdf import Project as _Project, Session as _Session
from lsl.common.sdf import UCF_USERNAME_RE, Observation, TBN, DRX, Solar, Jovian, Lunar, Stepped, BeamStep

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '1.2'
__all__ = ['UCF_USERNAME_RE', 'parse_time', 'Observer', 'ProjectOffice', 'Project', 'Session', 'Observation', 'TBF', 'TBN', 'DRX', 'Solar', 'Jovian', 'Lunar', 'Stepped', 'BeamStep', 'parse_sdf',  'get_observation_start_stop', 'is_valid']


_UTC = pytz.utc
_DRSUCapacityTB = 10
# Factors for computing the time it takes to read out a TBF from the number 
# of samples
_TBF_TIME_SCALE = 196000
_TBF_TIME_GAIN = 150


class Project(_Project):
    """
    Class to hold all the information about a specific session for a 
    project/proposal.
    
    .. versionchanged:: 1.2.1
        Added a new writeto() method to directly write the SDF to a file.
    """
    
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
        except (TypeError, IndexError):
            pos = None
        try:
            # Try to pull out the project office comments about the observations
            poo = self.project_office.observations[session]
        except (TypeError, IndexError):
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
        
        ## Session Information
        output += "SESSION_ID       %s\n" % (ses.id,)
        output += "SESSION_TITLE    %s\n" % ('None provided' if ses.name is None else ses.name,)
        output += "SESSION_REMPI    %s\n" % (ses.comments[:4090] if ses.comments else 'None provided',)
        output += "SESSION_REMPO    %s\n" % (f"Requested data return method is {ses.dataReturnMethod}" if pos == 'None' or pos is None else pos[:4090],)
        if ses.configuration_authority != 0:
            output += "SESSION_CRA      %i\n" % (ses.configuration_authority,)
        if ses.drx_beam != -1:
            output += "SESSION_DRX_BEAM %i\n" % (ses.drx_beam,)
        if ses.spcSetup[0] != 0 and ses.spcSetup[1] != 0:
            output += "SESSION_SPC      %i %i%s\n" % (ses.spcSetup[0], ses.spcSetup[1], '' if ses.spcMetatag is None else ses.spcMetatag)
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
            output += "OBS_REMPO        %s\n" % ("Estimated data volume for this observation is %s" % render_file_size(obs.dataVolume) if poo[i] == 'None' or poo[i] is None else poo[i],)
            output += "OBS_START_MJD    %i\n" % (obs.mjd,)
            output += "OBS_START_MPM    %i\n" % (obs.mpm,)
            output += "OBS_START        %s\n" % (obs.start.strftime("%Z %Y/%m/%d %H:%M:%S") if isinstance(obs.start, datetime) else obs.start,)
            output += "OBS_DUR          %i\n" % (obs.dur,)
            output += "OBS_DUR+         %s\n" % (obs.duration,)
            output += "OBS_MODE         %s\n" % (obs.mode,)
            if obs.beamDipole is not None:
                output += "OBS_BDM          %i %6.4f %6.4f %s\n" % (tuple(obs.beamDipole))
            if obs.mode == 'TBF':
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TBN':
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_RADEC':
                output += "OBS_RA           %.9f\n" % (obs.ra,)
                output += "OBS_DEC          %+.9f\n" % (obs.dec,)
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_SOL':
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_JOV':
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'TRK_LUN':
                output += "OBS_B            %s\n" % (obs.beam,)
                output += "OBS_FREQ1        %i\n" % (obs.freq1,)
                output += "OBS_FREQ1+       %.9f MHz\n" % (obs.frequency1/1e6,)
                output += "OBS_FREQ2        %i\n" % (obs.freq2,)
                output += "OBS_FREQ2+       %.9f MHz\n" % (obs.frequency2/1e6,)
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
            elif obs.mode == 'STEPPED':
                output += "OBS_BW           %i\n" % (obs.filter,)
                output += "OBS_BW+          %s\n" % (render_bandwidth(obs.filter, obs.filter_codes),)
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
            ## TBF settings
            if obs.mode == 'TBF':
                output += "OBS_TBF_SAMPLES  %i\n" % (obs.samples,)
            ## TBN gain
            if obs.mode == 'TBN':
                if obs.gain != -1:
                    output += "OBS_TBN_GAIN     %i\n" % (obs.gain,)
            ## DRX gain
            else:
                if obs.gain != -1:
                    output += "OBS_DRX_GAIN     %i\n" % (obs.gain,)
            output += "\n"
            
        return output


class TBF(Observation):
    """Sub-class of Observation specifically for TBF observations.  It features a
    reduced number of parameters needed to setup the observation and provides extra
    information about the number of number of samples.
    
    .. note::
        TBF read-out times in ms are calculated using (samples / 196000 + 1) * 150 per
        MCS
    
    Required Arguments:
     * observation name
     * observation target
     * observation start date/time (UTC YYYY/MM/DD HH:MM:SS.SSS string)
     * observation frequency (Hz) - 1
     * observation frequency (Hz) - 2
     * integer filter code
     * integer number of samples
    
    Optional Keywords:
     * comments - comments about the observation
    """
    
    filter_codes = DRXFilters
    
    def __init__(self, name, target, start, frequency1, frequency2, filter, samples, comments=None):
        self.samples = int(samples)
        
        duration = (self.samples / _TBF_TIME_SCALE + 1)*_TBF_TIME_GAIN*(2 if frequency2 != 0 else 1) + 5000
        durStr = '%02i:%02i:%06.3f' % (int(duration/1000.0)/3600, int(duration/1000.0)%3600/60, duration/1000.0%60)
        Observation.__init__(self, name, target, start, durStr, 'TBF', 0.0, 0.0, frequency1, frequency2, filter, comments=comments)
        
    def update(self):
        """Update the computed parameters from the string values."""
        
        self.beam = self.get_beam_type()
        
        # Update the duration based on the number of bits and samples used
        self.dur = (self.samples / _TBF_TIME_SCALE + 1)*_TBF_TIME_GAIN*(2 if self.freq2 != 0 else 1) + 5000
        
        self.dataVolume = self.estimate_bytes()
        
    def estimate_bytes(self):
        """Estimate the data volume for the specified type and duration of 
        observations.  For TBF:
        
            bytes = samples / samplesPerFrame * 1224 bytes
        """
        
        try:
            sample_rate = DRXFilters[self.filter]
        except KeyError:
            sample_rate = 0.0
        nFramesTime = self.samples / (196e6 / fC)
        nFramesChan = math.ceil(sample_rate / fC / 12)
        nBytes = nFramesTime * nFramesChan * TBFSize
        return nBytes
        
    def validate(self, verbose=False):
        """Evaluate the observation and return True if it is valid, False
        otherwise."""
        
        self.update()
        
        station = lwasv
        if self._parent is not None:
            station = self._parent.station
        backend = station.interface.get_module('backend')
        be_name = station.interface.backend.rsplit('.', 1)[1].upper()
        
        failures = 0
        # Basic - Sample size, frequency, and filter
        if self.samples > 5*196000000:
            if verbose:
                pid_print(f"Error: Invalid number of samples ({self.samples} > {5*196000000})")
            failures += 1
        if self.freq1 < backend.DRX_TUNING_WORD_MIN or self.freq1 > backend.DRX_TUNING_WORD_MAX:
            if verbose:
                pid_print(f"Error: Specified frequency for tuning 1 is outside of the {be_name} tuning range")
            failures += 1
        if (self.freq2 < backend.DRX_TUNING_WORD_MIN or self.freq2 > backend.DRX_TUNING_WORD_MAX) and self.freq2 != 0:
            if verbose:
                pid_print(f"Error: Specified frequency for tuning 2 is outside of the {be_name} tuning range")
            failures += 1
        if self.filter not in [1, 2, 3, 4, 5, 6, 7]:
            if verbose:
                pid_print(f"Error: Invalid filter code '{self.filter}'")
            failures += 1
            
        # Advanced - Data Volume
        if self.dataVolume >= (_DRSUCapacityTB*1024**4):
            if verbose:
                pid_print(f"Error: Data volume exceeds {_DRSUCapacityTB} TB DRSU limit")
            failures += 1
            
        # Advanced - ASP
        failures += not self._validate_asp(verbose=verbose)
        
        # Any failures indicates a bad observation
        if failures == 0:
            return True
        else:
            return False


class Session(_Session):
    """Class to hold all of the observations in a session."""
    
    _allowed_modes = (TBF, TBN, DRX, Stepped)
    
    def __init__(self, name, id, observations=None, data_return_method='DRSU', comments=None, station=lwasv):
        _Session.__init__(self, name, id, 
                          observations=observations, data_return_method=data_return_method, 
                          comments=comments, station=station)
        
    @property
    def station(self):
        return self._station
        
    @station.setter
    def station(self, value):
        """
        Update the station used by the project for source computations.
        
        .. versionadded:: 1.2.0
        """
        
        if value.interface.sdf != 'lsl.common.sdfADP':
            raise ValueError("Incompatible station: expected %s, got %s" % \
                             (value.interface.sdf, 'lsl.common.sdfADP'))
            
        self._station = value
        self.update()
        
    def validate(self, verbose=False):
        """Examine all of the observations associated with the session to check
        for validity.  If everything is valid, return True.  Otherwise, return
        False."""
        
        self.update()
        
        failures = 0
        is_valid = _Session.validate(self, verbose=verbose)
        if not is_valid:
            failures += 1
            
        # Validate beam number for TBF
        if len(self.observations) > 0:
            if self.observations[0].mode ==  'TBF':
                if self.drx_beam != 1:
                    if verbose:
                        pid_print("Error: TBF can only run on beam 1")
                    failures += 1
                    
        if failures == 0:
            return True
        else:
            return False


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
    
    # Build up a string representing the observation duration.
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
        pid_print(f"Obs {obs_temp['id']} is mode {mode}")
        
    if mode == 'TBF':
        obsOut = TBF(obs_temp['name'], obs_temp['target'], utcString, f1, f2, obs_temp['filter'], obs_temp['tbfSamples'], comments=obs_temp['comments'])
    elif mode == 'TBN':
        obsOut = TBN(obs_temp['name'], obs_temp['target'], utcString, durString, f1, obs_temp['filter'], gain=obs_temp['gain'], comments=obs_temp['comments'])
    elif mode == 'TRK_RADEC':
        obsOut = DRX(obs_temp['name'], obs_temp['target'], utcString, durString, obs_temp['ra'], obs_temp['dec'], f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'TRK_SOL':
        obsOut = Solar(obs_temp['name'], obs_temp['target'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'TRK_JOV':
        obsOut = Jovian(obs_temp['name'], obs_temp['target'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'TRK_LUN':
        obsOut = Lunar(obs_temp['name'], obs_temp['target'], utcString, durString, f1, f2, obs_temp['filter'], gain=obs_temp['gain'], max_snr=obs_temp['MaxSNR'], comments=obs_temp['comments'])
    elif mode == 'STEPPED':
        if verbose:
            pid_print(f"-> found {len(beam_temps)} steps")
            
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
        raise RuntimeError(f"Invalid mode encountered: {mode}")
        
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
                'stpRADec': True, 'tbwBits': 12, 'tbfSamples': 0, 'gain': -1, 
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
                    pid_print(f"Started obs {value}")
                    
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
            if keyword == 'OBS_TBF_SAMPLES':
                obs_temp['tbfSamples'] = int(value)
                continue
            if keyword == 'OBS_TBN_GAIN':
                obs_temp['gain'] = int(value)
                continue
            if keyword == 'OBS_DRX_GAIN':
                obs_temp['gain'] = int(value)
                continue
            
            # Keywords that might indicate this is for DP-based stations/actually an IDF
            if keyword in ('OBS_TBW_BITS', 'OBS_TBW_SAMPLES', 'RUN_ID'):
                raise RuntimeError(f"Invalid keyword encountered: {keyword}")
            
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
