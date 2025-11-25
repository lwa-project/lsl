"""
Module for working with an MCS meta-data tarball and extracting the useful bits out 
it and putting those bits into Python objects, e.g, :class:`lsl.common.stations.LWAStation` 
and :class:`lsl.common.sdm.SDM`.
"""

import os
import re
import copy
import glob
import logging
from functools import lru_cache
from datetime import datetime, timedelta

from lsl.common._metabundle_utils import *
from lsl.common import stations, sdmNDP, sdfNDP
from lsl.common.mcsNDP import *
from lsl.common.ndp import word_to_freq, fS
from lsl.common.color import colorfy

from lsl.logger import LSL_LOGGER

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.1'
__all__ = ['read_ses_file', 'read_obs_file', 'read_cs_file', 'get_sdm', 'get_beamformer_min_delay',
           'get_station', 'get_session_metadata', 'get_session_spec', 'get_observation_spec',
           'get_sdf', 'get_command_script', 'get_asp_configuration', 'get_asp_configuration_summary',
           'is_valid']


def read_ses_file(filename):
    """
    Read in a session specification file (MCS0030, Section 5) and return the data
    as a dictionary.
    """
    
    # Read the SES
    with open(filename, 'rb') as fh:
        bses = parse_c_struct(SSF_STRUCT, endianness='little')
        fh.readinto(bses)
        
        ## LWA-1 check
        #if bses.FORMAT_VERSION not in (6,):
        #	fh.close()
        #	raise RuntimeError("Version mis-match: File appears to be from LWA-1")
        
    record = {'ASP': bses.SESSION_MRP_ASP, 'NDP': bses.SESSION_MRP_DP_, 'SHL': bses.SESSION_MRP_SHL, 
              'MCS': bses.SESSION_MRP_MCS, 'DR1': bses.SESSION_MRP_DR1, 'DR2': bses.SESSION_MRP_DR2,
              'DR3': bses.SESSION_MRP_DR3, 'DR4': bses.SESSION_MRP_DR4}
    
    update = {'ASP': bses.SESSION_MUP_ASP, 'NDP': bses.SESSION_MUP_DP_, 'SHL': bses.SESSION_MUP_SHL, 
              'MCS': bses.SESSION_MUP_MCS, 'DR1': bses.SESSION_MUP_DR1, 'DR2': bses.SESSION_MUP_DR2,
              'DR3': bses.SESSION_MUP_DR3, 'DR4': bses.SESSION_MUP_DR4}
    
    return {'version': bses.FORMAT_VERSION,
            'project_id': bses.PROJECT_ID.lstrip().rstrip(), 'session_id': bses.SESSION_ID, 
            'configuration_authority': bses.SESSION_CRA,
            'drx_beam': bses.SESSION_DRX_BEAM, 'spcSetup': bses.SESSION_SPC, 
            'mjd': bses.SESSION_START_MJD, 'mpm': bses.SESSION_START_MPM, 'dur': bses.SESSION_DUR,
            'nobs': bses.SESSION_NOBS,
            'recordMIB': record, 'updateMIB': update, 
            'include_mcssch_log': bses.SESSION_LOG_SCH, 'include_mcsexe_log': bses.SESSION_LOG_EXE,
            'include_station_smib': bses.SESSION_INC_SMIB, 'include_station_design': bses.SESSION_INC_DES}


def read_obs_file(filename):
    """
    Read in a observation specification file (MCS0030, Section 6) and return the
    data as a dictionary.
    """
    
    # Read the OBS
    with open(filename, 'rb') as fh:
        bheader = parse_c_struct(OSF_STRUCT, endianness='little')
        bstep   = parse_c_struct(OSFS_STRUCT, endianness='little')
        bbeam   = parse_c_struct(BEAM_STRUCT, endianness='little')
        bfooter = parse_c_struct(OSF2_STRUCT, endianness='little')
        fh.readinto(bheader)
        
        ## LWA-1 check
        #if bheader.FORMAT_VERSION not in (6,):
        #	fh.close()
        #	raise RuntimeError("Version mis-match: File appears to be from LWA-1")
        
        steps = []
        for n in range(bheader.OBS_STP_N):
            fh.readinto(bstep)
            if bstep.OBS_STP_B == 3:
                fh.readinto(bbeam)
                bstep.delay = copy.deepcopy(bbeam.OBS_BEAM_DELAY)
                bstep.gain  = copy.deepcopy(flat_to_multi(bbeam.OBS_BEAM_GAIN, *bbeam.dims['OBS_BEAM_GAIN']))
            else:
                bstep.delay = []
                bstep.gain  = []
            
            steps.append(copy.deepcopy(bstep))
            
            alignment = parse_c_struct("""
            unsigned int block;
            """, endianness='little')
            fh.readinto(alignment)
            
            if alignment.block != (2**32 - 2):
                raise IOError("Byte alignment lost at byte %i" % fh.tell())
                
        fh.readinto(bfooter)
        
        if bfooter.alignment != (2**32 - 1):
            raise IOError("Byte alignment lost at byte %i" % fh.tell())
            
    output = {'version': bheader.FORMAT_VERSION,
              'project_id': bheader.PROJECT_ID.lstrip().rstrip(), 'session_id': bheader.SESSION_ID,
              'drx_beam': bheader.SESSION_DRX_BEAM, 'spcSetup': bheader.SESSION_SPC,
              'obs_id': bheader.OBS_ID,
              'mjd': bheader.OBS_START_MJD, 'mpm': bheader.OBS_START_MPM, 'dur': bheader.OBS_DUR,
              'mode': ObservingMode(bheader.OBS_MODE), 'beamdipole_mode': bheader.OBS_BDM, 
              'ra': bheader.OBS_RA, 'dec': bheader.OBS_DEC,
              'beam': bheader.OBS_B, 
              'freq1': word_to_freq(bheader.OBS_FREQ1), 'freq2': word_to_freq(bheader.OBS_FREQ2),
              'bw': bheader.OBS_BW,
              'nsteps': bheader.OBS_STP_N, 'is_radec': bheader.OBS_STP_RADEC,  'steps': steps, 
              'fee_power': flat_to_multi(bfooter.OBS_FEE, *bfooter.dims['OBS_FEE']), 
              'asp_filter': list(bfooter.OBS_ASP_FLT), 'asp_atten_1': list(bfooter.OBS_ASP_AT1), 
              'asp_atten_2': list(bfooter.OBS_ASP_AT2), 'asp_atten_split': list(bfooter.OBS_ASP_ATS)}
    output['tbf_samples'] = bfooter.OBS_TBF_SAMPLES
    output['tbf_gain'] = bfooter.OBS_TBF_GAIN
    output['drx_gain'] = bfooter.OBS_DRX_GAIN
    
    return output


def read_cs_file(filename):
    """
    Read in a command script file (MCS0030, currently undocumented) and return the
    data as a list of dictionaries.
    """
    
    # Read the CS file
    with open(filename, 'rb') as fh:
        commands = []
        while True:
            action = parse_c_struct("""
            long int tv[2];
            int bASAP;
            int sid;
            int cid;
            int len;
            """, endianness='little')
            
            try:
                fh.readinto(action)
                
                if action.tv[0] == 0:
                    break
                    
                if action.len > 0:
                    data = parse_c_struct("""
                    char data[%i];
                    """ % action.len, endianness='little')
                    
                    fh.readinto(data)
                    data = data.data
                else:
                    data = None
                
                actionPrime = {'time': action.tv[0] + action.tv[1]/1.0e6, 
                               'ignore_time': True if action.bASAP else False, 
                               'subsystem_id': SubsystemID(action.sid),
                               'command_id': CommandID(action.cid), 
                               'command_length': action.len, 'data': data}
                if actionPrime['subsystem_id'] == SubsystemID.DP:
                    raise RuntimeError("Command script references DP not NDP")
                if actionPrime['subsystem_id'] == SubsystemID.ADP:
                    raise RuntimeError("Command script references ADP not NDP")
                    
                commands.append( actionPrime )
            except IOError:
                break
                
    return commands


def get_sdm(tarname):
    """
    Given an MCS meta-data tarball, extract the information stored in the 
    dynamic/sdm.dat file and return a :class:`lsl.common.sdm.SDM` instance
    describing the dynamic condition of the station.
    
    If a sdm.dat file cannot be found in the tarball, None is returned.
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Extract the SDM file.  If the dynamic/sdm.dat file cannot be found, None
        # is returned via the try...except block.
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('dynamic/sdm.dat')
        except KeyError:
            return None
        tf.extractall(path=tempDir, members=[ti,])
        
        # Parse the SDM file and build the SDM instance
        dynamic = sdmNDP.parse_sdm(os.path.join(tempDir, 'dynamic', 'sdm.dat'))
        
    return dynamic


def get_station(tarname, apply_sdm=True):
    """
    Given an MCS meta-data tarball, extract the information stored in the ssmif.dat 
    file and return a :class:`lsl.common.stations.LWAStation` object.  Optionally, 
    update the :class:`lsl.common.stations.Antenna` instances associated whith the
    LWAStation object using the included SDM file.
    
    If a ssmif.dat file cannot be found in the tarball, None is returned.  
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Extract the SSMIF and SDM files.  If the ssmif.dat file cannot be found, None
        # is returned via the try...except block
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('ssmif.dat')
        except KeyError:
            return None
        tf.extractall(path=tempDir, members=[ti,])
        
        # Read in the SSMIF
        station = stations.parse_ssmif(os.path.join(tempDir, 'ssmif.dat'))
        
        # Get the beamformer minimum delay, if found
        mindelay = get_beamformer_min_delay(tarname)
        if mindelay is not None:
            station.beamformer_min_delay_samples = mindelay
            station.beamformer_min_delay = mindelay/fS
            
        # Get the SDM (if we need to)
        if apply_sdm:
            dynamic = get_sdm(tarname)
        else:
            dynamic = None
        
        # Update the SSMIF entries
        if dynamic is not None:
            newAnts = dynamic.update_antennas(station.antennas)
            station.antennas = newAnts
            
    # Return
    return station


def get_session_spec(tarname):
    """
    Given an MCS meta-data tarball, extract the session specification file (MCS0030, 
    Section 5) and return a dictionary of parameters.
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        path, basename = os.path.split(tarname)
        basename, ext = os.path.splitext(basename)
        
        # Extract the session specification file (.ses)
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('%s.ses' % basename)
        except KeyError:
            for ti in _get_members(tarname):
                if ti.name[-4:] == '.ses':
                    break
        tf.extractall(path=tempDir, members=[ti,])
        
        # Read in the SES
        ses = read_ses_file(os.path.join(tempDir, ti.name))
        
    # Return
    return ses


def get_observation_spec(tarname, obs_id=None):
    """
    Given an MCS meta-data tarball, extract one or more observation specification 
    file (MCS0030, Section 6) and return a list of dictionaries corresponding to
    each OBS file.  If the `obs_id` keyword is set to a list of observation
    numbers, only observations matching the numbers in `obs_id` are returned.
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Find all of the .obs files and extract them
        tf = _open_tarball(tarname)
        tis = []
        for ti in _get_members(tarname):
            if ti.name[-4:] == '.obs':
                tis.append(ti)
        tf.extractall(path=tempDir, members=tis)
        
        # Read in the OBS files
        obsList = []
        for of in sorted(glob.glob(os.path.join(tempDir, '*.obs'))):
            obsList.append( read_obs_file(of) )
            
        # Cull the list based on the observation ID selection
        if obs_id is None:
            outObs = obsList
        else:
            try:
                len(obs_id)
            except TypeError:
                obs_id = [obs_id]
                
            outObs = []
            for o in obsList:
                try:
                    if o['obs_id'] in obs_id:
                        outObs.append(o)
                except TypeError:
                    if o['obs_id'] == obs_id:
                        outObs.append(o)
                        
            if len(outObs) == 1:
                outObs = outObs[0]
                
    # Return
    return outObs


def get_sdf(tarname):
    """
    Given an MCS meta-data tarball, extract the session specification file, the 
    session meta-data file, and all observation specification files to build up
    a SDF-representation of the session.
    
    .. note::
        This function returns a full :class:`lsl.common.sdfNDP.Project` instance 
        with the session in question stored under `project.sessions[0]` and the 
        observations under `project.sessions[0].observations`.
    """
    
    # Find the SDF file contained in the tarball
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Find the right .txt file (not the metadata one) and extract it
        tf = _open_tarball(tarname)
        for ti in _get_members(tarname):
            if ti.name[-4:] == '.txt' and ti.name.find('metadata') == -1:
                break
        tf.extractall(path=tempDir, members=[ti,])
        
        # Parse it
        project = sdfNDP.parse_sdf(os.path.join(tempDir, ti.name))
        
    # Return the filled-in SDF instance
    return project


def get_command_script(tarname):
    """
    Given an MCS meta-data tarball, extract the command script and parse it.  The
    commands are returned as a list of dictionaries (one dictionary per command).
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Find the .cs file and extract it
        tf = _open_tarball(tarname)
        for ti in _get_members(tarname):
            if ti.name[-3:] == '.cs':
                break
        tf.extractall(path=tempDir, members=[ti,])
        
        # Read in the CS
        cs = read_cs_file(os.path.join(tempDir, ti.name))
        
    # Return
    return cs


def is_valid(tarname, verbose=False):
    """
    Given a filename, see if it is valid metadata tarball or not.
    
    .. versionadded:: 1.2.0
    """
    
    passes = 0
    failures = 0
    try:
        get_session_spec(tarname)
        passes += 1
        if verbose:
            print(colorfy("Session specification - {{%green OK}}"))
        LSL_LOGGER.info("Session specification - OK")
    except IOError as e:
        raise e
    except:
        failures += 1
        if verbose:
            print(colorfy("Session specification - {{%red {{%bold FAILED}}}}"))
        LSL_LOGGER.error("Session specification - FAILED")
        
    try:
        get_observation_spec(tarname)
        passes += 1
        if verbose:
            print(colorfy("Observation specification(s) - {{%green OK}}"))
        LSL_LOGGER.info("Observation specification(s) - OK")
    except:
        failures += 1
        if verbose:
            print(colorfy("Observation specification(s) - {{%red {{%bold FAILED}}}}"))
        LSL_LOGGER.error("Observation specification(s) - FAILED")
        
    try:
        get_command_script(tarname)
        passes += 1
        if verbose:
            print(colorfy("Command script - {{%green OK}}"))
        LSL_LOGGER.info("Command script - OK")
    except:
        failures += 1
        if verbose:
            print(colorfy("Command script - {{%red {{%bold FAILED}}}}"))
        LSL_LOGGER.error("Command script - FAILED")
        
    if verbose:
        print("---")
        print(f"{passes} passed / {failures} failed")
    LSL_LOGGER.info(f"{passes} passed / {failures} failed")
    
    return False if failures else True
