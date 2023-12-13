"""
Module for working with an MCS meta-data tarball and extracting the useful bits out 
it and putting those bits into Python objects, e.g, :class:`lsl.common.stations.LWAStation` 
and :class:`lsl.common.sdm.SDM`.
"""

import os
import re
import copy
import glob
import shutil
import tarfile
import tempfile
from functools import lru_cache
from datetime import datetime, timedelta

from lsl.common import stations, sdmADP, sdfADP
from lsl.common.mcsADP import *
from lsl.common.adp import word_to_freq, fS
from lsl.common.color import colorfy

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '1.1'
__all__ = ['read_ses_file', 'read_obs_file', 'read_cs_file', 'get_sdm', 'get_beamformer_min_delay',
           'get_station', 'get_session_metadata', 'get_session_spec', 'get_observation_spec',
           'get_sdf', 'get_command_script', 'get_asp_configuration', 'get_asp_configuration_summary',
           'is_valid']

# Regular expression for figuring out filenames
filenameRE = re.compile(r'(?P<projectID>[a-zA-Z0-9]{1,8})_(?P<sessionID>\d+)(_(?P<obsID>\d+)(_(?P<obsOutcome>\d+))?)?.*\..*')


@lru_cache(maxsize=1)
def _open_tarball(tarname):
    return tarfile.open(tarname, mode='r:*')


@lru_cache(maxsize=1)
def _get_members(tarname):
    tf = _open_tarball(tarname)
    return tf.getmembers()


class managed_mkdtemp(object):
    """
    Wrapper class around tempfile.mkdtemp to enable 'with' statements with 
    automatic cleanup.
    """
    
    def __init__(self, suffix='', prefix='tmp', dir=None):
        self._dir = tempfile.mkdtemp(suffix, prefix, dir)
        
    def __enter__(self):
        return self._dir
        
    def __exit__(self, type, value, tb):
        shutil.rmtree(self._dir, ignore_errors=True)
        
    @property
    def name(self):
        return self._dir


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
        
    record = {'ASP': bses.SESSION_MRP_ASP, 'ADP': bses.SESSION_MRP_DP_, 'SHL': bses.SESSION_MRP_SHL, 
              'MCS': bses.SESSION_MRP_MCS, 'DR1': bses.SESSION_MRP_DR1, 'DR2': bses.SESSION_MRP_DR2,
              'DR3': bses.SESSION_MRP_DR3, 'DR4': bses.SESSION_MRP_DR4}
    
    update = {'ASP': bses.SESSION_MUP_ASP, 'ADP': bses.SESSION_MUP_DP_, 'SHL': bses.SESSION_MUP_SHL, 
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
            
        if bheader.OBS_B > 2:
            ## Pre OBS_BDM
            fh.seek(0)
            
            newStruct = []
            for line in OSF_STRUCT.split('\n'):
                if line.find('OBS_BDM') != -1:
                    continue
                newStruct.append(line)
            newStruct = '\n'.join(newStruct)
            
            bheader = parse_c_struct(newStruct, endianness='little')
            fh.readinto(bheader)
            bheader.OBS_BDM = ''
            
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
                raise IOError(f"Byte alignment lost at byte {fh.tell()}")
                
        fh.readinto(bfooter)
        
        if bfooter.alignment != (2**32 - 1):
            raise IOError(f"Byte alignment lost at byte {fh.tell()}")
            
    output = {'version': bheader.FORMAT_VERSION,
              'project_id': bheader.PROJECT_ID.lstrip().rstrip(), 'session_id': bheader.SESSION_ID,
              'drx_beam': bheader.SESSION_DRX_BEAM, 'spcSetup': bheader.SESSION_SPC,
              'obs_id': bheader.OBS_ID,
              'mjd': bheader.OBS_START_MJD, 'mpm': bheader.OBS_START_MPM, 'dur': bheader.OBS_DUR,
              'mode': bheader.OBS_MODE, 'beamdipole_mode': bheader.OBS_BDM, 
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
    output['tbn_gain'] = bfooter.OBS_TBN_GAIN
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
                               'subsystem_id': sid_to_string(action.sid),
                               'command_id': cid_to_string(action.cid), 
                               'command_length': action.len, 'data': data}
                if actionPrime['subsystem_id'] == 'DP':
                    raise RuntimeError("Command script references DP not ADP")
                    
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
        dynamic = sdmADP.parse_sdm(os.path.join(tempDir, 'dynamic', 'sdm.dat'))
        
    return dynamic


def get_beamformer_min_delay(tarname):
    """
    Given an MCS meta-data tarball, extract the minimum beamformer delay in 
    samples and return it.  If no minimum delay can be found in the tarball, 
    None is returned.
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Extract the mindelay.txt file.  If mindelay.txt cannot be found, None
        # is returned via the try...except block.
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('mindelay.txt')
        except KeyError:
            return None
        tf.extractall(path=tempDir, members=[ti,])
        
        # Parse the SDM file and build the SDM instance
        with open(os.path.join(tempDir, 'mindelay.txt'), 'r') as fh:
            try:
                mindelay = int(fh.read(), 10)
            except ValueError:
                mindelay = None
                
    return mindelay


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


def get_session_metadata(tarname):
    """
    Given an MCS meta-data tarball, extract the session meta-data file (MCS0030, 
    Section 7) and return a dictionary of observations that contain dictionaries 
    of the OP_TAG (tag), DRSU Barcode (drsu), OBS_OUTCOME (outcome), and the 
    MSG (msg).
    
    .. versionchanged:: 0.6.5
        Update to the new _metadata.txt format
    """
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        path, basename = os.path.split(tarname)
        basename, ext = os.path.splitext(basename)
        
        # Extract the session meta-data file (_metadata.txt)
        tf = _open_tarball(tarname)
        try:
            ti = tf.getmember('%s_metadata.txt' % basename)
        except KeyError:
            for ti in _get_members(tarname):
                if ti.name[-13:] == '_metadata.txt':
                    break
        tf.extractall(path=tempDir, members=[ti,])
        
        # Read in the SMF
        filename = os.path.join(tempDir, ti.name)
        with open(filename, 'r') as fh:
            # Define a regular expresion to match the latest format
            lineRE = re.compile(r"\s*(?P<id>\d{1,}?)\s+\[(?P<tag>[\d_]+?)\]\s+\['?(?P<barcode>.+?)'?\]\s+(?P<outcome>\d)\s+\[(?P<msg>.*?)\]")
            
            result = {}
            for line in fh:
                line = line.replace('\n', '')
                if len(line) == 0:
                    continue
                    
                mtch = lineRE.search(line)
                if mtch is not None:
                    ## If it matches the new format
                    obsID = mtch.group('id')
                    opTag = mtch.group('tag')
                    drsuBarcode = mtch.group('barcode')
                    if drsuBarcode[:3] == 'Err':
                        try:
                            drsuBarcode = result[int(obsID)-1]['barcode']
                        except KeyError:
                            drsuBarcode = 'UNK'
                    obsOutcome = mtch.group('outcome')
                    msg = mtch.group('msg')
                    
                else:
                    ## Otherwise, I don't really know how the messages will look so we use this try...except
                    ## block should take care of the various situations.
                    try:
                        obsID, opTag, drsuBarcode, obsOutcome, msg = line.split(None, 4)
                        opTag = opTag.replace('[', '')
                        opTag = opTag.replace(']', '')
                        drsuBarcode = drsuBarcode.replace('[', '')
                        drsuBarcode = drsuBarcode.replace(']', '')
                        drsuBarcode = drsuBarcode.replace("'", '')
                    except ValueError:
                        try:
                            obsID, opTag, drsuBarcode, obsOutcome = line.split(None, 3)
                            msg = 'UNK'
                        except ValueError:
                            try:
                                obsID, opTag, obsOutcome = line.split(None, 2)
                                drsuBarcode = 'UNK'
                                obsOutcome = '-1'
                                msg = 'UNK'
                            except ValueError:
                                obsID, obsOutcome = line.split(None, 1)
                                drsuBarcode = 'UNK'
                                opTag = 'UNK'
                                msg = 'UNK'
                                
                obsID = int(obsID)
                obsOutcome = int(obsOutcome) if obsOutcome != 'Failed' else 1
                result[obsID] = {'tag': opTag, 'barcode': drsuBarcode, 
                                 'outcome': obsOutcome, 'msg': msg}
                
    # Return
    return result


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
        This function returns a full :class:`lsl.common.sdfADP.Project` instance 
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
        project = sdfADP.parse_sdf(os.path.join(tempDir, ti.name))
        
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


def get_asp_configuration(tarname, which='beginning'):
    """
    Given an MCS meta-data tarball, extract the ASP MIB contained in it and return 
    a dictionary of values for the filter, AT1, AT2, and ATSplit.  The 'which'
    keyword is used to specify whether or not the configuration returned is at the
    beginning (default) or end of the session.
    
    .. versionadded:: 0.6.5
    """
    
    which = which.lower()
    if which not in ('beginning', 'begin', 'end'):
        raise ValueError(f"Unknown configuration time '{which}'")
        
    # Stub ASP configuration
    aspConfig = {'asp_filter':      [-1 for i in range(264)],
                 'asp_atten_1':     [-1 for i in range(264)],
                 'asp_atten_2':     [-1 for i in range(264)],
                 'asp_atten_split': [-1 for i in range(264)]}
    
    with managed_mkdtemp(prefix='metadata-bundle-') as tempDir:
        # Find the .pag file and extract it
        tf = _open_tarball(tarname)
        mibs = []
        for ti in _get_members(tarname):
            if ti.name.find('_ASP_%s' % which[:5]) != -1:
                if ti.name[-4:] in ('.pag', '.gdb'):
                    mibs.append(ti)
                    
        if len(mibs) > 0:
            tf.extractall(path=tempDir, members=mibs)
            
            # Read in the right MIB
            aspMIB = {}
            for mib in mibs:
                if which[:5] == 'begin' and mib.name.find('_ASP_begin') == -1:
                    continue
                if which == 'end' and mib.name.find('_ASP_end') == -1:
                    continue
                    
                aspMIB = MIB()
                aspMIB.from_file(os.path.join(tempDir, mib.name))
                break
                
            # Extract the configuration
            for key in aspMIB.keys():
                values = [int(v) for v in key.split('.')]
                
                if values[0] == 3:
                    ## Filter
                    aspConfig['asp_filter'][values[1]-1] = int(aspMIB[key].value)
                    
                elif values[0] == 4:
                    ## Attenuators
                    if values[1] == 1:
                        ### AT1
                        aspConfig['asp_atten_1'][values[2]-1] = int(aspMIB[key].value)
                    elif values[1] == 2:
                        ### AT2
                        aspConfig['asp_atten_2'][values[2]-1] = int(aspMIB[key].value)
                    elif values[1] == 3:
                        ### ATSPLIT
                        aspConfig['asp_atten_split'][values[2]-1] = int(aspMIB[key].value)
                    else:
                        pass
                        
                else:
                    pass
                    
    return aspConfig


def get_asp_configuration_summary(tarname, which='beginning'):
    """
    Similar to get_asp_configuration, but returns only a single value for each
    of the four ASP paramters:  filter, AT, AT2, and ATSplit.  The values
    are based off the mode of the parameter.
    
    .. versionadded:: 0.6.5
    """
    
    # Get the full configuration
    aspConfig = get_asp_configuration(tarname, which=which)
    
    # Count
    count = {}
    for param in aspConfig.keys():
        count[param] = {}
        for ant in range(len(aspConfig[param])):
            value = aspConfig[param][ant]
            try:
                count[param][value] += 1
            except KeyError:
                count[param][value] = 1
                
    # Modes
    mode = {}
    for param in count.keys():
        best = 0
        mode[param] = 0
        
        for value in count[param].keys():
            num = count[param][value]
            if num > best:
                best = num
                mode[param] = value
                
    # Done
    return mode


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
    except IOError as e:
        raise e
    except:
        failures += 1
        if verbose:
            print(colorfy("Session specification - {{%red {{%bold FAILED}}}}"))
        
    try:
        get_observation_spec(tarname)
        passes += 1
        if verbose:
            print(colorfy("Observation specification(s) - {{%green OK}}"))
    except:
        failures += 1
        if verbose:
            print(colorfy("Observation specification(s) - {{%red {{%bold FAILED}}}}"))
            
    try:
        get_command_script(tarname)
        passes += 1
        if verbose:
            print(colorfy("Command script - {{%green OK}}"))
    except:
        failures += 1
        if verbose:
            print(colorfy("Command script - {{%red {{%bold FAILED}}}}"))
            
    if verbose:
        print("---")
        print("%i passed / %i failed" % (passes, failures))
        
    return False if failures else True
