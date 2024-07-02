"""
Module for working with any MCS meta-data tarball and extracting the useful bits
out it and putting those bits into Python objects, e.g, 
:class:`lsl.common.stations.LWAStation` and :class:`lsl.common.sdmDP.SDM`.
"""

import os

from lsl.common import metabundleDP, metabundleADP, metabundleNDP

from lsl.misc import telemetry
telemetry.track_module()

__version__ = '1.2'
__all__ = ['get_sdm', 'get_beamformer_min_delay', 'get_station',
           'get_session_metadata', 'get_session_spec', 'get_observation_spec',
           'get_sdf', 'get_command_script', 'get_asp_configuration',
           'get_asp_configuration_summary']


def get_sdm(tarname):
    """
    Given an MCS meta-data tarball, extract the information stored in the 
    dynamic/sdm.dat file and return a :class:`lsl.common.sdmDP.SDM` instance
    describing the dynamic condition of the station.
    
    If a sdm.dat file cannot be found in the tarball, None is returned.
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_sdm(tarname)
        except:
            pass
    raise RuntimeError("Failed to read SDM")


def get_beamformer_min_delay(tarname):
    """
    Given an MCS meta-data tarball, extract the minimum beamformer delay in 
    samples and return it.  If no minimum delay can be found in the tarball, 
    None is returned.
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_beamformer_min_delay(tarname)
        except:
            pass
    raise RuntimeError("Failed to read beamformer minimum delay")


def get_station(tarname, apply_sdm=True):
    """
    Given an MCS meta-data tarball, extract the information stored in the ssmif.dat 
    file and return a :class:`lsl.common.stations.LWAStation` object.  Optionally, 
    update the :class:`lsl.common.stations.Antenna` instances associated whith the
    LWAStation object using the included SDM file.
    
    If a ssmif.dat file cannot be found in the tarball, None is returned.  
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_station(tarname, apply_sdm=apply_sdm)
        except:
            pass
    raise RuntimeError("Failed to get station object")


def get_session_metadata(tarname):
    """
    Given an MCS meta-data tarball, extract the session meta-data file (MCS0030, 
    Section 7) and return a dictionary of observations that contain dictionaries 
    of the OP_TAG (tag), DRSU Barcode (drsu), OBS_OUTCOME (outcome), and the 
    MSG (msg).
    
    .. versionchanged:: 0.6.5
        Update to the new _metadata.txt format
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_session_metadata(tarname)
        except:
            pass
    raise RuntimeError("Failed to get session metadata")


def get_session_spec(tarname):
    """
    Given an MCS meta-data tarball, extract the session specification file (MCS0030, 
    Section 5) and return a dictionary of parameters.
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_session_spec(tarname)
        except:
            pass
    raise RuntimeError("Failed to get session specifications")


def get_observation_spec(tarname, obs_id=None):
    """
    Given an MCS meta-data tarball, extract one or more observation specification 
    file (MCS0030, Section 6) and return a list of dictionaries corresponding to
    each OBS file.  If the `obs_id` keyword is set to a list of observation
    numbers, only observations matching the numbers in `obs_id` are returned.
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_observation_spec(tarname, obs_id=obs_id)
        except:
            pass
    raise RuntimeError("Failed to get observation specifications")


def get_sdf(tarname):
    """
    Given an MCS meta-data tarball, extract the session specification file, the 
    session meta-data file, and all observation specification files to build up
    a SDF-representation of the session.
    
    .. note::
        This function returns a full :class:`lsl.common.sdf.Project` instance 
        with the session in question stored under `project.sessions[0]` and the 
        observations under `project.sessions[0].observations`.
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_sdf(tarname)
        except:
            pass
    raise RuntimeError("Failed to get SDF")


def get_command_script(tarname):
    """
    Given an MCS meta-data tarball, extract the command script and parse it.  The
    commands are returned as a list of dictionaries (one dictionary per command).
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_command_script(tarname)
        except:
            pass
    raise RuntimeError("Failed to get command script")


def get_asp_configuration(tarname, which='beginning'):
    """
    Given an MCS meta-data tarball, extract the ASP MIB contained in it and return 
    a dictionary of values for the filter, AT1, AT2, and ATSplit.  The 'which'
    keyword is used to specify whether or not the configuration returned is at the
    beginning (default) or end of the session.
    
    .. versionadded:: 0.6.5
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    which = which.lower()
    if which not in ('beginning', 'begin', 'end'):
        raise ValueError(f"Unknown configuration time '{which}'")
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_asp_configuration(tarname, which=which)
        except:
            pass
    raise RuntimeError("Failed to get ASP configurations for '%s'" % which)


def get_asp_configuration_summary(tarname, which='beginning'):
    """
    Similar to get_asp_configuration, but returns only a single value for each
    of the four ASP paramters:  filter, AT, AT2, and ATSplit.  The values
    are based off the mode of the parameter.
    
    .. versionadded:: 0.6.5
    """
    
    if not os.path.isfile(tarname) or not os.access(tarname, os.R_OK):
        raise OSError("%s does not exists or is not readable" % tarname)
        
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
            try:
                return backend.get_asp_configuration_summary(tarname, which=which)
            except:
                pass
    raise RuntimeError("Failed to get ASP configuration summary for '%s'" % which)
