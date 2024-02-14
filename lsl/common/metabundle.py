"""
Module for working with an MCS meta-data tarball and extracting the useful bits out 
it and putting those bits into Python objects, e.g, :class:`lsl.common.stations.LWAStation` 
and :class:`lsl.common.sdm.SDM`.
"""

__version__ = '1.2'
__all__ = ['get_sdm', 'get_beamformer_min_delay', 'get_station',
           'get_session_metadata', 'get_session_spec', 'get_observation_spec',
           'get_sdf', 'get_command_script', 'get_asp_configuration',
           'get_asp_configuration_summary']


from lsl.common import metabundleDP, metabundleADP, metabundleNDP


def get_sdm(tarname):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_sdm(tarname)
        except:
            pass
    raise RuntimeError("Failed to read SDM")


def get_beamformer_min_delay(tarname):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_beamformer_min_delay(tarname)
        except:
            pass
    raise RuntimeError("Failed to read beamformer minimum delay")


def get_station(tarname, apply_sdm=True):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_station(tarname, apply_sdm=apply_sdm)
        except:
            pass
    raise RuntimeError("Failed to get station object")


def get_session_metadata(tarname):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_session_metadata(tarname)
        except:
            pass
    raise RuntimeError("Failed to get session metadata")


def get_session_spec(tarname):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_session_spec(tarname)
        except:
            pass
    raise RuntimeError("Failed to get session specifications")


def get_observation_spec(tarname, obs_id=None):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_observation_spec(tarname, obs_id=obs_id)
        except:
            pass
    raise RuntimeError("Failed to get observation specifications")


def get_sdf(tarname):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_sdf(tarname)
        except:
            pass
    raise RuntimeError("Failed to get SDF")


def get_command_script(tarname):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_command_script(tarname)
        except:
            pass
    raise RuntimeError("Failed to get command script")


def get_asp_configuration(tarname, which='beginning'):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
        try:
            return backend.get_asp_configuration(tarname, which=which)
        except:
            pass
    raise RuntimeError("Failed to get ASP configurations for '%s'" % which)


def get_asp_configuration_summary(tarname, which='beginning'):
    for backend in (metabundleDP, metabundleADP, metabundleNDP):
            try:
                return backend.get_asp_configuration_summary(tarname, which=which)
            except:
                pass
    raise RuntimeError("Failed to get ASP configuration summary for '%s'" % which)
