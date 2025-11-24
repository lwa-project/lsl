from lsl.common.mcs import mjdmpm_to_datetime
from lsl.misc.ionosphere._utils import download_worker, load_mjd as base_load_mjd

from lsl.config import LSL_CONFIG
IONO_CONFIG = LSL_CONFIG.view('ionosphere')


__version__ = '0.1'
__all__ = ['FILENAME_TEMPLATE', 'FILENAME_TEMPLATE_ALT', 'load_mjd']


FILENAME_TEMPLATE = 'jplg%03i0.%02ii.Z'
FILENAME_TEMPLATE_ALT = 'jprg%03i0.%02ii.Z'


def _download(mjd, type='final', cache_dir=None):
    """
    Given an MJD value, download the corresponding JPL final data product 
    for that day.
    
    .. note::
        By default the "final" product is downloaded.  However, the "rapid" 
        data product may be downloaded if the 'type' keyword is set to 
        "rapid".
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    
    year = dt.year
    dayOfYear = int(dt.strftime('%j'), 10)
    
    # Figure out which file we need to download
    if type == 'final':
        ## Final
        filename = 'jplg%03i0.%02ii.Z' % (dayOfYear, year%100)
    elif type == 'rapid':
        ## Rapid
        filename = 'jprg%03i0.%02ii.Z' % (dayOfYear, year%100)
    else:
        ## ???
        raise ValueError(f"Unknown TEC file type '{type}'")
        
    # Attempt to download the data
    status = download_worker('%s/%04i/%03i/%s' % (IONO_CONFIG.get('jpl_url'), year, dayOfYear, filename), filename)
    if not status and IONO_CONFIG.get('jpl_mirror') is not None:
        status = download_worker('%s/%04i/%03i/%s' % (IONO_CONFIG.get('jpl_mirror'), year, dayOfYear, filename), filename)
    return status


def load_mjd(mjd):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    return base_load_mjd(mjd, FILENAME_TEMPLATE, FILENAME_TEMPLATE_ALT, _download)
