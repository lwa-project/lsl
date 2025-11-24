from lsl.common.mcs import mjdmpm_to_datetime
from lsl.misc.ionosphere._utils import download_worker, load_mjd as base_load_mjd

from lsl.config import LSL_CONFIG
IONO_CONFIG = LSL_CONFIG.view('ionosphere')


__version__ = '0.1'
__all__ = ['FILENAME_TEMPLATE', 'FILENAME_TEMPLATE_ALT', 'load_mjd']


FILENAME_TEMPLATE     = 'codg%03i0.%02ii.Z'
FILENAME_TEMPLATE_ALT = 'codg%03i0.%02ii.Z'


def _download(mjd, type='final'):
    """
    Given an MJD value, download the corresponding CODE final data product 
    for that day.
    
    .. note::
        The 'type' keyword is ignored in the call.  It is included for 
        compatiability with _download_igs().
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    
    year = dt.year
    dayOfYear = int(dt.strftime('%j'), 10)
    
    # Figure out which file we need to download
    filename = 'codg%03i0.%02ii.Z' % (dayOfYear, year%100)
    long_filename = 'COD0OPSFIN_%04i%03i0000_01D_01H_GIM.INX.gz' % (year, dayOfYear)
    
    # Attempt to download the data
    for fname in (long_filename, filename):
        status = download_worker('%s/%04i/%03i/%s' % (IONO_CONFIG.get('code_url'), year, dayOfYear, fname), filename)
        if not status:
            status = download_worker('%s/%04i/%03i/%s' % (IONO_CONFIG.get('code_mirror'), year, dayOfYear, fname), filename)
        if status:
            break
    return status


def load_mjd(mjd):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    return base_load_mjd(mjd, FILENAME_TEMPLATE, FILENAME_TEMPLATE_ALT, _download)
