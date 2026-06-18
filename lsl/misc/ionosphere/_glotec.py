import os
import h5py
import json
import tarfile
import numpy as np
from datetime import datetime, timedelta, timezone

from lsl.common.mcs import mjdmpm_to_datetime, datetime_to_mjdmpm
from lsl.misc.ionosphere._utils import get_cache_dir, download_worker, load_mjd as base_load_mjd

from lsl.config import LSL_CONFIG
IONO_CONFIG = LSL_CONFIG.view('ionosphere')


__version__ = '0.2'
__all__ = ['FILENAME_TEMPLATE', 'FILENAME_TEMPLATE_ALT', 'load_mjd']


FILENAME_TEMPLATE = 'GloTEC_TEC_%s.nc'
FILENAME_TEMPLATE_ALT = 'GloTEC_TEC_%s.nc'


_CACHE_DIR = get_cache_dir()


def _download(mjd, type='final'):
    """
    Given an MJD value, download the corresponding GloTEC data products for that
    day.
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    if dt.microsecond > 500000:
        dt = dt.replace(microsecond=0)
        dt += timedelta(seconds=1)
        
    year = dt.year
    month = dt.month
    dateStr = dt.strftime("%Y_%m_%d")
    # Build up the filename
    filename = f"glotec_{dateStr}.tar.gz"
    
    # Start check
    if dateStr < "2025_02_24":
        raise ValueError("Requested MJD is before the GloTEC start date of February 24, 2025")
        
    # Attempt to download the data
    status = download_worker('%s/%04i/%02i/%s' % (IONO_CONFIG.get('glotec_url'), year, month, filename), filename)
    if not status and IONO_CONFIG.get('glotec_mirror') is not None:
        status = download_worker('%s/%04i/%02i/%s' % (IONO_CONFIG.get('glotec_mirror'), year, month, filename), filename)
        
    if status:
        # Pull out the .nc file that we need
        with _CACHE_DIR.open(filename, 'rb') as fh:
            with tarfile.open(fileobj=fh, mode='r:*') as tf:
                h5name = FILENAME_TEMPLATE % dateStr
                h5name = os.path.join("glotec_%s" % dateStr, h5name)
                ti = tf.extractfile(h5name)
                if ti:
                    with _CACHE_DIR.open(os.path.basename(h5name), 'wb') as hh:
                        hh.write(ti.read())
                else:
                    status = False
                    
        # Cleanup the full tarball we've downloaded
        _CACHE_DIR.remove(filename)
        
    return status


def _parse_glotec_map(filename_or_fh):
    """
    Parse an individual TEC map from the GloTEC project.  This returns a four-
    element tuple of:
     * 2-D array of latitude values for the maps in degrees
     * 2-D array of longitude values for the maps in degrees
     * 2-D array of TEC values in TECU.  The dimensions are latitude by 
        longitude
     * 2-D array of TEC RMS values in TECU.  The dimensions are latitude 
        by longitude.
    """
    
    # Open the TEC file for reading
    try:
        fh = open(filename_or_fh, 'r')
        do_close = True
    except TypeError:
        fh = filename_or_fh
        do_close = False
        
    try:
        with h5py.File(fh, 'r') as f:
            dates = f['time'][...]
            lats = f['latitude'][...]
            lngs = f['longitude'][...]
            data = f['TEC'][...]
            quality = f['quality_flag'][...]
            
    finally:
        if do_close:
            fh.close()
            
    # Date conversion
    mjds = []
    for ts in dates:
        dt = datetime.fromtimestamp(ts, tz=timezone.utc)
        mjd, mpm = datetime_to_mjdmpm(dt)
        mjds.append(mjd + mpm/1000.0/86400.)
        
    # Bring it into NumPy
    lats = np.array(lats)
    lngs = np.array(lngs)
    data = np.array(data)

    # GloTEC does not provide an uncertainty/RMS map, so derive one from the
    # per-pixel quality flag (0=no observations ... 5=5 or more observations)
    # scaled by the TEC value
    rdata = (6 - np.array(quality)) / 100.0 * data
    
    # Build the lat/lng grip
    lngs, lats = np.meshgrid(lngs, lats)
    
    # Reverse
    lats = lats[::-1,:]
    lngs = lngs[::-1,:]
    data = data[:,::-1,:]
    rdata = rdata[:,::-1,:]
    
    # Done
    return mjds, lats, lngs, data, rdata


def load_mjd(mjd):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    dateStr = dt.strftime("%Y_%m_%d")
    
    ## Figure out the filename
    filename = FILENAME_TEMPLATE % (dateStr)
    
    ## Is the primary file in the disk cache?
    if filename not in _CACHE_DIR:
        ### Can we download it?
        status = _download(mjd)
        
    else:
        ## Good we have the primary file
        pass
        
    ## Parse it and populate the MJD, TEC, and RMS lists lists (timestamp per item)
    with _CACHE_DIR.open(filename, 'rb') as fh:
        daily_mjds, daily_lats, daily_lngs, daily_tec, daily_rms = _parse_glotec_map(fh)
        
    # list -> np.array
    daily_mjds = np.array(daily_mjds)
    daily_tec = np.array(daily_tec)
    daily_rms = np.array(daily_rms)
    
    # Build the output
    daily_output = {'dates': daily_mjds,
                    'lats': daily_lats, 'lngs': daily_lngs,
                    'tec': daily_tec, 'rms': daily_rms}
    
    return daily_output
