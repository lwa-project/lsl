import json
import numpy as np
from datetime import timedelta

from lsl.common.mcs import mjdmpm_to_datetime, datetime_to_mjdmpm
from lsl.misc.ionosphere._utils import get_cache_dir, download_worker, load_mjd as base_load_mjd

from lsl.config import LSL_CONFIG
IONO_CONFIG = LSL_CONFIG.view('ionosphere')


__version__ = '0.1'
__all__ = ['FILENAME_TEMPLATE', 'FILENAME_TEMPLATE_ALT', 'load_mjd']


FILENAME_TEMPLATE = 'glotec_icao_%sZ.geojson'
FILENAME_TEMPLATE_ALT = 'glotec_icao_%sZ.geojson'


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
    dateStr = dt.strftime("%Y%m%dT%H%M%S")
    # Build up the filename
    filename = f"glotec_icao_{dateStr}Z.geojson"
    
    # Attempt to download the data
    return download_worker('%s/%s' % (IONO_CONFIG.get('glotec_url'), filename), filename)


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
        fh = gzip.GzipFile(filename_or_fh, 'rb')
        do_close = True
    except TypeError:
        fh = gzip.GzipFile(fileobj=filename_or_fh, mode='rb')
        do_close = False
        
    try:
        jdata = json.load(fh)
        
        lngs = [[]]
        lats = [[]]
        data = [[]]
        rdata = [[]]
        for feature in jdata['features']:
            lng, lat = feature['geometry']['coordinates']
            tec, rms = feature['properties']['tec'], feature['properties']['anomaly']
            
            if len(lngs[-1]) and lng < lngs[-1][-1]:
                lats.append([])
                lngs.append([])
                data.append([])
                rdata.append([])
            lats[-1].append(lat)
            lngs[-1].append(lng)
            data[-1].append(tec)
            rdata[-1].append(rms)
            
    finally:
        if do_close:
            fh.close()
            
    # Bring it into NumPy
    lats = np.array(lats)
    lngs = np.array(lngs)
    data = np.array(data)
    rdata = np.array(rdata)
    
    # Reverse
    lats = lats[::-1,:]
    lngs = lngs[::-1,:]
    data = data[::-1,:]
    rdata = rdata[::-1,:]
    
    # Done
    return lats, lngs, data, rdata


def load_mjd(mjd):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    
    # Loop over 10 min intervals in the day
    daily_mjds = []
    daily_tec = []
    daily_rms = []
    
    idt = dt.replace(minute=0, second=0, microsecond=0)
    idt -= timedelta(seconds=5*60)
    for i in range(-5, 1*60+10, 10):
        ## Get the YMDHMS string
        imjd, impm = datetime_to_mjdmpm(idt)
        imjd = imjd + impm/1000./86400
        dateStr = idt.strftime("%Y%m%dT%H%M%S")
        
        ## Figure out the filename
        filename = FILENAME_TEMPLATE % (dateStr)
        filenameComp = filename+'.gz'
        
        ## Is the primary file in the disk cache?
        if filename not in _CACHE_DIR:
            ### Can we download it?
            status = _download(imjd)
            
            ### Compress it
            with _CACHE_DIR.open(filename, 'rb') as fh:
                with _CACHE_DIR.open(filenameComp, 'wb') as gh:
                    with gzip.GzipFile(fileobj=gh, mode='wb') as gh:
                        gh.write(fh.read())
            _CACHE_DIR.remove(filename)
            
        else:
            ## Good we have the primary file
            pass
            
        ## Parse it and add it to the list
        daily_mjds.append(imjd)
        with _CACHE_DIR.open(filenameComp, 'rb') as fh:
            daily_lats, daily_lngs, tec, rms = _parse_glotec_map(fh)
        daily_tec.append(tec)
        daily_rms.append(rms)
        
        ## Update the time
        idt += timedelta(seconds=10*60)
        
    # list -> np.array
    daily_mjds = np.array(daily_mjds)
    daily_tec = np.array(daily_tec)
    daily_rms = np.array(daily_rms)
    
    # Build the output
    daily_output = {'dates': daily_mjds,
                    'lats': daily_lats, 'lngs': daily_lngs,
                    'tec': daily_tec, 'rms': daily_rms}
    
    return daily_output
