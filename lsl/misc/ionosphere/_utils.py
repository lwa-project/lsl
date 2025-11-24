import os
import sys
import gzip
import numpy as np
import socket
import warnings
import subprocess
from io import StringIO
from urllib.request import Request, urlopen
from datetime import datetime, timedelta
from ftplib import FTP_TLS, error_perm as FTP_ERROR_PERM, error_temp as FTP_ERROR_TEMP

from lsl.common.data_access import download_file
from lsl.common.mcs import mjdmpm_to_datetime, datetime_to_mjdmpm
from lsl.common.progress import DownloadBar
from lsl.common.color import colorfy
from lsl.misc.file_cache import FileCache, MemoryCache

from lsl.config import LSL_CONFIG
IONO_CONFIG = LSL_CONFIG.view('ionosphere')
DOWN_CONFIG = LSL_CONFIG.view('download')


__version__ = '0.1'
__all__ = ['get_cache_dir', 'download_worker', 'load_mjd']


# Create the cache directory
try:
    _CACHE_DIR = FileCache(os.path.join(LSL_CONFIG.dirname, 'ionospheric_cache'),
                           max_size=lambda: IONO_CONFIG.get('max_cache_size'))
except OSError:
    _CACHE_DIR = MemoryCache(max_size=lambda: IONO_CONFIG.get('max_cache_size'))
    warnings.warn(colorfy("{{%yellow Cannot create or write to on-disk data cache, using in-memory data cache}}"), RuntimeWarning)


def get_cache_dir():
    """
    Return the current instance of FileCache/MemoryCache used for data file
    caching.
    """
    
    return _CACHE_DIR


def _convert_to_gzip(filename):
    """
    Given a unix compressed .Z file, convert it to a gzip .gz file and update
    the cache.
    """
    
    # Load in the file
    with _CACHE_DIR.open(filename, 'rb') as fh:
        cached_filename = fh.name
        uncompressed = subprocess.check_output(['gzip', '-d', '-c', cached_filename])
        
    # Write it back out
    with _CACHE_DIR.open(filename, 'wb') as fh:
        with gzip.GzipFile(fileobj=fh, mode='wb') as gh:
            gh.write(uncompressed)


def _download_worker_cddis(url, filename):
    """
    Download the URL from gdc.cddis.eosdis.nasa.gov via FTP-SSL and save it to a file.
    """
    
    is_interactive = sys.__stdin__.isatty()
    
    # Attempt to download the data
    print(f"Downloading {url}")
    ## Login
    ftps = FTP_TLS("gdc.cddis.eosdis.nasa.gov", timeout=DOWN_CONFIG.get('timeout'))
    status = ftps.login("anonymous", "lwa@unm.edu")
    if not status.startswith("230"):
        ftps.close()
        return False
        
    ## Secure
    status = ftps.prot_p()
    if not status.startswith("200"):
        ftps.close()
        return False
        
    ## Download
    remote_path = url.split("gdc.cddis.eosdis.nasa.gov", 1)[1]
    try:
        remote_size = ftps.size(remote_path)
    except (FTP_ERROR_TEMP, FTP_ERROR_PERM):
        ftps.close()
        return False
        
    with _CACHE_DIR.open(filename, 'wb') as fh:
        pbar = DownloadBar(max=remote_size)
        def write(data):
            fh.write(data)
            pbar.inc(len(data))
            if is_interactive:
                sys.stdout.write(pbar.show()+'\r')
                sys.stdout.flush()
                
        try:
            status = ftps.retrbinary(f"RETR {remote_path}", write, blocksize=DOWN_CONFIG.get('block_size'))
            if is_interactive:
                sys.stdout.write(pbar.show()+'\n')
                sys.stdout.flush()
        except (FTP_ERROR_TEMP, FTP_ERROR_PERM):
            status = 'FAILED'
            
    if not status.startswith("226"):
        _CACHE_DIR.remove(filename)
        ftps.close()
        return False
        
    ## Further processing, if needed
    if os.path.splitext(url)[1] == '.Z':
        ## Save it to a regular gzip'd file after uncompressing it.
        _convert_to_gzip(filename)
        
    # Done
    ftps.close()
    return True


def _download_worker_standard(url, filename):
    """
    Download the URL and save it to a file.
    """
    
    status = False
    try:
        # Attempt to download the data
        with _CACHE_DIR.open(filename, 'wb') as fh:
            _, received, mtime = download_file(url, fh)
            
        # Did we get anything or, at least, enough of something like it looks like 
        # a real file?
        if received < 3:
            raise RuntimeError("Received too few bytes")
            
        # Further processing, if needed
        if os.path.splitext(url)[1] == '.Z':
            ## Save it to a regular gzip'd file after uncompressing it.
            _convert_to_gzip(filename)
            
        status = True
        
    except RuntimeError:
        # Failed to download, remove the empty file
        _CACHE_DIR.remove(filename)
        
    return status


def download_worker(url, filename):
    """
    Download the URL and save it to a file.
    """
    
    # Attempt to download the data
    if url.find('gdc.cddis.eosdis.nasa.gov') != -1:
        status = _download_worker_cddis(url, filename)
    else:
        status = _download_worker_standard(url, filename)
        
    return status


def _parse_tec_map(filename_or_fh):
    """
    Given the name of a file containing a TEC map from the IGC, parse it 
    and return a dictionary containing the files data.
    
    The dictionary keys are:
     * dates - array of MJD values for each time step in the map
     * lats - 2-D array of latitude values for the maps in degrees
     * lngs - 2-D array of longitude values for the maps in degrees
     * height - height for the ionospheric pierce point in km
     * tec - 3-D array of TEC values in TECU.  The dimensions are time by
             latitude by longitude.
     * rms - 3-D array of TEC RMS values in TECU.  The dimensions are time
             by latitude by longitude.
    """
    
    # Variables to hold the map sequences
    dates = []
    tecMaps = []
    rmsMaps = []
    
    # State control variables to help keep up with where we are
    inMap = False
    inBlock = False
    
    try:
        fh = gzip.GzipFile(filename_or_fh, 'rb')
    except TypeError:
        fh = gzip.GzipFile(fileobj=filename_or_fh, mode='rb')
        
    try:
        for line in fh:
            try:
                line = line.decode('ascii', errors='ignore')
            except AttributeError:
                pass
                
            ## Are we beginning a map?
            line = line.replace('\n', '')
            if line.find('START OF TEC MAP') != -1 or line.find('START OF RMS MAP') != -1:
                inMap = True
                continue
                
            ## Have we just ended a map?
            if line.find('END OF TEC MAP') != -1 or line.find('END OF RMS MAP') != -1:
                if line.find('TEC') != -1:
                    tecMaps.append( cmap )
                else:
                    rmsMaps.append( cmap )
                
                inMap = False
                continue
                
            ## Are we in a map?
            if inMap:
                ## Is this part of the preamble? 
                if line.find('EPOCH OF CURRENT MAP') != -1:
                    ### Parse the date/time string
                    year, month, day, hour, minute, second = line.split(None, 6)[:6]
                    year = int(year)
                    month = int(month)
                    day = int(day)
                    hour = int(hour)
                    minute = int(minute)
                    second = int(second)
                    
                    ### Figure out the MJD
                    try:
                        dt = datetime(year, month, day, hour, minute, second, 0)
                    except ValueError:
                        if hour >= 24:
                            dt = datetime(year, month, day, hour-24, minute, second, 0)
                            dt += timedelta(days=1)
                        else:
                            continue
                    mjd, mpm = datetime_to_mjdmpm(dt)
                    mjd = mjd + mpm/1000.0/3600.0/24.0
                    if mjd not in dates:
                        dates.append( mjd )
                        
                    ### Initialize the map and the coorindates
                    cmap = []
                    lats = []
                    lngs = []
                    
                    continue
                    
                ## Is this a different part of the preamble? 
                elif line.find('LAT/LON1/LON2/DLON/H') != -1:
                    lat = float(line[3:8])
                    lng1 = float(line[8:14])
                    lng2 = float(line[14:20])
                    dlng = float(line[20:26])
                    height = float(line[26:32])
                    
                    cmap.append( [] )
                    lats.append( lat )
                    lngs = list(np.arange(lng1, lng2+dlng, dlng))
                    
                    inBlock = True
                    continue
                    
                ## Process the data block keeping in mind that missing values are stored 
                ## as 9999
                if inBlock:
                    fields = np.array([float(v)/10.0 for v in line.split(None)])
                    fields[np.where( fields == 999.9 )] = np.nan
                    
                    cmap[-1].extend( fields )
                    continue
                    
    finally:
        fh.close()
        
    # Combine everything together
    dates = np.array(dates, dtype=np.float64)
    tec = np.array(tecMaps, dtype=np.float32)
    rms = np.array(rmsMaps, dtype=np.float32)
    lats = np.array(lats, dtype=np.float32)
    lngs = np.array(lngs, dtype=np.float32)
    
    # Do we have a valid RMS map?  If not, make one.
    if rms.size != tec.size:
        rms = tec*0.05
        
    # Make lats and lngs 2-D to match the data
    lngs, lats = np.meshgrid(lngs, lats)
    
    # Build up the output
    output = {'dates': dates, 'lats': lats, 'lngs': lngs, 'height': height, 'tec': tec, 'rms': rms}
    
    # Done
    return output


def load_mjd(mjd, filenameTemplate, filenameAltTemplate, downloader):
    """
    Given an MJD value, filenaming templates, and a function to download
    missing files, load the corresponding TEC map.  If the map is not
    already avaliable on disk, download it with the downloader.
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    
    # Pull out the year and the day-of-year
    year = dt.year
    dayOfYear = int(dt.strftime('%j'), 10)
    
    # Figure out the filenames in order of preference.  We'd rather have
    # final values than rapid values
    filename = filenameTemplate % (dayOfYear, year%100)
    filenameAlt = filenameAltTemplate % (dayOfYear, year%100)
    
    # Is the primary file in the disk cache?
    if filename not in _CACHE_DIR:
        ## Can we download it?
        status = downloader(mjd, type='final')
        if not status:
            ## Nope, now check for the secondary file on disk
            if filenameAlt not in _CACHE_DIR:
                ## Can we download it?
                status = downloader(mjd, type='rapid')
                if status:
                    ### Good, we have the secondary file
                    filename = filenameAlt
            else:
                ### Good, we have the secondary file
                filename = filenameAlt
    else:
        ## Good we have the primary file
        pass
        
    # Parse it
    with _CACHE_DIR.open(filename, 'rb') as fh:
        data_set = _parse_tec_map(fh)
        
    return data_set
