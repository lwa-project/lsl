import os
import numpy as np
import tarfile
from io import StringIO
from datetime import datetime

from scipy.interpolate import RectBivariateSpline

from lsl.common.mcs import mjdmpm_to_datetime, datetime_to_mjdmpm
from lsl.misc.ionosphere._utils import get_cache_dir, download_worker, load_mjd as base_load_mjd

from lsl.config import LSL_CONFIG
IONO_CONFIG = LSL_CONFIG.view('ionosphere')


__version__ = '0.1'
__all__ = ['FILENAME_TEMPLATE', 'FILENAME_TEMPLATE_ALT', 'load_mjd']


FILENAME_TEMPLATE = '%s_ustec.tar.gz'
FILENAME_TEMPLATE_ALT = '%s_ustec.tar.gz'


_CACHE_DIR = get_cache_dir()


def _download(mjd, type='final'):
    """
    Given an MJD value, download the corresponding USTEC data product for that
    day.
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    
    year = dt.year
    month = dt.month
    dateStr = dt.strftime("%Y%m%d")
    altDateStr = dt.strftime("%Y_%m_%d")
    # Build up the filename
    filename = f"{dateStr}_ustec.tar.gz"
    # New style filename (this should probably be the default)
    alt_filename = f"ustec_{altDateStr}.tar.gz"
    
    # Attempt to download the data
    for fname in (alt_filename, filename):
        status = download_worker('%s/%04i/%02i/%s' % (IONO_CONFIG.get('ustec_url'), year, month, fname), filename)
        if not status and IONO_CONFIG.get('ustec_url') is not None:
            status = download_worker('%s/%04i/%02i/%s' % (IONO_CONFIG.get('ustec_mirror'), year, month, fname), filename)
        if status:
            break
    return status


def _parse_ustec_individual(filename_or_fh, rmsname_or_fh=None):
    """
    Parse an individual TEC map from the USTEC project.  This returns a four-
    element tuple of:
     * 2-D array of latitude values for the maps in degrees
     * 2-D array of longitude values for the maps in degrees
     * 2-D array of TEC values in TECU.  The dimensions are latitude by 
        longitude
     * 2-D array of TEC RMS values in TECU.  The dimensions are latitude 
        by longitude.
    
    Format Reference:
    https://www.ngdc.noaa.gov/stp/iono/ustec/README.html
    """
    
    # Open the TEC file for reading
    try:
        fh = open(filename_or_fh, 'r')
        do_close = True
    except TypeError:
        fh = filename_or_fh
        do_close = False
        
    try:
        # Go!
        inBlock = False
        lats = []
        data = []
        for line in fh:
            if line[0] in ('#', ':'):
                ## Comments
                continue
            elif len(line) < 5:
                ## Blank lines
                continue
                
            ## Start the parsing
            if not inBlock:
                ### The first row consists of a list of longitudes
                fields = line.split()
                lngs = [float(f)/10.0 for f in fields[1:]]
                inBlock = True
                continue
            else:
                ### Slant TEC values are stored at the end of the file
                if line[:3] == '999':
                    inBlock = False
                    break
                    
                ### Before the satellite we have the TEC values, one for each latitude
                fields = line.split()
                lats.append( float(fields[0])/10 )
                data.append( [float(f)/10 for f in fields[1:]] )
                
    finally:
        if do_close:
            fh.close()
            
    # Bring it into NumPy
    lats = np.array(lats)
    lngs = np.array(lngs)
    lngs,lats = np.meshgrid(lngs,lats)
    data = np.array(data)
    
    # Check for an associated RMS file
    if rmsname_or_fh is not None:
        ## Oh good, we have one
        try:
            fh = open(rmsname_or_fh, 'r')
            do_close = True
        except TypeError:
            fh = rmsname_or_fh
            do_close = False
            
        try:
            ## Go! (again)
            inBlock = False
            rlats = []
            rdata = []
            for line in fh:
                if line[0] in ('#', ':'):
                    ## Comments
                    continue
                elif len(line) < 3:
                    ## Blank lines
                    continue
                    
                ## Start the parsing
                if not inBlock:
                    ### The first row consists of a list of longitudes
                    fields = line.split()
                    rlngs = [float(f)/10.0 for f in fields[1:]]
                    inBlock = True
                    continue
                else:
                    ### Slant TEC values are stored at the end of the file
                    if line[:3] == '999':
                        inBlock = False
                        break
                        
                    ### Before the satellite we have the TEC RMS values, one for each latitude
                    fields = line.split()
                    rlats.append( float(fields[0])/10 )
                    rdata.append( [float(f)/10 for f in fields[1:]] )
                    
        finally:
            if do_close:
                fh.close()
                
        # Bring it into NumPy
        rlats = np.array(rlats)
        rlngs = np.array(rlngs)
        rdata = np.array(rdata)
        
        # For some reason the RMS map has a lower resolution than the actual TEC map.
        # Interpolate up the resolution of the actual TEC map so that we have an
        # uncertainty at each point
        interpFunction = RectBivariateSpline(rlats, rlngs, rdata, kx=1, ky=1)
        rms = data*0.0
        for i in range(lats.shape[0]):
            for j in range(lats.shape[0]):
                rms[i,j] = interpFunction(lats[i,j], lngs[i,j])[0,0]
    else:
        ## Sadness, no RMS file found...
        rms = data*0.05
        
    # Reverse
    lats = lats[::-1,:]
    lngs = lngs[::-1,:]
    data = data[::-1,:]
    rms  = rms[::-1,:]
    
    # Done
    return lats, lngs, data, rms


def _parse_ustec_height(filename_or_fh):
    """
    Parse emperical orthonormal functions to come up with an effective 
    height for the ionosphere.
    
    Format Reference:
    https://www.ngdc.noaa.gov/stp/iono/ustec/README.html
    """
    
    # Open the EOF file for reading
    try:
        fh = open(filename_or_fh, 'r')
        do_close = True
    except TypeError:
        fh = filename_or_fh
        do_close = False
        
    try:
        # Go!
        inBlock = False
        heights = []
        data = []
        for line in fh:
            if line[0] in ('#', ':'):
                ## Comments
                continue
            elif len(line) < 3:
                ## Blank lines
                continue
                
            ## Start the parsing
            if not inBlock:
                fields = line.split()
                height = float(fields[2])
                step = float(fields[3])
                inBlock = True
            else:
                fields = line.split()
                heights.append( height - 6371)
                height += step
                data.append( float(fields[0]) )
                
    finally:
        if do_close:
            fh.close()
            
    # Bring it into Numpy and get an average height
    heights = np.array(heights)
    data = np.array(data)
    height = (heights*data).sum() / data.sum()
    
    # Done
    return height


def _parse_ustec_map(filename_or_fh):
    """
    Given the name of a file containing a TEC map from the USTEC project, 
    parse it and return a dictionary containing the files data.
    
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
    
    try:
        tf = tarfile.open(filename_or_fh, 'r:*')
        do_close = True
    except TypeError:
        tf = tarfile.open(fileobj=filename_or_fh, mode='r:*')
        do_close = False
        
    valid_filename = lambda x: ((x.find('_TEC.txt') != -1) \
                                or (x.find('_ERR.txt') != -1) \
                                or (x.find('_EOF.txt') != -1))
    try:
        tecFiles = {}
        errFiles = {}
        eofFiles = {}
        for entry in tf:
            if not valid_filename(entry.name):
                continue
                
            contents = tf.extractfile(entry.name).read()
            try:
                contents = contents.decode()
            except AttributeError:
                # Already a string
                pass
            contents = StringIO(contents)
            
            if entry.name.find('_TEC.txt') != -1:
                tecFiles[entry.name] = contents
            elif entry.name.find('_ERR.txt') != -1:
                errFiles[entry.name] = contents
            elif entry.name.find('_EOF.txt') != -1:
                eofFiles[entry.name] = contents
                
        # Variables to hold the map sequences
        dates = []
        tecMaps = []
        rmsMaps = []
        
        # Get all of the TEC map files and load them in
        tecfilenames = list(tecFiles.keys())
        tecfilenames.sort()
        for tecfilename in tecfilenames:
            tecfh = tecFiles[tecfilename]
            try:
                rmsfilename = tecfilename.replace('_TEC', '_ERR')
                rmsfh = errFiles[rmsfilename]
            except KeyError:
                rmsfh = None
            lats, lngs, tec, rms = _parse_ustec_individual(tecfh, rmsname_or_fh=rmsfh)
            
            ### Figure out the MJD
            dt = os.path.basename(tecfilename).split('_')[0]
            dt = datetime.strptime(dt, "%Y%m%d%H%M")
            mjd, mpm = datetime_to_mjdmpm(dt)
            mjd = mjd + mpm/1000.0/3600.0/24.0
            if mjd not in dates:
                dates.append( mjd )
                
            # Stack on the new TEC and RMS maps
            tecMaps.append( tec )
            rmsMaps.append( rms )
                
            #except Exception as e:
                #pass
                
        # Get the mean ionospheric height
        eoffilename = list(eofFiles.keys())[0]
        #try:
        height = _parse_ustec_height(eofFiles[eoffilename])
        #except:
        #	height = 450
        
    finally:
        if do_close:
            tf.close()
            
    # Combine everything together
    dates = np.array(dates, dtype=np.float64)
    tec = np.array(tecMaps, dtype=np.float32)
    rms = np.array(rmsMaps, dtype=np.float32)
    
    # Build up the output
    output = {'dates': dates, 'lats': lats, 'lngs': lngs, 'height': height, 'tec': tec, 'rms': rms}
    
    # Done
    return output


def load_mjd(mjd):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm_to_datetime(int(mjd), mpm)
    
    # Pull out a YMD string
    dateStr = dt.strftime("%Y%m%d")
    
    # Figure out the filenames in order of preference.  We'd rather have
    # final values than rapid values
    filename = FILENAME_TEMPLATE % (dateStr)
    
    # Is the primary file in the disk cache?
    if filename not in _CACHE_DIR:
        ## Can we download it?
        status = _download(mjd)
        
    else:
        ## Good we have the primary file
        pass
        
    # Parse it
    with _CACHE_DIR.open(filename, 'rb') as fh:
        data_set = _parse_ustec_map(fh)
        
    return data_set
