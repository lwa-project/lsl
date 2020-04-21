"""
A collection of utilities for retrieving parameters that may be relevant 
for ionospheric corrections.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import glob
import gzip
import numpy
import socket
import shutil
import tarfile
import tempfile
import warnings
import subprocess
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
from datetime import datetime, timedelta

from astropy.constants import R_earth

from scipy.special import lpmv
try:
    from scipy.misc import factorial
except ImportError:
    from scipy.special import factorial
from scipy.optimize import fmin
from scipy.interpolate import RectBivariateSpline

from lsl.common.stations import geo2ecef
from lsl.common.paths import data as dataPath
from lsl.common.mcs import mjdmpm2datetime, datetime2mjdmpm

from lsl.misc import telemetry
telemetry.track_module()


__version__ = "0.5"
__all__ = ['get_magnetic_field', 'compute_magnetic_declination', 'compute_magnetic_inclination', 
           'get_tec_value', 'get_ionospheric_pierce_point']


# Create the cache directory
if not os.path.exists(os.path.join(os.path.expanduser('~'), '.lsl')):
    os.mkdir(os.path.join(os.path.expanduser('~'), '.lsl'))
_CACHE_DIR = os.path.join(os.path.expanduser('~'), '.lsl', 'ionospheric_cache')
if not os.path.exists(_CACHE_DIR):
    os.mkdir(_CACHE_DIR)


# Create the on-line cache
_CACHE = {}

# Radius of the Earth in meters for the IGRF
_RADIUS_EARTH = R_earth.to('m').value


def _load_igrf(filename):
    """
    Given a filename pointing to a list of IGRF coefficients, load in the 
    data and return a dictionary containing the raw coefficients.
    
    The dictionary keys are:
     * years - list of years for each of the models
     * g - dictionary of cosine term coefficients
     * h - dictionary of sine term coefficients
    
    The g and h dictionaries are keyed off the degree of the Legendre 
    function and each stores a list of n orders.  Each order is composed
    of len(years)+1 values, one for each each plus a secular evolution
    term.
    """
    
    # Open the file
    with open(filename, 'r') as fh:
        # Go!
        dataCos = {}
        dataSin = {}
        for line in fh:
            ## Is this line a comment?
            line = line.replace('\n', '')
            if line.find('#') != -1:
                continue
                
            ## Is it a header?
            if line.find('IGRF') != -1:
                continue
            if line.find('g/h') != -1:
                fields = line.split(None)
                years = [float(value) for value in fields[3:-1]]
                continue
                
            ## Must be data...  parse it
            fields = line.split(None)
            t, n, m = fields[0], int(fields[1]), int(fields[2])
            c = numpy.array([float(v) for v in fields[3:]])
            
            ## Sort out cosine (g) vs. sine (h)
            if t == 'g':
                try:
                    dataCos[n][m] = c
                except KeyError:
                    dataCos[n] = [numpy.zeros(len(years)+1) for i in range(n+1)]
                    dataCos[n][m] = c
            else:
                try:
                    dataSin[n][m] = c
                except KeyError:
                    dataSin[n] = [numpy.zeros(len(years)+1) for i in range(n+1)]
                    dataSin[n][m] = c
                    
    # Build the output
    output = {'years': years, 'g': dataCos, 'h': dataSin}
    
    # Done
    return output


def _compute_igrf_coefficents(year, coeffs):
    """
    Given a decimal year and a coefficient dictionary from _parseIGRFModel(),
    compute the actual coefficients for that epoch and return a dictionary
    containing the cosine and sine coefficients for that year.
    
    The dictionary keys are:
    * year - year the coefficients are computed for
    * g - dictionary of cosine term coefficients
    * h - dictionary of sine term coefficients
    
    The g and h dictionaries are keyed off the degree of the Legendre 
    function and each stores a list of n orders.
    """
    
    # Figure out the closest model point(s) to the requested year taking into
    # account that a new model comes out every five years
    best = numpy.where( numpy.abs(year - numpy.array(coeffs['years'])) < 5 )[0]
    
    if year < min(coeffs['years']):
        # If the requested year is before 1900 we can't do anything
        raise RuntimeError("Invalid year %i" % year)
    else:
        # Otherwise, figure out the coefficients using a simple linear interpolation 
        # or extrapolation using the secular changes in the model
        coeffsCos = {}
        coeffsSin = {}
        
        # Loop over the degrees
        for n in coeffs['g'].keys():
            ## Loop over the orders
            for m in range(0, n+1):
                if year > max(coeffs['years']):
                    ### If we are beyond the last year in the model, use the secular changes
                    slope = coeffs['g'][n][m][-1]
                else:
                    ### Otherwise, do a linear interpolation between the two nearest models
                    slope = coeffs['g'][n][m][best[-1]] - coeffs['g'][n][m][best[0]]
                    if best[-1] == best[0]:
                        slope = 0.0
                    else:
                        slope /= coeffs['years'][best[-1]] - coeffs['years'][best[0]]
                        
                ### Compute the cosine terms
                try:
                    coeffsCos[n][m] = slope*(year - coeffs['years'][best[0]]) + coeffs['g'][n][m][best[0]]
                except:
                    coeffsCos[n] = [0.0 for i in range(n+1)]
                    coeffsCos[n][m] = slope*(year - coeffs['years'][best[0]]) + coeffs['g'][n][m][best[0]]
                    
                if year > max(coeffs['years']):
                    ### If we are beyond the last year in the model, use the secular changes
                    slope = coeffs['h'][n][m][-1]
                else:
                    ### Otherwise, do a linear interpolation between the two nearest models
                    slope = coeffs['h'][n][m][best[-1]] - coeffs['h'][n][m][best[0]]
                    if best[-1] == best[0]:
                        slope = 0.0
                    else:
                        slope /= coeffs['years'][best[-1]] - coeffs['years'][best[0]]
                        
                ### Compute the sine terms
                try:
                    coeffsSin[n][m] = slope*(year - coeffs['years'][best[0]]) + coeffs['h'][n][m][best[0]]
                except:
                    coeffsSin[n] = [0.0 for i in range(n+1)]
                    coeffsSin[n][m] = slope*(year - coeffs['years'][best[0]]) + coeffs['h'][n][m][best[0]]
                    
    # Build the output
    output = {'year': year, 'g': coeffsCos, 'h': coeffsSin}
    
    # Done
    return output


def _Snm(n, m):
    """
    Compute the factor needed to convert an unnormalized associated Legendre
    function of degree n and order m into a Schmidt quasi-normalized
    associated Legendre function.
    """
    
    if m == 0:
        return numpy.sqrt(factorial(n-m)/factorial(n+m))
    else:
        return numpy.sqrt(2.0*factorial(n-m)/factorial(n+m))


def _Pnm(n, m, mu):
    """
    Compute the value of an unnormalized associated Legendre function of
    degree n and order m at mu following the convention of Abramowitz 
    and Stegun (1972).
    """
    
    return (-1)**m*lpmv(m, n, mu)



def _dPnm(n, m, mu):
    """
    Compute d Pnm(cos(theta)) / d theta for an unnormalized associated 
    Legendre function of degree n and order m at a cos(theta) of mu.
    """
    
    o = n*mu*_Pnm(n, m, mu) - (n+m)*_Pnm(n-1, m, mu)
    o /= numpy.sqrt(1.0 - mu**2)
    return o


def get_magnetic_field(lat, lng, elev, mjd=None, ecef=False):
    """
    Given a geodetic location described by a latitude in degrees (North 
    positive), a longitude in degrees (West negative), an elevation 
    in meters and an MJD value, compute the Earth's magnetic field in that 
    location and return a three-element tuple of the magnetic field's 
    components in nT.  By default these are in topocentric coordinates of
    (North, East, Up).  To return values in ECEF, set the 'ecef' keyword to
    True.  If the MJD file is None, the current time is used.
    
    .. note::
        The convention used for the topocentric coordinates differs
        from what the IGRF uses in the sense that the zenith direction
        points up rather than down.
    """
    
    # Get the current time if mjd is None
    if mjd is None:
        mjd, mpm = datetime2mjdmpm( datetime.utcnow() )
        mjd = mjd + mpm/1000.0/3600.0/24.0
        
    # Convert the MJD to a decimal year.  This is a bit tricky
    ## Break the MJD into an integer MJD and an MPM in order to build a datetime instance
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000.0)
    mjd0 = mjdmpm2datetime(int(mjd), mpm)
    ## Convert the datetime instance to January 1
    mjd0 = mjd0.replace(month=1, day=1, hour=0, second=0, microsecond=0)
    ## Figure out January 1 for the following year
    mjd1 = mjd0.replace(year=mjd0.year+1)
    ## Figure out how long the year is in days
    diffDays = mjd1-mjd0
    diffDays = diffDays.days + diffDays.seconds/86400.0 + diffDays.microseconds/1e6/86400.0
    ## Convert the January 1 date back to an MJD
    mjd0, mpm0 = datetime2mjdmpm(mjd0)
    mjd0 = mjd0 + mpm/1000.0/3600.0/24.0
    year = (mjd1.year - 1) + (mjd - mjd0) / diffDays
    
    # Convert the geodetic position provided to a geocentric one for calculation
    ## Deal with the poles
    if 90.0 - lat < 0.001:
        xyz = numpy.array(geo2ecef(89.999*numpy.pi/180, lng*numpy.pi/180, elev))
    elif 90.0 + lat < 0.001:
        xyz = numpy.array(geo2ecef(-89.999*numpy.pi/180, lng*numpy.pi/180, elev))
    else:
        xyz = numpy.array(geo2ecef(lat*numpy.pi/180, lng*numpy.pi/180, elev))
    ## To geocentric
    r = numpy.sqrt( (xyz**2).sum() )
    lt = numpy.arcsin(xyz[2]/r)
    ln = numpy.arctan2(xyz[1], xyz[0])
    
    # Load in the coefficients
    try:
        coeffs = _CACHE['IGRF']
    except KeyError:
        filename = os.path.join(dataPath, 'igrf13coeffs.txt')
        _CACHE['IGRF'] = _load_igrf(filename)
        
        coeffs = _CACHE['IGRF']
        
    # Compute the coefficients for the epoch
    coeffs = _compute_igrf_coefficents(year, coeffs)
    
    # Compute the field strength in spherical coordinates
    Br, Bth, Bph = 0.0, 0.0, 0.0
    for n in coeffs['g'].keys():
        for m in range(0, n+1):
            Br  += (n+1.0)*(_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['g'][n][m]*numpy.cos(m*ln) * _Pnm(n, m, numpy.sin(lt))
            Br  += (n+1.0)*(_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['h'][n][m]*numpy.sin(m*ln) * _Pnm(n, m, numpy.sin(lt))
            
            Bth -= (_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['g'][n][m]*numpy.cos(m*ln) * _dPnm(n, m, numpy.sin(lt))
            Bth -= (_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['h'][n][m]*numpy.sin(m*ln) * _dPnm(n, m, numpy.sin(lt))
            
            Bph += (_RADIUS_EARTH/r)**(n+2)/numpy.cos(lt) * _Snm(n,m)*coeffs['g'][n][m]*m*numpy.sin(m*ln) * _Pnm(n, m, numpy.sin(lt))
            Bph -= (_RADIUS_EARTH/r)**(n+2)/numpy.cos(lt) * _Snm(n,m)*coeffs['h'][n][m]*m*numpy.cos(m*ln) * _Pnm(n, m, numpy.sin(lt))
    ## And deal with NaNs
    if Br != Br:
        Br = 0.0
    if Bth != Bth:
        Bth = 0.0
    if Bph != Bph:
        Bph = 0.0
        
    # Convert from spherical to ECEF
    Bx = Br*numpy.cos(lt)*numpy.cos(ln) + Bth*numpy.sin(lt)*numpy.cos(ln) - Bph*numpy.sin(ln)
    By = Br*numpy.cos(lt)*numpy.sin(ln) + Bth*numpy.sin(lt)*numpy.sin(ln) + Bph*numpy.cos(ln)
    Bz = Br*numpy.sin(lt) - Bth*numpy.cos(lt)
    
    # Are we done?
    if ecef:
        # For ECEF we don't need to do anything else
        outputField = Bx, By, Bz
        
    else:
        # Convert from ECEF to topocentric (geodetic)
        ## Update the coordinates for geodetic
        lt = lat*numpy.pi/180.0
        if 90.0 - lat < 0.001:
            lt = 89.999*numpy.pi/180.0
        elif 90.0 + lat < 0.001:
            lt = -89.999*numpy.pi/180.0
        else:
            lt = lat*numpy.pi/180.0
        ln = lng*numpy.pi/180.0
        
        ## Build the rotation matrix for ECEF to SEZ
        rot = numpy.array([[ numpy.sin(lt)*numpy.cos(ln), numpy.sin(lt)*numpy.sin(ln), -numpy.cos(lt)], 
                    [-numpy.sin(ln),               numpy.cos(ln),                0            ],
                    [ numpy.cos(lt)*numpy.cos(ln), numpy.cos(lt)*numpy.sin(ln),  numpy.sin(lt)]])
                    
        ## Apply and extract
        sez = numpy.dot(rot, numpy.array([Bx,By,Bz]))
        Bn, Be, Bz = -sez[0], sez[1], sez[2]
        
        outputField = Bn, Be, Bz
    
    # Done
    return outputField


def compute_magnetic_declination(Bn, Be, Bz):
    """
    Given the topocentric output of get_magnetic_field(), compute and return 
    the magnetic declination (deviation between magnetic north and true 
    north) in degrees.
    
    .. note::
        The declination is poorly defined (NaN) around the magnetic poles
        where the horizontal field strength is less than 100 nT.
    """
    
    # Compute the horizontal field strength
    Bh = numpy.sqrt(Bn**2+Be**2)
    
    # Compute the declination
    decl = 2.0*numpy.arctan2(Be, Bh+Bn)
    
    # Check for bounds
    if Bh < 100.0:
        decl = numpy.nan
        
    # Convert to degrees and done
    return decl*180.0/numpy.pi


def compute_magnetic_inclination(Bn, Be, Bz):
    """
    Given the topocentric output of get_magnetic_field(), compute and return 
    the magnetic inclination or magnetic dip (angle between the magnetic 
    field and the horizontal) in degrees.
    """
    
    # Compute the horizontal field strength
    Bh = numpy.sqrt(Bn**2+Be**2)
    
    # Compute the inclination.  This has an extra negative sign because of the
    # convention used in get_magnetic_field().
    incl = numpy.arctan2(-Bz, Bh)
    
    # Convert to degrees and done
    return incl*180.0/numpy.pi


def _download_igs(mjd, base_url='ftp://cddis.gsfc.nasa.gov/gps/products/ionex/', timeout=120, type='final'):
    """
    Given an MJD value, download the corresponding IGS final data product 
    for that day.
    
    .. note::
        By default the "final" product is downloaded.  However, the "rapid" 
        data product may be downloaded if the 'type' keyword is set to 
        "rapid".
    """
    
    # Convert the MJD to a datetime instance so that we can pull out the year
    # and the day-of-year
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
    dt = mjdmpm2datetime(int(mjd), mpm)
    
    year = dt.year
    dayOfYear = int(dt.strftime('%j'), 10)
    
    # Figure out which file we need to download
    if type == 'final':
        ## Final
        filename = 'igsg%03i0.%02ii.Z' % (dayOfYear, year%100)
    elif type == 'rapid':
        ## Rapid
        filename = 'igrg%03i0.%02ii.Z' % (dayOfYear, year%100)
    else:
        ## ???
        raise ValueError("Unknown TEC file type '%s'" % type)
        
    # Attempt to download the data
    try:
        tecFH = urlopen('%s/%04i/%03i/%s' % (base_url, year, dayOfYear, filename), timeout=timeout)
        data = tecFH.read()
        tecFH.close()
    except IOError as e:
        warnings.warn('Error downloading file for %i, %i: %s' % (dayOfYear, year, str(e)), RuntimeWarning)
        data = ''
    except socket.timeout:
        data = ''
        
    # Did we get anything?
    if len(data) == 0:
        ## Fail
        return False
    else:
        ## Success!  Save it to a file and then decompress it with 'gunzip' 
        ## since I can't figure out how to decompress this in Python.
        fh = open(os.path.join(_CACHE_DIR, filename), 'wb')
        fh.write(data)
        fh.close()
        
        subprocess.check_call(['gunzip', '-f', os.path.join(_CACHE_DIR, filename)])
        subprocess.check_call(['gzip', os.path.join(_CACHE_DIR, os.path.splitext(filename)[0])])
        
        return True


def _download_jpl(mjd, base_url='ftp://cddis.gsfc.nasa.gov/gps/products/ionex/', timeout=120, type='final'):
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
    dt = mjdmpm2datetime(int(mjd), mpm)
    
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
        raise ValueError("Unknown TEC file type '%s'" % type)
        
    # Attempt to download the data
    try:
        tecFH = urlopen('%s/%04i/%03i/%s' % (base_url, year, dayOfYear, filename), timeout=timeout)
        data = tecFH.read()
        tecFH.close()
    except IOError as e:
        warnings.warn('Error downloading file for %i, %i: %s' % (dayOfYear, year, str(e)), RuntimeWarning)
        data = ''
    except socket.timeout:
        data = ''
        
    # Did we get anything?
    if len(data) == 0:
        ## Fail
        return False
    else:
        ## Success!  Save it to a file and then decompress it with 'gunzip' 
        ## since I can't figure out how to decompress this in Python.
        fh = open(os.path.join(_CACHE_DIR, filename), 'wb')
        fh.write(data)
        fh.close()
        
        subprocess.check_call(['gunzip', '-f', os.path.join(_CACHE_DIR, filename)])
        subprocess.check_call(['gzip', os.path.join(_CACHE_DIR, os.path.splitext(filename)[0])])
        
        return True


def _download_uqr(mjd, base_url='ftp://cddis.gsfc.nasa.gov/gps/products/ionex/', timeout=120, type='final'):
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
    dt = mjdmpm2datetime(int(mjd), mpm)
    
    year = dt.year
    dayOfYear = int(dt.strftime('%j'), 10)
    
    # Figure out which file we need to download
    if type == 'final':
        ## Final
        filename = 'uqrg%03i0.%02ii.Z' % (dayOfYear, year%100)
    elif type == 'rapid':
        ## Rapid
        filename = 'uqrg%03i0.%02ii.Z' % (dayOfYear, year%100)
    else:
        ## ???
        raise ValueError("Unknown TEC file type '%s'" % type)
        
    # Attempt to download the data
    try:
        tecFH = urlopen('%s/%04i/%03i/%s' % (base_url, year, dayOfYear, filename), timeout=timeout)
        data = tecFH.read()
        tecFH.close()
    except IOError as e:
        warnings.warn('Error downloading file for %i, %i: %s' % (dayOfYear, year, str(e)), RuntimeWarning)
        data = ''
    except socket.timeout:
        data = ''
        
    # Did we get anything?
    if len(data) == 0:
        ## Fail
        return False
    else:
        ## Success!  Save it to a file and then decompress it with 'gunzip' 
        ## since I can't figure out how to decompress this in Python.
        fh = open(os.path.join(_CACHE_DIR, filename), 'wb')
        fh.write(data)
        fh.close()
        
        subprocess.check_call(['gunzip', '-f', os.path.join(_CACHE_DIR, filename)])
        subprocess.check_call(['gzip', os.path.join(_CACHE_DIR, os.path.splitext(filename)[0])])
        
        return True


def _download_code(mjd, base_url='ftp://ftp.aiub.unibe.ch/CODE/', timeout=120, type='final'):
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
    dt = mjdmpm2datetime(int(mjd), mpm)
    
    year = dt.year
    dayOfYear = int(dt.strftime('%j'), 10)
    
    # Figure out which file we need to download
    filename = 'CODG%03i0.%02iI.Z' % (dayOfYear, year%100)
    
    # Attempt to download the data
    try:
        tecFH = urlopen('%s/%04i/%s' % (base_url, year, filename), timeout=timeout)
        data = tecFH.read()
        tecFH.close()
    except IOError as e:
        warnings.warn('Error downloading file for %i, %i: %s' % (dayOfYear, year, str(e)), RuntimeWarning)
        data = ''
    except socket.timeout:
        data = ''
        
    # Did we get anything?
    if len(data) == 0:
        ## Fail
        return False
    else:
        ## Success!  Save it to a file and then decompress it with 'gunzip' 
        ## since I can't figure out how to decompress this in Python.
        fh = open(os.path.join(_CACHE_DIR, filename), 'wb')
        fh.write(data)
        fh.close()
        
        subprocess.check_call(['gunzip', '-f', os.path.join(_CACHE_DIR, filename)])
        subprocess.check_call(['gzip', os.path.join(_CACHE_DIR, os.path.splitext(filename)[0])])
        
        return True


def _download_ustec(mjd, base_url='http://www.ngdc.noaa.gov/stp/iono/ustec/products/', timeout=120):
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
    dt = mjdmpm2datetime(int(mjd), mpm)
    
    year = dt.year
    month = dt.month
    dateStr = dt.strftime("%Y%m%d")
    # Build up the filename
    filename = '%s_ustec.tar.gz' % dateStr
    
    # Attempt to download the data
    try:
        tecFH = urlopen('%s/%04i/%02i/%s' % (base_url, year, month, filename), timeout=timeout)
        data = tecFH.read()
        tecFH.close()
    except IOError as e:
        warnings.warn('Error downloading file for %s: %s' % (dateStr, str(e)), RuntimeWarning)
        data = ''
    except socket.timeout:
        data = ''
        
    # Did we get anything?
    if len(data) == 0:
        ## Fail
        return False
    else:
        ## Success!  Save it to a file and then decompress it with 'gunzip' 
        ## since I can't figure out how to decompress this in Python.
        fh = open(os.path.join(_CACHE_DIR, filename), 'wb')
        fh.write(data)
        fh.close()
        
        return True


def _parse_tec_map(filename):
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
    
    # Go
    with gzip.open(filename, 'r') as fh:
        for line in fh:
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
                    mjd, mpm = datetime2mjdmpm(dt)
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
                    lngs = list(numpy.arange(lng1, lng2+dlng, dlng))
                    
                    inBlock = True
                    continue
                    
                ## Process the data block keeping in mind that missing values are stored 
                ## as 9999
                if inBlock:
                    fields = numpy.array([float(v)/10.0 for v in line.split(None)])
                    fields[numpy.where( fields == 999.9 )] = numpy.nan
                    
                    cmap[-1].extend( fields )
                    continue
                    
    # Combine everything together
    dates = numpy.array(dates, dtype=numpy.float64)
    tec = numpy.array(tecMaps, dtype=numpy.float32)
    rms = numpy.array(rmsMaps, dtype=numpy.float32)
    lats = numpy.array(lats, dtype=numpy.float32)
    lngs = numpy.array(lngs, dtype=numpy.float32)
    
    # Do we have a valid RMS map?  If not, make one.
    if rms.size != tec.size:
        rms = tec*0.05
        
    # Make lats and lngs 2-D to match the data
    lngs, lats = numpy.meshgrid(lngs, lats)
    
    # Build up the output
    output = {'dates': dates, 'lats': lats, 'lngs': lngs, 'height': height, 'tec': tec, 'rms': rms}
    
    # Done
    return output


def _parse_ustec_individual(filename):
    """
    Parse an individual TEC map from the USTEC project.  This returns a five-
    element tuple of:
     * datetime for the start of the map
     * 2-D array of latitude values for the maps in degrees
     * 2-D array of longitude values for the maps in degrees
     * 2-D array of TEC values in TECU.  The dimensions are latitude by 
        longitude
     * 2-D array of TEC RMS values in TECU.  The dimensions are latitude 
        by longitude.
    
    Format Reference:
    https://www.ngdc.noaa.gov/stp/iono/ustec/README.html
    """
    
    # Get the date/time from the filename
    dt = os.path.basename(filename).split('_')[0]
    dt = datetime.strptime(dt, "%Y%m%d%H%M")
    
    # Open the TEC file for reading
    with open(filename, 'r') as fh:
        # Go!
        inBlock = False
        lats = []
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
                
    # Bring it into NumPy
    lats = numpy.array(lats)
    lngs = numpy.array(lngs)
    lngs,lats = numpy.meshgrid(lngs,lats)
    data = numpy.array(data)
    
    # Check for an associated RMS file
    rmsname = filename.replace('_TEC', '_ERR')
    if os.path.exists(rmsname):
        ## Oh good, we have one
        with open(rmsname, 'r') as fh:
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
                    
        # Bring it into NumPy
        rlats = numpy.array(rlats)
        rlngs = numpy.array(rlngs)
        rdata = numpy.array(rdata)
        
        # For some reason the RMS map has a lower resolution than the actual TEC map.
        # Interpolate up the the resolution of the actual TEC map so that we have an
        # uncertainty at each point
        interpFunction = RectBivariateSpline(rlats, rlngs, rdata, kx=1, ky=1)
        rms = data*0.0
        for i in range(lats.shape[0]):
            for j in range(lats.shape[0]):
                rms[i,j] = interpFunction(lats[i,j], lngs[i,j])
    else:
        ## Sadness, no RMS file found...
        rms = data*0.05
        
    # Reverse
    lats = lats[::-1,:]
    lngs = lngs[::-1,:]
    data = data[::-1,:]
    rms  = rms[::-1,:]
    
    # Done
    return dt, lats, lngs, data, rms


def _parse_ustec_height(filename):
    """
    Parse emperical orthonormal functions to come up with an effective 
    height for the ionosphere.
    
    Format Reference:
    https://www.ngdc.noaa.gov/stp/iono/ustec/README.html
    """
    
    # Open the EOF file for reading
    with open(filename, 'r') as fh:
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
                
    # Bring it into Numpy and get an average height
    heights = numpy.array(heights)
    data = numpy.array(data)
    height = (heights*data).sum() / data.sum()
    
    # Done
    return height


def _parse_ustec_map(filename):
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
    
    tempDir = tempfile.mkdtemp(prefix='ionosphere-')
    
    tf = tarfile.open(filename, 'r:*')
    tecFiles = [tio for tio in tf.getmembers() if tio.name.find('_TEC.txt') != -1]
    errFiles = [tio for tio in tf.getmembers() if tio.name.find('_ERR.txt') != -1]
    eofFiles = [tio for tio in tf.getmembers() if tio.name.find('_EOF.txt') != -1]
    tf.extractall(path=tempDir, members=tecFiles)
    tf.extractall(path=tempDir, members=errFiles)
    tf.extractall(path=tempDir, members=eofFiles)
    
    # Variables to hold the map sequences
    dates = []
    tecMaps = []
    rmsMaps = []
    
    # Get all of the TEC map files and load them in
    tecfilenames = glob.glob(os.path.join(tempDir, '*_TEC.txt'))
    tecfilenames.sort()
    for tecfilename in tecfilenames:
        #try:
        dt, lats, lngs, tec, rms = _parse_ustec_individual(tecfilename)
        
        ### Figure out the MJD
        mjd, mpm = datetime2mjdmpm(dt)
        mjd = mjd + mpm/1000.0/3600.0/24.0
        if mjd not in dates:
            dates.append( mjd )
            
        # Stack on the new TEC and RMS maps
        tecMaps.append( tec )
        rmsMaps.append( rms )
            
        #except Exception as e:
            #pass
            
    # Get the mean ionospheric height
    eoffilename = glob.glob(os.path.join(tempDir, '*_EOF.txt'))[0]
    #try:
    height = _parse_ustec_height(eoffilename)
    #except:
    #	height = 450
        
    # Cleanup
    tf.close()
    shutil.rmtree(tempDir, ignore_errors=True)
    
    # Combine everything together
    dates = numpy.array(dates, dtype=numpy.float64)
    tec = numpy.array(tecMaps, dtype=numpy.float32)
    rms = numpy.array(rmsMaps, dtype=numpy.float32)
    
    # Build up the output
    output = {'dates': dates, 'lats': lats, 'lngs': lngs, 'height': height, 'tec': tec, 'rms': rms}
    
    # Done
    return output


def _load_map(mjd, timeout=120, type='IGS'):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    # Figure out which map to use
    if type == 'IGS':
        ## Cache entry name
        cacheName = 'TEC-IGS-%i' % mjd
        
        ## Download helper
        downloader = _download_igs
        
        ## Filename templates
        filenameTemplate = 'igsg%03i0.%02ii.gz'
        filenameAltTemplate = 'igrg%03i0.%02ii.gz'
        
    elif type == 'JPL':
        ## Cache entry name
        cacheName = 'TEC-JPL-%i' % mjd
        
        ## Download helper
        downloader = _download_jpl
        
        ## Filename templates
        filenameTemplate = 'jplg%03i0.%02ii.gz'
        filenameAltTemplate = 'jprg%03i0.%02ii.gz'
        
    elif type == 'UQR':
        ## Cache entry name
        cacheName = 'TEC-UQR-%i' % mjd
        
        ## Download helper
        downloader = _download_uqr
        
        ## Filename templates
        filenameTemplate = 'uqrg%03i0.%02ii.gz'
        filenameAltTemplate = 'uqrg%03i0.%02ii.gz'
        
    elif type == 'CODE':
        ## Cache entry name
        cacheName = 'TEC-CODE-%i' % mjd
        
        ## Download helper
        downloader = _download_code
        
        ## Filename templates
        filenameTemplate = 'CODG%03i0.%02iI.gz'
        filenameAltTemplate = 'CODG%03i0.%02iI.gz'
        
    elif type == 'USTEC':
        ## Cache entry name
        cacheName = 'TEC-USTEC-%i' % mjd
        
        ## Download helper
        downloader = _download_ustec
        
        ## Filename templates
        filenameTemplate = '%s_ustec.tar.gz'
        filenameAltTemplate = '%s_ustec.tar.gz'
        
    else:
        raise ValueError("Unknown data source '%s'" % type)
        
    try:
        # Is it already in the on-line cache?
        tecMap = _CACHE[cacheName]
    except KeyError:
        # Nope, we need to fetch it
        
        # Convert the MJD to a datetime instance so that we can pull out the year
        # and the day-of-year
        mpm = int((mjd - int(mjd))*24.0*3600.0*1000)
        dt = mjdmpm2datetime(int(mjd), mpm)
        
        if type == 'USTEC':
            # Pull out a YMD string
            dateStr = dt.strftime("%Y%m%d")
            
            # Figure out the filenames in order of preference.  We'd rather have
            # final values than rapid values
            filename = filenameTemplate % (dateStr)
            
            # Is the primary file in the disk cache?
            if not os.path.exists(os.path.join(_CACHE_DIR, filename)):
                ## Can we download it?
                status = downloader(mjd, timeout=timeout)
                
            else:
                ## Good we have the primary file
                pass
                
            # Parse it
            _CACHE[cacheName] = _parse_ustec_map(os.path.join(_CACHE_DIR, filename))
            
        else:
            
            # Pull out the year and the day-of-year
            year = dt.year
            dayOfYear = int(dt.strftime('%j'), 10)
            
            # Figure out the filenames in order of preference.  We'd rather have
            # final values than rapid values
            filename = filenameTemplate % (dayOfYear, year%100)
            filenameAlt = filenameAltTemplate % (dayOfYear, year%100)
            
            # Is the primary file in the disk cache?
            if not os.path.exists(os.path.join(_CACHE_DIR, filename)):
                ## Can we download it?
                status = downloader(mjd, timeout=timeout, type='final')
                if not status:
                    ## Nope, now check for the secondary file on disk
                    if not os.path.exists(os.path.join(_CACHE_DIR, filenameAlt)):
                        ## Can we download it?
                        status = downloader(mjd, timeout=timeout, type='rapid')
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
            _CACHE[cacheName] = _parse_tec_map(os.path.join(_CACHE_DIR, filename))
            
        tecMap = _CACHE[cacheName]
        
    # Done
    return tecMap


def get_tec_value(mjd, lat=None, lng=None, include_rms=False, timeout=120, type='IGS'):
    """
    Given an MJD value and, optionally, a latitude and longitude in degrees, 
    compute the TEC value in TECU above that location using data from the 
    IGS or CODE (depending on the value of the 'type' keyword).  If the 
    'include_rms' keyword is set to True, also return the RMS in the TEC 
    value.
    """
    
    # Load in the right map
    tecMap = _load_map(mjd, timeout=timeout, type=type)
    
    if type == 'USTEC':
        # Figure out the closest model point(s) to the requested MJD taking into
        # account that a new model is generated every fifteen minutes
        best = numpy.where( numpy.abs((tecMap['dates']-mjd)) < 15/60./24.0 )[0]
    else:
        # Figure out the closest model point(s) to the requested MJD taking into
        # account that a new model is generated every two hours
        best = numpy.where( numpy.abs((tecMap['dates']-mjd)) < 2/24.0 )[0]
        
    # Interpolate in time
    ## TEC
    slope = tecMap['tec'][best[-1],:,:] - tecMap['tec'][best[0],:,:]
    if best[-1] == best[0]:
        slope = 0.0
    else:
        slope /= tecMap['dates'][best[-1]] - tecMap['dates'][best[0]]
    tec = slope*(mjd - tecMap['dates'][best[0]]) + tecMap['tec'][best[0],:,:]
    
    ## RMS
    if include_rms:
        slope = tecMap['rms'][best[-1],:,:] - tecMap['rms'][best[0],:,:]
        if best[-1] == best[0]:
            slope = 0.0
        else:
            slope /= tecMap['dates'][best[-1]] - tecMap['dates'][best[0]]
        rms = slope*(mjd - tecMap['dates'][best[0]]) + tecMap['rms'][best[0],:,:]
        
    # Interpolate in location, if desired
    if lat is not None and lng is not None:
        ## TEC
        interpFunction = RectBivariateSpline(tecMap['lats'][::-1,0], tecMap['lngs'][0,:], tec[::-1,:], kx=1, ky=1)
        tec = interpFunction(lat, lng)
        
        ## RMS
        if include_rms:
            interpFunction = RectBivariateSpline(tecMap['lats'][::-1,0], tecMap['lngs'][0,:], rms[::-1,:], kx=1, ky=1)
            rms = interpFunction(lat, lng)
            
    # Done
    if include_rms:
        return tec, rms
    else:
        return tec


def get_ionospheric_pierce_point(site, az, el, height=450e3, verbose=False):
    """
    Given a site and a pointing direction (azimuth and elevation in degrees),
    compute the location of the ionospheric pierce  point.  Since the height
    assumed for the ionosphere is model-dependent the 'height' keyword sets 
    the elevation to use in meters.  Returns a three-element tuple of 
    latitude (degrees, North is positive), longitude (degrees, East is 
    positive), and elevation (meters).
    """
    
    # Create the various functions we need to figure out this optimization 
    # problem.  Maybe there is a better way to do this.
    
    ## Function that takes in a two-element tuple of latitude and longitude
    ## and returns the azimuth and elevation relative to the observer assuming
    ## a particular height.
    def func(params, xdata, site=site, elev=height):
        lat,lon = params
        
        az,el,d = site.getPointingAndDistance((lat, lon, elev))
        az %= (2*numpy.pi)
        
        az *= 180/numpy.pi
        el *= 180/numpy.pi
        
        return numpy.array([az, el])
        
    ## Error function that computes the difference between the input and the 
    ## model used in func().
    def err(params, ydata, xdata):
        return ydata - func(params, xdata)
        
    ## Secondary error function that computes a sum(err**2) that can be used
    ## with scipy.optimize.fmin()
    def err2(params, ydata, xdata):
        return (err(params, ydata, xdata)**2).sum()
        
    # Initial conditions for the optimization - we start directly overhead
    lat = site.lat * 180/numpy.pi
    lon = site.lon * 180/numpy.pi
    elev = site.elev + height
    x0 = (lat, lon)
    
    # Optimize
    output = fmin(err2, x0, args=(numpy.array([az, el]), []), disp=verbose)
    
    # Done
    return output[0], output[1], height
