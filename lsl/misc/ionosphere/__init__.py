"""
A collection of utilities for retrieving parameters that may be relevant 
for ionospheric corrections.
"""

import os
import sys
import ephem
import numpy as np
import warnings
from functools import lru_cache
from datetime import datetime, timedelta

from scipy.special import lpmv
try:
    from scipy.misc import factorial
except ImportError:
    from scipy.special import factorial
from scipy.optimize import fmin
from scipy.interpolate import RectBivariateSpline

from astropy.coordinates import Angle as AstroAngle

from lsl.common.stations import geo_to_ecef
from lsl.common.data_access import DataAccess
from lsl.common.mcs import mjdmpm_to_datetime, datetime_to_mjdmpm
from lsl.common.color import colorfy

from lsl.misc.ionosphere import _igs, _jpl, _emr, _uqr, _code, _ustec, _glotec

from lsl.misc import telemetry
telemetry.track_module()


__version__ = "0.8"
__all__ = ['get_magnetic_field', 'compute_magnetic_declination', 'compute_magnetic_inclination', 
           'get_tec_value', 'get_ionospheric_pierce_point']


# Radius of the Earth in meters for the IGRF
_RADIUS_EARTH = 6371.2*1e3


@lru_cache(maxsize=8)
def _load_igrf():
    """
    Load in the list of IGRF coefficient data and return a dictionary
    containing the raw coefficients.
    
    The dictionary keys are:
     * years - list of years for each of the models
     * g - dictionary of cosine term coefficients
     * h - dictionary of sine term coefficients
    
    The g and h dictionaries are keyed off the degree of the Legendre 
    function and each stores a list of n orders.  Each order is composed
    of len(years)+1 values, one for each each plus a secular evolution
    term.
    """
    
    # Go!
    dataCos = {}
    dataSin = {}
    with DataAccess.open('geo/igrf13coeffs.txt', 'r') as fh:
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
            c = np.array([float(v) for v in fields[3:]])
            
            ## Sort out cosine (g) vs. sine (h)
            if t == 'g':
                try:
                    dataCos[n][m] = c
                except KeyError:
                    dataCos[n] = [np.zeros(len(years)+1) for i in range(n+1)]
                    dataCos[n][m] = c
            else:
                try:
                    dataSin[n][m] = c
                except KeyError:
                    dataSin[n] = [np.zeros(len(years)+1) for i in range(n+1)]
                    dataSin[n][m] = c
                    
    # Build the output
    output = {'years': years, 'g': dataCos, 'h': dataSin}
    
    # Done
    return output


def _compute_igrf_coefficents(year, coeffs):
    """
    Given a decimal year and a coefficient dictionary from _load_igrf(),
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
    best = np.where( np.abs(year - np.array(coeffs['years'])) < 5 )[0]
    
    if year < min(coeffs['years']):
        # If the requested year is before 1900 we can't do anything
        raise RuntimeError(f"Invalid year {year}")
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
        return np.sqrt(factorial(n-m)/factorial(n+m))
    else:
        return np.sqrt(2.0*factorial(n-m)/factorial(n+m))


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
    o /= np.sqrt(1.0 - mu**2)
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
        mjd, mpm = datetime_to_mjdmpm( datetime.utcnow() )
        mjd = mjd + mpm/1000.0/3600.0/24.0
        
    # Convert the MJD to a decimal year.  This is a bit tricky
    ## Break the MJD into an integer MJD and an MPM in order to build a datetime instance
    mpm = int((mjd - int(mjd))*24.0*3600.0*1000.0)
    mjd0 = mjdmpm_to_datetime(int(mjd), mpm)
    ## Convert the datetime instance to January 1
    mjd0 = mjd0.replace(month=1, day=1, hour=0, second=0, microsecond=0)
    ## Figure out January 1 for the following year
    mjd1 = mjd0.replace(year=mjd0.year+1)
    ## Figure out how long the year is in days
    diffDays = mjd1-mjd0
    diffDays = diffDays.days + diffDays.seconds/86400.0 + diffDays.microseconds/1e6/86400.0
    ## Convert the January 1 date back to an MJD
    mjd0, mpm0 = datetime_to_mjdmpm(mjd0)
    mjd0 = mjd0 + mpm/1000.0/3600.0/24.0
    year = (mjd1.year - 1) + (mjd - mjd0) / diffDays
    
    # Convert the geodetic position provided to a geocentric one for calculation
    ## Deal with the poles
    if 90.0 - lat < 0.001:
        xyz = np.array(geo_to_ecef(89.999*np.pi/180, lng*np.pi/180, elev))
    elif 90.0 + lat < 0.001:
        xyz = np.array(geo_to_ecef(-89.999*np.pi/180, lng*np.pi/180, elev))
    else:
        xyz = np.array(geo_to_ecef(lat*np.pi/180, lng*np.pi/180, elev))
    ## To geocentric
    r = np.sqrt( (xyz**2).sum() )
    lt = np.arcsin(xyz[2]/r)
    ln = np.arctan2(xyz[1], xyz[0])
    
    # Load in the coefficients
    coeffs = _load_igrf()
    
    # Compute the coefficients for the epoch
    coeffs = _compute_igrf_coefficents(year, coeffs)
    
    # Compute the field strength in spherical coordinates
    Br, Bth, Bph = 0.0, 0.0, 0.0
    for n in coeffs['g'].keys():
        for m in range(0, n+1):
            Br  += (n+1.0)*(_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['g'][n][m]*np.cos(m*ln) * _Pnm(n, m, np.sin(lt))
            Br  += (n+1.0)*(_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['h'][n][m]*np.sin(m*ln) * _Pnm(n, m, np.sin(lt))
            
            Bth -= (_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['g'][n][m]*np.cos(m*ln) * _dPnm(n, m, np.sin(lt))
            Bth -= (_RADIUS_EARTH/r)**(n+2) * _Snm(n,m)*coeffs['h'][n][m]*np.sin(m*ln) * _dPnm(n, m, np.sin(lt))
            
            Bph += (_RADIUS_EARTH/r)**(n+2)/np.cos(lt) * _Snm(n,m)*coeffs['g'][n][m]*m*np.sin(m*ln) * _Pnm(n, m, np.sin(lt))
            Bph -= (_RADIUS_EARTH/r)**(n+2)/np.cos(lt) * _Snm(n,m)*coeffs['h'][n][m]*m*np.cos(m*ln) * _Pnm(n, m, np.sin(lt))
    ## And deal with NaNs
    if np.isnan(Br):
        Br = 0.0
    if np.isnan(Bth):
        Bth = 0.0
    if np.isnan(Bph):
        Bph = 0.0
        
    # Convert from spherical to ECEF
    Bx = Br*np.cos(lt)*np.cos(ln) + Bth*np.sin(lt)*np.cos(ln) - Bph*np.sin(ln)
    By = Br*np.cos(lt)*np.sin(ln) + Bth*np.sin(lt)*np.sin(ln) + Bph*np.cos(ln)
    Bz = Br*np.sin(lt) - Bth*np.cos(lt)
    
    # Are we done?
    if ecef:
        # For ECEF we don't need to do anything else
        outputField = Bx, By, Bz
        
    else:
        # Convert from ECEF to topocentric (geodetic)
        ## Update the coordinates for geodetic
        lt = lat*np.pi/180.0
        if 90.0 - lat < 0.001:
            lt = 89.999*np.pi/180.0
        elif 90.0 + lat < 0.001:
            lt = -89.999*np.pi/180.0
        else:
            lt = lat*np.pi/180.0
        ln = lng*np.pi/180.0
        
        ## Build the rotation matrix for ECEF to SEZ
        rot = np.array([[ np.sin(lt)*np.cos(ln), np.sin(lt)*np.sin(ln), -np.cos(lt)], 
                        [-np.sin(ln),            np.cos(ln),             0         ],
                        [ np.cos(lt)*np.cos(ln), np.cos(lt)*np.sin(ln),  np.sin(lt)]])
                    
        ## Apply and extract
        sez = np.dot(rot, np.array([Bx,By,Bz]))
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
    Bh = np.sqrt(Bn**2+Be**2)
    
    # Compute the declination
    decl = 2.0*np.arctan2(Be, Bh+Bn)
    
    # Check for bounds
    if Bh < 100.0:
        decl = np.nan
        
    # Convert to degrees and done
    return decl*180.0/np.pi


def compute_magnetic_inclination(Bn, Be, Bz):
    """
    Given the topocentric output of get_magnetic_field(), compute and return 
    the magnetic inclination or magnetic dip (angle between the magnetic 
    field and the horizontal) in degrees.
    """
    
    # Compute the horizontal field strength
    Bh = np.sqrt(Bn**2+Be**2)
    
    # Compute the inclination.  This has an extra negative sign because of the
    # convention used in get_magnetic_field().
    incl = np.arctan2(-Bz, Bh)
    
    # Convert to degrees and done
    return incl*180.0/np.pi


@lru_cache(maxsize=64)
def _load_map(mjd, type='IGS'):
    """
    Given an MJD value, load the corresponding TEC map.  If the map is not
    avaliable on disk, download it.
    """
    
    # Figure out which map to use
    if type.upper() == 'IGS':
        ## Cache entry name
        cacheName = 'TEC-IGS-%i' % mjd
        
        ## Loading helper
        loader = _igs.load_mjd
        
    elif type.upper() == 'JPL':
        ## Cache entry name
        cacheName = 'TEC-JPL-%i' % mjd
        
        ## Loading helper
        loader = _jpl.load_mjd
        
    elif type.upper() == 'EMR':
        ## Cache entry name
        cacheName = 'TEC-EMR-%i' % mjd
        
        ## Loading helper
        loader = _emr.load_mjd
        
    elif type.upper() == 'UQR':
        ## Cache entry name
        cacheName = 'TEC-UQR-%i' % mjd
        
        ## Loading helper
        loader = _uqr.load_mjd
        
    elif type.upper() == 'CODE':
        ## Cache entry name
        cacheName = 'TEC-CODE-%i' % mjd
        
        ## Loading helper
        loader = _code.load_mjd
        
    elif type.upper() == 'USTEC':
        ## Cache entry name
        cacheName = 'TEC-USTEC-%i' % mjd
        
        ## Loading helper
        loader = _ustec.load_mjd
        
    elif type.upper() == 'GLOTEC':
        # Cache entry name
        hour = (mjd - int(mjd)) * 24
        cacheName = 'TEC-GLOTEC-%i-%i' % (mjd, hour)
        
        ## Loading helper
        loader = _glotec.load_mjd
        
    else:
        raise ValueError(f"Unknown data source '{type}'")
        
    tecMap = loader(mjd)    
    
    # Done
    return tecMap


def get_tec_value(mjd, lat=None, lng=None, include_rms=False, type='IGS'):
    """
    Given an MJD value and, optionally, a latitude and longitude in degrees, 
    compute the TEC value in TECU above that location using data from the 
    IGS or CODE (depending on the value of the 'type' keyword).  If the 
    'include_rms' keyword is set to True, also return the RMS in the TEC 
    value.
    """
    
    # Load in the right map
    tecMap = _load_map(mjd, type=type)
    
    if type.upper() == 'GLOTEC':
        # Figure out the closest model point(s) to the requested MJD taking into
        # account that a new model is generated every ten minutes
        best = np.where( np.abs((tecMap['dates']-mjd)) < 10/60./24.0 )[0]
    elif type.upper() == 'USTEC':
        # Figure out the closest model point(s) to the requested MJD taking into
        # account that a new model is generated every fifteen minutes
        best = np.where( np.abs((tecMap['dates']-mjd)) < 15/60./24.0 )[0]
    else:
        # Figure out the closest model point(s) to the requested MJD taking into
        # account that a new model is generated every two hours
        best = np.where( np.abs((tecMap['dates']-mjd)) < 2/24.0 )[0]
        
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
        
        az,el,d = site.get_pointing_and_distance((lat, lon, elev))
        az %= (2*np.pi)
        
        az *= 180/np.pi
        el *= 180/np.pi
        
        return np.array([az, el])
        
    ## Error function that computes the difference between the input and the 
    ## model used in func().
    def err(params, ydata, xdata):
        return ydata - func(params, xdata)
        
    ## Secondary error function that computes a sum(err**2) that can be used
    ## with scipy.optimize.fmin()
    def err2(params, ydata, xdata):
        return (err(params, ydata, xdata)**2).sum()
        
    # Initial conditions for the optimization - we start directly overhead
    lat = site.lat * 180/np.pi
    lon = site.lon * 180/np.pi
    elev = site.elev + height
    x0 = (lat, lon)
    
    # Convert
    if isinstance(az, AstroAngle):
        az = az.deg
    elif isinstance(az, ephem.Angle):
        az = az * 180/np.pi
        
    if isinstance(el, AstroAngle):
        el = el.deg
    elif isinstance(el, ephem.Angle):
        el = el * 180/np.pi
        
    # Optimize
    output = fmin(err2, x0, args=(np.array([az, el]), []), disp=verbose)
    
    # Done
    return output[0], output[1], height
