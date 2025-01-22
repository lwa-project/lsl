"""
Module to support imaging correlated data.  This module provides utilities to 
read FITS IDI files into :class:`lsl.imaging.data.VisibilityDataSet` or 
:class:`lsl.imaging.data.VisibilityData` object (as described in 
:mod:`lsl.imaging.data`) and build AIPY ImgW instances from the data.

.. versionadded:: 0.5.0

.. versionchanged:: 1.0.0
    Added support for UVFITS files and CASA measurement sets

.. versionchanged:: 1.0.1
    Added the plot_gridded_image() function
    
.. versionchanged:: 1.1.0
    Added the get_image_radec() and get_image_azalt() functions to complement
    plot_gridded_image() and make it easier to work with phase centers 
    that are not at zenith.  Added in the ImgWPlus class to add support
    for imaging weighting and tapering.
    
.. versionchanged:: 1.1.2
    Added support for CASA MeasurementSets that are stored as a tarball
    
.. versionchanged:: 1.1.4
    Fixed a conjugation problem in the visibilities read from a FITS-IDI file
"""

t = {}
for c in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ':
    t[ord(c)] = None
def strip_letters(s, t=t):
    return s.translate(t)
    
import os
import re
import sys
import aipy
import ephem
import numpy as np
import atexit
import shutil
import tarfile
import tempfile
import warnings
from calendar import timegm
from datetime import datetime

from astropy import units as astrounits
from astropy.constants import c as vLight
from astropy.io import fits as astrofits
from astropy.time import Time as AstroTime
from astropy.coordinates import EarthLocation, SkyCoord, ICRS, FK5, AltAz, CartesianRepresentation
from astropy.wcs import WCS

from lsl import astro
from lsl.statistics import robust
from lsl.common import stations
from lsl.sim import vis as simVis
from lsl.writer.fitsidi import NUMERIC_STOKES
from lsl.writer.measurementset import NUMERIC_STOKES as NUMERIC_STOKESMS
from lsl.common.color import colorfy
from lsl.testing import SilentVerbose

from lsl.imaging._gridder import WProjection
from lsl.imaging.data import PolarizationDataSet, VisibilityDataSet, VisibilityData

from scipy import fftpack
fft2Function = fftpack.fft2
ifft2Function = fftpack.ifft2

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '1.0'
__all__ = ['CorrelatedData', 'CorrelatedDataIDI', 'CorrelatedDataUV', 'CorrelatedDataMS',
           'convert_to_stokes', 'convert_to_linear', 'ImgWPlus', 'build_gridded_image',
           'plot_gridded_image', 'get_image_radec', 'get_image_azalt']



vLight = vLight.to('m/s').value


# Regular expression for trying to get the stand number out of an antenna
# name
_annameRE = re.compile(r'^.*?(?P<id>\d{1,3})$')


def CorrelatedData(filename, verbose=False):
    """
    Read in and work with FITS IDI and UVFITS files.  Returns either a 
    CorrelateDataIDI or CorrelatedDataUV instance.
    """
    
    valid = False
    
    # Basic filesystem validation
    if not os.path.exists(filename):
        raise IOError(f"File '{filename}' does not exists")
    if not os.access(filename, os.R_OK):
        raise IOError(f"File '{filename}' cannot be read")
        
    if os.path.isdir(filename):
        # Only a MS can be a directory
        try:
            return CorrelatedDataMS(filename)
        except Exception as e:
            if verbose:
                print("MS - ERROR: %s" % str(e))
            raise RuntimeError(f"Directory '{filename}' does not appear to be a MeasurmentSet")
            
    else:
        # Standard files
        ## FITS IDI
        try:
            return CorrelatedDataIDI(filename)
        except Exception as e:
            if verbose:
                print("FITSIDI - ERROR: %s" % str(e))
                
        ## UVFITS
        try:
            return CorrelatedDataUV(filename)
        except Exception as e:
            if verbose:
                print("UVFITS - ERROR: %s" % str(e))
                
        ## Measurment Set as a compressed entity
        try:
            return CorrelatedDataMS(filename)
        except Exception as e:
            if verbose:
                print("MS - ERROR: %s" % str(e))
                
    if not valid:
        raise RuntimeError(f"File '{filename}' does not appear to be either a FITS IDI file, UV FITS file, or MeasurmentSet")


class CorrelatedDataBase(object):
    """
    Base class for acessing visibility data stored in a file.
    """
    
    def __init__(self, filename):
        self.filename = filename
        
    def __str__(self):
        return f"{self.__name__} @ {self.filename}"
        
    def __enter__(self):
        return self
        
    def __exit__(self, type, value, tb):
        self.close()
        
    def get_antennaarray(self):
        """
        Return an AIPY AntennaArray instance for the array that made the 
        observations contained here.
        """
        
        # Get the date of observations
        refJD = astro.unix_to_utcjd(timegm(self.date_obs.timetuple()))
        
        # Return
        return simVis.build_sim_array(self.station, self.antennas, self.freq/1e9, jd=refJD)
        
    def get_observer(self):
        """
        Return a ephem.Observer instances for the array described in the file.
        """
        
        return self.station
        
    def get_data_set(self, sets, include_auto=False, sort=True, min_uv=0, max_uv=np.inf):
        """
        Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
        :class:`lsl.imaging.data.VisibilityData` object for all 
        baselines for a given set of observations for the specified data set.  
        By default this excludes the autocorrelations.  To include 
        autocorrelations set the value of 'include_auto' to True.  Setting the
        'sort' keyword to False will disable the baseline sorting.  Optionally,
        baselines with lengths between min_uv and max_uv can only be returned.

        .. note::
            min_uv and max_uv should be specified in lambda
        
        .. versionchanged:: 1.1.0
            'set' can now be either an integer or a list to pull back multiple 
            integrations.
        """
        raise NotImplementedError
        
    def data_set_sequence(self, include_auto=False, sort=True, min_uv=0, max_uv=np.inf):
        """
        Return a generator that yields :class:`lsl.imaging.data.VisibilityDataSet` 
        objects for each integration contained in the file one at a time.
        
        .. note::
            min_uv and max_uv should be specified in lambda
        """
        
        for i in range(self.integration_count):
            yield self.get_data_set(i+1, include_auto=include_auto, sort=sort, 
                                    min_uv=min_uv, max_uv=max_uv)
            
    def close(self):
        """
        Close out the object.
        """
        
        pass


class CorrelatedDataIDI(CorrelatedDataBase):
    """
    Class to make accessing information about a FITS IDI easy.  This wraps 
    all of the "messy" machinery needed to extract both the metadata and data 
    from the file and return them as common LSL objects.
    
    This class has three main attributes to interact with:
     * get_antennaarray - Return a :class:`lsl.sim.vim.AntennaArray` instance
                          that represents the array where the data was obtained.
                          This is useful for simulation proposes and computing 
                          source positions.
     * get_observer - Return a ephem.Observer instance representing the array
     * get_data_set - Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
                      :class:`lsl.imaging.data.VisibilityData` object for all 
                      baselines for a given set of observations
        
    The class also includes a variety of useful metadata attributes:
     * pols - Numpy array of polarization product codes
     * freq - Numpy array of frequency channels in Hz
     * station - LSL :class:`lsl.common.stations.LWAStation` instance for the
                 array
     * date_obs - Datetime object for the reference date of the FIT IDI file
     * antennas - List of :class:`lsl.common.stations.Antenna` instances
    
    .. note::
        The CorrelatedData.antennas attribute should be used over 
        CorrelatedData.station.antennas since the mapping in the FITS IDI
        file may not be the same as the digitizer order.
    """
    
    def __init__(self, filename):
        """
        Initialize a new CorrelatedDataIDI instance from a FITS IDI file and 
        fill in the metadata.
        """
        
        CorrelatedDataBase.__init__(self, filename)
        
        # Open the file, check if it looks like FITS IDI, and pull out the UV_DATA table
        with astrofits.open(self.filename) as hdulist:
            tbls = [i.header['EXTNAME'] for i in hdulist[1:]]
            for tbl in ('ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'SOURCE', 'UV_DATA'):
                if tbl not in tbls:
                    raise RuntimeError(f"Cannot find table '{tbl}' in '{self.filename}'")
            
            self.extended = False
            self.conjugate = True
            try:
                if hdulist[0].header['LWATYPE'] == 'IDI-EXTENDED-ZA':
                    self.extended = True
            except KeyError:
                ## Catch for LEDA64-NM data
                pass
            try:
                if hdulist[0].header['LWAMAJV'] < 3:
                    self.conjugate = False
            except KeyError:
                pass
            ag = hdulist['ARRAY_GEOMETRY']
            fq = hdulist['FREQUENCY']
            uvData = hdulist['UV_DATA']
            
            # Antennas
            try:
                mapper = hdulist['NOSTA_MAPPER']
                
                nosta = mapper.data.field('NOSTA')
                noact = mapper.data.field('NOACT')
                anname = mapper.data.field('ANNAME')
            except KeyError:
                nosta = ag.data.field('NOSTA')
                noact = ag.data.field('NOSTA')
                anname = ag.data.field('ANNAME')
                
            # Station/telescope information
            try:
                self.telescope = hdulist[0].header['TELESCOP']
                self.date_obs = datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
            except ValueError:
                ## Catch for DiFX FITS-IDI files
                self.date_obs = datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%d")
            except KeyError:
                ## Catch for LEDA64-NM data
                self.telescope = uvData.header['TELESCOP']
                self.date_obs = datetime.strptime(uvData.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
                
            ## Extract the site position
            geo = np.array([ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ']])
            site = stations.ecef_to_geo(*geo)
            
            ## Try to back out the "real" stand names
            noact2 = []
            for nam in anname:
                try:
                    mtch =  _annameRE.match(nam)
                    id = int(mtch.group('id'))
                    noact2.append(id)
                except (ValueError, AttributeError):
                    break
            if len(noact2) == len(noact):
                noact = np.array(noact2)
                
            ## Build up the station
            self.station = stations.LWAStation(ag.header['ARRNAM'], site[0]*180/np.pi, site[1]*180/np.pi, site[2])
            self.station.date = astro.unix_to_utcjd(timegm(self.date_obs.timetuple())) \
                                - astro.DJD_OFFSET
            
            ## Build up the list of antennas
            antennas = []
            for line,act in zip(ag.data, noact):
                el = EarthLocation.from_geocentric(*(line['STABXYZ']*astrounits.m))
                enz = self.station.get_enz_offset(el)
                
                stand = stations.Stand(act, *enz)
                antennas.append(stations.Antenna(2*(stand.id-1)+1, stand=stand, pol=0))
                
            ## Add in the antennas
            self.station.antennas = antennas
            
            self.stand_map = {}
            self.stands = []
            for sta, act in zip(nosta, noact):
                self.stand_map[sta] = act
                self.stands.append(act)
                
            self.antenna_map = {}
            self.antennas = []
            for stand in self.stands:
                for ant in self.station.antennas:
                    if ant.stand.id == stand and ant.pol == 0:
                        self.antennas.append(ant)
                        self.antenna_map[ant.stand.id] = ant
                        break
                        
            # Polarization and frequency
            self.pols  = np.arange(1, uvData.header['MAXIS2']+1) - uvData.header['CRPIX2']
            self.pols *= uvData.header['CDELT2'] 
            self.pols += uvData.header['CRVAL2']
            self.freq  = np.array([], dtype=np.float64)
            for i in range(uvData.header['NO_BAND']):
                width = fq.data['CH_WIDTH'][0][i] if uvData.header['NO_BAND'] > 1 else fq.data['CH_WIDTH'][0]
                offset = fq.data['BANDFREQ'][0][i] if uvData.header['NO_BAND'] > 1 else fq.data['BANDFREQ'][0]
                
                freqIF = np.arange(1, uvData.header['NO_CHAN']+1, dtype=np.float64) - uvData.header['REF_PIXL']
                freqIF *= width
                freqIF += uvData.header['REF_FREQ'] + offset
                self.freq = np.concatenate([self.freq, freqIF])
                
            # Total baseline count
            self.total_baseline_count = len(uvData.data['BASELINE'])
            
            # Data set times and integration count
            jd = uvData.data['DATE'] + uvData.data['TIME']
            self._times = np.unique(jd)
            self.integration_count = len(self._times)
            
    def get_data_set(self, sets, include_auto=False, sort=True, min_uv=0, max_uv=np.inf):
        """
        Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
        :class:`lsl.imaging.data.VisibilityData` object for all 
        baselines for a given set of observations for the specified data set.  
        By default this excludes the autocorrelations.  To include 
        autocorrelations set the value of 'include_auto' to True.  Setting the
        'sort' keyword to False will disable the baseline sorting.  Optionally,
        baselines with lengths between min_uv and max_uv can only be returned.

        .. note::
            min_uv and max_uv should be specified in lambda
        
        .. versionchanged:: 1.1.0
            'set' can now be either an integer or a list to pull back multiple 
            integrations.
        """
        
        # Open the file
        with astrofits.open(self.filename) as hdulist:
            uvData = hdulist['UV_DATA']
            
            dataSets = VisibilityData()
            try:
                len(sets)
            except TypeError:
                sets = range(sets, sets+1)
            for set in sets:
                # Set the time to look for
                targetTime = self._times[set-1]
                targetJD = targetTime
                
                # Figure out what rows we need
                selection = np.where( (uvData.data['DATE'] + uvData.data['TIME']) == targetTime )[0]
                
                # Figure out the source we are working on and create a phase center
                # if there is only a single source
                phase_center = None
                src_id = uvData.data['SOURCE'][selection]
                src_id = np.unique(src_id)
                if len(src_id) == 1:
                    src_id = src_id[0]
                    srcData = hdulist['SOURCE']
                    for row in srcData.data:
                        if row['SOURCE_ID'] == src_id:
                            phase_time = AstroTime(row['EPOCH'], format='jyear', scale='utc')
                            phase_center = aipy.amp.RadioFixedBody(row['RAEPO'] * np.pi/180, 
                                                                   row['DECEPO'] * np.pi/180, 
                                                                   name=row['SOURCE'], 
                                                                   epoch=phase_time.jd - astro.DJD_OFFSET)
                            
                # Figure out if we have seperate WEIGHT data or not
                seperateWeights = False
                for col in uvData.data.columns:
                    if col.name == 'WEIGHT':
                        seperateWeights = True
                        break
                        
                # Pull out the raw data from the table
                bl = uvData.data['BASELINE'][selection]
                jd = uvData.data['DATE'][selection] + uvData.data['TIME'][selection]
                try:
                    u, v, w = uvData.data['UU'][selection], uvData.data['VV'][selection], uvData.data['WW'][selection]
                except KeyError:
                    u, v, w = uvData.data['UU---SIN'][selection], uvData.data['VV---SIN'][selection], uvData.data['WW---SIN'][selection]
                vis = np.ascontiguousarray(uvData.data['FLUX'][selection], dtype=np.float32)
                if seperateWeights:
                    wgt = np.ascontiguousarray(uvData.data['WEIGHT'][selection])
                else:
                    wgt = None
                    
                # Re-work the data into something more useful
                ## Axis sizes
                nFreq = len(self.freq)
                nStk = len(self.pols)
                nCmp = 2 if seperateWeights else 3
                if vis.size//nFreq//nStk//nCmp != len(bl):
                    ### Catch for FITS-IDI files generate by interfits
                    nCmp = 2 if nCmp == 3 else 3
                ## Frequency for converting the u, v, and w coordinates
                freq = self.freq*1.0
                freq.shape += (1,)
                ## Convert u, v, and w from seconds to wavelengths and then into one massive array
                u = (u*freq).T
                v = (v*freq).T
                w = (w*freq).T
                uvw = np.array([u,v,w], dtype=np.float32)
                ## Reshape the visibilities and weights
                vis.shape = (vis.size//nFreq//nStk//nCmp, nFreq, nStk, nCmp)
                if seperateWeights:
                    if wgt.shape != nFreq*nStk:
                        ## Catch for some old stuff
                        wgt = np.concatenate([wgt,]*len(self.pols))
                    wgt.shape = (wgt.size//nFreq//nStk, nFreq, nStk)
                else:
                    try:
                        wgt = vis[:,:,:,2]
                        vis = vis[:,:,:,:2]
                    except IndexError:
                        ### Catch for FITS-IDI files generate by interfits
                        wgt = np.ones([vis.shape[i] for i in range(3)], dtype=np.float32)
                ## Back to complex
                vis = vis.view(np.complex64)
                vis = vis[...,0]
                if self.conjugate:
                    ## NOTE: This is this conjugate since there seems to be a convention mis-match
                    ##       between LSL and AIPS/the FITS-IDI convention.
                    vis = vis.conj()
                ## Scale
                try:
                    scl = uvData.header['VIS_SCAL']
                    vis /= scl
                except KeyError:
                    pass
                    
                # Setup the output data
                baselines = []
                select = []
                for b in range(bl.size):
                    if not self.extended:
                        i = self.stand_map[(bl[b] >> 8) & 255]
                        j = self.stand_map[bl[b] & 255]
                    else:
                        i = self.stand_map[(bl[b] >> 16) & 65535]
                        j = self.stand_map[bl[b] & 65535]
                    if i == j and not include_auto:
                        ## Skip auto-correlations
                        continue
                    ri = np.where(self.stands == i)[0][0]
                    rj = np.where(self.stands == j)[0][0]
                    baselines.append( (ri,rj) )
                    select.append( b )
                    
                # Update the phase center information
                aa = self.get_antennaarray()
                aa.date = jd[0] - astro.DJD_OFFSET
                phase_center.compute(aa)
                
                # Build the output data set
                dataSet = VisibilityDataSet(jd[0], self.freq*1.0, baselines=baselines, 
                                            uvw=uvw[:,select,:].transpose(1,0,2), 
                                            antennaarray=aa, 
                                            phase_center=phase_center)
                for p,l in enumerate(self.pols):
                    name = NUMERIC_STOKES[l]
                    polDataSet = PolarizationDataSet(name, data=vis[select,:,p], weight=wgt[select,:,p])
                    dataSet.append(polDataSet)
                dataSets.append( dataSet )
                
        # Sort
        if sort:
            dataSets.sort()
            
        # Prune
        if min_uv != 0 or max_uv != np.inf:
            dataSets = dataSets.get_uv_range(min_uv=min_uv, max_uv=max_uv)
            
        # Prune a different way
        if len(dataSets) == 1:
            dataSets = dataSets.pop()
            
        # Return
        return dataSets


class CorrelatedDataUV(CorrelatedDataBase):
    """
    Class to make accessing information about a UVFITS file easy.  This wraps 
    all of the "messy" machinery needed to extract both the metadata and data 
    from the file and return them as common LSL objects.
    
    This class has three main attributes to interact with:
     * get_antennaarray - Return a :class:`lsl.sim.vim.AntennaArray` instance
                          that represents the array where the data was obtained.
                          This is useful for simulation proposes and computing 
                          source positions.
     * get_observer - Return a ephem.Observer instance representing the array
     * get_data_set - Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
                      :class:`lsl.imaging.data.VisibilityData` object of all 
                      baselines for a given set of observations
        
    The class also includes a variety of useful metadata attributes:
     * pols - Numpy array of polarization product codes
     * freq - Numpy array of frequency channels in Hz
     * station - LSL :class:`lsl.common.stations.LWAStation` instance for the
                 array
     * date_obs - Datetime object for the reference date of the FIT IDI file
     * antennas - List of :class:`lsl.common.stations.Antenna` instances
    
    .. note::
        The CorrelatedDataUV.antennas attribute should be used over 
        CorrelatedDataUV.station.antennas since the mapping in the UVFITS
        file may not be the same as the digitizer order.
    """
    
    def __init__(self, filename):
        """
        Initialize a new CorrelatedDataUV instance from a UVFITS file and 
        fill in the metadata.
        """
        
        CorrelatedDataBase.__init__(self, filename)
        
        # Open the various tables that we need
        with astrofits.open(filename) as hdulist:
            uvData = hdulist[0]
            ag = hdulist['AIPS AN']
            
            # Antennas
            nosta = ag.data.field('NOSTA')
            noact = ag.data.field('NOSTA')
            anname = ag.data.field('ANNAME')
            
            # Station/telescope information
            self.telescope = hdulist[0].header['TELESCOP']
            dt = hdulist[0].header['DATE-OBS']
            dt = dt.rsplit('.', 1)[0]
            try:
                self.date_obs = datetime.strptime(dt, "%Y-%m-%dT%H:%M:%S")
            except ValueError:
                ## Catch for AIPS UVFITS files which only have a date set
                self.date_obs = datetime.strptime(dt, "%Y-%m-%d")
                
            ## Extract the site position
            geo = np.array([ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ']])
            site = stations.ecef_to_geo(*geo)
            
            ## Try to back out the "real" stand names
            noact2 = []
            for nam in anname:
                try:
                    mtch =  _annameRE.match(nam)
                    id = int(mtch.group('id'))
                    noact2.append(id)
                except (ValueError, AttributeError):
                    break
            if len(noact2) == len(noact):
                noact = np.array(noact2)
                
            ## Build up the station
            self.station = stations.LWAStation(ag.header['ARRNAM'], site[0]*180/np.pi, site[1]*180/np.pi, site[2])
            self.station.date = astro.unix_to_utcjd(timegm(self.date_obs.timetuple())) \
                                - astro.DJD_OFFSET
            
            ## Build up the list of antennas
            antennas = []
            for line,act in zip(ag.data, noact):
                el = EarthLocation.from_geocentric(*(line['STABXYZ']*astrounits.m))
                enz = self.station.get_enz_offset(el)
                
                stand = stations.Stand(act, *enz)
                antennas.append(stations.Antenna(2*(stand.id-1)+1, stand=stand, pol=0))
                
            ## Add in the antennas
            self.station.antennas = antennas
            
            self.stand_map = {}
            self.stands = []
            for nosta, noact in zip(nosta, noact):
                self.stand_map[nosta] = noact
                self.stands.append(noact)
                
            self.antenna_map = {}
            self.antennas = []
            for stand in self.stands:
                for ant in self.station.antennas:
                    if ant.stand.id == stand and ant.pol == 0:
                        self.antennas.append(ant)
                        self.antenna_map[ant.stand.id] = ant
                        break
                        
            # Polarization and frequency
            self.pols  = np.arange(1, uvData.header['NAXIS3']+1) - uvData.header['CRPIX3']
            self.pols *= uvData.header['CDELT3'] 
            self.pols += uvData.header['CRVAL3']
            nChan = uvData.header['NAXIS4']
            if uvData.header['CTYPE5'] == 'IF':
                ## Merge the frequency and IF columns
                nChan *= uvData.header['NAXIS5']
            self.freq  = np.arange(1, nChan+1, dtype=np.float64) - uvData.header['CRPIX4']
            self.freq *= uvData.header['CDELT4']
            self.freq += uvData.header['CRVAL4']
            
            # Total baseline count
            self.total_baseline_count = len(hdulist[0].data['BASELINE'])
            
            # Data set times and integration count
            try:
                jd = hdulist[0].data['DATE'] + hdulist[0].data['_DATE']
            except KeyError:
                jd = hdulist[0].data['DATE']
            self._times = np.unique(jd)
            self.integration_count = len(self._times)
            
    def get_data_set(self, sets, include_auto=False, sort=True, min_uv=0, max_uv=np.inf):
        """
        Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
        :class:`lsl.imaging.data.VisibilityData` object for the specified data set.  
        By default this excludes the autocorrelations.  To include 
        autocorrelations set the value of 'include_auto' to True.  Setting the
        'sort' keyword to False will disable the baseline sorting.  Optionally,
        baselines with lengths between min_uv and max_uv can only be returned.

        .. note::
            min_uv and max_uv should be specified in lambda
        
        .. versionchanged:: 1.1.0
            'set' can now be either an integer or a list to pull back multiple 
            integrations.
        """
        
        # Open the file
        with astrofits.open(self.filename) as hdulist:
            uvData = hdulist[0]
            
            dataSets = VisibilityData()
            try:
                len(sets)
            except TypeError:
                sets = range(sets, sets+1)
            for set in sets:
                # Set the time to look for
                targetTime = self._times[set-1]
                targetJD = targetTime
                
                # Figure out what rows we need
                try:
                    selection = np.where( uvData.data['DATE']+uvData.data['_DATE'] == targetTime )[0]
                except KeyError:
                    selection = np.where( uvData.data['DATE'] == targetTime )[0]
                    
                # Figure out the source we are working on and create a phase center
                # if there is only a single source
                phase_center = None
                src_id = uvData.data['SOURCE'][selection]
                src_id = np.unique(src_id)
                if len(src_id) == 1:
                    src_id = src_id[0]
                    srcData = hdulist['AIPS SU']
                    for row in srcData.data:
                        if row['ID. NO.'] == src_id:
                            phase_time = AstroTime(row['EPOCH'], format='jyear', scale='utc')
                            phase_center = aipy.amp.RadioFixedBody(row['RAEPO'] * np.pi/180, 
                                                                   row['DECEPO'] * np.pi/180, 
                                                                   name=row['SOURCE'], 
                                                                   epoch=phase_time.jd - astro.DJD_OFFSET)
                            
                            
                # Pull out the raw data from the table
                bl = uvData.data['BASELINE'][selection]
                try:
                    jd = uvData.data['DATE'][selection] + uvData.data['_DATE'][selection]
                except KeyError:
                    jd = uvData.data['DATE'][selection]
                try:
                    u, v, w = uvData.data['UU'][selection], uvData.data['VV'][selection], uvData.data['WW'][selection]
                except KeyError:
                    u, v, w = uvData.data['UU---SIN'][selection], uvData.data['VV---SIN'][selection], uvData.data['WW---SIN'][selection]
                vis = uvData.data['DATA'][selection]
                wgt = None
                
                # Re-work the data into something more useful
                ## Axis sizes
                nFreq = len(self.freq)
                nStk = len(self.pols)
                nCmp = vis.shape[-1]
                ## Frequency for converting the u, v, and w coordinates
                freq = self.freq*1.0
                freq.shape += (1,)
                ## Convert u, v, and w from seconds to wavelengths and then into one massive array
                u = (u*freq).T
                v = (v*freq).T
                w = (w*freq).T
                uvw = np.array([u,v,w], dtype=np.float32)
                ## Reshape the visibilities and weights
                if len(vis.shape) == 7:
                    ### Merge the frequency and IF columns
                    vis = vis[:,0,0,:,:,:,:]
                    vis.shape = (vis.shape[0], vis.shape[1]*vis.shape[2], vis.shape[3], vis.shape[4])
                else:
                    vis = vis[:,0,0,:,:,:]
                if vis.shape[-1] == 3:
                    wgt = vis[:,:,:,2]
                    vis = vis[:,:,:,:2]
                else:
                    wgt = np.ones((vis.shape[0], vis.shape[1], vis.shape[2]), dtype=np.float32)
                ## Back to complex
                vis = vis[:,:,:,0] + 1j*vis[:,:,:,1]
                
                # Setup the output data
                baselines = []
                select = []
                for b in range(bl.size):
                    if bl[b] >= 65536:
                        i = self.stand_map[int((bl[b] - 65536) / 2048)]
                        j = self.stand_map[int((bl[b] - 65536) % 2048)]
                    else:
                        i = self.stand_map[int(bl[b] / 256)]
                        j = self.stand_map[int(bl[b] % 256)]
                    if i == j and not include_auto:
                        ## Skip auto-correlations
                        continue
                    ri = np.where(self.stands == i)[0][0]
                    rj = np.where(self.stands == j)[0][0]
                    baselines.append( (ri,rj) )
                    select.append( b )
                    
                # Update the phase center information
                aa = self.get_antennaarray()
                aa.date = jd[0] - astro.DJD_OFFSET
                phase_center.compute(aa)
                
                # Build the output data set
                dataSet = VisibilityDataSet(jd[0], self.freq*1.0, baselines=baselines, 
                                            uvw=uvw[:,select,:].transpose(1,0,2), 
                                            antennaarray=aa, 
                                            phase_center=phase_center)
                for p,l in enumerate(self.pols):
                    name = NUMERIC_STOKES[l]
                    polDataSet = PolarizationDataSet(name, data=vis[select,:,p], weight=wgt[select,:,p])
                    dataSet.append(polDataSet)
                dataSets.append( dataSet )
                
        # Sort
        if sort:
            dataSets.sort()
            
        # Prune
        if min_uv != 0 or max_uv != np.inf:
            dataSets = dataSets.get_uv_range(min_uv=min_uv, max_uv=max_uv)
            
        # Prune a different way
        if len(dataSets) == 1:
            dataSets = dataSets.pop()
            
        # Return
        return dataSets


try:
    from casacore.tables import table
    
    # Stokes codes for CASA Measurement Sets
    NUMERIC_STOKESMS = {1:'I',   2:'Q',   3:'U',   4:'V', 
                        5:'RR',  6:'RL',  7:'LR',  8:'LL',
                        9:'XX', 10:'XY', 11:'YX', 12:'YY'}
    
    class CorrelatedDataMS(CorrelatedDataBase):
        """
        Class to make accessing information about a MS easy.  This wraps 
        all of the "messy" machinery needed to extract both the metadata and data 
        from the file and return them as common LSL objects.
        
        This class has three main attributes to interact with:
         * get_antennaarray - Return a :class:`lsl.sim.vim.AntennaArray` instance
                              that represents the array where the data was obtained.
                              This is useful for simulation proposes and computing 
                              source positions.
         * get_observer - Return a ephem.Observer instance representing the array
         * get_data_set - Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
                      :class:`lsl.imaging.data.VisibilityData` object of all 
                      baselines for a given set of observations
        
        The class also includes a variety of useful metadata attributes:
         * pols - Numpy array of polarization product codes
         * freq - Numpy array of frequency channels in Hz
         * station - LSL :class:`lsl.common.stations.LWAStation` instance for the
                     array
         * date_obs - Datetime object for the reference date of the FIT IDI file
         * antennas - List of :class:`lsl.common.stations.Antenna` instances
        
        .. note::
            The CorrelatedDataMS.antennas attribute should be used over 
            CorrelatedDataMS.station.antennas since the mapping in the MS
            may not be the same as the digitizer order.
        """
        
        def __init__(self, filename):
            """
            Initialize a new CorrelatedData instance from a MS and fill 
            in the metadata.
            """
            
            CorrelatedDataBase.__init__(self, filename)
            
            if not os.path.isdir(self.filename) and tarfile.is_tarfile(self.filename):
                # LASI generate compressed tarballs that contain the MS.  Deal with 
                # those in a transparent manner by automatically unpacking them
                tempdir = tempfile.mkdtemp(prefix='CorrelatedMS-')
                tf = tarfile.open(self.filename, mode='r:*')
                tf.extractall(tempdir)
                tf.close()
                
                # Find the first directory that could be a MS
                self.filename = None
                for path in os.listdir(tempdir):
                    path = os.path.join(tempdir, path)
                    if os.path.isdir(path):
                        self.filename = path
                        break
                        
                # Clean up the temporary directory when the script exists
                atexit.register(lambda: shutil.rmtree(tempdir))
                
            # Open the various tables that we need
            data = table(self.filename, ack=False)
            try:
                ants = table(os.path.join(self.filename, 'ANTENNA'), ack=False)
            except:
                raise RuntimeError(f"Cannot find table 'ANTENNA' in '{self.filename}'")
            try:
                pols = table(os.path.join(self.filename, 'POLARIZATION'), ack=False)
            except:
                raise RuntimeError(f"Cannot find table 'POLARIZATION' in '{self.filename}'")
            try:
                obs = table(os.path.join(self.filename, 'OBSERVATION'), ack=False)
            except:
                raise RuntimeError(f"Cannot find table 'OBSERVATION' in '{self.filename}'")
            try:
                src = table(os.path.join(self.filename, 'SOURCE'), ack=False)
            except:
                src = None
                warnings.warn(colorfy("{{%%yellow Cannot find table 'SOURCE' in '%s', assuming zenith is the only source}}" % self.filename), RuntimeWarning)
            try:
                fld = table(os.path.join(self.filename, 'FIELD'), ack=False)
            except:
                raise RuntimeError(f"Cannot find table 'FIELD' in '{self.filename}'")
            try:
                spw = table(os.path.join(self.filename, 'SPECTRAL_WINDOW'), ack=False)
            except:
                raise RuntimeError(f"Cannot find table 'SPECTRAL_WINDOW' in '{self.filename}'")
            try:
                dsc = table(os.path.join(self.filename, 'DATA_DESCRIPTION'), ack=False)
            except:
                raise RuntimeError(f"Cannot find table 'DATA_DESCRIPTION' in '{self.filename}'")
                
            # Station/telescope information
            self.telescope = obs.col('TELESCOPE_NAME')[0]
            
            ## Get latitude and longitude for all antennas
            lat = np.array([], dtype=np.float64)
            lng = np.array([], dtype=np.float64)
            elv = np.array([], dtype=np.float64)
            for row in ants.col('POSITION'):
                la,ln,el = stations.ecef_to_geo(*row)
                lat = np.append(lat, la*180/np.pi)
                lng = np.append(lng, ln*180/np.pi)
                elv = np.append(elv, el)
                
            ## Estimate the center of the station
            sLat = robust.mean(lat)
            sLng = robust.mean(lng)
            sElv = robust.mean(elv)
            
            ## Build a preliminayr represenation of the station
            self.station = stations.LWAStation(ants.col('STATION')[0], sLat, sLng, sElv)
            
            ## Fill in the antennas instances
            antennas = []
            for i in range(lat.size):
                enz = self.station.get_enz_offset((lat[i], lng[i], elv[i]))
                sid = int(strip_letters(ants.col('NAME')[i]))
                
                stand = stations.Stand(sid, *enz)
                antennas.append( stations.Antenna(2*(stand.id-1)-1, stand=stand) )
            self.station.antennas = antennas
            
            # Antennas
            self.stand_map = {}
            self.stands = []
            for i,noact in enumerate(ants.col('NAME')):
                noact = int(noact[3:])
                self.stand_map[i] = noact
                self.stands.append(noact)
                
            self.antenna_map = {}
            self.antennas = []
            for stand in self.stands:
                for ant in self.station.antennas:
                    if ant.stand.id == stand and ant.pol == 0:
                        self.antennas.append(ant)
                        self.antenna_map[ant.stand.id] = ant
                        break
                        
            # Polarization and frequency
            self.pols = pols.col('CORR_TYPE')[0]
            self.freq = np.array([], dtype=np.float64)
            for freqIF in spw.col('CHAN_FREQ'):
                self.freq  = np.concatenate([self.freq, freqIF])
                
            # Total baseline count
            self.total_baseline_count = data.nrows()
            
            # Data set times
            self._times = np.unique(data.getcol('TIME'))
            jd = self._times[0] / 86400.0 + astro.MJD_OFFSET
            self.date_obs = datetime.utcfromtimestamp(astro.utcjd_to_unix(jd))
            self.station.date = astro.unix_to_utcjd(timegm(self.date_obs.timetuple())) \
                                - astro.DJD_OFFSET
            
            # Data set sources
            self._sources = []
            if src is not None:
                for sdir,name in zip(src.col('DIRECTION'), src.col('NAME')):
                    sdir = sdir.flatten()
                    self._sources.append( aipy.amp.RadioFixedBody(*sdir, name=name, epoch=ephem.J2000) )
            else:
                self._sources.append('z')
                
            # Data set fields
            self._fields = []
            for s in fld.col('SOURCE_ID'):
                self._fields.append( s )
                
            # Data set spectral windows
            self._windows = []
            for d in dsc.col('SPECTRAL_WINDOW_ID'):
                self._windows.append( d )
                
            # Integration count
            jd = np.array(self._times) / 3600.0 / 24.0 + astro.MJD_OFFSET
            self.integration_count = len(np.unique(jd))
            
            # Close
            data.close()
            ants.close()
            pols.close()
            obs.close()
            if src is not None:
                src.close()
            fld.close()
            spw.close()
            dsc.close()
            
        def get_data_set(self, sets, include_auto=False, sort=True, min_uv=0, max_uv=np.inf):
            """
            Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
            :class:`lsl.imaging.data.VisibilityData` object for the specified data set.  
            By default this excludes the autocorrelations.  To include 
            autocorrelations set the value of 'include_auto' to True.  Setting the
            'sort' keyword to False will disable the baseline sorting.  Optionally,
            baselines with lengths between min_uv and max_uv can only be returned.

            .. note::
                min_uv and max_uv should be specified in lambda
            """
            
            # Open the data table
            data = table(self.filename, ack=False)
            
            dataSets = VisibilityData()
            try:
                len(sets)
            except TypeError:
                sets = range(sets, sets+1)
            for set in sets:
                # Set the time to look for
                targetTime = self._times[set-1]
                targetJD = targetTime / 3600.0 / 24.0 + astro.MJD_OFFSET
                
                # Pull out the data
                targetData = data.query('TIME == %.16f AND FLAG_ROW == false' % targetTime, sortlist='DATA_DESC_ID,ANTENNA1,ANTENNA2')
                uvw  = targetData.getcol('UVW')
                try:
                    wgt  = None
                    wgtS = targetData.getcol('WEIGHT SPECTRUM')
                except RuntimeError:
                    wgt  = targetData.getcol('WEIGHT')
                    wgtS = None
                ant1 = targetData.getcol('ANTENNA1')
                ant2 = targetData.getcol('ANTENNA2')
                try:
                    vis = targetData.getcol('CORRECTED_DATA')
                except RuntimeError:
                    vis  = targetData.getcol('DATA')
                fld  = targetData.getcol('FIELD_ID')
                dsc  = targetData.getcol('DATA_DESC_ID')
                
                # Figure out the source we are working on and create a phase center
                # if there is only a single source
                phase_center = None
                fld = np.unique(fld)
                if len(fld) == 1:
                    phase_center = self._sources[self._fields[fld[0]]]
                    
                # Compute the full uvw coordinates
                uvw.shape += (1,)
                uscl = self.freq / vLight
                uscl.shape = (1,1)+uscl.shape
                uvw = uvw*uscl
                
                # Expand the weights, if needed
                if wgtS is None:
                    wgtS = np.ones(vis.shape, dtype=wgt.dtype)
                    for i in range(wgtS.shape[1]):
                        wgtS[:,i,:] = wgt
                    wgt = wgtS
                    
                # Update the phase center information
                aa = self.get_antennaarray()
                aa.date = targetJD - astro.DJD_OFFSET
                phase_center.compute(aa)
                
                # Setup the output data
                # NOTE: This assumes that the data are stored in an array that 
                #       is window, basline per the lsl.writer.measurementset 
                #       module
                baselines = []
                select = []
                for b,(a1,a2) in enumerate(zip(ant1,ant2)):
                    if a1 == a2 and not include_auto:
                        ## Skip auto-correlations
                        continue
                    baselines.append( (a1,a2) )
                    select.append( b )
                if len(self._windows) > 1:
                    baselines = baselines[:len(baselines)//2]
                    selectU = select[:len(baselines)]
                else:
                    selectU = select
                # Build the output data set
                dataSet = VisibilityDataSet(targetJD, self.freq*1.0, baselines=baselines, 
                                            uvw=uvw[selectU,:,:], 
                                            antennaarray=aa, 
                                            phase_center=phase_center)
                for p,l in enumerate(self.pols):
                    name = NUMERIC_STOKESMS[l]
                    subvis = vis[select,:,p]
                    subwgt = wgt[select,:,p]
                    
                    # Deal with multiple spectral windows
                    # NOTE: This assumes that the data are stored in an array that 
                    #       is window, baseline per the lsl.writer.measurementset 
                    #       module
                    if len(self._windows) > 1:
                        nBand = len(self._windows)
                        subvis.shape = (nBand,subvis.shape[0]//nBand,subvis.shape[1])
                        subwgt.shape = subvis.shape
                        
                        subvis = subvis.transpose(1,0,2)
                        subwgt = subwgt.transpose(1,0,2)
                        
                        subvis = subvis.reshape(subvis.shape[0],subvis.shape[1]*subvis.shape[2])
                        subwgt = subwgt.reshape(*subvis.shape)
                        
                    polDataSet = PolarizationDataSet(name, data=subvis, weight=subwgt)
                    dataSet.append(polDataSet)
                dataSets.append( dataSet )
            # Close
            data.close()
            
            # Sort
            if sort:
                dataSets.sort()
                
            # Prune
            if min_uv != 0 or max_uv != np.inf:
                dataSets = dataSets.get_uv_range(min_uv=min_uv, max_uv=max_uv)
                
            # Prune a different way
            if len(dataSets) == 1:
                dataSets = dataSets.pop()
                
            # Return
            return dataSets
            
except ImportError:
    warnings.warn(colorfy('{{%yellow Cannot import casacore.tables, MS support disabled}}'), RuntimeWarning)
    
    class CorrelatedDataMS(object):
        """
        Class to make accessing information about a MS easy.  This wraps 
        all of the "messy" machinery needed to extract both the metadata and data 
        from the file and return them as common LSL objects.
        """
        
        def __init__(self, filename):
            raise RuntimeError("Cannot import casacore.tables, MS support disabled")


def convert_to_stokes(data_set):
    """
    Given a :class:`lsl.imaging.data.VisibilityDataSet` object that contains
    linear polarization products, convert to Stokes parameters and return a
    new :class:`lsl.imaging.data.VisibilityDataSet` containing them.
    """
    
    # Catch VisibilityData objects before we go any further so we can iterate 
    # over them
    if isinstance(data_set, VisibilityData):
        new_data = VisibilityData()
        for ds in data_set:
            new_data_set = convert_to_stokes(ds)
            new_data.append(new_data_set)
        return new_data
        
    if 'I' in data_set.pols or 'Q' in data_set.pols or 'U' in data_set.pols or 'V' in data_set.pols:
        raise RuntimeError("Data already appear to be represented as Stokes parameters")
        
    pairs = 0
    if 'XX' in data_set.pols and 'YY' in data_set.pols:
        pairs += 1
    if 'XY' in data_set.pols and 'YX' in data_set.pols:
        pairs += 1
    if pairs == 0:
        raise RuntimeError("Too few linear polarization products to form any Stokes parameters")
        
    try:
        XX = data_set.XX
        YY = data_set.YY
        
        Id = XX.data + YY.data
        Iw = (XX.weight + YY.weight)/2.0
        Im = XX.mask | YY.mask
        
        Qd = XX.data - YY.data
        Qw = (XX.weight + YY.weight)/2.0
        Qm = XX.mask | YY.mask
    except AttributeError:
        Id = Iw = Im = None
        Qd = Qw = Qm = None
        
    try:
        XY = data_set.XY
        YX = data_set.YX
        
        Ud = XY.data + YX.data
        Uw = (XY.weight + YX.weight)/2.0
        Um = XY.mask | YX.mask
        
        Vd = 1j*(XY.data - YX.data)
        Vw = (XY.weight + YX.weight)/2.0
        Vm = XY.mask | YX.mask
    except AttributeError:
        Ud = Uw = Um = None
        Vd = Vw = Vm = None
        
    new_data_set = data_set.copy(include_pols=False)
    for p,d,w,m in zip(('I','Q','U','V'), (Id,Qd,Ud,Vd), (Iw,Qw,Uw,Vw), (Im,Qm,Um,Vm)):
        if d is None:
            continue
            
        pol = PolarizationDataSet(p, data=d, weight=w, mask=m)
        new_data_set.append(pol)
        
    return new_data_set


def convert_to_linear(data_set):
    """
    Given a :class:`lsl.imaging.data.VisibilityDataSet` object that contains
    Stokes parameters, convert to linear polarization products and return a
    new :class:`lsl.imaging.data.VisibilityDataSet` containing them.
    """
    
    # Catch VisibilityData objects before we go any further so we can iterate 
    # over them
    if isinstance(data_set, VisibilityData):
        new_data = VisibilityData()
        for ds in data_set:
            new_data_set = convert_to_linear(ds)
            new_data.append(new_data_set)
        return new_data
        
    if 'XX' in data_set.pols or 'YY' in data_set.pols or 'XY' in data_set.pols or 'YX' in data_set.pols:
        raise RuntimeError("Data already appear to be represented as linear polarization products")
        
    pairs = 0
    if 'I' in data_set.pols and 'Q' in data_set.pols:
        pairs += 1
    if 'U' in data_set.pols and 'V' in data_set.pols:
        pairs += 1
    if pairs == 0:
        raise RuntimeError("Too few Stokes parameters to form any linear polarization products")
        
    try:
        I = data_set.I
        Q = data_set.Q
        
        XXd = (I.data + Q.data)/2.0
        XXw = (I.weight + Q.weight)/2.0
        XXm = I.mask | Q.mask
        
        YYd = (I.data - Q.data)/2.0
        YYw = (I.weight + Q.weight)/2.0
        YYm = I.mask | Q.mask
    except AttributeError:
        XXd = XXw = XXm = None
        YYd = YYw = YYm = None
    
    try:
        U = data_set.U
        V = data_set.V
        
        XYd = (U.data - 1j*V.data)/2.0
        XYw = (U.weight + V.weight)/2.0
        XYm = U.mask | V.mask
        
        YXd = (U.data + 1j*V.data)/2.0
        YXw = (U.weight + V.weight)/2.0
        YXm = U.mask | V.mask
    except AttributeError:
        XYd = XYw = XYm = None
        YXd = YXw = YXm = None
        
    new_data_set = data_set.copy(include_pols=False)
    for p,d,w,m in zip(('XX','YY','XY','YX'), (XXd,YYd,XYd,YXd), (XXw,YYw,XYw,YXw), (XXm,YYm,XYm,YXm)):
        if d is None:
            continue
            
        pol = PolarizationDataSet(p, data=d, weight=w, mask=m)
        new_data_set.append(pol)
        
    return new_data_set


class ImgWPlus(aipy.img.ImgW):
    """
    Sub-class of the aipy.img.ImgW class that adds support for different 
    visibility weighting scheme and uv plane tapering.  This class also
    adds in a couple of additional methods that help determine the size of
    the field of view and the pixels near the phase center.
    """
    
    _mjd = None
    _phase_center = None
    _antennaarray = None
    
    def __iadd__(self, uv, bm=None):
        if isinstance(uv, ImgWPlus):
            if self.shape != uv.shape:
                raise ValueError(f"Shape mis-match: {self.shape} != {uv.shape}")
            if len(self.bm) != len(uv.bm):
                raise ValueError(f"Order mis-match: {len(self.bm)} != {len(uv.bm)}")
                
            self.uv += uv.uv
            for i in range(len(self.bm)):
                self.bm[i] += uv.bm[i]
                
        elif isinstance(uv, np.ndarray):
            if not isinstance(bm, (list, tuple)):
                raise ValueError("Expected bm to be a list or tuple")
            if self.uv.shape != uv.shape:
                raise ValueError(f"Shape mis-match: {self.uv.shape} != {uv.shape}")
            if len(self.bm) != len(bm):
                raise ValueError(f"Order mis-match: {len(self.bm)} != {len(bm)}")
                
            self.uv += uv
            for i in range(len(self.bm)):
                self.bm[i] += bm[i]
                
        else:
            raise TypeError("Expected an ImgWPlus or np.ndarray instance")
            
    @property
    def field_of_view(self):
        """
        Return the approximate size of the field of view in radians.  The 
        field of view calculate is based off the maximum and minimum values
        of L found for the inverted uv matrix.
        """
        
        im_fov = 360/np.pi * np.arcsin(np.min([1.0, 1/self.res/2]))
        return im_fov * np.pi/180
        
    @property
    def pixel_size(self):
        """
        Return the approximate size of pixels at the phase center in radians.
        """
        
        im_res = 360/np.pi * np.arcsin(1/self.size/2)
        return im_res * np.pi/180
        
    def _gen_img(self, data, center=(0,0), weighting='natural', robust=0.0, taper=(0.0, 0.0)):
        """
        Return the inverse FFT of the provided data, with the 0,0 point 
        moved to 'center'.  In the images return north is up and east is 
        to the left.
        
        There are a few keywords that control how the image is formed.  
        There are:
          * weighting - The weighting scheme ('natural', 'uniform', or 
                        'briggs') used on the data;
          * robust - The value for the weighting robustness under the 
                     'briggs' method; and
          * taper - The size of u and v Gaussian tapers at the 30% level.
        """
        
        # Make sure that we have a valid weighting scheme to use
        if weighting not in ('natural', 'uniform', 'briggs'):
            raise ValueError(f"Unknown weighting scheme '{weighting}'")
            
        # Apply the weighting
        if weighting == 'natural':
            ## Natural weighting - we already have it
            pass
        
        elif weighting == 'uniform':
            ## Uniform weighting
            wgt = self.bm[0].real
            data *= 1.0/wgt
            data[np.where(wgt == 0)] = 0.0
            
        elif weighting == 'briggs':
            ## Robust weighting - we need to calculate it
            wgt = self.bm[0].real
            wavg = wgt.sum()
            w2avg = (wgt**2).sum()
            S2 = 25 * 10**(-2*robust) / (w2avg / wavg)
            data *= 1.0 / (1.0 + wgt*S2)
            data[np.where(wgt == 0)] = 0.0
            
        # Make sure that we have the right type to taper with
        try:
            taper1 = taper[0]
            taper2 = taper[1]
            taper = (taper1, taper2)
        except TypeError:
            taper = (taper, taper)
            
        # Apply the taper
        if taper[0] > 0.0 or taper[1] > 0.0:
            u,v = self.get_uv()
            
            taper1 = 1.0
            if taper[0] > 0.0:
                cu = np.log(0.3) / taper[0]**2
                taper1 = np.exp(cu*u**2)
                
            taper2 = 1.0
            if taper[1] > 0.0:
                cv = np.log(0.3) / taper[1]**2
                taper2 = np.exp(cv*v**2)
                
            data = data*taper1*taper2
            
        return aipy.img.recenter(ifft2Function(data).real.astype(np.float32)/self.kern_corr, center)
        
    def image(self, center=(0,0), weighting='natural', robust=0.0, taper=(0.0, 0.0)):
        """Return the inverse FFT of the UV matrix, with the 0,0 point moved
        to 'center'.  In the images return north is up and east is 
        to the left.
        
        There are a few keywords that control how the image is formed.  
        There are:
          * weighting - The weighting scheme ('natural', 'uniform', or 
                        'briggs') used on the data;
          * robust - The value for the weighting robustness under the 
                     'briggs' method; and
          * taper - The size of u and v Gaussian tapers at the 30% level.
        """
        
        return self._gen_img(self.uv, center=center, weighting=weighting, robust=robust, taper=taper)
        
    def bm_image(self, center=(0,0), term=None, weighting='natural', robust=0.0, taper=(0.0, 0.0)):
        """Return the inverse FFT of the sample weightings (for all mf_order
        terms, or the specified term if supplied), with the 0,0 point
        moved to 'center'.  In the images return north is up and east is 
        to the left.
        
        There are a few keywords that control how the image is formed.  
        There are:
          * weighting - The weighting scheme ('natural', 'uniform', or 
                       'briggs') used on the data;
          * robust - The value for the weighting robustness under the 
                     'briggs' method; and
          * taper - The size of u and v Gaussian tapers at the 30% level.
        """
        
        if not term is None:
            return self._gen_img(self.bm[term], center=center, weighting=weighting, robust=robust, taper=taper)
        else:
            return [self._gen_img(b, center=center, weighting=weighting, robust=robust, taper=taper) for b in self.bm]
            
    @property
    def mjd(self):
        """
        Observation/image MJD.
        """
        
        return self._mjd
        
    @mjd.setter
    def mjd(self, value):
        if isinstance(value, AstroTime):
            self._mjd = value.utc.mjd
        elif isinstance(value, (int, np.integer, float, np.floating)):
            self._mjd = float(value)
        else:
            raise TypeError("Expected an object of type astropy.time.Time, int, or float")
            
    @property
    def phase_center(self):
        """
        Observation/image phase center.
        """
        
        return self._phase_center
        
    @phase_center.setter
    def phase_center(self, value):
        if isinstance(value, (SkyCoord, ICRS, FK5, ephem.Body, astro.equ_posn)):
            self._phase_center = value
        else:
            raise TypeError("Expected an object of type SkyCoord, ICRS, FK5, ephem.Body, or astro.equ_posn")
            
    @property
    def antennaarray(self):
        """
        Observation/image antenna array.
        """
        
        return self._antennaarray
        
    @antennaarray.setter
    def antennaarray(self, value):
        if isinstance(value, simVis.AntennaArray):
            self._antennaarray = value
        else:
            raise TypeError("Expected an object of type sim.vis.AntennaArray")
            
    @property
    def wcs(self):
        """
        World Coordinate System (WCS) associated with the image in equatorial
        (RA/dec.) coordinates.
        """
        
        if self._mjd is None:
            raise RuntimeError("Cannot create WCS: unknown observation MJD")
        if self._phase_center is None:
            raise RuntimeError("Cannot create WCS: unknown phase center")
            
        # Extract what we need from the phase center
        obs_t = AstroTime(self._mjd, format='mjd', scale='utc')
        if isinstance(self._phase_center, astro.equ_posn):
            pc_ra = self._phase_center.ra
            pc_dec = self._phase_center.dec
            pc_equinox = 'J2000'
            pc = FK5(ra=pc_ra*astrounits.deg, dec=pc_dec*astrounits.deg, equinox=pc_equinox)
            pc = self._phase_center.transform_to(FK5(equinox=obs_t))
        else:
            try:
                pc = self._phase_center.transform_to(FK5(equinox=obs_t))
                pc_ra = pc.ra.deg
                pc_dec = pc.dec.deg
            except AttributeError:
                pc_ra = np.degrees(self._phase_center.ra)
                pc_dec = np.degrees(self._phase_center.dec)
                
        return WCS(header={'WCSAXES': 2,
                           'RADESYS': 'FK5',
                           'EQUINOX': obs_t.jyear,
                           'CTYPE1':  'RA---SIN',
                           'CRPIX1':  np.round(self.size/self.res) // 2 + 1,
                           'CDELT1':  -1 * np.degrees(self.pixel_size),
                           'CRVAL1':  pc_ra,
                           'CUNIT1':  'deg',
                           'CTYPE2':  'DEC--SIN',
                           'CRPIX2':  np.round(self.size/self.res) // 2 + 1,
                           'CDELT2':  np.degrees(self.pixel_size),
                           'CRVAL2':  pc_dec,
                           'CUNIT2':  'deg',
                           'LONPOLE': 180.0,
                           'LATPOLE': 90.0})


def build_gridded_image(data_set, size=80, res=0.50, im_size=None, im_res=None, wres=0.10, pol='XX',chan=None, im=None, verbose=True):
    """
    Given a :class:`lsl.imaging.data.VisibilityDataSet` object, build an aipy.img.ImgW 
    object of gridded uv data which can be used for imaging.  The ImgW object 
    itself is returned by this function to make it more versatile.
      * size and res are the aperture plane size and resolution.
      
      * im_size and im_res are the image plane size and resolution.
    
    """
            
    # Catch VisibilityData objects before we go any further so we can iterate 
    # over them
    if isinstance(data_set, VisibilityData):
        for ds in data_set:
            im = build_gridded_image(ds, size=size, res=res, im_size=im_size, im_res=im_res, wres=wres, 
                                     pol=pol, chan=chan, im=im, verbose=verbose)
        return im

    # Catch Input Shapes
    if (im_size is not None) and (im_res is not None): # Case: User provides Image plane inputs
        res = (2 * im_size * np.sin(np.pi * im_res / 360))**-1
        size = im_size * res
        
    # Make sure we have the right kinds of objects
    ## im
    if im is None:
        im = ImgWPlus(size=size, res=res, wres=wres)
    elif not isinstance(im, ImgWPlus):
        raise TypeError("Expected 'im' to be a ImgWPlus instance")
    if im.size != size:
        raise ValueError("Size of 'im' does not match 'size' keyword")
    if im.res != res:
        raise ValueError("Resolution of 'im' does not match 'res' keyword")
    if im.wres != wres:
        raise ValueError("w Resolution of 'im' does not match 'wres' keyword")
        
    ## data_set
    if not isinstance(data_set, VisibilityDataSet):
        raise TypeError("Expected data to be stored in a VisibiltyData or VisibilityDataSet object")
        
    # Make sure we have the right polarization
    if pol not in data_set.pols:
        raise RuntimeError(f"Data dictionary does not have data for polarization '{pol}'")
    pds = getattr(data_set, pol)
        
    if chan is not None:
        # Make sure that `chan' is an array by trying to find its length
        try:
            len(chan)
        except TypeError:
            chan = [chan]
            
        # Build up the data using only the specified channels
        uvw = data_set.uvw[:,:,chan]
        vis = pds.data[:,chan]
        wgt = pds.weight[:,chan]
        
    else:
        uvw = data_set.uvw
        vis = pds.data
        wgt = pds.weight
        
    uvw = np.concatenate(uvw, axis=1)
    vis = np.concatenate(vis)
    wgt = np.concatenate(wgt)
    
    with SilentVerbose(stdout=not verbose):
        uvw, vis, wgt = im.append_hermitian(uvw, vis, wgts=wgt)
        u,v,w = uvw
        order = np.argsort(w)
        u,v,w = u.take(order), v.take(order), w.take(order)
        vis,wgt = vis.take(order), np.array([wg.take(order) for wg in wgt]).squeeze()
        if wgt.dtype != np.complex64:
            wgt = wgt.astype(np.complex64)
            
        im.uv, im.bm[0], im.kern_corr = WProjection(u, v, w, vis, wgt, size, np.float64(res), np.float64(wres))
        
    # Update MJD, phase center, and antenna array information
    im.mjd = data_set.mjd
    im.phase_center = data_set.phase_center
    im.antennaarray = data_set.antennaarray
    
    return im


def plot_gridded_image(ax, gimg, shifted=True, origin='lower', interpolation='nearest', **kwargs):
    """
    Given a matplotlib axes instance with projection set to a astropy.wcs.WCS
    and a gridded image generated by the build_gridded_image() function, plot
    the image on the axes.  This function returns the matplotlib object added to
    the plot
    
    .. versionchanged:: 3.0.0
       This function now assumes that gimg.wcs has been used to initalize the
       axes being used.
    
    .. versionchanged:: 1.2.1
        Changed the function to return the matplotlib object plotted so
        that colorbars can be added
        
    .. versionchanged:: 1.1.0
        Added a 'shifted' keyword to control whether or not the image
        is centered or not.
        
    .. versionadded:: 1.0.1
    """
    
    # Build the unshifted image
    img = gimg.image()
    
    # Shift the image so that it is centered in the frame
    if shifted:
        imgSize = img.shape[0]	# should be square
        img = np.roll(img, imgSize//2, axis=0)
        img = np.roll(img, imgSize//2, axis=1)
        
    # Plot
    return ax.imshow(img, origin=origin, interpolation=interpolation, **kwargs)


def get_image_radec(gimg, shifted=True):
    """
    Given a gridded image generated by the build_gridded_image() function
    and an AntennaArray instance, return a two-element tuple containing
    the RA and dec. values (in radians) for each pixel in the image.  
    
    .. versionchanged:: 3.0.0
        Dropped the `aa` argument and the `phase_center` keyword.
        
    .. versionadded: 1.1.0
    """
    
    top = gimg.get_top()
    im_shape = top[0].shape
    x, y = np.meshgrid(np.arange(im_shape[1]), np.arange(im_shape[0]))
    
    ra, dec = gimg.wcs.all_pix2world(x.ravel(), y.ravel(), 0)
    ra = ra.reshape(im_shape) * np.pi/180
    dec = dec.reshape(im_shape) * np.pi/180
    
    # Shift, if needed
    if not shifted:
        raSize = ra.shape[0]	# should be square
        ra = np.roll(ra, raSize//2, axis=0)
        ra = np.roll(ra, raSize//2, axis=1)
        
        decSize = dec.shape[0]	# should be square
        dec = np.roll(dec, decSize//2, axis=0)
        dec = np.roll(dec, decSize//2, axis=1)
        
    # Done
    return ra, dec


def get_image_azalt(gimg, shifted=True):
    """
    Given a gridded image generated by the build_gridded_image() function
    and an AntennaArray instance, return a two-element tuple containing
    the azimuth and altitude, both in radians, for each pixel in the image.
    
    .. versionchanged:: 3.0.0
        Dropped the `phase_center` keyword.
        
    .. versionadded: 1.1.0
    """
    
    # Get the RA and dec. coordinates for each pixel
    ot = AstroTime(gimg.mjd, format='mjd', scale='utc')
    ra, dec = get_image_radec(gimg, shifted=shifted)
    
    # Convert to azimuth and altitude
    eq = FK5(ra*astrounits.rad, dec*astrounits.rad, equinox=ot)
    tc = eq.transform_to(AltAz(location=gimg.antennaarray.earth_location, obstime=ot))
    az, alt = tc.az.rad, tc.alt.rad
    
    # Done
    return az, alt
