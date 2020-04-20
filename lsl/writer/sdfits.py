"""
Module for writing spectrometer output to a SDFITS file.  The SDFITS created by this 
modulefiles closely follow the Parkes variant of the SDFITS convention 
(see [http://fits.gsfc.nasa.gov/registry/sdfits.html]).  The main differences between 
the two is that the LWA SDFITS do not contain calbration or weather information.  This,
however, should not stop the files from being loaded into CASA with the ATNF Spectral 
Analysis Package (ASAP).

.. versionchanged:: 0.5.0
    The classes and functions defined in this module are based heavily off 
    the :mod:`lsl.writer.fitsidi` writer.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import gc
import re
import numpy
import warnings
from astropy.io import fits as astrofits
from datetime import datetime
from functools import cmp_to_key

from lsl import astro
from lsl.common.stations import lwa1
from lsl.writer.fitsidi import WriterBase

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.5'
__all__ = ['Sd', 'STOKES_CODES', 'NUMERIC_STOKES']


SD_VERSION = (1, 6)

STOKES_CODES = { 'I':  1,  'Q': 2,   'U':  3,  'V':  4, 
                'RR': -1, 'LL': -2, 'RL': -3, 'LR': -4, 
                'XX': -5, 'YY': -6, 'XY': -7, 'YX': -8}

NUMERIC_STOKES = { 1: 'I',   2: 'Q',   3: 'U',   4: 'V', 
                  -1: 'RR', -2: 'LL', -3: 'RL', -4: 'RL', 
                  -5: 'XX', -6: 'YY', -7: 'XY', -8: 'YX'}


class Sd(WriterBase):
    """
    Class for storing spectrometer data and writing the data, along with array
    frequency setup, etc., to a SDFITS file that can be read into CASA via the
    sd.scantable() function.
    """
   
    _STOKES_CODES = STOKES_CODES
    
    class _SpectrometerData(object):
        """
        Represents one spectrum for a given observation time.
        """
        
        def __init__(self, obsTime, intTime, dataDict, pol=STOKES_CODES['XX']):
            self.obsTime = obsTime
            self.intTime = intTime
            self.dataDict = dataDict
            self.pol = pol
        
        def time(self):
            return self.obsTime
            
    def __init__(self, filename, ref_time=0.0, verbose=False, memmap=None, clobber=False):
        """
        Initialize a new SDFITS object using a filename and a reference time 
        given in seconds since the UNIX 1970 ephem, a python datetime object, or a 
        string in the format of 'YYYY-MM-DDTHH:MM:SS'.
        
        .. versionchanged:: 1.1.2
            Added the 'memmap' and 'clobber' keywords to control if the file
            is memory mapped and whether or not to overwrite an existing file, 
            respectively.
        """

        # File-specific information
        WriterBase.__init__(self, filename, ref_time=ref_time, verbose=verbose)
        
        # Observation-specific information
        self.site = lwa1
        
        # Misc.
        self.tSys = 250
        self.observer = 'UKNOWN'
        self.project = 'UNKNOWN'
        self.mode = 'UNKNOWN'
        

        # Open the file and get going
        if os.path.exists(filename):
            if clobber:
                os.unlink(filename)
            else:
                raise IOError("File '%s' already exists" % filename)
        self.FITS = astrofits.open(filename, mode='append', memmap=memmap)

    def set_site(self, site):
        """
        Set the TELESCOP keyword in the primary HDU using an :class:`lsl.common.stations.LWAStation`
        object.
        """

        self.site = site
        
    def set_observer(self, observer, project='UNKNOWN', mode='UNKNOWN'):
        """
        Set the observer name, project, and observation mode (if given) to the 
        self.observer, self.project, and self.mode attributes, respectively.
        """
        
        self.observer = observer
        self.project = project
        self.mode = mode
        
    def add_comment(self, comment):
        """
        Add a comment to data.
        
        .. versionadded:: 1.2.4
        """
        
        try:
            self._comments.append( comment )
        except AttributeError:
            self._comments = [comment,]
            
    def add_history(self, history):
        """
        Add a history entry to the data.
        
        .. versionadded:: 1.2.4
        """
        
        try:
            self._history.append( history )
        except AttributeError:
            self._history = [history,]
            
    def add_data_set(self, obsTime, intTime, beam, data, pol='XX'):
        """
        Create a SpectrometerData object to store a collection of spectra.
        """
        
        if type(pol) == str:
            numericPol = STOKES_CODES[pol.upper()]
        else:
            numericPol = pol
        
        dataDict = {}
        for i,b in enumerate(beam):
            dataDict[b] = data[i,:]

        self.data.append( self._SpectrometerData(obsTime, intTime, dataDict, pol=numericPol) )

    def write(self):
        """
        Fill in the SDFITS file will all of the tables in the correct order.
        """
        
        def __sortData(x, y):
            """
            Function to sort the self.data list in order of time and then 
            polarization code.
            """
            
            xID = x.obsTime*10000000 + abs(x.pol)
            yID = y.obsTime*10000000 + abs(y.pol)
            
            if xID > yID:
                return 1
            elif xID < yID:
                return -1
            else:
                return 0
                
        # Sort the data set
        self.data.sort(key=cmp_to_key(__sortData))
        
        self._write_primary_hdu()
        self._write_singledish_hdu()
        
        # Clear out the data section
        del(self.data[:])
        gc.collect()

    def close(self):
        """
        Close out the file.
        """

        self.FITS.flush()
        self.FITS.close()

    def _write_primary_hdu(self):
        """
        Write the primary HDU to file.
        """

        primary = astrofits.PrimaryHDU()
        
        primary.header['NAXIS'] = (0, 'indicates SD file')
        primary.header['EXTEND'] = (True, 'indicates SD file')
        ts = str(astro.get_date_from_sys())
        primary.header['DATE'] = (ts.split()[0], 'IDI file creation date')
        primary.header['ORIGIN'] = 'LSL SDFITS writer'
        primary.header['TELESCOP'] = (self.site.name, 'Telescope name')
        
        # Write the comments and history
        try:
            for comment in self._comments:
                primary.header['COMMENT'] = comment
            del self._comments
        except AttributeError:
            pass
        primary.header['COMMENT'] = " FITS (Flexible Image Transport System) format is defined in 'Astronomy and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"
        try:
            for hist in self._history:
                primary.header['HISTORY'] = hist
            del self._history
        except AttributeError:
            pass
            
        self.FITS.append(primary)
        self.FITS.flush()
        
    def _write_singledish_hdu(self):
        """
        Define the SINGLE DISH table.
        """
        
        scanList = []
        dateList = []
        timeList = []
        intTimeList = []
        beamList = []
        mList = []
        rawList = []
        scanCount = 1
        for i,dataSet in enumerate(self.data):
            if dataSet.pol == self.stokes[0]:
                tempMList = {}
                for stokes in self.stokes:
                    tempMList[stokes] = {}
        
            beams = list(dataSet.dataDict.keys())
            beams.sort()
            for b in beams:
                specData = dataSet.dataDict[b]
                
                # Load the data into a matrix
                tempMList[dataSet.pol][b] = specData.ravel()
                
                if dataSet.pol == self.stokes[0]:
                    # Observation date and time
                    utc = astro.taimjd_to_utcjd(dataSet.obsTime)
                    date = astro.get_date(utc)
                    date.hours = 0
                    date.minutes = 0
                    date.seconds = 0
                    utc0 = date.to_jd()
                        
                    scanList.append(scanCount)
                    dateList.append('%4i-%02i-%02i' % (date.years, date.months, date.days))
                    timeList.append((utc - utc0)*24*3600)
                    intTimeList.append(dataSet.intTime)
                    beamList.append(b.id)
                    rawList.append(b)
            
            if dataSet.pol == self.stokes[-1]:
                for b in rawList:
                    matrix = numpy.zeros((self.nStokes,self.nChan), dtype=numpy.float32)
                    for p in range(self.nStokes):
                        try:
                            matrix[p,:] = tempMList[self.stokes[p]][b]
                        except KeyError:
                            warnings.warn('Key mis-match %s %s' % (str(b), str(tempMList[self.stokes[p]].keys())), RuntimeWarning)
                            
                    mList.append(matrix.ravel())
                scanCount += 1
                rawList = []
        
        # Scan number
        c1  = astrofits.Column(name='SCAN', format='1I', 
                        array=numpy.array(scanList))
        ## Cycle
        #c2 = astrofits.Column(name='CYCLE', format='1J', 
                        #array=numpy.array([1 for s in scanList]))
        # DATE-OBS
        c3  = astrofits.Column(name='DATE-OBS', format='10A', 
                        array = numpy.array(dateList))
        # Time elapsed since 0h
        c4  = astrofits.Column(name='TIME', format='1D', unit = 's', 
                        array = numpy.array(timeList))
        # Integration time (seconds)
        c5  = astrofits.Column(name='EXPOSURE', format='1E', unit='s', 
                        array=numpy.array(intTimeList, dtype=numpy.float32))
        # Object name
        c6  = astrofits.Column(name='OBJECT', format='16A', 
                        array=numpy.array(['LWA_OBS' for s in scanList]))
        # Object position (deg and deg)
        c7  = astrofits.Column(name='OBJ-RA', format='1D', unit='deg', 
                        array=numpy.array([0.0 for s in scanList]))
        c8  = astrofits.Column(name='OBJ-DEC', format='1D', unit='deg', 
                        array=numpy.array([0.0 for s in scanList]))
        # Rest frequency (Hz)
        c9  = astrofits.Column(name='RESTFRQ', format='1D', unit='Hz', 
                        array=numpy.array([0.0 for s in scanList]))
        # Observation mode
        c10 = astrofits.Column(name='OBSMODE', format='16A', 
                        array=numpy.array([self.mode for s in scanList]))
        # Beam (tuning)
        c11 = astrofits.Column(name='BEAM', format='1I', 
                        array=numpy.array(beamList))
        # IF
        c12 = astrofits.Column(name='IF', format='1I', 
                        array=numpy.array([self.freq[0].id for s in scanList]))
        # Frequency resolution (Hz)
        c13 = astrofits.Column(name='FREQRES', format='1D', unit='Hz', 
                        array=numpy.array([self.freq[0].chWidth for s in scanList]))
        # Bandwidth of the system (Hz)
        c14 = astrofits.Column(name='BANDWID', format='1D', unit='Hz', 
                        array=numpy.array([self.freq[0].totalBW for s in scanList]))
        # Frequency axis - 1
        c15 = astrofits.Column(name='CRPIX1', format='1E',
                        array=numpy.array([self.refPix for s in scanList]))
        c16 = astrofits.Column(name='CRVAL1', format='1D', unit='Hz', 
                        array=numpy.array([self.refVal for s in scanList]))
        c17 = astrofits.Column(name='CDELT1', format='1D', unit='Hz', 
                        array=numpy.array([self.freq[0].chWidth for s in scanList]))
        c18 = astrofits.Column(name='CRVAL3', format='1D', unit='deg', 
                        array=numpy.array([0.0 for s in scanList]))
        # Dec. axis - 4
        c19 = astrofits.Column(name='CRVAL4', format='1D', unit='deg', 
                        array=numpy.array([0.0 for s in scanList]))
        ## Scan rate
        #c20 = astrofits.Column(name='SCANRATE', format='2E', unit='deg/s', 
                        #array=numpy.array([[0,0] for s in scanList]))
                        
        #
        # Calibration information (currently not implemented)
        #
        ## System temperature  *** UNKNOWN ***
        #c21 =  astrofits.Column(name='TSYS', format='2E', unit='K', 
                        #array=numpy.array([[self.tSys,self.tSys] for s in scanList]))
        ## CALFCTR *** UNKNOWN ***
        #c22 =  astrofits.Column(name='CALFCTR', format='2E', unit='K', 
                        #array=numpy.array([[1,1] for s in scanList]))
        
        # Data
        c23 = astrofits.Column(name='DATA', format='%iE' % (self.nStokes*self.nChan), unit='UNCALIB', 
                        array=numpy.array(mList))
                        
        #
        # Data masking table (currently not implemented)
        #
        # Flag table
        #c24 = astrofits.Column(name='FLAGGED', format='%iB' % (self.nStokes*self.nChan), 
                        #array=numpy.array([[0,]*self.nStokes*self.nChan for s in scanList]))
        
        #
        # Calibration information (currently not implemented)
        #
        ## TCAL *** UNKNOWN ***
        #c25 = astrofits.Column(name='TCAL', format='2E', unit='Jy', 
                        #array=numpy.array([[1,1] for s in scanList]))
        ## TCALTIME *** UNKNOWN ***
        #c26 = astrofits.Column(name='TCALTIME', format='16A', 
                        #array=numpy.array(['UNKNOWN' for s in scanList]))
        
        #
        # Pointing information (currently not implemented)
        #
        ## Azimuth *** UNKNOWN ***
        #c27 = astrofits.Column(name='AZIMUTH', format='1E', unit='deg', 
                        #array=numpy.array([0 for s in scanList]))
        ## Elevation *** UNKNOWN ***
        #c28 = astrofits.Column(name='ELEVATIO', format='1E', unit='deg', 
                        #array=numpy.array([0 for s in scanList]))
        ## Parallactic angle *** UNKNOWN ***
        #c29 = astrofits.Column(name='PARANGLE', format='1E', unit='deg', 
                        #array=numpy.array([0 for s in scanList]))
        
        #
        # Focusing information (currently not implemented and probably never will be)
        #
        ## FOCUSAXI *** NOT NEEDED ***
        #c30 = astrofits.Column(name='FOCUSAXI', format='1E', unit='m', 
                        #array=numpy.array([0 for s in scanList]))
        ## FOCUSTAN *** NOT NEEDED ***
        #c31 = astrofits.Column(name='FOCUSTAN', format='1E', unit='m', 
                        #array=numpy.array([0 for s in scanList]))
        ## FOCUSROT *** NOT NEEDED ***
        #c32 = astrofits.Column(name='FOCUSROT', format='1E', unit='deg', 
                        #array=numpy.array([0 for s in scanList]))
        
        #
        # Weather information (currently not implemented)
        #
        ## Ambient temperature *** UNKNOWN ***
        #c33 = astrofits.Column(name='TAMBIENT', format='1E', unit='C', 
                        #array=numpy.array([0 for s in scanList]))
        ## Air pressure *** UNKNOWN ***
        #c34 = astrofits.Column(name='PRESSURE', format='1E', unit='Pa', 
                        #array=numpy.array([0 for s in scanList]))
        ## Humidity *** UNKNOWN ***
        #c35 = astrofits.Column(name='HUMIDITY', format='1E', unit='%', 
                        #array=numpy.array([0 for s in scanList]))
        ## Wind speed *** UNKNOWN ***
        #c36 = astrofits.Column(name='WINDSPEE', format='1E', unit='m/s', 
                        #array=numpy.array([0 for s in scanList]))
        ## Wind direction *** UNKNOWN ***
        #c37 = astrofits.Column(name='WINDDIRE', format='1E', unit='deg', 
                        #array=numpy.array([0 for s in scanList]))
        
        # Gather together all of the needed columns and figure out which ones
        # store the data and flag tables.  This information is needed later to
        # set the appropriate TDIM keywords.
        cs = []
        dataIndex = 0
        #flagIndex = 0
        n = 1
        for i in range(1, 38):
            try:
                cs.append(eval('c%i' % i))
                if eval('c%i.name' %i) == 'DATA':
                    dataIndex = n
                #if eval('c%i.name' %i) == 'FLAGGED':
                    #flagIndex = n
                n += 1
            except NameError:
                pass
        colDefs = astrofits.ColDefs(cs)
        
        # Create the SINGLE DISH table and update its header
        sd = astrofits.BinTableHDU.from_columns(colDefs)
        
        ## Single disk keywords - order seems to matter
        sd.header['EXTNAME'] = ('SINGLE DISH', 'SDFITS table name')
        sd.header['NMATRIX'] = 1
        sd.header['OBSERVER'] = (self.observer, 'Observer name(s)')
        sd.header['PROJID']   = (self.project, 'Project name')
        sd.header['TELESCOP'] = (self.site.name, 'Telescope name')
        x,y,z = self.site.get_geocentric_location()
        sd.header['OBSGEO-X'] = (x, '[m] Antenna ECEF X-coordinate')
        sd.header['OBSGEO-Y'] = (y, '[m] Antenna ECEF Y-coordinate')
        sd.header['OBSGEO-Z'] = (z, '[m] Antenna ECEF Z-coordinate')
        
        sd.header['SPECSYS'] = ('LSRK', 'Doppler reference frame (transformed)')
        sd.header['SSYSOBS'] = ('TOPOCENT', 'Doppler reference frame of observation')
        sd.header['EQUINOX'] = (2000.0, 'Equinox of equatorial coordinates')
        sd.header['RADESYS'] = ('FK5', 'Equatorial coordinate system frame')
        
        ## Data and flag table dimensionality
        sd.header['TDIM%i' % dataIndex] = ('(%i,%i,1,1)' % (self.nChan, self.nStokes))
        #sd.header.set('TDIM%i' % flagIndex, '(%i,%i,1,1)' % (self.nChan, self.nStokes), after='TFORM%i' % flagIndex)
        
        ## Data and flag table axis descriptions
        ### Frequency
        sd.header['CTYPE1'] = ('FREQ', 'axis 1 is FREQ (frequency)')
        sd.header['CDELT1'] = self.freq[0].chWidth
        sd.header['CRPIX1'] = self.refPix
        sd.header['CRVAL1'] = self.refVal
        ### Stokes
        sd.header['CTYPE2'] = ('STOKES', 'axis 2 is STOKES axis (polarization)')
        if self.stokes[0] < 0:
            sd.header['CDELT2'] = -1.0
        else:
            sd.header['CDELT2'] = 1.0
        sd.header['CRPIX2'] = 1.0
        sd.header['CRVAL2'] = float(self.stokes[0])
        ### RA
        sd.header['CTYPE3'] = ('RA', 'axis 3 is RA axis (pointing)')
        sd.header['CRPIX3'] = 1.0
        sd.header['CDELT3'] = -1.0
        ### Dec
        sd.header['CTYPE4'] = ('DEC', 'axis 4 is Dec. axis (pointing)')
        sd.header['CRPIX4'] = 1.0
        sd.header['CDELT4'] = 1.0
        
        self.FITS.append(sd)
        self.FITS.flush()
