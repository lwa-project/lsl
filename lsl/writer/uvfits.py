# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    from functools import cmp_to_key
    
"""
Module for writing correlator output to a UVFITS file.  The classes and 
functions defined in this module are based heavily off the lwda_fits library.

.. note::
    For arrays with between 256 and 2048 antennas, the baseline packing
    follows the MIRIAD convention.
"""

import os
import gc
import re
import math
import ephem
import numpy
from astropy.io import fits as astrofits
from datetime import datetime

from lsl import astro
from lsl.misc import geodesy
from lsl.common import constants
from lsl.writer.fitsidi import StokesCodes
from lsl.misc.total_sorting import cmp_to_total

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['UV', 'StokesCodes', 'NumericStokes', '__version__', '__revision__', '__all__']


UVVersion = (1, 0)


def merge_baseline(ant1, ant2):
    """
    Merge two stand ID numbers into a single baseline.
    """
    
    if ant1 > 255 or ant2 > 255:
        baseline = ant1*2048 + ant2 + 65536
    else:
        baseline = ant1*256 + ant2
        
    return baseline


def split_baseline(baseline):
    """
    Given a baseline, split it into it consistent stand ID numbers.
    """
    
    if baseline >= 65536:
        ant1 = int((baseline - 65536) // 2048)
        ant2 = int((baseline - 65536) % 2048)
    else:
        ant1 = int(baseline // 256)
        ant2 = int(baseline % 256)
        
    return ant1,ant2


class UV(object):
    """
    Class for storing visibility data and writing the data, along with array
    geometry, frequency setup, etc., to a UVFITS file that can be read into 
    AIPS via the UVLOD task.
    """
    
    class _Antenna(object):
        """
        Holds information describing the location and properties of an antenna.
        """
        
        def __init__(self, id, x, y, z, bits=8):
            self.id = id
            self.x = x
            self.y = y
            self.z = z
            self.levels = bits
            self.polA = {'Type': 'X', 'Angle': 0.0, 'Cal': [0.0, 0.0]}
            self.polB = {'Type': 'Y', 'Angle': 90.0, 'Cal': [0.0, 0.0]}
            
        def getName(self):
            return "LWA%03i" % self.id
            
    class _Frequency:
        """
        Holds information about the frequency setup used in the file.
        """
        
        def __init__(self, offset, channelWidth, bandwidth):
            self.id = 1
            self.bandFreq = offset
            self.chWidth = channelWidth
            self.totalBW = bandwidth
            self.sideBand = 1
            self.baseBand = 0
            
    @cmp_to_total
    class _UVData(object):
        """
        Represents one UV visibility data set for a given observation time.
        """
    
        def __init__(self, obsTime, intTime, baselines, visibilities, pol=StokesCodes['XX'], source='z'):
            self.obsTime = obsTime
            self.intTime = intTime
            self.baselines = baselines
            self.visibilities = visibilities
            self.pol = pol
            self.source = source
            
        def __cmp__(self, y):
            """
            Function to sort the self.data list in order of time and then 
            polarization code.
            """
            
            sID = self.obsTime*10000000 + abs(self.pol)
            yID =    y.obsTime*10000000 + abs(   y.pol)
            
            if sID > yID:
                return 1
            elif sID < yID:
                return -1
            else:
                return 0
                
        def time(self):
            return self.obsTime
            
        def get_uvw(self, HA, dec, obs):
            Nbase = len(self.baselines)
            uvw = numpy.zeros((Nbase,3), dtype=numpy.float32)
            
            # Phase center coordinates
            # Convert numbers to radians and, for HA, hours to degrees
            HA2 = HA * 15.0 * numpy.pi/180
            dec2 = dec * numpy.pi/180
            lat2 = obs.lat
            
            # Coordinate transformation matrices
            trans1 = numpy.matrix([[0, -numpy.sin(lat2), numpy.cos(lat2)],
                            [1,  0,               0],
                            [0,  numpy.cos(lat2), numpy.sin(lat2)]])
            trans2 = numpy.matrix([[ numpy.sin(HA2),                  numpy.cos(HA2),                 0],
                            [-numpy.sin(dec2)*numpy.cos(HA2),  numpy.sin(dec2)*numpy.sin(HA2), numpy.cos(dec2)],
                            [ numpy.cos(dec2)*numpy.cos(HA2), -numpy.cos(dec2)*numpy.sin(HA2), numpy.sin(dec2)]])
                    
            for i,(a1,a2) in enumerate(self.baselines):
                # Go from a east, north, up coordinate system to a celestial equation, 
                # east, north celestial pole system
                xyzPrime = a1.stand - a2.stand
                xyz = trans1*numpy.matrix([[xyzPrime[0]],[xyzPrime[1]],[xyzPrime[2]]])
                
                # Go from CE, east, NCP to u, v, w
                temp = trans2*xyz
                uvw[i,:] = numpy.squeeze(temp) / constants.c
                
            return uvw
                
        def argsort(self, mapper=None, shift=16):
            packed = []
            for a1,a2 in self.baselines:
                if mapper is None:
                    s1, s2 = a1.stand.id, a2.stand.id
                else:
                    s1, s2 = mapper[a1.stand.id], mapper[a2.stand.id]
                packed.append( merge_baseline(s1, s2) )
            packed = numpy.array(packed, dtype=numpy.int32)
            
            return numpy.argsort(packed)
            
    def parse_time(self, ref_time):
        """
        Given a time as either a integer, float, string, or datetime object, 
        convert it to a string in the formation 'YYYY-MM-DDTHH:MM:SS'.
        """
        
        # Valid time string (modulo the 'T')
        timeRE = re.compile(r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?')
        
        if type(ref_time) in (int, long, float):
            refDateTime = datetime.utcfromtimestamp(ref_time)
            ref_time = refDateTime.strftime("%Y-%m-%dT%H:%M:%S")
        elif type(ref_time) == datetime:
            ref_time = ref_time.strftime("%Y-%m-%dT%H:%M:%S")
        elif type(ref_time) == str:
            # Make sure that the string times are of the correct format
            if re.match(timeRE, ref_time) is None:
                raise RuntimeError("Malformed date/time provided: %s" % ref_time)
            else:
                ref_time = ref_time.replace(' ', 'T', 1)
        else:
            raise RuntimeError("Unknown time format provided.")
            
        return ref_time
        
    def ref_time2AstroDate(self):
        """
        Convert a reference time string to an :class:`lsl.astro.date` object.
        """
        
        dateStr = self.ref_time.replace('T', '-').replace(':', '-').split('-')
        return astro.date(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), int(dateStr[3]), int(dateStr[4]), float(dateStr[5]))
        
    def __init__(self, filename, ref_time=0.0, verbose=False, memmap=None, clobber=False):
        """
        Initialize a new UVFITS object using a filename and a reference time 
        given in seconds since the UNIX 1970 ephem, a python datetime object, or a 
        string in the format of 'YYYY-MM-DDTHH:MM:SS'.
        
        .. versionchanged:: 1.1.2
            Added the 'memmap' and 'clobber' keywords to control if the file
            is memory mapped and whether or not to overwrite an existing file, 
            respectively.
        """
        
        # File-specific information
        self.filename = filename
        self.verbose = verbose
        
        # Observatory-specific information
        self.siteName = 'Unknown'
        
        # Observation-specific information
        self.ref_time = self.parse_time(ref_time)
        self.nAnt = 0
        self.nChan = 0
        self.nStokes = 0
        self.refVal = 0
        self.refPix = 0
        self.channelWidth = 0
        
        # Parameters that store the meta-data and data
        self.array = []
        self.freq = []
        self.stokes = []
        self.data = []
        
        # Open the file and get going
        if os.path.exists(filename):
            if clobber:
                os.unlink(filename)
            else:
                raise IOError("File '%s' already exists" % filename)
        self.FITS = astrofits.open(filename, mode='append', memmap=memmap)
        
    def set_stokes(self, polList):
        """
        Given a list of Stokes parameters, update the object's parameters.
        """
        
        for pol in polList:
            if type(pol) == str:
                numericPol = StokesCodes[pol.upper()]
            else:
                numericPol = pol
                
            if numericPol not in self.stokes:
                self.stokes.append(numericPol)
                
        # Sort into order of 'XX', 'YY', 'XY', and 'YX' or 'I', 'Q', 'U', and 'V'
        self.stokes.sort()
        if self.stokes[0] < 0:
            self.stokes.reverse()
            
        self.nStokes = len(self.stokes)
        
    def set_frequency(self, freq):
        """
        Given a numpy array of frequencies, set the relevant common observation
        parameters and add an entry to the self.freq list.
        """
        
        self.nChan = len(freq)
        self.refVal = freq[0]
        self.refPix = 1
        self.channelWidth = numpy.abs(freq[1] - freq[0])
        totalWidth = numpy.abs(freq[-1] - freq[0])
        
        freqSetup = self._Frequency(0.0, self.channelWidth, totalWidth)
        self.freq.append(freqSetup)
        
    def set_geometry(self, site, antennas, bits=8):
        """
        Given a station and an array of stands, set the relevant common observation
        parameters and add entries to the self.array list.
        
        .. versionchanged:: 0.4.0
            Switched over to passing in Antenna instances generated by the
            :mod:`lsl.common.stations` module instead of a list of stand ID
            numbers.
        """
        
        # Make sure that we have been passed 2047 or fewer stands
        if len(antennas) > 2047:
            raise RuntimeError("UVFITS supports up to 2047 antennas only, given %i" % len(antennas))
            
        # Update the observatory-specific information
        self.siteName = site.name
        
        stands = []
        for ant in antennas:
            stands.append(ant.stand.id)
        stands = numpy.array(stands)
        
        arrayX, arrayY, arrayZ = site.getGeocentricLocation()
        
        xyz = numpy.zeros((len(stands),3))
        i = 0
        for ant in antennas:
            xyz[i,0] = ant.stand.x
            xyz[i,1] = ant.stand.y
            xyz[i,2] = ant.stand.z
            i += 1
            
        # Create the stand mapper
        mapper = {}
        if stands.max() > 2047:
            enableMapper = True
        else:
            enableMapper = False
            
        ants = []
        topo2eci = site.getECITransform()
        for i in xrange(len(stands)):
            eci = numpy.dot(topo2eci, xyz[i,:])
            ants.append( self._Antenna(stands[i], eci[0], eci[1], eci[2], bits=bits) )
            if enableMapper:
                mapper[stands[i]] = i+1
            else:
                mapper[stands[i]] = stands[i]
                
        # If the mapper has been enabled, tell the user about it
        if enableMapper and self.verbose:
            print("UVFITS: stand ID mapping enabled")
            for key, value in mapper.iteritems():
                print("UVFITS:  stand #%i -> mapped #%i" % (key, value))
                
        self.nAnt = len(ants)
        self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'enableMapper': enableMapper, 'inputAnts': antennas} )
        
    def add_data_set(self, obsTime, intTime, baselines, visibilities, pol='XX', source='z'):
        """
        Create a UVData object to store a collection of visibilities.
        
        .. versionchanged:: 0.4.0
            Switched over to passing in Antenna instances generated by the
            :mod:`lsl.common.stations` module instead of a list of stand ID
            as part of the baselines.
        """
        
        if type(pol) == str:
            numericPol = StokesCodes[pol.upper()]
        else:
            numericPol = pol
            
        self.data.append( self._UVData(obsTime, intTime, baselines, visibilities, pol=numericPol, source=source) )
        
    def write(self):
        """
        Fill in the FITS-IDI file will all of the tables in the 
        correct order.
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
        try:
            self.data.sort(cmp=__sortData)
        except TypeError:
            self.data.sort(key=cmp_to_key(__sortData))
            
        self._write_primary_hdu()
        self._write_aipsan_hdu()
        
        # Clear out the data section
        del(self.data[:])
        gc.collect()
        
    def close(self):
        """
        Close out the file.
        """
        
        self.FITS.flush()
        self.FITS.close()
        
    def _add_common_keywords(self, hdr, name, revision):
        """
        Added keywords common to all table headers.
        """
        
        hdr['EXTNAME'] = (name, 'UVFITS table name')
        hdr['EXTVER'] = (1, 'table instance number') 
        hdr['TABREV'] = (revision, 'table format revision number')
        
        date = self.ref_time.split('-')
        name = "ZA%s%s%s" % (date[0][2:], date[1], date[2])
        hdr['OBSCODE'] = (name, 'zenith all-sky image')
        
        hdr['ARRNAM'] = self.siteName
        hdr['RDATE'] = (self.ref_time, 'file data reference date')
        
    def _write_primary_hdu(self):
        """
        Write the primary HDU to file.
        """
        
        self._write_aipsan_hdu(dummy=True)
        hrz = astro.hrz_posn(0, 90)
        (arrPos, ag) = self.read_array_geometry(dummy=True)
        (mapper, inverseMapper) = self.read_array_mapper(dummy=True)
        ids = ag.keys()
        
        obs = ephem.Observer()
        obs.lat = arrPos.lat * numpy.pi/180
        obs.lon = arrPos.lng * numpy.pi/180
        obs.elev = arrPos.elv * numpy.pi/180
        obs.pressure = 0
        
        first = True
        mList = []
        uList = []
        vList = []
        wList = []
        dateList = []
        blineList = []
        for dataSet in self.data:
            # Sort the data by packed baseline
            try:
                order
            except NameError:
                order = dataSet.argsort(mapper=mapper)
                
            # Deal with defininig the values of the new data set
            if dataSet.pol == self.stokes[0]:
                ## Figure out the new date/time for the observation
                utc = astro.taimjd_to_utcjd(dataSet.obsTime)
                date = astro.get_date(utc)
                date.hours = 0
                date.minutes = 0
                date.seconds = 0
                utc0 = date.to_jd()
                
                ## Update the observer so we can figure out where the source is
                obs.date = utc - astro.DJD_OFFSET
                if dataSet.source == 'z':
                    ### Zenith pointings
                    equ = astro.equ_posn( obs.sidereal_time()*180/numpy.pi, obs.lat*180/numpy.pi )
                    
                    ### format 'source' name based on local sidereal time
                    raHms = astro.deg_to_hms(equ.ra)
                    (tsecs, secs) = math.modf(raHms.seconds)
                    name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
                else:
                    ### Real-live sources (ephem.Body instances)
                    name = dataSet.source.name
                    
                ## Compute the uvw coordinates of all baselines
                if dataSet.source == 'z':
                    RA = obs.sidereal_time()
                    HA = 0.0
                    dec = equ.dec
                else:
                    RA = dataSet.source.ra * 180/numpy.pi
                    HA = obs.sidereal_time() - dataSet.source.ra
                    dec = dataSet.source.dec * 180/numpy.pi
                    
                if first is True:
                    sourceRA, sourceDec = RA, dec
                    first = False
                    
                uvwCoords = dataSet.get_uvw(HA, dec, obs)
                
                ## Populate the metadata
                ### Add in the new baselines
                try:
                    blineList.extend( baselineMapped )
                except NameError:
                    baselineMapped = []
                    for o in order:
                        antenna1, antenna2 = dataSet.baselines[o]
                        if mapper is None:
                            stand1, stand2 = antenna1.stand.id, antenna2.stand.id
                        else:
                            stand1, stand2 = mapper[antenna1.stand.id], mapper[antenna2.stand.id]
                        baselineMapped.append( merge_baseline(stand1, stand2) ) 
                    blineList.extend( baselineMapped )
                    
                ### Add in the new u, v, and w coordinates
                uList.extend( uvwCoords[order,0] )
                vList.extend( uvwCoords[order,1] )
                wList.extend( uvwCoords[order,2] )
                
                ### Add in the new date/time
                dateList.extend( [utc for bl in dataSet.baselines] )
                
                ### Zero out the visibility data
                try:
                    matrix *= 0.0
                except NameError:
                    matrix = numpy.zeros((len(order), 1, 1, self.nChan, self.nStokes, 2), dtype=numpy.float32)
                    
            # Save the visibility data in the right order
            matrix[:,0,0,:,self.stokes.index(dataSet.pol),0] = dataSet.visibilities[order,:].real
            matrix[:,0,0,:,self.stokes.index(dataSet.pol),1] = dataSet.visibilities[order,:].imag
            
            # Deal with saving the data once all of the polarizations have been added to 'matrix'
            if dataSet.pol == self.stokes[-1]:
                mList.append( matrix*1.0 )
                
        nBaseline = len(blineList)
        
        # Create the UV Data table and update its header
        uv = astrofits.GroupData(numpy.concatenate(mList), parnames=['UU', 'VV', 'WW', 'BASELINE', 'DATE'], 
                            pardata=[numpy.array(uList, dtype=numpy.float32), numpy.array(vList, dtype=numpy.float32), 
                                    numpy.array(wList, dtype=numpy.float32), numpy.array(blineList), 
                                    numpy.array(dateList)], bitpix=-32)
        primary = astrofits.GroupsHDU(uv)
        
        primary.header['EXTEND'] = (True, 'indicates UVFITS file')
        primary.header['GROUPS'] = (True, 'indicates UVFITS file')
        primary.header['OBJECT'] = 'BINARYTB'
        primary.header['TELESCOP'] = self.siteName
        primary.header['INSTRUME'] = self.siteName
        primary.header['OBSERVER'] = ('ZASKY', 'zenith all-sky image')
        primary.header['ORIGIN'] = 'LSL'
        primary.header['CORRELAT'] = ('LWASWC', 'Correlator used')
        primary.header['FXCORVER'] = ('1', 'Correlator version')
        primary.header['LWATYPE'] = ('UV-ZA', 'LWA FITS file type')
        primary.header['LWAMAJV'] = (UVVersion[0], 'LWA UVFITS file format major version')
        primary.header['LWAMINV'] = (UVVersion[1], 'LWA UVFITS file format minor version')
        primary.header['DATE-OBS'] = (self.ref_time, 'UVFITS file data collection date')
        ts = str(astro.get_date_from_sys())
        primary.header['DATE-MAP'] = (ts.split()[0], 'UVFITS file creation date')
        
        primary.header['CTYPE2'] = ('COMPLEX', 'axis 2 is COMPLEX axis')
        primary.header['CDELT2'] = 1.0
        primary.header['CRPIX2'] = 1.0
        primary.header['CRVAL2'] = 1.0
        
        primary.header['CTYPE3'] = ('STOKES', 'axis 3 is STOKES axis (polarization)')
        if self.stokes[0] < 0:
            primary.header['CDELT3'] = -1.0
        else:
            primary.header['CDELT3'] = 1.0
        primary.header['CRPIX3'] = 1.0
        primary.header['CRVAL3'] = float(self.stokes[0])
        
        primary.header['CTYPE4'] = ('FREQ', 'axis 4 is FREQ axis (frequency)')
        primary.header['CDELT4'] = self.freq[0].chWidth
        primary.header['CRPIX4'] = self.refPix
        primary.header['CRVAL4'] = self.refVal
        
        primary.header['CTYPE5'] = ('RA', 'axis 5 is RA axis (position of phase center)')
        primary.header['CDELT5'] = 0.0
        primary.header['CRPIX5'] = 1.0
        primary.header['CRVAL5'] = sourceRA
        
        primary.header['CTYPE6'] = ('DEC', 'axis 6 is DEC axis (position of phase center)')
        primary.header['CDELT6'] = 0.0
        primary.header['CRPIX6'] = 1.0
        primary.header['CRVAL6'] = sourceDec
        
        primary.header['TELESCOP'] = self.siteName
        primary.header['OBSERVER'] = 'ZASKY'
        primary.header['SORT'] = ('TB', 'data is sorted in [time,baseline] order')
        
        primary.header['VISSCALE'] = (1.0, 'UV data scale factor')
        
        self.FITS.append(primary)
        self.FITS.flush()
        
    def _write_aipsan_hdu(self, dummy=False):
        """
        Define the 'AIPS AN' table .
        """
        
        i = 0
        names = []
        xyz = numpy.zeros((self.nAnt,3), dtype=numpy.float64)
        for ant in self.array[0]['ants']:
            xyz[i,:] = numpy.array([ant.x, ant.y, ant.z])
            names.append(ant.getName())
            i = i + 1
            
        # Antenna name
        c1 = astrofits.Column(name='ANNAME', format='A8', 
                        array=numpy.array([ant.getName() for ant in self.array[0]['ants']]))
        # Station coordinates in meters
        c2 = astrofits.Column(name='STABXYZ', unit='METERS', format='3D', 
                        array=xyz)
        # Station number
        c3 = astrofits.Column(name='NOSTA', format='1J', 
                        array=numpy.array([self.array[0]['mapper'][ant.id] for ant in self.array[0]['ants']]))
        # Mount type (0 == alt-azimuth)
        c4 = astrofits.Column(name='MNTSTA', format='1J', 
                        array=numpy.zeros((self.nAnt,), dtype=numpy.int32))
        # Axis offset in meters
        c5 = astrofits.Column(name='STAXOF', unit='METERS', format='3E', 
                        array=numpy.zeros((self.nAnt,3), dtype=numpy.float32))
        # Feed A polarization label
        c6 = astrofits.Column(name='POLTYA', format='A1', 
                        array=numpy.array([ant.polA['Type'] for ant in self.array[0]['ants']]))
        # Feed A orientation in degrees
        c7 = astrofits.Column(name='POLAA', format='1E', 
                        array=numpy.array([ant.polA['Angle'] for ant in self.array[0]['ants']], dtype=numpy.float32))
        # Feed A polarization parameters
        c8 = astrofits.Column(name='POLCALA', format='2E', 
                        array=numpy.array([ant.polA['Cal'] for ant in self.array[0]['ants']], dtype=numpy.float32))
        # Feed B polarization label
        c9 = astrofits.Column(name='POLTYB', format='A1', 
                        array=numpy.array([ant.polB['Type'] for ant in self.array[0]['ants']]))
        # Feed B orientation in degrees
        c10 = astrofits.Column(name='POLAB', format='1E', 
                        array=numpy.array([ant.polB['Angle'] for ant in self.array[0]['ants']], dtype=numpy.float32))
        # Feed B polarization parameters
        c11 = astrofits.Column(name='POLCALB', format='2E', 
                        array=numpy.array([ant.polB['Cal'] for ant in self.array[0]['ants']], dtype=numpy.float32))
                        
        # Define the collection of columns
        colDefs = astrofits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11])
        
        # Create the table and fill in the header
        ag = astrofits.BinTableHDU.from_columns(colDefs)
        self._add_common_keywords(ag.header, 'AIPS AN', 1)
        
        ag.header['EXTVER'] = (1, 'array ID')
        ag.header['ARRNAM'] = self.siteName
        ag.header['FRAME'] = ('GEOCENTRIC', 'coordinate system')
        ag.header['NUMORB'] = (0, 'number of orbital parameters')
        ag.header['FREQ'] = (self.refVal, 'reference frequency (Hz)')
        ag.header['TIMSYS'] = ('UTC', 'time coordinate system')

        date = self.ref_time2AstroDate()
        utc0 = date.to_jd()
        gst0 = astro.get_apparent_sidereal_time(utc0)
        ag.header['GSTIA0'] = (gst0 * 15, 'GAST (deg) at RDATE 0 hours')
        
        utc1 = utc0 + 1
        gst1 = astro.get_apparent_sidereal_time(utc1)
        if gst1 < gst0:
            gst1 += 24.0
        ds = gst1 - gst0
        deg = ds * 15.0      
        ag.header['DEGPDY'] = (360.0 + deg, 'rotation rate of the earth (deg/day)')
        
        refDate = self.ref_time2AstroDate()
        refMJD = refDate.to_jd() - astro.MJD_OFFSET
        eop = geodesy.getEOP(refMJD)
        if eop is None:
            eop = geodesy.EOP(mjd=refMJD)
            
        ag.header['UT1UTC'] = (eop.utDiff, 'difference UT1 - UTC for reference date')
        ag.header['IATUTC'] = (astro.leap_secs(utc0), 'TAI - UTC for reference date')
        ag.header['POLARX'] = eop.x
        ag.header['POLARY'] = eop.y
        
        ag.header['ARRAYX'] = (self.array[0]['center'][0], 'array ECI X coordinate (m)')
        ag.header['ARRAYY'] = (self.array[0]['center'][1], 'array ECI Y coordinate (m)')
        ag.header['ARRAYZ'] = (self.array[0]['center'][2], 'array ECI Z coordinate (m)')
        
        ag.header['NOSTAMAP'] = (int(self.array[0]['enableMapper']), 'Mapping enabled for stand numbers')
        
        if dummy:
            self.an = ag
            if self.array[0]['enableMapper']:
                self._write_mapper_hdu(dummy=True)
                
        else:
            ag.name = 'AIPS AN'
            self.FITS.append(ag)
            self.FITS.flush()
            
            if self.array[0]['enableMapper']:
                self._write_mapper_hdu()
                
    def _write_mapper_hdu(self, dummy=False):
        """
        Write a fits table that contains information about mapping stations 
        numbers to actual antenna numbers.  This information can be backed out of
        the names, but this makes the extraction more programmatic.
        """
        
        c1 = astrofits.Column(name='ANNAME', format='A8', 
                        array=numpy.array([ant.getName() for ant in self.array[0]['ants']]))
        c2 = astrofits.Column(name='NOSTA', format='1J', 
                        array=numpy.array([self.array[0]['mapper'][ant.id] for ant in self.array[0]['ants']]))
        c3 = astrofits.Column(name='NOACT', format='1J', 
                        array=numpy.array([ant.id for ant in self.array[0]['ants']]))
                        
        colDefs = astrofits.ColDefs([c1, c2, c3])
        
        # Create the ID mapping table and update its header
        nsm = astrofits.BinTableHDU.from_columns(colDefs)
        self._add_common_keywords(nsm.header, 'NOSTA_MAPPER', 1)
        
        if dummy:
            self.am = nsm
        else:
            nsm.name = 'NOSTA_MAPPER'
            self.FITS.append(nsm)
            self.FITS.flush()
            
    def read_array_geometry(self, dummy=False):
        """
        Return a tuple with the array geodetic position and the local 
        positions for all antennas defined in the AIPS AN table.
        """
        
        if dummy:
            try:
                ag = self.an
            except AttributeError:
                raise RuntimeError("Temporary 'AIPS AN' table not found.")
                
        else:
            try:
                ag = self.FITS['AIPS AN']
            except IndexError:
                raise RuntimeError("File does not have an 'AIPS AN' table.")
                
        # Array position
        arrayGeo = astro.rect_posn(ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ'])
        arrayGeo = astro.get_geo_from_rect(arrayGeo)
        
        # Antenna positions
        antennaGeo = {}
        antenna = ag.data
        antennaID = antenna.field('NOSTA')
        antennaPos = iter(antenna.field('STABXYZ'))
        for id in antennaID:
            antennaGeo[id] = next(antennaPos)
            
        # Return
        return (arrayGeo, antennaGeo)
        
    def read_array_mapper(self, dummy=False):
        """
        Return a tuple with the array NOSTA mapper and inverse mapper (both
        dictionaries.  If the stand IDs have not been mapped, return None for
        both.
        """
        
        if dummy:
            try:
                nsm = self.am
            except AttributeError:
                return (None, None)
                
        else:
            try:
                nsm = self.FITS['NOSTA_MAPPER']
            except KeyError:
                return (None, None)
                
        # Build the mapper and inverseMapper
        mapper = {}
        inverseMapper = {}
        nosta = nsm.data.field('NOSTA')
        noact = nsm.data.field('NOACT')
        for idMapped, idActual in zip(nosta, noact):
            mapper[idActual] = idMapped
            inverseMapper[idMapped] = idActual
            
        # Return
        return (mapper, inverseMapper)
