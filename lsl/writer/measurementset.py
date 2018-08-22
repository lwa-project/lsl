# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range

"""
Module for writing correlator output to a CASA measurement set.

.. versionadded:: 1.2.1
"""

import os
import gc
import re
import glob
import math
import ephem
import numpy
import shutil
from datetime import datetime

from lsl import astro
from lsl.misc import geodesy
from lsl.common import constants
from lsl.misc.total_sorting import cmp_to_total

try:
    from collections import OrderedDict
except ImportError:
    from lsl.misc.OrderedDict import OrderedDict

from casacore.tables import table, tableutil


__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['MS', 'StokesCodes', 'NumericStokes', 
           '__version__', '__revision__', '__all__']


StokesCodes = {'I':  1, 'Q':  2, 'U':  3, 'V':  4, 
               'RR': 5, 'RL': 6, 'LR': 7, 'LL': 8,
               'XX': 9, 'XY':10, 'YX':11, 'YY':12}
               

NumericStokes = { 1:'I',   2:'Q',   3:'U',   4:'V', 
                  5:'RR',  6:'RL',  7:'LR',  8:'LL',
                  9:'XX', 10:'XY', 11:'YX', 12:'YY'}

def mergeBaseline(ant1, ant2, shift=16):
    """
    Merge two stand ID numbers into a single baseline using the specified bit 
    shift size.
    """
    
    return (ant1 << shift) | ant2

def splitBaseline(baseline, shift=16):
    """
    Given a baseline, split it into it consistent stand ID numbers.
    """
    
    part = 2**shift - 1
    return (baseline >> shift) & part, baseline & part


class MS(object):
    """
    Class for storing visibility data and writing the data, along with array
    geometry, frequency setup, etc., to a CASA measurement set.
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
            
        def getUVW(self, HA, dec, obs):
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
                uvw[i,:] = numpy.squeeze(temp)
                
            return uvw
                
        def argsort(self, mapper=None, shift=16):
            packed = []
            for a1,a2 in self.baselines:
                if mapper is None:
                    s1, s2 = a1.stand.id, a2.stand.id
                else:
                    s1, s2 = mapper.index(a1.stand.id), mapper.index(a2.stand.id)
                packed.append( mergeBaseline(s1, s2, shift=shift) )
            packed = numpy.array(packed, dtype=numpy.int32)
            
            return numpy.argsort(packed)
            
    def parseRefTime(self, refTime):
        """
        Given a time as either a integer, float, string, or datetime object, 
        convert it to a string in the formation 'YYYY-MM-DDTHH:MM:SS'.
        """
        
        # Valid time string (modulo the 'T')
        timeRE = re.compile(r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?')
        
        if type(refTime) in (int, long, float):
            refDateTime = datetime.utcfromtimestamp(refTime)
            refTime = refDateTime.strftime("%Y-%m-%dT%H:%M:%S")
        elif type(refTime) == datetime:
            refTime = refTime.strftime("%Y-%m-%dT%H:%M:%S")
        elif type(refTime) == str:
            # Make sure that the string times are of the correct format
            if re.match(timeRE, refTime) is None:
                raise RuntimeError("Malformed date/time provided: %s" % refTime)
            else:
                refTime = refTime.replace(' ', 'T', 1)
        else:
            raise RuntimeError("Unknown time format provided.")
            
        return refTime
        
    def refTime2AstroDate(self):
        """
        Convert a reference time string to an :class:`lsl.astro.date` object.
        """
        
        dateStr = self.refTime.replace('T', '-').replace(':', '-').split('-')
        return astro.date(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), int(dateStr[3]), int(dateStr[4]), float(dateStr[5]))
        
    def __init__(self, filename, refTime=0.0, verbose=False, memmap=None, clobber=False):
        """
        Initialize a new FITS IDI object using a filename and a reference time 
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
        self.refTime = self.parseRefTime(refTime)
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
                shutil.rmtree(filename, ignore_errors=False)
            else:
                raise IOError("File '%s' already exists" % filename)
        self.basename = filename
        
    def setStokes(self, polList):
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
        
    def setFrequency(self, freq):
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
        
    def setGeometry(self, site, antennas, bits=8):
        """
        Given a station and an array of stands, set the relevant common observation
        parameters and add entries to the self.array list.
        
        .. versionchanged:: 0.4.0
            Switched over to passing in Antenna instances generated by the
            :mod:`lsl.common.stations` module instead of a list of stand ID
            numbers.
        """
        
        # Update the observatory-specific information
        self.siteName = site.name
        
        stands = []
        for ant in antennas:
            stands.append(ant.stand.id)
        stands = numpy.array(stands)
        
        arrayX, arrayY, arrayZ = site.getGeocentricLocation()
        
        xyz = numpy.zeros((len(stands),3))
        for i,ant in enumerate(antennas):
            xyz[i,0] = ant.stand.x
            xyz[i,1] = ant.stand.y
            xyz[i,2] = ant.stand.z
            
        # Create the stand mapper
        mapper = []
        ants = []
        topo2eci = site.getECITransform()
        for i in xrange(len(stands)):
            eci = numpy.dot(topo2eci, xyz[i,:])
            ants.append( self._Antenna(stands[i], eci[0], eci[1], eci[2], bits=bits) )
            mapper.append( stands[i] )
            
        self.nAnt = len(ants)
        self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'inputAnts': antennas} )
        
    def addDataSet(self, obsTime, intTime, baselines, visibilities, pol='XX', source='z'):
        """
        Create a UVData object to store a collection of visibilities.
        
        .. versionchanged:: 0.4.0
            Switched over to passing in Antenna instances generated by the
            :mod:`lsl.common.stations` module instead of a list of stand ID
            as part of the baselines.
            
        .. versionchanged:: 1.1.0
            Added a new 'source' keyword to set the phase center for the data.
            This can either by 'z' for zenith or a ephem.Body instances for a
            point on the sky.
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
        
        # Sort the data set
        self.data.sort()
        
        # Write the tables
        self._writeMain()
        self._writeAntenna()
        self._writePolarization()
        self._writeObservation()
        self._writeSpectralWindow()
        self._writeMisc()
        
        # Fixup the keywords
        tb = table("%s" % self.basename, readonly=False, ack=False)
        tb.putkeyword('MS_VERSION', numpy.float32(2.0))
        for filename in sorted(glob.glob('%s/*' % self.basename)):
            if os.path.isdir(filename):
                tname = os.path.basename(filename)
                stb = table("%s/%s" % (self.basename, tname), ack=False)
                tb.putkeyword(tname, stb)
                stb.close()
        tb.close()
            
        # Clear out the data section
        del(self.data[:])
        gc.collect()
        
    def close(self):
        """
        Close out the file.
        """
        
        pass
        
    def _writeAntenna(self):
        """
        Write the antenna table.
        """
        
        col1 = tableutil.makearrcoldesc('OFFSET', 0.0, 1, 
                                        keywords={'QuantumUnits':['m','m','m'], 
                                                  'MEASINFO':{'type':'position', 'Ref':'ITRF'}
                                        })
        col2 = tableutil.makearrcoldesc('POSITION', 0.0, 1,
                                        keywords={'QuantumUnits':['m','m','m'], 
                                                  'MEASINFO':{'type':'position', 'Ref':'ITRF'}
                                                  })
        col3 = tableutil.makescacoldesc('TYPE', "ground-based")
        col4 = tableutil.makescacoldesc('DISH_DIAMETER', 2.0, 
                                        keywords={'QuantumUnits':['m',]})
        col5 = tableutil.makescacoldesc('FLAG_ROW', False)
        col6 = tableutil.makescacoldesc('MOUNT', "alt-az")
        col7 = tableutil.makescacoldesc('NAME', "none")
        col8 = tableutil.makescacoldesc('STATION', self.siteName)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8])
        tb = table("%s/ANTENNA" % self.basename, desc, nrow=self.nAnt, ack=False)
        
        for i,ant in enumerate(self.array[0]['ants']):
            tb.putcell('OFFSET', i, [0.0, 0.0, 0.0])
            tb.putcell('POSITION', i, [ant.x + self.array[0]['center'][0],
                                       ant.y + self.array[0]['center'][1], 
                                       ant.z + self.array[0]['center'][2]])
            tb.putcell('TYPE', i, 'GROUND-BASED')
            tb.putcell('DISH_DIAMETER', i, 2.0)
            tb.putcell('FLAG_ROW', i, False)
            tb.putcell('MOUNT', i, 'ALT-AZ')
            tb.putcell('NAME', i, ant.getName())
            tb.putcell('STATION', i, self.siteName)
            
        tb.done()
        
    def _writePolarization(self):
        """
        Write the polarization table.
        """
        
        # Polarization
        
        stks = numpy.array(self.stokes)
        prds = numpy.zeros((2,self.nStokes), dtype=numpy.int32)
        for i,stk in enumerate(self.stokes):
            stks[i] = stk
            if stk > 4:
                prds[0,i] = ((stk-1) % 4) / 2
                prds[1,i] = ((stk-1) % 4) % 2
            else:
                prds[0,i] = 1
                prds[1,i] = 1
                
        col1 = tableutil.makearrcoldesc('CORR_TYPE', 0, 1)
        col2 = tableutil.makearrcoldesc('CORR_PRODUCT', 0, 2)
        col3 = tableutil.makescacoldesc('FLAG_ROW', False)
        col4 = tableutil.makescacoldesc('NUM_CORR', self.nStokes)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4])
        tb = table("%s/POLARIZATION" % self.basename, desc, nrow=1, ack=False)
        
        tb.putcell('CORR_TYPE', 0, self.stokes)
        tb.putcell('CORR_PRODUCT', 0, prds.T)
        tb.putcell('FLAG_ROW', 0, False)
        tb.putcell('NUM_CORR', 0, self.nStokes)
        
        tb.done()
        
        # Feed
        
        col1  = tableutil.makearrcoldesc('POSITION', 0.0, 1, 
                                         keywords={'QuantumUnits':['m','m','m'], 
                                                   'MEASINFO':{'type':'position', 'Ref':'ITRF'}
                                                   })
        col2  = tableutil.makearrcoldesc('BEAM_OFFSET', 0.0, 2, 
                                         keywords={'QuantumUnits':['rad','rad'], 
                                                   'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                   })
        col3  = tableutil.makearrcoldesc('POLARIZATION_TYPE', 'X', 1)
        col4  = tableutil.makearrcoldesc('POL_RESPONSE', 1j, 2, valuetype='complex')
        col5  = tableutil.makearrcoldesc('RECEPTOR_ANGLE', 0.0, 1,  
                                         keywords={'QuantumUnits':['rad',]})
        col6  = tableutil.makescacoldesc('ANTENNA_ID', 0)
        col7  = tableutil.makescacoldesc('BEAM_ID', -1)
        col8  = tableutil.makescacoldesc('FEED_ID', 0)
        col9  = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         keywords={'QuantumUnits':['s',]})
        col10 = tableutil.makescacoldesc('NUM_RECEPTORS', 2)
        col11 = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', -1)
        col12 = tableutil.makescacoldesc('TIME', 0.0, 
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, 
                                        col9, col10, col11, col12])
        tb = table("%s/FEED" % self.basename, desc, nrow=self.nAnt, ack=False)
        
        ptype = [None, None]
        presp = numpy.zeros((2,2), dtype=numpy.complex64)
        if self.stokes[0] > 8:
            ptype = ['X', 'Y']
            presp[0,0] = 1.0
            presp[0,1] = 0.0
            presp[1,0] = 0.0
            presp[1,1] = 1.0
        elif self.stokes[0] > 4:
            ptype = ['R', 'L']
            presp[0,0] = 1.0 - 1.0j
            presp[0,1] = 0.0
            presp[1,0] = 0.0
            presp[1,1] = 1.0 + 1.0j
        else:
            ptype = ['X', 'Y']
            presp[0,0] = 1.0
            presp[0,1] = 0.0
            presp[1,0] = 0.0
            presp[1,1] = 1.0
            
        for i,ant in enumerate(self.array[0]['ants']):
            tb.putcell('POSITION', i, numpy.zeros(3))
            tb.putcell('BEAM_OFFSET', i, numpy.zeros((2,2)))
            tb.putcell('POLARIZATION_TYPE', i, ptype)
            tb.putcell('POL_RESPONSE', i, presp)
            tb.putcell('RECEPTOR_ANGLE', i, numpy.zeros(2))
            tb.putcell('ANTENNA_ID', i, i)
            tb.putcell('BEAM_ID', i, -1)
            tb.putcell('FEED_ID', i, 0)
            tb.putcell('INTERVAL', i, 0.0)
            tb.putcell('NUM_RECEPTORS', i, 2)
            tb.putcell('SPECTRAL_WINDOW_ID', i, -1)
            tb.putcell('TIME', i, 0.0)
        
        tb.done()
        
    def _writeObservation(self):
        """
        Write the observation table.
        """
        
        # Observation
        
        col1 = tableutil.makearrcoldesc('TIME_RANGE', 0.0, 1, 
                                        keywords={'QuantumUnits':['s',], 
                                                  'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                  })
        col2 = tableutil.makearrcoldesc('LOG', 'none', 1)
        col3 = tableutil.makearrcoldesc('SCHEDULE', 'none', 1)
        col4 = tableutil.makescacoldesc('FLAG_ROW', False)
        col5 = tableutil.makescacoldesc('OBSERVER', 'ZASKY')
        col6 = tableutil.makescacoldesc('PROJECT', 'ZASKY')
        col7 = tableutil.makescacoldesc('RELEASE_DATE', 0.0, 
                                        keywords={'QuantumUnits':['s',], 
                                                  'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                  })
        col8 = tableutil.makescacoldesc('SCHEDULE_TYPE', 'none')
        col9 = tableutil.makescacoldesc('TELESCOPE_NAME', self.siteName)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tb = table("%s/OBSERVATION" % self.basename, desc, nrow=1, ack=False)
        
        tStart = astro.taimjd_to_utcjd(self.data[ 0].obsTime) - astro.MJD_OFFSET
        tStop  = astro.taimjd_to_utcjd(self.data[-1].obsTime) - astro.MJD_OFFSET
        
        tb.putcell('TIME_RANGE', 0, [tStart*86400, tStop*86400])
        tb.putcell('LOG', 0, 'Not provided')
        tb.putcell('SCHEDULE', 0, 'Not provided')
        tb.putcell('FLAG_ROW', 0, False)
        tb.putcell('OBSERVER', 0, 'ZASKY')
        tb.putcell('PROJECT', 0, 'ZASKY')
        tb.putcell('RELEASE_DATE', 0, tStop*86400)
        tb.putcell('SCHEDULE_TYPE', 0, 'None')
        tb.putcell('TELESCOPE_NAME', 0, self.siteName)
        
        tb.done()
        
        # Source
        
        arrayGeo = astro.rect_posn(*self.array[0]['center'])
        arrayGeo = astro.get_geo_from_rect(arrayGeo)
        
        obs = ephem.Observer()
        obs.lat = arrayGeo.lat * numpy.pi/180
        obs.lon = arrayGeo.lng * numpy.pi/180
        obs.elev = arrayGeo.elv * numpy.pi/180
        obs.pressure = 0
        
        nameList = []
        posList = []
        sourceID = 0
        lastSourceName = None
        for dataSet in self.data:
            if dataSet.pol == self.stokes[0]:
                utc = astro.taimjd_to_utcjd(dataSet.obsTime)
                date = astro.get_date(utc)
                date.hours = 0
                date.minutes = 0
                date.seconds = 0
                utc0 = date.to_jd()
                
                obs.date = utc - astro.DJD_OFFSET
                
                try:
                    currSourceName = dataSet.source.name
                except AttributeError:
                    currSourceName = dataSet.source
                
                if currSourceName != lastSourceName:
                    sourceID += 1
                    
                    if dataSet.source == 'z':
                        ## Zenith pointings
                        equ = astro.equ_posn( obs.sidereal_time()*180/numpy.pi, obs.lat*180/numpy.pi )
                        
                        # format 'source' name based on local sidereal time
                        raHms = astro.deg_to_hms(equ.ra)
                        (tsecs, secs) = math.modf(raHms.seconds)
                        
                        name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
                        equPo = astro.get_equ_prec2(equ, utc, astro.J2000_UTC_JD)
                        
                    else:
                        ## Real-live sources (ephem.Body instances)
                        name = dataSet.source.name
                        equ = astro.equ_posn(dataSet.source.ra*180/numpy.pi, dataSet.source.dec*180/numpy.pi)
                        equPo = astro.equ_posn(dataSet.source.a_ra*180/numpy.pi, dataSet.source.a_dec*180/numpy.pi)
                        
                    # J2000 zenith equatorial coordinates
                    posList.append( [equPo.ra*numpy.pi/180, equPo.dec*numpy.pi/180] )
                    
                    # name
                    nameList.append(name)
                    
                    # Update
                    lastSourceName = name
                    
        nSource = len(nameList)
        
        # Save these for later since we might need them
        self._sourceTable = nameList
        
        col1  = tableutil.makearrcoldesc('DIRECTION', 0.0, 1, 
                                         keywords={'QuantumUnits':['rad','rad'], 
                                                   'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                   })
        col2  = tableutil.makearrcoldesc('PROPER_MOTION', 0.0, 1, 
                                         keywords={'QuantumUnits':['rad/s',]})
        col3  = tableutil.makescacoldesc('CALIBRATION_GROUP', 0)
        col4  = tableutil.makescacoldesc('CODE', "none")
        col5  = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         keywords={'QuantumUnits':['s',]})
        col6  = tableutil.makescacoldesc('NAME', "none")
        col7  = tableutil.makescacoldesc('NUM_LINES', 0)
        col8  = tableutil.makescacoldesc('SOURCE_ID', 0)
        col9  = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', -1)
        col10 = tableutil.makescacoldesc('TIME', 0.0,
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        #col11 = tableutil.makearrcoldesc('TRANSITION', 'none', 1)
        #col12 = tableutil.makearrcoldesc('REST_FREQUENCY', 0.0, 1, 
                                         #keywords={'QuantumUnits':['Hz',], 
                                                   #'MEASINFO':{'type':'frequency', 'Ref':'LSRK'}
                                                   #})
        #col13 = tableutil.makearrcoldesc('SYSVEL', 0.0, 1, 
                                         #keywords={'QuantumUnits':['m/s',], 
                                                   #'MEASINFO':{'type':'radialvelocity', 'Ref':'LSRK'}
                                                   #})
                                        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9, 
                                      col10])#, col11, col12, col13])
        tb = table("%s/SOURCE" % self.basename, desc, nrow=nSource, ack=False)
        
        for i in xrange(nSource):
            tb.putcell('DIRECTION', i, posList[i])
            tb.putcell('PROPER_MOTION', i, [0.0, 0.0])
            tb.putcell('CALIBRATION_GROUP', i, 0)
            tb.putcell('CODE', i, 'none')
            tb.putcell('INTERVAL', i, 0.0)
            tb.putcell('NAME', i, nameList[i])
            tb.putcell('NUM_LINES', i, 0)
            tb.putcell('SOURCE_ID', i, i)
            tb.putcell('SPECTRAL_WINDOW_ID', i, -1)
            tb.putcell('TIME', i, (tStart+tStop)/2*86400)
            #tb.putcell('TRANSITION', i, [])
            #tb.putcell('REST_FREQUENCY', i, [])
            #tb.putcell('SYSVEL', i, [])
            
        tb.close()
        
        # Field
        
        col1 = tableutil.makearrcoldesc('DELAY_DIR', 0.0, 2, 
                                        keywords={'QuantumUnits':['rad','rad'], 
                                                  'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                  })
        col2 = tableutil.makearrcoldesc('PHASE_DIR', 0.0, 2, 
                                        keywords={'QuantumUnits':['rad','rad'], 
                                                  'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                  })
        col3 = tableutil.makearrcoldesc('REFERENCE_DIR', 0.0, 2, 
                                        keywords={'QuantumUnits':['rad','rad'], 
                                                  'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                  })
        col4 = tableutil.makescacoldesc('CODE', "none")
        col5 = tableutil.makescacoldesc('FLAG_ROW', False)
        col6 = tableutil.makescacoldesc('NAME', "none")
        col7 = tableutil.makescacoldesc('NUM_POLY', 0)
        col8 = tableutil.makescacoldesc('SOURCE_ID', 0)
        col9 = tableutil.makescacoldesc('TIME', (tStart+tStop)/2, 
                                        keywords={'QuantumUnits':['s',],
                                                  'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                  })
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tb = table("%s/FIELD" % self.basename, desc, nrow=nSource, ack=False)
        
        for i in xrange(nSource):
            tb.putcell('DELAY_DIR', i, numpy.array([posList[i],]))
            tb.putcell('PHASE_DIR', i, numpy.array([posList[i],]))
            tb.putcell('REFERENCE_DIR', i, numpy.array([posList[i],]))
            tb.putcell('CODE', i, 'None')
            tb.putcell('FLAG_ROW', i, False)
            tb.putcell('NAME', i, nameList[i])
            tb.putcell('NUM_POLY', i, 0)
            tb.putcell('SOURCE_ID', i, i)
            tb.putcell('TIME', i, (tStart+tStop)/2*86400)
            
        tb.close()
        
    def _writeSpectralWindow(self):
        """
        Write the spectral window table.
        """
        
        # Spectral Window
        
        col1  = tableutil.makescacoldesc('MEAS_FREQ_REF', 0)
        col2  = tableutil.makearrcoldesc('CHAN_FREQ', 0.0, 1, 
                                         keywords={'QuantumUnits':['Hz',], 
                                                   'MEASINFO':{'type':'frequency', 
                                                               'VarRefCol':'MEAS_FREQ_REF', 
                                                               'TabRefTypes':['REST','LSRK','LSRD','BARY','GEO','TOPO','GALACTO','LGROUP','CMB','Undefined'],
                                                               'TabRefCodes':[0,1,2,3,4,5,6,7,8,64]}
                                                   })
        col3  = tableutil.makescacoldesc('REF_FREQUENCY', self.refVal, 
                                         keywords={'QuantumUnits':['Hz',], 
                                                   'MEASINFO':{'type':'frequency', 
                                                               'VarRefCol':'MEAS_FREQ_REF', 
                                                               'TabRefTypes':['REST','LSRK','LSRD','BARY','GEO','TOPO','GALACTO','LGROUP','CMB','Undefined'],
                                                               'TabRefCodes':[0,1,2,3,4,5,6,7,8,64]}
                                                   })
        col4  = tableutil.makearrcoldesc('CHAN_WIDTH', 0.0, 1, 
                                         keywords={'QuantumUnits':['Hz',]})
        col5  = tableutil.makearrcoldesc('EFFECTIVE_BW', 0.0, 1, 
                                         keywords={'QuantumUnits':['Hz',]})
        col6  = tableutil.makearrcoldesc('RESOLUTION', 0.0, 1, 
                                         keywords={'QuantumUnits':['Hz',]})
        col7  = tableutil.makescacoldesc('FLAG_ROW', False)
        col8  = tableutil.makescacoldesc('FREQ_GROUP', 1)
        col9  = tableutil.makescacoldesc('FREQ_GROUP_NAME', "group1")
        col10 = tableutil.makescacoldesc('IF_CONV_CHAIN', 0)
        col11 = tableutil.makescacoldesc('NAME', "%i channels" % self.nChan)
        col12 = tableutil.makescacoldesc('NET_SIDEBAND', 0)
        col13 = tableutil.makescacoldesc('NUM_CHAN', 0)
        col14 = tableutil.makescacoldesc('TOTAL_BANDWIDTH', 0.0, 
                                         keywords={'QuantumUnits':['Hz',]})
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9, 
                                        col10, col11, col12, col13, col14])
        tb = table("%s/SPECTRAL_WINDOW" % self.basename, desc, nrow=1, ack=False)
        
        tb.putcell('MEAS_FREQ_REF', 0, 0)
        tb.putcell('CHAN_FREQ', 0, self.refVal + numpy.arange(self.nChan)*self.channelWidth)
        tb.putcell('REF_FREQUENCY', 0, self.refVal)
        tb.putcell('CHAN_WIDTH', 0, [self.channelWidth for i in xrange(self.nChan)])
        tb.putcell('EFFECTIVE_BW', 0, [self.channelWidth for i in xrange(self.nChan)])
        tb.putcell('RESOLUTION', 0, [self.channelWidth for i in xrange(self.nChan)])
        tb.putcell('FLAG_ROW', 0, False)
        tb.putcell('FREQ_GROUP', 0, 1)
        tb.putcell('FREQ_GROUP_NAME', 0, 'group1')
        tb.putcell('IF_CONV_CHAIN', 0, 0)
        tb.putcell('NAME', 0, "%i channels" % self.nChan)
        tb.putcell('NET_SIDEBAND', 0, 0)
        tb.putcell('NUM_CHAN', 0, self.nChan)
        tb.putcell('TOTAL_BANDWIDTH', 0, self.channelWidth*self.nChan)
        
        tb.done
        
    def _writeMain(self):
        """
        Write the main table.
        """
        
        # Main
        
        arrayGeo = astro.rect_posn(*self.array[0]['center'])
        arrayGeo = astro.get_geo_from_rect(arrayGeo)
        
        obs = ephem.Observer()
        obs.lat = arrayGeo.lat * numpy.pi/180
        obs.lon = arrayGeo.lng * numpy.pi/180
        obs.elev = arrayGeo.elv * numpy.pi/180
        obs.pressure = 0
        
        mapper = self.array[0]['mapper']
        
        col1  = tableutil.makearrcoldesc('UVW', 0.0, 1, 
                                         keywords={'QuantumUnits':['m','m','m'], 
                                                   'MEASINFO':{'type':'uvw', 'Ref':'ITRF'}
                                                   })
        col2  = tableutil.makearrcoldesc('FLAG', False, 2)
        col3  = tableutil.makearrcoldesc('FLAG_CATEGORY', False, 3,  
                                         keywords={'CATEGORY':['',]})
        col4  = tableutil.makearrcoldesc('WEIGHT', 1.0, 1, valuetype='float')
        col5  = tableutil.makearrcoldesc('SIGMA', 9999., 1, valuetype='float')
        col6  = tableutil.makescacoldesc('ANTENNA1', 0)
        col7  = tableutil.makescacoldesc('ANTENNA2', 0)
        col8  = tableutil.makescacoldesc('ARRAY_ID', 0)
        col9  = tableutil.makescacoldesc('DATA_DESC_ID', 0)
        col10 = tableutil.makescacoldesc('EXPOSURE', 0.0, 
                                         keywords={'QuantumUnits':['s',]})
        col11 = tableutil.makescacoldesc('FEED1', 0)
        col12 = tableutil.makescacoldesc('FEED2', 0)
        col13 = tableutil.makescacoldesc('FIELD_ID', 0)
        col14 = tableutil.makescacoldesc('FLAG_ROW', False)
        col15 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         keywords={'QuantumUnits':['s',]})
        col16 = tableutil.makescacoldesc('OBSERVATION_ID', 0)
        col17 = tableutil.makescacoldesc('PROCESSOR_ID', -1)
        col18 = tableutil.makescacoldesc('SCAN_NUMBER', 1)
        col19 = tableutil.makescacoldesc('STATE_ID', -1)
        col20 = tableutil.makescacoldesc('TIME', 0.0, 
                                         keywords={'QuantumUnits':['s',],
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        col21 = tableutil.makescacoldesc('TIME_CENTROID', 0.0, 
                                         keywords={'QuantumUnits':['s',],
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        col22 = tableutil.makearrcoldesc("DATA", 0j, 2, valuetype='complex')
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9, 
                                        col10, col11, col12, col13, col14, col15, col16, 
                                        col17, col18, col19, col20, col21, col22])
        tb = table("%s" % self.basename, desc, nrow=0, ack=False)
        
        i = 0
        s = 1
        _sourceTable = []
        for dataSet in self.data:
            # Sort the data by packed baseline
            try:
                order
            except NameError:
                order = dataSet.argsort(mapper=mapper, shift=16)
                
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
                    
                ## Update the source ID
                try:
                    sourceID = _sourceTable.index(name)
                except ValueError:
                    _sourceTable.append(name)
                    sourceID = _sourceTable.index(name)
                    
                ## Compute the uvw coordinates of all baselines
                if dataSet.source == 'z':
                    HA = 0.0
                    dec = equ.dec
                else:
                    HA = (obs.sidereal_time() - dataSet.source.ra) * 12/numpy.pi
                    dec = dataSet.source.dec * 180/numpy.pi
                uvwCoords = dataSet.getUVW(HA, dec, obs)
                
                ## Populate the metadata
                ### Add in the baselines
                try:
                    ant1List
                    ant2List
                except NameError:
                    a1List, a2List = [], []
                    for o in order:
                        antenna1, antenna2 = dataSet.baselines[o]
                        a1 = mapper.index(antenna1.stand.id)
                        a2 = mapper.index(antenna2.stand.id)
                        a1List.append( a1 )
                        a2List.append( a2 )
                    ant1List = a1List
                    ant2List = a2List
                    
                ### Add in the new u, v, and w coordinates
                uvwList = uvwCoords[order,:]
                
                ### Add in the new date/time and integration time
                timeList = [utc - astro.MJD_OFFSET for bl in dataSet.baselines]
                intTimeList = [dataSet.intTime for bl in dataSet.baselines]
                
                ### Add in the new new source ID and name
                sourceList = [sourceID for bl in dataSet.baselines]
                 
                ### Zero out the visibility data
                try:
                    matrix *= 0.0
                except NameError:
                    matrix = numpy.zeros((len(order), self.nStokes, self.nChan), dtype=numpy.complex64)
                    
            # Save the visibility data in the right order
            matrix[:,self.stokes.index(dataSet.pol),:] = dataSet.visibilities[order,:]
            
            # Deal with saving the data once all of the polarizations have been added to 'matrix'
            if dataSet.pol == self.stokes[-1]:
                mList = matrix
                
                tb.addrows(uvwList.shape[0])
                for j in xrange(uvwList.shape[0]):
                    tb.putcell('UVW', i, uvwList[j,:])
                    tb.putcell('FLAG', i, numpy.zeros((self.nStokes,self.nChan), dtype=numpy.bool).T)
                    tb.putcell('FLAG_CATEGORY', i, numpy.zeros((self.nStokes,self.nChan,1), dtype=numpy.bool).T)
                    tb.putcell('WEIGHT', i, numpy.ones(self.nStokes))
                    tb.putcell('SIGMA', i, numpy.ones(self.nStokes)*9999)
                    tb.putcell('ANTENNA1', i, ant1List[j])
                    tb.putcell('ANTENNA2', i, ant2List[j])
                    tb.putcell('ARRAY_ID', i, 0)
                    tb.putcell('DATA_DESC_ID', i, 0)
                    tb.putcell('EXPOSURE', i, intTimeList[j])
                    tb.putcell('FEED1', i, 0)
                    tb.putcell('FEED2', i, 0)
                    tb.putcell('FIELD_ID', i, sourceList[j])
                    tb.putcell('FLAG_ROW', i, False)
                    tb.putcell('INTERVAL', i, intTimeList[j])
                    tb.putcell('OBSERVATION_ID', i, 0)
                    tb.putcell('PROCESSOR_ID', i, -1)
                    tb.putcell('SCAN_NUMBER', i, s)
                    tb.putcell('STATE_ID', i, -1)
                    tb.putcell('TIME', i, (timeList[j] + intTimeList[j]/2)*86400)
                    tb.putcell('TIME_CENTROID', i, (timeList[j] + intTimeList[j]/2)*86400)
                    tb.putcell('DATA', i, mList[j,:,:].T)
                    i += 1
                s += 1
            
        tb.close()
        
        # Data description
        
        col1 = tableutil.makescacoldesc('FLAG_ROW', False)
        col2 = tableutil.makescacoldesc('POLARIZATION_ID', 0)
        col3 = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', 0)
        
        desc = tableutil.maketabdesc([col1, col2, col3])
        tb = table("%s/DATA_DESCRIPTION" % self.basename, desc, nrow=1, ack=False)
        
        tb.putcell('FLAG_ROW', 0, False)
        tb.putcell('POLARIZATION_ID', 0, 0)
        tb.putcell('SPECTRAL_WINDOW_ID', 0, 0)
        
        tb.close()
        
    def _writeMisc(self):
        """
        Write the other tables that are part of the measurement set but 
        don't contain anything by default.
        """
        
        # Flag command
        
        col1 = tableutil.makescacoldesc('TIME', 0.0, 
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        col2 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         keywords={'QuantumUnits':['s',]})
        col3 = tableutil.makescacoldesc('TYPE', 'flag')
        col4 = tableutil.makescacoldesc('REASON', 'reason')
        col5 = tableutil.makescacoldesc('LEVEL', 0)
        col6 = tableutil.makescacoldesc('SEVERITY', 0)
        col7 = tableutil.makescacoldesc('APPLIED', False)
        col8 = tableutil.makescacoldesc('COMMAND', 'command')
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8])
        tb = table("%s/FLAG_CMD" % self.basename, desc, nrow=0, ack=False)
        
        tb.close()
        
        # History
        
        col1 = tableutil.makescacoldesc('TIME', 0.0, 
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        col2 = tableutil.makescacoldesc('OBSERVATION_ID', 0)
        col3 = tableutil.makescacoldesc('MESSAGE', 'message')
        col4 = tableutil.makescacoldesc('PRIORITY', 'NORMAL')
        col5 = tableutil.makescacoldesc('ORIGIN', 'origin')
        col6 = tableutil.makescacoldesc('OBJECT_ID', 0)
        col7 = tableutil.makescacoldesc('APPLICATION', 'application')
        col8 = tableutil.makearrcoldesc('CLI_COMMAND', 'command', 1)
        col9 = tableutil.makearrcoldesc('APP_PARAMS', 'params', 1)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tb = table("%s/HISTORY" % self.basename, desc, nrow=0, ack=False)
        
        tb.close()
        
        # POINTING
        
        col1 = tableutil.makescacoldesc('ANTENNA_ID', 0)
        col2 = tableutil.makescacoldesc('TIME', 0.0, 
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        col3 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         keywords={'QuantumUnits':['s',]})
        col4 = tableutil.makescacoldesc('NAME', 'name')
        col5 = tableutil.makescacoldesc('NUM_POLY', 0)
        col6 = tableutil.makescacoldesc('TIME_ORIGIN', 0.0, 
                                         keywords={'QuantumUnits':['s',], 
                                                   'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   })
        col7 = tableutil.makearrcoldesc('DIRECTION', 0.0, 2, 
                                         keywords={'QuantumUnits':['rad','rad'], 
                                                   'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                   })
        col8 = tableutil.makearrcoldesc('TARGET', 0.0, 2, 
                                         keywords={'QuantumUnits':['rad','rad'], 
                                                   'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                   })
        col9 = tableutil.makescacoldesc('TRACKING', True)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tb = table("%s/POINTING" % self.basename, desc, nrow=0, ack=False)
        
        tb.close()
        
        # Processor
        
        col1 = tableutil.makescacoldesc('TYPE', 'type')
        col2 = tableutil.makescacoldesc('SUB_TYPE', 'subtype')
        col3 = tableutil.makescacoldesc('TYPE_ID', 0)
        col4 = tableutil.makescacoldesc('MODE_ID', 0)
        col5 = tableutil.makescacoldesc('FLAG_ROW', False)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5])
        tb = table("%s/PROCESSOR" % self.basename, desc, nrow=0, ack=False)
        
        tb.close()
        
        # State
        
        col1 = tableutil.makescacoldesc('SIG', True)
        col2 = tableutil.makescacoldesc('REF', True)
        col3 = tableutil.makescacoldesc('CAL', 0.0, 
                                        keywords={'QuantumUnits':['K',]})
        col4 = tableutil.makescacoldesc('LOAD', 0.0, 
                                        keywords={'QuantumUnits':['K',]})
        col5 = tableutil.makescacoldesc('SUB_SCAN', 0)
        col6 = tableutil.makescacoldesc('OBS_MODE', 'mode')
        col7 = tableutil.makescacoldesc('FLAG_ROW', False)
        
        desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7])
        tb = table("%s/STATE" % self.basename, desc, nrow=0, ack=False)
        
        tb.close()
        
        ## Syscal
        
        #col1 = tableutil.makescacoldesc('ANTENNA_ID', 0)
        #col2 = tableutil.makescacoldesc('FEED_ID', 0)
        #col3 = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', 0)
        #col4 = tableutil.makescacoldesc('TIME', 0.0, 
                                         #keywords={'QuantumUnits':['s',], 
                                                   #'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                   #})
        #col5 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                         #keywords={'QuantumUnits':['s',]})
        
        #desc = tableutil.maketabdesc([col1, col2, col3, col4, col5])
        #tb = table("%s/SYSCAL" % self.basename, desc, nrow=0, ack=False)
        
        #tb.close()
