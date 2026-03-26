"""
Module for writing correlator output to a CASA measurement set.

.. versionadded:: 1.2.1
"""

import os
import gc
import glob
import math
import numpy as np
import shutil
import warnings
from datetime import datetime

from astropy import units as astrounits
from astropy.time import Time as AstroTime
from astropy.coordinates import EarthLocation, AltAz, ITRS, FK5

from lsl import astro
from lsl.reader.base import FrameTimestamp
from lsl.writer.fitsidi import WriterBase
from lsl.common.color import colorfy

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['Ms', 'STOKES_CODES', 'NUMERIC_STOKES']


STOKES_CODES = {'I': 1,  'Q': 2,  'U': 3,  'V': 4, 
               'RR': 5, 'RL': 6, 'LR': 7, 'LL': 8,
               'XX': 9, 'XY':10, 'YX':11, 'YY':12}
               

NUMERIC_STOKES = { 1:'I',   2:'Q',   3:'U',   4:'V', 
                   5:'RR',  6:'RL',  7:'LR',  8:'LL',
                   9:'XX', 10:'XY', 11:'YX', 12:'YY'}


def merge_baseline(ant1, ant2, shift=16):
    """
    Merge two stand ID numbers into a single baseline using the specified bit 
    shift size.
    """
    
    return (ant1 << shift) | ant2


def split_baseline(baseline, shift=16):
    """
    Given a baseline, split it into it consistent stand ID numbers.
    """
    
    part = 2**shift - 1
    return (baseline >> shift) & part, baseline & part


try:
    from casacore.tables import table, tableutil
    
    class Ms(WriterBase):
        """
        Class for storing visibility data and writing the data, along with array
        geometry, frequency setup, etc., to a CASA measurement set.
        """
        
        _STOKES_CODES = STOKES_CODES
        
        class _MS_UVData(WriterBase._UVData):
            """
            Represents one MS UV visibility data set for a given observation time.
            """
        
            def get_uvw(self, HA, dec, el):
                Nbase = len(self.baselines)
                uvw = np.zeros((Nbase,3), dtype=np.float32)
                
                # Phase center coordinates
                # Convert numbers to radians and, for HA, hours to degrees
                HA2 = HA * 15.0 * np.pi/180
                dec2 = dec * np.pi/180
                lat2 = el.lat.rad
                
                # Coordinate transformation matrices
                trans1 = np.array([[0, -np.sin(lat2), np.cos(lat2)],
                                   [1,  0,            0           ],
                                   [0,  np.cos(lat2), np.sin(lat2)]])
                trans2 = np.array([[ np.sin(HA2),               np.cos(HA2),              0           ],
                                   [-np.sin(dec2)*np.cos(HA2),  np.sin(dec2)*np.sin(HA2), np.cos(dec2)],
                                   [ np.cos(dec2)*np.cos(HA2), -np.cos(dec2)*np.sin(HA2), np.sin(dec2)]])
                        
                for i,(a1,a2) in enumerate(self.baselines):
                    # Go from a east, north, up coordinate system to a celestial equator, 
                    # east, north celestial pole system
                    xyzPrime = a1.stand - a2.stand
                    xyz = np.dot(trans1, np.array([[xyzPrime[0]],[xyzPrime[1]],[xyzPrime[2]]]))
                    
                    # Go from CE, east, NCP to u, v, w
                    temp = np.dot(trans2, xyz)
                    uvw[i,:] = np.squeeze(temp)
                    
                return uvw
                    
            def argsort(self, mapper=None, shift=16):
                packed = []
                for a1,a2 in self.baselines:
                    if mapper is None:
                        s1, s2 = a1.stand.id, a2.stand.id
                    else:
                        s1, s2 = mapper.index(a1.stand.id), mapper.index(a2.stand.id)
                    packed.append( merge_baseline(s1, s2, shift=shift) )
                packed = np.array(packed, dtype=np.int32)
                
                return np.argsort(packed)
                
        def __init__(self, filename, ref_time=0.0, verbose=False, memmap=None, overwrite=False):
            """
            Initialize a new Measurment Set object using a filename and a reference time
            given in seconds since the UNIX 1970 epoch, a python datetime object, or a
            string in the format of 'YYYY-MM-DDTHH:MM:SS'.
            """
            
            # File-specific information
            WriterBase. __init__(self, filename, ref_time=ref_time, verbose=verbose)
            
            # Open the file and get going
            if os.path.exists(filename):
                if overwrite:
                    shutil.rmtree(filename, ignore_errors=False)
                else:
                    raise IOError(f"File '{filename}' already exists")
            self.basename = filename
            
        def set_geometry(self, site, antennas, bits=8):
            """
            Given a station and an array of stands, set the relevant common observation
            parameters and add entries to the self.array list.
            """
            
            # Update the observatory-specific information
            self.siteName = site.name
            
            stands = []
            for ant in antennas:
                stands.append(ant.stand.id)
            stands = np.array(stands)
            
            arrayX, arrayY, arrayZ = site.geocentric_location
            
            xyz = np.zeros((len(stands),3))
            for i,ant in enumerate(antennas):
                ecef = ant.stand.earth_location.itrs
                xyz[i,:] = ecef.cartesian.xyz.to('m').value
                
            # Create the stand mapper
            mapper = []
            ants = []
            for i in range(len(stands)):
                ants.append( self._Antenna(stands[i], xyz[i,0], xyz[i,1], xyz[i,2], bits=bits) )
                mapper.append( stands[i] )
                
            self.nAnt = len(ants)
            self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'inputAnts': antennas} )
            
        def add_header_keyword(self, name, value, comment=None):
            """
            Add an additional entry to the header of the primary HDU.
            
            .. versionadded:: 2.0.0
            """
            
            raise NotImplementedError
            
        def add_data_set(self, obsTime, intTime, baselines, visibilities, pol='XX', source='z'):
            """
            Create a UVData object to store a collection of visibilities for the
            specified TAI MJD time.
            """
            
            if isinstance(pol, str):
                numericPol = self._STOKES_CODES[pol.upper()]
            else:
                numericPol = pol
                
            if isinstance(obsTime, FrameTimestamp):
                obsTime = obsTime.tai_mjd
            elif isinstance(obsTime, AstroTime):
                obsTime = obsTime.tai.mjd
                
            self.data.append( self._MS_UVData(obsTime, intTime, baselines, visibilities, pol=numericPol, source=source) )
            
        def write(self):
            """
            Fill in the Measurement Set will all of the tables in the
            correct order.
            """
            
            # Validate
            if self.nStokes == 0:
                raise RuntimeError("No polarization setups defined")
            if len(self.freq) == 0:
                raise RuntimeError("No frequency setups defined")
            if self.nAnt == 0:
                raise RuntimeError("No array geometry defined")
            if len(self.data) == 0:
                raise RuntimeError("No visibility data defined")
                
            # Sort the data set
            self.data.sort()
            
            # Write the tables
            self._write_main_table()
            self._write_antenna_table()
            self._write_polarization_table()
            self._write_observation_table()
            self._write_spectralwindow_table()
            self._write_misc_required_tables()
            
            # Fixup the info and keywords for the main table
            tb = table("%s" % self.basename, readonly=False, ack=False)
            tb.putinfo({'type':'Measurement Set', 
                        'readme':'This is a MeasurementSet Table holding measurements from a Telescope'})
            tb.putkeyword('MS_VERSION', np.float32(2.0))
            for filename in sorted(glob.glob('%s/*' % self.basename)):
                if os.path.isdir(filename):
                    tname = os.path.basename(filename)
                    stb = table("%s/%s" % (self.basename, tname), ack=False)
                    tb.putkeyword(tname, stb)
                    stb.close()
            tb.flush()
            tb.close()
                
            # Clear out the data section
            del(self.data[:])
            gc.collect()
            
        def close(self):
            """
            Close out the file.
            """
            
            pass
            
        def _write_antenna_table(self):
            """
            Write the antenna table.
            """
            
            col1 = tableutil.makearrcoldesc('OFFSET', 0.0, 1, 
                                            comment='Axes offset of mount to FEED REFERENCE point', 
                                            keywords={'QuantumUnits':['m','m','m'], 
                                                      'MEASINFO':{'type':'position', 'Ref':'ITRF'}
                                            })
            col2 = tableutil.makearrcoldesc('POSITION', 0.0, 1,
                                            comment='Antenna X,Y,Z phase reference position', 
                                            keywords={'QuantumUnits':['m','m','m'], 
                                                      'MEASINFO':{'type':'position', 'Ref':'ITRF'}
                                                      })
            col3 = tableutil.makescacoldesc('TYPE', "ground-based", 
                                            comment='Antenna type (e.g. SPACE-BASED)')
            col4 = tableutil.makescacoldesc('DISH_DIAMETER', 2.0, 
                                            comment='Physical diameter of dish', 
                                            keywords={'QuantumUnits':['m',]})
            col5 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='Flag for this row')
            col6 = tableutil.makescacoldesc('MOUNT', "alt-az", 
                                            comment='Mount type e.g. alt-az, equatorial, etc.')
            col7 = tableutil.makescacoldesc('NAME', "none", 
                                            comment='Antenna name, e.g. VLA22, CA03')
            col8 = tableutil.makescacoldesc('STATION', self.siteName, 
                                            comment='Station (antenna pad) name')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8])
            tb = table("%s/ANTENNA" % self.basename, desc, nrow=self.nAnt, ack=False)
            
            tb.putcol('OFFSET', np.zeros((self.nAnt,3)), 0, self.nAnt)
            tb.putcol('TYPE', ['GROUND-BASED,']*self.nAnt, 0, self.nAnt)
            tb.putcol('DISH_DIAMETER', [2.0,]*self.nAnt, 0, self.nAnt)
            tb.putcol('FLAG_ROW', [False,]*self.nAnt, 0, self.nAnt)
            tb.putcol('MOUNT', ['ALT-AZ',]*self.nAnt, 0, self.nAnt)
            tb.putcol('NAME', [ant.get_name() for ant in self.array[0]['ants']], 0, self.nAnt)
            tb.putcol('STATION', [self.siteName,]*self.nAnt, 0, self.nAnt)
            
            for i,ant in enumerate(self.array[0]['ants']):
                #tb.putcell('OFFSET', i, [0.0, 0.0, 0.0])
                tb.putcell('POSITION', i, [ant.x + self.array[0]['center'][0],
                                           ant.y + self.array[0]['center'][1], 
                                           ant.z + self.array[0]['center'][2]])
                #tb.putcell('TYPE', i, 'GROUND-BASED')
                #tb.putcell('DISH_DIAMETER', i, 2.0)
                #tb.putcell('FLAG_ROW', i, False)
                #tb.putcell('MOUNT', i, 'ALT-AZ')
                #tb.putcell('NAME', i, ant.get_name())
                #tb.putcell('STATION', i, self.siteName)
                
            tb.flush()
            tb.close()
            
        def _write_polarization_table(self):
            """
            Write the polarization table.
            """
            
            # Polarization
            
            stks = np.array(self.stokes)
            prds = np.zeros((2,self.nStokes), dtype=np.int32)
            for i,stk in enumerate(self.stokes):
                stks[i] = stk
                if stk > 4:
                    prds[0,i] = ((stk-1) % 4) / 2
                    prds[1,i] = ((stk-1) % 4) % 2
                else:
                    prds[0,i] = 1
                    prds[1,i] = 1
                    
            col1 = tableutil.makearrcoldesc('CORR_TYPE', 0, 1, 
                                            comment='The polarization type for each correlation product, as a Stokes enum.')
            col2 = tableutil.makearrcoldesc('CORR_PRODUCT', 0, 2, 
                                            comment='Indices describing receptors of feed going into correlation')
            col3 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='flag')
            col4 = tableutil.makescacoldesc('NUM_CORR', self.nStokes, 
                                            comment='Number of correlation products')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4])
            tb = table("%s/POLARIZATION" % self.basename, desc, nrow=1, ack=False)
            
            tb.putcell('CORR_TYPE', 0, self.stokes)
            tb.putcell('CORR_PRODUCT', 0, prds.T)
            tb.putcell('FLAG_ROW', 0, False)
            tb.putcell('NUM_CORR', 0, self.nStokes)
            
            tb.flush()
            tb.close()
            
            # Feed
            
            col1  = tableutil.makearrcoldesc('POSITION', 0.0, 1, 
                                             comment='Position of feed relative to feed reference position', 
                                             keywords={'QuantumUnits':['m','m','m'], 
                                                       'MEASINFO':{'type':'position', 'Ref':'ITRF'}
                                                       })
            col2  = tableutil.makearrcoldesc('BEAM_OFFSET', 0.0, 2, 
                                             comment='Beam position offset (on sky but in antennareference frame)', 
                                             keywords={'QuantumUnits':['rad','rad'], 
                                                       'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                       })
            col3  = tableutil.makearrcoldesc('POLARIZATION_TYPE', 'X', 1, 
                                             comment='Type of polarization to which a given RECEPTOR responds')
            col4  = tableutil.makearrcoldesc('POL_RESPONSE', 1j, 2,
                                             valuetype='complex',
                                             comment='D-matrix i.e. leakage between two receptors')
            col5  = tableutil.makearrcoldesc('RECEPTOR_ANGLE', 0.0, 1,  
                                             comment='The reference angle for polarization', 
                                             keywords={'QuantumUnits':['rad',]})
            col6  = tableutil.makescacoldesc('ANTENNA_ID', 0, 
                                             comment='ID of antenna in this array')
            col7  = tableutil.makescacoldesc('BEAM_ID', -1, 
                                             comment='Id for BEAM model')
            col8  = tableutil.makescacoldesc('FEED_ID', 0, 
                                             comment='Feed id')
            col9  = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                             comment='Interval for which this set of parameters is accurate', 
                                             keywords={'QuantumUnits':['s',]})
            col10 = tableutil.makescacoldesc('NUM_RECEPTORS', 2, 
                                             comment='Number of receptors on this feed (probably 1 or 2)')
            col11 = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', -1, 
                                             comment='ID for this spectral window setup')
            col12 = tableutil.makescacoldesc('TIME', 0.0, 
                                             comment='Midpoint of time for which this set of parameters is accurate', 
                                             keywords={'QuantumUnits':['s',], 
                                                       'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                       })
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, 
                                          col9, col10, col11, col12])
            tb = table("%s/FEED" % self.basename, desc, nrow=self.nAnt, ack=False)
            
            presp = np.zeros((self.nAnt,2,2), dtype=np.complex64)
            if self.stokes[0] > 8:
                ptype = [['X', 'Y'] for i in range(self.nAnt)]
                presp[:,0,0] = 1.0
                presp[:,0,1] = 0.0
                presp[:,1,0] = 0.0
                presp[:,1,1] = 1.0
            elif self.stokes[0] > 4:
                ptype = [['R', 'L'] for i in range(self.nAnt)]
                presp[:,0,0] = 1.0
                presp[:,0,1] = -1.0j
                presp[:,1,0] = 1.0j
                presp[:,1,1] = 1.0
            else:
                ptype = [['X', 'Y'] for i in range(self.nAnt)]
                presp[:,0,0] = 1.0
                presp[:,0,1] = 0.0
                presp[:,1,0] = 0.0
                presp[:,1,1] = 1.0
                
            tb.putcol('POSITION', np.zeros((self.nAnt,3)), 0, self.nAnt)
            tb.putcol('BEAM_OFFSET', np.zeros((self.nAnt,2,2)), 0, self.nAnt)
            tb.putcol('POLARIZATION_TYPE', np.array(ptype, dtype='S'), 0, self.nAnt)
            tb.putcol('POL_RESPONSE', presp, 0, self.nAnt)
            tb.putcol('RECEPTOR_ANGLE', np.zeros((self.nAnt,2)), 0, self.nAnt)
            tb.putcol('ANTENNA_ID', list(range(self.nAnt)), 0, self.nAnt)
            tb.putcol('BEAM_ID', [-1,]*self.nAnt, 0, self.nAnt)
            tb.putcol('FEED_ID', [0,]*self.nAnt, 0, self.nAnt)
            tb.putcol('INTERVAL', [0.0,]*self.nAnt, 0, self.nAnt)
            tb.putcol('NUM_RECEPTORS', [2,]*self.nAnt, 0, self.nAnt)
            tb.putcol('SPECTRAL_WINDOW_ID', [-1,]*self.nAnt, 0, self.nAnt)
            tb.putcol('TIME', [0.0,]*self.nAnt, 0, self.nAnt)
            
            tb.flush()
            tb.close()
            
        def _write_observation_table(self):
            """
            Write the observation table.
            """
            
            # Observation
            
            col1 = tableutil.makearrcoldesc('TIME_RANGE', 0.0, 1, 
                                            comment='Start and end of observation', 
                                            keywords={'QuantumUnits':['s',], 
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            col2 = tableutil.makearrcoldesc('LOG', 'none', 1,
                                            comment='Observing log')
            col3 = tableutil.makearrcoldesc('SCHEDULE', 'none', 1,
                                            comment='Observing schedule')
            col4 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='Row flag')
            col5 = tableutil.makescacoldesc('OBSERVER', self.observer, 
                                            comment='Name of observer(s)')
            col6 = tableutil.makescacoldesc('PROJECT', self.project, 
                                            comment='Project identification string')
            col7 = tableutil.makescacoldesc('RELEASE_DATE', 0.0, 
                                            comment='Release date when data becomes public', 
                                            keywords={'QuantumUnits':['s',], 
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            col8 = tableutil.makescacoldesc('SCHEDULE_TYPE', self.mode, 
                                            comment='Observing schedule type')
            col9 = tableutil.makescacoldesc('TELESCOPE_NAME', self.siteName, 
                                            comment='Telescope Name (e.g. WSRT, VLBA)')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
            tb = table("%s/OBSERVATION" % self.basename, desc, nrow=1, ack=False)
            
            tStart = astro.taimjd_to_utcjd(self.data[ 0].obsTime) - astro.MJD_OFFSET
            tStop  = astro.taimjd_to_utcjd(self.data[-1].obsTime) - astro.MJD_OFFSET
            
            tb.putcell('TIME_RANGE', 0, [tStart*86400, tStop*86400])
            tb.putcell('LOG', 0, ['Not provided',])
            tb.putcell('SCHEDULE', 0, ['Not provided',])
            tb.putcell('FLAG_ROW', 0, False)
            tb.putcell('OBSERVER', 0, self.observer)
            tb.putcell('PROJECT', 0, self.project)
            tb.putcell('RELEASE_DATE', 0, tStop*86400)
            tb.putcell('SCHEDULE_TYPE', 0, self.mode)
            tb.putcell('TELESCOPE_NAME', 0, self.siteName)
            
            tb.flush()
            tb.close()
            
            # Source
            
            el = EarthLocation.from_geocentric(self.array[0]['center'][0]*astrounits.m,
                                               self.array[0]['center'][1]*astrounits.m,
                                               self.array[0]['center'][2]*astrounits.m,)
            
            nameList = []
            posList = []
            sourceID = 0
            for dataSet in self.data:
                if dataSet.pol == self.stokes[0]:
                    date = AstroTime(dataSet.obsTime, format='mjd', scale='tai')
                    
                    try:
                        currSourceName = dataSet.source.name
                    except AttributeError:
                        currSourceName = dataSet.source
                    
                    if currSourceName not in nameList:
                        sourceID += 1
                        
                        if dataSet.source == 'z':
                            ## Zenith pointings
                            tc = AltAz(0.0*astrounits.deg, 90.0*astrounits.deg,
                                       location=el, obstime=date)
                            equ = tc.transform_to(FK5(equinox=date))
                            
                            # format 'source' name based on local sidereal time
                            raHms = astro.deg_to_hms(equ.ra.deg)
                            (tsecs, secs) = math.modf(raHms.seconds)
                            
                            name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
                            equPo = equ.transform_to(FK5(equinox='J2000'))
                            
                            equ = astro.equ_posn.from_astropy(equ)
                            equPo = astro.equ_posn.from_astropy(equPo)
                            
                        else:
                            ## Real-live sources (ephem.Body instances)
                            name = dataSet.source.name
                            equ = astro.equ_posn(dataSet.source.ra*180/np.pi, dataSet.source.dec*180/np.pi)
                            equPo = astro.equ_posn(dataSet.source.a_ra*180/np.pi, dataSet.source.a_dec*180/np.pi)
                            
                        # J2000 zenith equatorial coordinates
                        posList.append( [[equPo.ra*np.pi/180, equPo.dec*np.pi/180],] )
                        
                        # name
                        nameList.append(name)
                        
            nSource = len(nameList)
            
            # Save these for later since we might need them
            self._sourceTable = nameList
            
            col1  = tableutil.makearrcoldesc('DIRECTION', 0.0, 1, 
                                             comment='Direction (e.g. RA, DEC).', 
                                             keywords={'QuantumUnits':['rad','rad'], 
                                                       'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                       })
            col2  = tableutil.makearrcoldesc('PROPER_MOTION', 0.0, 1, 
                                             comment='Proper motion', 
                                             keywords={'QuantumUnits':['rad/s',]})
            col3  = tableutil.makescacoldesc('CALIBRATION_GROUP', 0, 
                                             comment='Number of grouping for calibration purpose.')
            col4  = tableutil.makescacoldesc('CODE', "none", 
                                             comment='Special characteristics of source, e.g. Bandpass calibrator')
            col5  = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                             comment='Interval of time for which this set of parameters is accurate', 
                                             keywords={'QuantumUnits':['s',]})
            col6  = tableutil.makescacoldesc('NAME', "none", 
                                             comment='Name of source as given during observations')
            col7  = tableutil.makescacoldesc('NUM_LINES', 0, 
                                             comment='Number of spectral lines')
            col8  = tableutil.makescacoldesc('SOURCE_ID', 0, 
                                             comment='Source id')
            col9  = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', -1, 
                                             comment='ID for this spectral window setup')
            col10 = tableutil.makescacoldesc('TIME', 0.0,
                                             comment='Midpoint of time for which this set of parameters is accurate.', 
                                             keywords={'QuantumUnits':['s',], 
                                                       'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                       })
            col11 = tableutil.makearrcoldesc('TRANSITION', 'none', 1, 
                                             comment='Line Transition name')
            col12 = tableutil.makearrcoldesc('REST_FREQUENCY', 1.0, 1, 
                                             comment='Line rest frequency', 
                                             keywords={'QuantumUnits':['Hz',], 
                                                       'MEASINFO':{'type':'frequency', 
                                                                   'Ref':'LSRK'}
                                                       })
            col13 = tableutil.makearrcoldesc('SYSVEL', 1.0, 1, 
                                             comment='Systemic velocity at reference', 
                                             keywords={'QuantumUnits':['m/s',], 
                                                       'MEASINFO':{'type':'radialvelocity', 
                                                                   'Ref':'LSRK'}
                                                       })
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9, 
                                          col10, col11, col12, col13])
            tb = table("%s/SOURCE" % self.basename, desc, nrow=nSource, ack=False)
            
            for i in range(nSource):
                tb.putcell('DIRECTION', i, np.array(posList[i])[0])
                tb.putcell('PROPER_MOTION', i, np.array([[0.0, 0.0],])[0])
                tb.putcell('CALIBRATION_GROUP', i, 0)
                tb.putcell('CODE', i, 'none')
                tb.putcell('INTERVAL', i, 0.0)
                tb.putcell('NAME', i, str(nameList[i]))
                tb.putcell('NUM_LINES', i, 0)
                tb.putcell('SOURCE_ID', i, i)
                tb.putcell('SPECTRAL_WINDOW_ID', i, -1)
                tb.putcell('TIME', i, (tStart+tStop)/2*86400)
                #tb.putcell('TRANSITION', i, [])
                #tb.putcell('REST_FREQUENCY', i, [])
                #tb.putcell('SYSVEL', i, [])
                
            tb.flush()
            tb.close()
            
            # Field
            
            col1 = tableutil.makearrcoldesc('DELAY_DIR', 0.0, 2, 
                                            comment='Direction of delay center (e.g. RA, DEC)as polynomial in time.', 
                                            keywords={'QuantumUnits':['rad','rad'], 
                                                      'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                      })
            col2 = tableutil.makearrcoldesc('PHASE_DIR', 0.0, 2, 
                                            comment='Direction of phase center (e.g. RA, DEC).', 
                                            keywords={'QuantumUnits':['rad','rad'], 
                                                      'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                      })
            col3 = tableutil.makearrcoldesc('REFERENCE_DIR', 0.0, 2, 
                                            comment='Direction of REFERENCE center (e.g. RA, DEC).as polynomial in time.', 
                                            keywords={'QuantumUnits':['rad','rad'], 
                                                      'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                      })
            col4 = tableutil.makescacoldesc('CODE', "none", 
                                            comment='Special characteristics of field, e.g. Bandpass calibrator')
            col5 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='Row Flag')
            col6 = tableutil.makescacoldesc('NAME', "none", 
                                            comment='Name of this field')
            col7 = tableutil.makescacoldesc('NUM_POLY', 0, 
                                            comment='Polynomial order of _DIR columns')
            col8 = tableutil.makescacoldesc('SOURCE_ID', 0, 
                                            comment='Source id')
            col9 = tableutil.makescacoldesc('TIME', (tStart+tStop)/2, 
                                            comment='Time origin for direction and rate', 
                                            keywords={'QuantumUnits':['s',],
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
            tb = table("%s/FIELD" % self.basename, desc, nrow=nSource, ack=False)
            
            for i in range(nSource):
                tb.putcell('DELAY_DIR', i, np.array(posList[i]))
                tb.putcell('PHASE_DIR', i, np.array(posList[i]))
                tb.putcell('REFERENCE_DIR', i, np.array(posList[i]))
                tb.putcell('CODE', i, 'None')
                tb.putcell('FLAG_ROW', i, False)
                tb.putcell('NAME', i, nameList[i])
                tb.putcell('NUM_POLY', i, 0)
                tb.putcell('SOURCE_ID', i, i)
                tb.putcell('TIME', i, (tStart+tStop)/2*86400)
                
            tb.flush()
            tb.close()
            
        def _write_spectralwindow_table(self):
            """
            Write the spectral window table.
            """
            
            # Spectral Window
            
            nBand = len(self.freq)
            
            col1  = tableutil.makescacoldesc('MEAS_FREQ_REF', 0, 
                                             comment='Frequency Measure reference')
            col2  = tableutil.makearrcoldesc('CHAN_FREQ', 0.0, 1, 
                                             comment='Center frequencies for each channel in the data matrix', 
                                             keywords={'QuantumUnits':['Hz',], 
                                                       'MEASINFO':{'type':'frequency', 
                                                                   'VarRefCol':'MEAS_FREQ_REF', 
                                                                   'TabRefTypes':['REST','LSRK','LSRD','BARY','GEO','TOPO','GALACTO','LGROUP','CMB','Undefined'],
                                                                   'TabRefCodes':np.array([0,1,2,3,4,5,6,7,8,64], dtype=np.uint32)}
                                                       })
            col3  = tableutil.makescacoldesc('REF_FREQUENCY', self.refVal, 
                                             comment='The reference frequency', 
                                             keywords={'QuantumUnits':['Hz',], 
                                                       'MEASINFO':{'type':'frequency', 
                                                                   'VarRefCol':'MEAS_FREQ_REF', 
                                                                   'TabRefTypes':['REST','LSRK','LSRD','BARY','GEO','TOPO','GALACTO','LGROUP','CMB','Undefined'],
                                                                   'TabRefCodes':np.array([0,1,2,3,4,5,6,7,8,64], dtype=np.uint32)}
                                                       })
            col4  = tableutil.makearrcoldesc('CHAN_WIDTH', 0.0, 1, 
                                             comment='Channel width for each channel', 
                                             keywords={'QuantumUnits':['Hz',]})
            col5  = tableutil.makearrcoldesc('EFFECTIVE_BW', 0.0, 1, 
                                             comment='Effective noise bandwidth of each channel', 
                                             keywords={'QuantumUnits':['Hz',]})
            col6  = tableutil.makearrcoldesc('RESOLUTION', 0.0, 1, 
                                             comment='The effective noise bandwidth for each channel', 
                                             keywords={'QuantumUnits':['Hz',]})
            col7  = tableutil.makescacoldesc('FLAG_ROW', False, 
                                             comment='flag')
            col8  = tableutil.makescacoldesc('FREQ_GROUP', 1, 
                                             comment='Frequency group')
            col9  = tableutil.makescacoldesc('FREQ_GROUP_NAME', "group1", 
                                             comment='Frequency group name')
            col10 = tableutil.makescacoldesc('IF_CONV_CHAIN', 0, 
                                             comment='The IF conversion chain number')
            col11 = tableutil.makescacoldesc('NAME', "%i channels" % self.nChan, 
                                             comment='Spectral window name')
            col12 = tableutil.makescacoldesc('NET_SIDEBAND', 0, 
                                             comment='Net sideband')
            col13 = tableutil.makescacoldesc('NUM_CHAN', 0, 
                                             comment='Number of spectral channels')
            col14 = tableutil.makescacoldesc('TOTAL_BANDWIDTH', 0.0, 
                                             comment='The total bandwidth for this window', 
                                             keywords={'QuantumUnits':['Hz',]})
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9, 
                                          col10, col11, col12, col13, col14])
            tb = table("%s/SPECTRAL_WINDOW" % self.basename, desc, nrow=nBand, ack=False)
            
            for i,freq in enumerate(self.freq):
                tb.putcell('MEAS_FREQ_REF', i, 0)
                tb.putcell('CHAN_FREQ', i, self.refVal + freq.bandFreq + np.arange(self.nChan)*self.channelWidth)
                tb.putcell('REF_FREQUENCY', i, self.refVal)
                tb.putcell('CHAN_WIDTH', i, [freq.chWidth for j in range(self.nChan)])
                tb.putcell('EFFECTIVE_BW', i, [freq.chWidth for j in range(self.nChan)])
                tb.putcell('RESOLUTION', i, [freq.chWidth for j in range(self.nChan)])
                tb.putcell('FLAG_ROW', i, False)
                tb.putcell('FREQ_GROUP', i, i+1)
                tb.putcell('FREQ_GROUP_NAME', i, 'group%i' % (i+1))
                tb.putcell('IF_CONV_CHAIN', i, i)
                tb.putcell('NAME', i, "IF %i, %i channels" % (i+1, self.nChan))
                tb.putcell('NET_SIDEBAND', i, 0)
                tb.putcell('NUM_CHAN', i, self.nChan)
                tb.putcell('TOTAL_BANDWIDTH', i, freq.totalBW)
                
            tb.flush()
            tb.close()
            
        def _write_main_table(self):
            """
            Write the main table.
            """
            
            # Main
            
            nBand = len(self.freq)
            
            el = EarthLocation.from_geocentric(self.array[0]['center'][0]*astrounits.m,
                                               self.array[0]['center'][1]*astrounits.m,
                                               self.array[0]['center'][2]*astrounits.m,)
            
            mapper = self.array[0]['mapper']
            
            col1  = tableutil.makearrcoldesc('UVW', 0.0, 1, 
                                             comment='Vector with uvw coordinates (in meters)', 
                                             keywords={'QuantumUnits':['m','m','m'], 
                                                       'MEASINFO':{'type':'uvw', 'Ref':'ITRF'}
                                                       })
            col2  = tableutil.makearrcoldesc('FLAG', False, 2, 
                                             comment='The data flags, array of bools with same shape as data')
            col3  = tableutil.makearrcoldesc('FLAG_CATEGORY', False, 3,  
                                             comment='The flag category, NUM_CAT flags for each datum', 
                                             keywords={'CATEGORY':['',]})
            col4  = tableutil.makearrcoldesc('WEIGHT', 1.0, 1, 
                                             valuetype='float', 
                                             comment='Weight for each polarization spectrum')
            col5  = tableutil.makearrcoldesc('SIGMA', 9999., 1, 
                                             valuetype='float', 
                                             comment='Estimated rms noise for channel with unity bandpass response')
            col6  = tableutil.makescacoldesc('ANTENNA1', 0, 
                                             comment='ID of first antenna in interferometer')
            col7  = tableutil.makescacoldesc('ANTENNA2', 0, 
                                             comment='ID of second antenna in interferometer')
            col8  = tableutil.makescacoldesc('ARRAY_ID', 0, 
                                             comment='ID of array or subarray')
            col9  = tableutil.makescacoldesc('DATA_DESC_ID', 0, 
                                             comment='The data description table index')
            col10 = tableutil.makescacoldesc('EXPOSURE', 0.0, 
                                             comment='he effective integration time', 
                                             keywords={'QuantumUnits':['s',]})
            col11 = tableutil.makescacoldesc('FEED1', 0, 
                                             comment='The feed index for ANTENNA1')
            col12 = tableutil.makescacoldesc('FEED2', 0, 
                                             comment='The feed index for ANTENNA2')
            col13 = tableutil.makescacoldesc('FIELD_ID', 0, 
                                             comment='Unique id for this pointing')
            col14 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                             comment='Row flag - flag all data in this row if True')
            col15 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                             comment='The sampling interval', 
                                             keywords={'QuantumUnits':['s',]})
            col16 = tableutil.makescacoldesc('OBSERVATION_ID', 0, 
                                             comment='ID for this observation, index in OBSERVATION table')
            col17 = tableutil.makescacoldesc('PROCESSOR_ID', -1, 
                                             comment='Id for backend processor, index in PROCESSOR table')
            col18 = tableutil.makescacoldesc('SCAN_NUMBER', 1, 
                                             comment='Sequential scan number from on-line system')
            col19 = tableutil.makescacoldesc('STATE_ID', -1, 
                                             comment='ID for this observing state')
            col20 = tableutil.makescacoldesc('TIME', 0.0, 
                                             comment='Modified Julian Day', 
                                             keywords={'QuantumUnits':['s',],
                                                       'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                       })
            col21 = tableutil.makescacoldesc('TIME_CENTROID', 0.0, 
                                             comment='Modified Julian Day', 
                                             keywords={'QuantumUnits':['s',],
                                                       'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                       })
            col22 = tableutil.makearrcoldesc("DATA", 0j, 2, 
                                             valuetype='complex',
                                             comment='The data column')
            
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
                    if len(dataSet.visibilities) != len(order):
                        raise NameError
                except NameError:
                    order = dataSet.argsort(mapper=mapper, shift=16)
                    try:
                        del baselineMapped
                    except NameError:
                        pass
                        
                # Deal with defininig the values of the new data set
                if dataSet.pol == self.stokes[0]:
                    ## Figure out the new date/time for the observation
                    utc = astro.taimjd_to_utcjd(dataSet.obsTime)
                    date = AstroTime(utc, format='jd', scale='utc')
                    utc0 = AstroTime(f"{date.ymdhms[0]}-{date.ymdhms[1]}-{date.ymdhms[2]} 00:00:00", format='iso', scale='utc')
                    utc0 = utc0.jd
                    
                    ## Update the observer so we can figure out where the source is
                    if dataSet.source == 'z':
                        ### Zenith pointings
                        tc = AltAz(0.0*astrounits.deg, 90.0*astrounits.deg,
                                   location=el, obstime=date)
                        equ = tc.transform_to(FK5(equinox=date))
                        
                        ### format 'source' name based on local sidereal time
                        raHms = astro.deg_to_hms(equ.ra.deg)
                        (tsecs, secs) = math.modf(raHms.seconds)
                        name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
                    else:
                        ### Real-live sources (ephem.Body instances)
                        equ = FK5(dataSet.source.a_ra*astrounits.rad, dataSet.source.a_dec*astrounits.rad,
                                  equinox=date)
                        
                        name = dataSet.source.name
                        
                    ## Update the source ID
                    try:
                        sourceID = _sourceTable.index(name)
                    except ValueError:
                        _sourceTable.append(name)
                        sourceID = _sourceTable.index(name)
                        
                    ## Compute the uvw coordinates of all baselines
                    it = equ.transform_to(ITRS(location=el, obstime=date))
                    HA = ((el.lon - it.spherical.lon).wrap_at('180deg')).hourangle
                    dec = it.spherical.lat.deg
                    uvwCoords = dataSet.get_uvw(HA, dec, el)
                    
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
                    #timeList = [utc - astro.MJD_OFFSET,]*len(dataSet.baselines)
                    intTimeList = [dataSet.intTime,]*len(dataSet.baselines)
                    timeList = [(utc - astro.MJD_OFFSET)*86400 + dataSet.intTime/2.0,]*len(dataSet.baselines)
                    
                    ### Add in the new new source ID and name
                    sourceList = [sourceID,]*len(dataSet.baselines)
                     
                    ### Zero out the visibility data
                    try:
                        if matrix.shape[0] != len(order):
                            raise NameError
                        matrix.shape = (len(order), self.nStokes, nBand*self.nChan)
                        matrix *= 0.0
                    except NameError:
                        matrix = np.zeros((len(order), self.nStokes, self.nChan*nBand), dtype=np.complex64)
                        
                # Save the visibility data in the right order
                matrix[:,self.stokes.index(dataSet.pol),:] = dataSet.visibilities[order,:]
                
                # Deal with saving the data once all of the polarizations have been added to 'matrix'
                if dataSet.pol == self.stokes[-1]:
                    nBL = uvwList.shape[0]
                    tb.addrows(nBand*nBL)
                    
                    matrix.shape = (len(order), self.nStokes, nBand, self.nChan)
                    
                    for j in range(nBand):
                        fg = np.zeros((nBL,self.nStokes,self.nChan), dtype=bool)
                        fc = np.zeros((nBL,self.nStokes,self.nChan,1), dtype=bool)
                        wg = np.ones((nBL,self.nStokes))
                        sg = np.ones((nBL,self.nStokes))*9999
                        
                        tb.putcol('UVW', uvwList, i, nBL)
                        tb.putcol('FLAG', fg.transpose(0,2,1), i, nBL)
                        tb.putcol('FLAG_CATEGORY', fc.transpose(0,3,2,1), i, nBL)
                        tb.putcol('WEIGHT', wg, i, nBL)
                        tb.putcol('SIGMA', sg, i, nBL)
                        tb.putcol('ANTENNA1', ant1List, i, nBL)
                        tb.putcol('ANTENNA2', ant2List, i, nBL)
                        tb.putcol('ARRAY_ID', [0,]*nBL, i, nBL)
                        tb.putcol('DATA_DESC_ID', [j,]*nBL, i, nBL)
                        tb.putcol('EXPOSURE', intTimeList, i, nBL)
                        tb.putcol('FEED1', [0,]*nBL, i, nBL)
                        tb.putcol('FEED2', [0,]*nBL, i, nBL)
                        tb.putcol('FIELD_ID', sourceList, i, nBL)
                        tb.putcol('FLAG_ROW', [False,]*nBL, i, nBL)
                        tb.putcol('INTERVAL', intTimeList, i, nBL)
                        tb.putcol('OBSERVATION_ID', [0,]*nBL, i, nBL)
                        tb.putcol('PROCESSOR_ID', [-1,]*nBL, i, nBL)
                        tb.putcol('SCAN_NUMBER', [s,]*nBL, i, nBL)
                        tb.putcol('STATE_ID', [-1,]*nBL, i, nBL)
                        tb.putcol('TIME', timeList, i, nBL)
                        tb.putcol('TIME_CENTROID', timeList, i, nBL)
                        tb.putcol('DATA', matrix[...,j,:].transpose(0,2,1), i, nBL)
                        i += nBL
                    s += 1
                    
            tb.flush()
            tb.close()
            
            # Data description
            
            col1 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='Flag this row')
            col2 = tableutil.makescacoldesc('POLARIZATION_ID', 0, 
                                            comment='Pointer to polarization table')
            col3 = tableutil.makescacoldesc('SPECTRAL_WINDOW_ID', 0, 
                                            comment='Pointer to spectralwindow table')
            
            desc = tableutil.maketabdesc([col1, col2, col3])
            tb = table("%s/DATA_DESCRIPTION" % self.basename, desc, nrow=nBand, ack=False)
            
            for i in range(nBand):
                tb.putcell('FLAG_ROW', i, False)
                tb.putcell('POLARIZATION_ID', i, 0)
                tb.putcell('SPECTRAL_WINDOW_ID', i, i)
                
            tb.flush()
            tb.close()
            
        def _write_misc_required_tables(self):
            """
            Write the other tables that are part of the measurement set but 
            don't contain anything by default.
            """
            
            # Flag command
            
            col1 = tableutil.makescacoldesc('TIME', 0.0, 
                                            comment='Midpoint of interval for which this flag is valid', 
                                            keywords={'QuantumUnits':['s',], 
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            col2 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                            comment='Time interval for which this flag is valid', 
                                            keywords={'QuantumUnits':['s',]})
            col3 = tableutil.makescacoldesc('TYPE', 'flag', 
                                            comment='Type of flag (FLAG or UNFLAG)')
            col4 = tableutil.makescacoldesc('REASON', 'reason', 
                                            comment='Flag reason')
            col5 = tableutil.makescacoldesc('LEVEL', 0, 
                                            comment='Flag level - revision level')
            col6 = tableutil.makescacoldesc('SEVERITY', 0, 
                                            comment='Severity code (0-10)')
            col7 = tableutil.makescacoldesc('APPLIED', False, 
                                            comment='True if flag has been applied to main table')
            col8 = tableutil.makescacoldesc('COMMAND', 'command', 
                                            comment='Flagging command')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8])
            tb = table("%s/FLAG_CMD" % self.basename, desc, nrow=0, ack=False)
            
            tb.flush()
            tb.close()
            
            # History
            
            col1 = tableutil.makescacoldesc('TIME', 0.0, 
                                            comment='Timestamp of message', 
                                            keywords={'QuantumUnits':['s',], 
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            col2 = tableutil.makescacoldesc('OBSERVATION_ID', 0, 
                                            comment='Observation id (index in OBSERVATION table)')
            col3 = tableutil.makescacoldesc('MESSAGE', 'message', 
                                            comment='Log message')
            col4 = tableutil.makescacoldesc('PRIORITY', 'NORMAL', 
                                            comment='Message priority')
            col5 = tableutil.makescacoldesc('ORIGIN', 'origin', 
                                            comment='(Source code) origin from which message originated')
            col6 = tableutil.makescacoldesc('OBJECT_ID', 0, 
                                            comment='Originating ObjectID')
            col7 = tableutil.makescacoldesc('APPLICATION', 'application', 
                                            comment='Application name')
            col8 = tableutil.makearrcoldesc('CLI_COMMAND', 'command', 1, 
                                            comment='CLI command sequence')
            col9 = tableutil.makearrcoldesc('APP_PARAMS', 'params', 1, 
                                            comment='Application parameters')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
            tb = table("%s/HISTORY" % self.basename, desc, nrow=0, ack=False)
            
            tb.flush()
            tb.close()
            
            # POINTING
            
            col1 = tableutil.makescacoldesc('ANTENNA_ID', 0, 
                                            comment='Antenna Id')
            col2 = tableutil.makescacoldesc('TIME', 0.0, 
                                            comment='Time interval midpoint', 
                                            keywords={'QuantumUnits':['s',], 
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            col3 = tableutil.makescacoldesc('INTERVAL', 0.0, 
                                            comment='Time interval', 
                                            keywords={'QuantumUnits':['s',]})
            col4 = tableutil.makescacoldesc('NAME', 'name', 
                                            comment='Pointing position name')
            col5 = tableutil.makescacoldesc('NUM_POLY', 0, 
                                            comment='Series order')
            col6 = tableutil.makescacoldesc('TIME_ORIGIN', 0.0, 
                                            comment='Time origin for direction', 
                                            keywords={'QuantumUnits':['s',], 
                                                      'MEASINFO':{'type':'epoch', 'Ref':'UTC'}
                                                      })
            col7 = tableutil.makearrcoldesc('DIRECTION', 0.0, 2, 
                                            comment='Antenna pointing direction as polynomial in time', 
                                            keywords={'QuantumUnits':['rad','rad'], 
                                                      'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                      })
            col8 = tableutil.makearrcoldesc('TARGET', 0.0, 2, 
                                            comment='target direction as polynomial in time',
                                            keywords={'QuantumUnits':['rad','rad'], 
                                                      'MEASINFO':{'type':'direction', 'Ref':'J2000'}
                                                      })
            col9 = tableutil.makescacoldesc('TRACKING', True, 
                                            comment='Tracking flag - True if on position')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7, col8, col9])
            tb = table("%s/POINTING" % self.basename, desc, nrow=0, ack=False)
            
            tb.flush()
            tb.close()
            
            # Processor
            
            col1 = tableutil.makescacoldesc('TYPE', 'type', 
                                            comment='Processor type')
            col2 = tableutil.makescacoldesc('SUB_TYPE', 'subtype', 
                                            comment='Processor sub type')
            col3 = tableutil.makescacoldesc('TYPE_ID', 0, 
                                            comment='Processor type id')
            col4 = tableutil.makescacoldesc('MODE_ID', 0, 
                                            comment='Processor mode id')
            col5 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='flag')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5])
            tb = table("%s/PROCESSOR" % self.basename, desc, nrow=0, ack=False)
            
            tb.flush()
            tb.close()
            
            # State
            
            col1 = tableutil.makescacoldesc('SIG', True, 
                                            comment='True for a source observation')
            col2 = tableutil.makescacoldesc('REF', False, 
                                            comment='True for a reference observation')
            col3 = tableutil.makescacoldesc('CAL', 0.0, 
                                            comment='Noise calibration temperature', 
                                            keywords={'QuantumUnits':['K',]})
            col4 = tableutil.makescacoldesc('LOAD', 0.0, 
                                            comment='Load temperature', 
                                            keywords={'QuantumUnits':['K',]})
            col5 = tableutil.makescacoldesc('SUB_SCAN', 0, 
                                            comment='Sub scan number, relative to scan number')
            col6 = tableutil.makescacoldesc('OBS_MODE', 'mode', 
                                            comment='Observing mode, e.g., OFF_SPECTRUM')
            col7 = tableutil.makescacoldesc('FLAG_ROW', False, 
                                            comment='Row flag')
            
            desc = tableutil.maketabdesc([col1, col2, col3, col4, col5, col6, col7])
            tb = table("%s/STATE" % self.basename, desc, nrow=0, ack=False)
            
            tb.flush()
            tb.close()
            
except ImportError:
    warnings.warn(colorfy('{{%yellow Cannot import casacore.tables, MS support disabled}}'), RuntimeWarning)
    
    class Ms(WriterBase):
        """
        Class for storing visibility data and writing the data, along with array
        geometry, frequency setup, etc., to a CASA measurement set.
        """
        
        _STOKES_CODES = STOKES_CODES
        
        def __init__(self, filename, ref_time=0.0, verbose=False, memmap=None, overwrite=False):
            """
            Initialize a new Measurement Set object using a filename and a reference time
            given in seconds since the UNIX 1970 epoch, a python datetime object, or a 
            string in the format of 'YYYY-MM-DDTHH:MM:SS'.
            """
            
            raise RuntimeError("Cannot import casacore.tables, MS support disabled")
