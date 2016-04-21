# -*- coding: utf-8 -*-

"""
Module for writing correlator output to a FITS IDI file.  The classes and 
functions defined in this module are based heavily off the lwda_fits library.

.. note::
	Some software that deal with FITS IDI files, e.g. AIPS, does not support
	FITS IDI files with more than 99 antennas.  See the :class:`lsl.writer.fitsidi.aips`
	class for a FITS IDI writer that respects this limit.
"""

import os
import gc
import re
import sys
import math
import ephem
import numpy
import pyfits
from datetime import datetime

from lsl import astro
from lsl.common import constants
from lsl.common import dp as dp_common
from lsl.common.stations import geo2ecef
from lsl.correlator import uvUtils
from lsl.misc import mathutil
from lsl.misc import geodesy

__version__ = '0.7'
__revision__ = '$Rev$'
__all__ = ['IDI', 'AIPS', 'ExtendedIDI', 'StokesCodes', 'NumericStokes', '__version__', '__revision__', '__all__']


IDIVersion = (2, 0)

StokesCodes = { 'I':  1,  'Q': 2,   'U':  3,  'V':  4, 
			'RR': -1, 'LL': -2, 'RL': -3, 'LR': -4, 
			'XX': -5, 'YY': -6, 'XY': -7, 'YX': -8}

NumericStokes = { 1: 'I',   2: 'Q',   3: 'U',   4: 'V', 
			  -1: 'RR', -2: 'LL', -3: 'RL', -4: 'LR', 
			  -5: 'XX', -6: 'YY', -7: 'XY', -8: 'YX'}


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



class IDI(object):
	"""
	Class for storing visibility data and writing the data, along with array
	geometry, frequency setup, etc., to a FITS IDI file that can be read into 
	AIPS via the FITLD task.
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
			
	class _UVData(object):
		"""
		Represents one UV visibility data set for a given observation time.
		"""
    
		def __init__(self, obsTime, intTime, dataDict, pol=StokesCodes['XX'], source='z'):
			self.obsTime = obsTime
			self.intTime = intTime
			self.dataDict = dataDict
			self.pol = pol
			self.source = source
			
		def time(self):
			return self.obsTime
			
	def parseRefTime(self, refTime):
		"""
		Given a time as either a integer, float, string, or datetime object, 
		convert it to a string in the formation 'YYYY-MM-DDTHH:MM:SS'.
		"""
		
		# Valid time string (modulo the 'T')
		timeRE = re.compile(r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?')
		
		if type(refTime).__name__ in ['int', 'long', 'float']:
			refDateTime = datetime.utcfromtimestamp(refTime)
			refTime = refDateTime.strftime("%Y-%m-%dT%H:%M:%S")
		elif type(refTime).__name__ in ['datetime']:
			refTime = refTime.strftime("%Y-%m-%dT%H:%M:%S")
		elif type(refTime).__name__ in ['str']:
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
		
	def __init__(self, filename, refTime=0.0, verbose=False):
		"""
		Initialize a new FITS IDI object using a filename and a reference time 
		given in seconds since the UNIX 1970 ephem, a python datetime object, or a 
		string in the format of 'YYYY-MM-DDTHH:MM:SS'.
		"""
		
		# File-specific information
		self.filename = filename
		self.verbose = verbose
		
		# Observator-specific information
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
		self.FITS = pyfits.open(filename, mode='append')
		
	def setStokes(self, polList):
		"""
		Given a list of Stokes parameters, update the object's parameters.
		"""
		
		for pol in polList:
			if type(pol).__name__ == 'str':
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
		
		# Make sure that we have been passed 255 or fewer stands
		if len(antennas) > 255:
			raise RuntimeError("FITS IDI supports up to 255 antennas only, given %i" % len(antennas))
			
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
			
		# Create the stand mapper to deal with the fact that stands range from 
		# 1 to 258, not 1 to 255
		mapper = {}
		if stands.max() > 255:
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
			print "FITS IDI: stand ID mapping enabled"
			for key, value in mapper.iteritems():
				print "FITS IDI:  stand #%i -> mapped #%i" % (key, value)
				
		self.nAnt = len(ants)
		self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'enableMapper': enableMapper, 'inputAnts': antennas} )
		
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
		
		if type(pol).__name__ == 'str':
			numericPol = StokesCodes[pol.upper()]
		else:
			numericPol = pol
			
		dataDict = {}
		for bl, (ant1,ant2) in enumerate(baselines):
			baseline = mergeBaseline(ant1.stand.id, ant2.stand.id, shift=16)
			dataDict[baseline] = visibilities[bl,:].squeeze()
			
		self.data.append( self._UVData(obsTime, intTime, dataDict, pol=numericPol, source=source) )
		
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
		self.data.sort(cmp=__sortData)
		
		self._writePrimary()
		self._writeGeometry()
		self._writeFrequency()
		self._writeAntenna()
		self._writeBandpass()
		self._writeSource()
		self._writeData()
		
		# Clear out the data section
		del(self.data[:])
		gc.collect()
		
	def close(self):
		"""
		Close out the file.
		"""
		
		self.FITS.flush()
		self.FITS.close()
		
	def _addCommonKeywords(self, hdr, name, revision):
		"""
		Added keywords common to all table headers.
		"""
		
		hdr['EXTNAME'] = (name, 'FITS-IDI table name')
		hdr['EXTVER'] = (1, 'table instance number') 
		hdr['TABREV'] = (revision, 'table format revision number')
		hdr['NO_STKD'] = (self.nStokes, 'number of Stokes parameters')
		hdr['STK_1'] = (self.stokes[0], 'first Stokes parameter')
		hdr['NO_BAND'] = (1, 'number of frequency bands')
		hdr['NO_CHAN'] = (self.nChan, 'number of frequency channels')
		hdr['REF_FREQ'] = (self.refVal, 'reference frequency (Hz)')
		hdr['CHAN_BW'] = (self.channelWidth, 'channel bandwidth (Hz)')
		hdr['REF_PIXL'] = (float(self.refPix), 'reference frequency bin')
		
		date = self.refTime.split('-')
		name = "ZA%s%s%s" % (date[0][2:], date[1], date[2])
		hdr['OBSCODE'] = (name, 'zenith all-sky image')
		
		hdr['ARRNAM'] = self.siteName      
		hdr['RDATE'] = (self.refTime, 'file data reference date')
		
	def _makeAppendTable(self, extension, AddRows=1):
		"""
		Private function to make a temporary table for appending data.
		"""
		
		nrows = self.hdulist[extension].data.shape[0]
		tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+AddRows)
		for key in list(self.hdulist[extension].header.keys()):
			tempHDU.header[key] = self.hdulist[extension].header[key]
			
		return tempHDU
		
	def _applyAppendTable(self, extension, tempHDU):
		"""
		Private function to replace the given extension with the temporary
		table.
		"""
		
		self.hdulist[extension] = tempHDU
		self.flush()
		
	def _writePrimary(self):
		"""
		Write the primary HDU to file.
		"""
		
		primary = pyfits.PrimaryHDU()
		
		primary.header['NAXIS'] = (0, 'indicates IDI file')
		primary.header['EXTEND'] = (True, 'indicates IDI file')
		primary.header['GROUPS'] = (True, 'indicates IDI file')
		primary.header['GCOUNT'] = 0
		primary.header['PCOUNT'] = 0
		primary.header['OBJECT'] = 'BINARYTB'
		primary.header['TELESCOP'] = self.siteName
		primary.header['INSTRUME'] = self.siteName
		primary.header['OBSERVER'] = ('ZASKY', 'zenith all-sky image')
		primary.header['ORIGIN'] = 'LSL'
		primary.header['CORRELAT'] = ('LWASWC', 'Correlator used')
		primary.header['FXCORVER'] = ('1', 'Correlator version')
		primary.header['LWATYPE'] = ('IDI-ZA', 'LWA FITS file type')
		primary.header['LWAMAJV'] = (IDIVersion[0], 'LWA FITS file format major version')
		primary.header['LWAMINV'] = (IDIVersion[1], 'LWA FITS file format minor version')
		primary.header['DATE-OBS'] = (self.refTime, 'IDI file data collection date')
		ts = str(astro.get_date_from_sys())
		primary.header['DATE-MAP'] = (ts.split()[0], 'IDI file creation date')
		
		self.FITS.append(primary)
		self.FITS.flush()
		
	def _writeGeometry(self):
		"""
		Define the Array_Geometry table (group 1, table 1).
		"""
		
		i = 0
		names = []
		xyz = numpy.zeros((self.nAnt,3), dtype=numpy.float64)
		for ant in self.array[0]['ants']:
			xyz[i,:] = numpy.array([ant.x, ant.y, ant.z])
			names.append(ant.getName())
			i = i + 1
			
		# Antenna name
		c1 = pyfits.Column(name='ANNAME', format='A8', 
						array=numpy.array([ant.getName() for ant in self.array[0]['ants']]))
		# Station coordinates in meters
		c2 = pyfits.Column(name='STABXYZ', unit='METERS', format='3D', 
						array=xyz)
		# First order derivative of station coordinates in m/s
		c3 = pyfits.Column(name='DERXYZ', unit='METERS/S', format='3E', 
						array=numpy.zeros((self.nAnt,3), dtype=numpy.float32))
		# Orbital elements
		c4 = pyfits.Column(name='ORBPARM', format='1D', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.float64))
		# Station number
		c5 = pyfits.Column(name='NOSTA', format='1J', 
						array=numpy.array([self.array[0]['mapper'][ant.id] for ant in self.array[0]['ants']]))
		# Mount type (0 == alt-azimuth)
		c6 = pyfits.Column(name='MNTSTA', format='1J', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.int32))
		# Axis offset in meters
		c7 = pyfits.Column(name='STAXOF', unit='METERS', format='3E', 
						array=numpy.zeros((self.nAnt,3), dtype=numpy.float32))
						
		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
		
		# Create the table and fill in the header
		ag = pyfits.new_table(colDefs)
		self._addCommonKeywords(ag.header, 'ARRAY_GEOMETRY', 1)
		
		ag.header['EXTVER'] = (1, 'array ID')
		ag.header['ARRNAM'] = self.siteName
		ag.header['FRAME'] = ('GEOCENTRIC', 'coordinate system')
		ag.header['NUMORB'] = (0, 'number of orbital parameters')
		ag.header['FREQ'] = (self.refVal, 'reference frequency (Hz)')
		ag.header['TIMSYS'] = ('UTC', 'time coordinate system')
		
		date = self.refTime2AstroDate()
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
		
		refDate = self.refTime2AstroDate()
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
		
		ag.name = 'ARRAY_GEOMETRY'
		self.FITS.append(ag)
		self.FITS.flush()
		
		if self.array[0]['enableMapper']:
			self._writeMapper()
			
	def _writeFrequency(self):
		"""
		Define the Frequency table (group 1, table 3).
		"""
		
		# Frequency setup number
		c1 = pyfits.Column(name='FREQID', format='1J', 
						array=numpy.array([self.freq[0].id], dtype=numpy.int32))
		# Frequency offsets in Hz
		c2 = pyfits.Column(name='BANDFREQ', format='1D', unit='HZ', 
						array=numpy.array([self.freq[0].bandFreq], dtype=numpy.float64))
		# Channel width in Hz
		c3 = pyfits.Column(name='CH_WIDTH', format='1E', unit='HZ', 
						array=numpy.array([self.freq[0].chWidth], dtype=numpy.float32))
		# Total bandwidths of bands
		c4 = pyfits.Column(name='TOTAL_BANDWIDTH', format='1E', unit='HZ', 
						array=numpy.array([self.freq[0].totalBW], dtype=numpy.float32))
		# Sideband flag
		c5 = pyfits.Column(name='SIDEBAND', format='1J', 
						array=numpy.array([self.freq[0].sideBand], dtype=numpy.int32))
		# Baseband channel
		c6 = pyfits.Column(name='BB_CHAN', format='1J', 
						array=numpy.array([self.freq[0].baseBand], dtype=numpy.int32))
						
		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6])
		
		# Create the table and header
		fq = pyfits.new_table(colDefs)
		self._addCommonKeywords(fq.header, 'FREQUENCY', 1)
		
		# Add the table to the file
		fq.name = 'FREQUENCY'
		self.FITS.append(fq)
		self.FITS.flush()
		
	def _writeAntenna(self):
		"""
		Define the Antenna table (group 2, table 1).
		"""
		
		# Central time of period covered by record in days
		c1 = pyfits.Column(name='TIME', unit='DAYS', format='1D', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.float64))
		# Duration of period covered by record in days
		c2 = pyfits.Column(name='TIME_INTERVAL', unit='DAYS', format='1E', 
						array=(2*numpy.ones((self.nAnt,), dtype=numpy.float32)))
		# Antenna name
		c3 = pyfits.Column(name='ANNAME', format='A8', 
						array=self.FITS['ARRAY_GEOMETRY'].data.field('ANNAME'))
		# Antenna number
		c4 = pyfits.Column(name='ANTENNA_NO', format='1J', 
						array=self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA'))
		# Array number
		c5 = pyfits.Column(name='ARRAY', format='1J', 
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Frequency setup number
		c6 = pyfits.Column(name='FREQID', format='1J', 
						array=(numpy.zeros((self.nAnt,), dtype=numpy.int32) + self.freq[0].id))
		# Number of digitizer levels
		c7 = pyfits.Column(name='NO_LEVELS', format='1J', 
						array=numpy.array([ant.levels for ant in self.array[0]['ants']]))
		# Feed A polarization label
		c8 = pyfits.Column(name='POLTYA', format='A1', 
						array=numpy.array([ant.polA['Type'] for ant in self.array[0]['ants']]))
		# Feed A orientation in degrees
		c9 = pyfits.Column(name='POLAA', format='1E', 
						array=numpy.array([ant.polA['Angle'] for ant in self.array[0]['ants']], dtype=numpy.float32))
		# Feed A polarization parameters
		c10 = pyfits.Column(name='POLCALA', format='2E', 
						array=numpy.array([ant.polA['Cal'] for ant in self.array[0]['ants']], dtype=numpy.float32))
		# Feed B polarization label
		c11 = pyfits.Column(name='POLTYB', format='A1', 
						array=numpy.array([ant.polB['Type'] for ant in self.array[0]['ants']]))
		# Feed B orientation in degrees
		c12 = pyfits.Column(name='POLAB', format='1E', 
						array=numpy.array([ant.polB['Angle'] for ant in self.array[0]['ants']], dtype=numpy.float32))
		# Feed B polarization parameters
		c13 = pyfits.Column(name='POLCALB', format='2E', 
						array=numpy.array([ant.polB['Cal'] for ant in self.array[0]['ants']], dtype=numpy.float32))
						
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13])
							
		# Create the Antenna table and update it's header
		an = pyfits.new_table(colDefs)
		self._addCommonKeywords(an.header, 'ANTENNA', 1)
		
		an.header['NOPCAL'] = (2, 'number of polarization parameters')
		an.header['POLTYPE'] = ('X-Y LIN', 'polarization parameterization')
		
		an.name = 'ANTENNA'
		self.FITS.append(an)
		self.FITS.flush()
		
	def _writeBandpass(self):
		"""
		Define the Bandpass table (group 2, table 3).
		"""
		
		# Central time of period covered by record in days
		c1 = pyfits.Column(name='TIME', unit='DAYS', format='1D', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.float64))
		# Duration of period covered by record in days
		c2 = pyfits.Column(name='TIME_INTERVAL', unit='DAYS', format='1E',
						array=(2*numpy.ones((self.nAnt,), dtype=numpy.float32)))
		# Source ID
		c3 = pyfits.Column(name='SOURCE_ID', format='1J', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.int32))
		# Antenna number
		c4 = pyfits.Column(name='ANTENNA_NO', format='1J', 
						array=self.FITS['ANTENNA'].data.field('ANTENNA_NO'))
		# Array number
		c5 = pyfits.Column(name='ARRAY', format='1J', 
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Frequency setup number
		c6 = pyfits.Column(name='FREQID', format='1J',
						array=(numpy.zeros((self.nAnt,), dtype=numpy.int32) + self.freq[0].id))
		# Bandwidth in Hz
		c7 = pyfits.Column(name='BANDWIDTH', unit='HZ', format='1E',
						array=(numpy.zeros((self.nAnt,), dtype=numpy.float32)+self.freq[0].totalBW))
		# Band frequency in Hz
		c8 = pyfits.Column(name='BAND_FREQ', unit='HZ', format='1D',
						array=(numpy.zeros((self.nAnt,), dtype=numpy.float64)+self.freq[0].bandFreq))
		# Reference antenna number (pol. 1)
		c9 = pyfits.Column(name='REFANT_1', format='1J',
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Real part of the bandpass (pol. 1)
		c10 = pyfits.Column(name='BREAL_1', format='%dE' % self.nChan,
						array=numpy.ones((self.nAnt,self.nChan), dtype=numpy.float32))
		# Imaginary part of the bandpass (pol. 1)
		c11 = pyfits.Column(name='BIMAG_1', format='%dE' % self.nChan,
						array=numpy.zeros((self.nAnt,self.nChan), dtype=numpy.float32))
		# Reference antenna number (pol. 2)
		c12 = pyfits.Column(name='REFANT_2', format='1J',
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Real part of the bandpass (pol. 2)
		c13 = pyfits.Column(name='BREAL_2', format='%dE' % self.nChan,
						array=numpy.ones((self.nAnt,self.nChan), dtype=numpy.float32))
		# Imaginary part of the bandpass (pol. 2)
		c14 = pyfits.Column(name='BIMAG_2', format='%dE' % self.nChan,
						array=numpy.zeros((self.nAnt,self.nChan), dtype=numpy.float32))
						
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13, c14])
							
		# Create the Bandpass table and update its header
		bp = pyfits.new_table(colDefs)
		self._addCommonKeywords(bp.header, 'BANDPASS', 1)
		
		bp.header['NO_ANT'] = self.nAnt
		bp.header['NO_POL'] = 2
		bp.header['NO_BACH'] = self.nChan
		bp.header['STRT_CHN'] = self.refPix
		
		bp.name = 'BANDPASS'
		self.FITS.append(bp)
		self.FITS.flush()
		
	def _writeSource(self):
		"""
		Define the Source table (group 1, table 2).
		"""
		
		(arrPos, ag) = self.readArrayGeometry()
		ids = ag.keys()
		
		obs = ephem.Observer()
		obs.lat = arrPos.lat * numpy.pi/180
		obs.lon = arrPos.lng * numpy.pi/180
		obs.elev = arrPos.elv * numpy.pi/180
		obs.pressure = 0
		
		nameList = []
		raList = []
		decList = []
		raPoList = []
		decPoList = []
		sourceID = 0
		lastSource = None
		for dataSet in self.data:
			if dataSet.pol == self.stokes[0]:
				utc = astro.taimjd_to_utcjd(dataSet.obsTime)
				date = astro.get_date(utc)
				date.hours = 0
				date.minutes = 0
				date.seconds = 0
				utc0 = date.to_jd()
				
				obs.date = utc - astro.DJD_OFFSET
				
				if dataSet.source != lastSource:
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
						lastSource = dataSet.source
						
						name = dataSet.source.name
						equ = astro.equ_posn(dataSet.source.ra*180/numpy.pi, dataSet.source.dec*180/numpy.pi)
						equPo = astro.equ_posn(dataSet.source.a_ra*180/numpy.pi, dataSet.source.a_dec*180/numpy.pi)
						
					# current apparent zenith equatorial coordinates
					raList.append(equ.ra)
					decList.append(equ.dec)
					
					# J2000 zenith equatorial coordinates
					raPoList.append(equPo.ra)
					decPoList.append(equPo.dec)
					
					# name
					nameList.append(name)
					
		nSource = len(nameList)
		
		# Save these for later since we might need them
		self._sourceTable = nameList
		
		# Source ID number
		c1 = pyfits.Column(name='SOURCE_ID', format='1J', 
						array=numpy.arange(1, nSource+1, dtype=numpy.int32))
		# Source name
		c2 = pyfits.Column(name='SOURCE', format='A16', 
						array=numpy.array(nameList))
		# Source qualifier
		c3 = pyfits.Column(name='QUAL', format='1J', 
						array=numpy.zeros((nSource,), dtype=numpy.int32))
		# Calibrator code
		c4 = pyfits.Column(name='CALCODE', format='A4', 
						array=numpy.array(('   ',)).repeat(nSource))
		# Frequency group ID
		c5 = pyfits.Column(name='FREQID', format='1J', 
						array=(numpy.zeros((nSource,), dtype=numpy.int32)+self.freq[0].id))
		# Stokes I flux density in Jy
		c6 = pyfits.Column(name='IFLUX', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
		# Stokes I flux density in Jy
		c7 = pyfits.Column(name='QFLUX', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
		# Stokes I flux density in Jy
		c8 = pyfits.Column(name='UFLUX', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
		# Stokes I flux density in Jy
		c9 = pyfits.Column(name='VFLUX', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
		# Spectral index
		c10 = pyfits.Column(name='ALPHA', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
		# Frequency offset in Hz
		c11 = pyfits.Column(name='FREQOFF', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
		# Mean equinox and epoch
		c12 = pyfits.Column(name='EQUINOX', format='A8',
						array=numpy.array(('J2000',)).repeat(nSource))
		c13 = pyfits.Column(name='EPOCH', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64) + 2000.0)
		# Apparent right ascension in degrees
		c14 = pyfits.Column(name='RAAPP', format='1D', 
						array=numpy.array(raList))
		# Apparent declination in degrees
		c15 = pyfits.Column(name='DECAPP', format='1D', 
						array=numpy.array(decList))
		# Right ascension at mean equinox in degrees
		c16 = pyfits.Column(name='RAEPO', format='1D', 
						array=numpy.array(raPoList))
		# Declination at mean equinox in degrees
		c17 = pyfits.Column(name='DECEPO', format='1D', 
						array=numpy.array(decPoList))
		# Systemic velocity in m/s
		c18 = pyfits.Column(name='SYSVEL', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64))
		# Velocity type
		c19 = pyfits.Column(name='VELTYP', format='A8', 
						array=numpy.array(('GEOCENTR',)).repeat(nSource))
		# Velocity definition
		c20 = pyfits.Column(name='VELDEF', format='A8', 
						array=numpy.array(('OPTICAL',)).repeat(nSource))
		# Line rest frequency in Hz
		c21 = pyfits.Column(name='RESTFREQ', format='1D', 
						array=(numpy.zeros((nSource,), dtype=numpy.float64) + self.refVal))
		# Proper motion in RA in degrees/day
		c22 = pyfits.Column(name='PMRA', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64))
		# Proper motion in Dec in degrees/day
		c23 = pyfits.Column(name='PMDEC', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64))
		# Parallax of source in arc sec.
		c24 = pyfits.Column(name='PARALLAX', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))
						
		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, 
							c21, c22, c23, c24])
							
		# Create the Source table and update its header
		sr = pyfits.new_table(colDefs)
		self._addCommonKeywords(sr.header, 'SOURCE', 1)
		
		sr.name = 'SOURCE'
		self.FITS.append(sr)
		self.FITS.flush()
		
	def _writeData(self):
		"""
		Define the UV_Data table (group 3, table 1).
		"""
		
		(arrPos, ag) = self.readArrayGeometry()
		(mapper, inverseMapper) = self.readArrayMapper()
		ids = ag.keys()
		
		obs = ephem.Observer()
		obs.lat = arrPos.lat * numpy.pi/180
		obs.lon = arrPos.lng * numpy.pi/180
		obs.elev = arrPos.elv * numpy.pi/180
		obs.pressure = 0
		
		# Retrieve the original list of Antenna objects and convert them to
		# a dictionary index by the stand ID number
		inputAnts = {}
		for ant in self.array[0]['inputAnts']:
			inputAnts[ant.stand.id] = ant
			
		mList = []
		uList = []
		vList = []
		wList = []
		timeList = []
		dateList = []
		intTimeList = []
		blineList = []
		rawList = []
		nameList = []
		sourceList = []
		for dataSet in self.data:
			if dataSet.pol == self.stokes[0]:
				utc = astro.taimjd_to_utcjd(dataSet.obsTime)
				date = astro.get_date(utc)
				date.hours = 0
				date.minutes = 0
				date.seconds = 0
				utc0 = date.to_jd()
				
				obs.date = utc - astro.DJD_OFFSET
				
				if dataSet.source == 'z':
					## Zenith pointings
					equ = astro.equ_posn( obs.sidereal_time()*180/numpy.pi, obs.lat*180/numpy.pi )
					
					# format 'source' name based on local sidereal time
					raHms = astro.deg_to_hms(equ.ra)
					(tsecs, secs) = math.modf(raHms.seconds)
					name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
				else:
					## Real-live sources (ephem.Body instances)
					name = dataSet.source.name
					
				sourceID = self._sourceTable.index(name) + 1
				
				# Compute the uvw coordinates of all baselines
				if inverseMapper is not None:
					standIDs = []
					for stand in self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA'):
						standIDs.append(inverseMapper[stand])
					standIDs = numpy.array(standIDs)
				else:
					standIDs = self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA')
				antennaStands = []
				for standID in standIDs:
					antennaStands.append( inputAnts[standID] )
				
				uvwBaselines = numpy.array([mergeBaseline(a1.stand.id, a2.stand.id, shift=16) for a1,a2 in uvUtils.getBaselines(antennaStands)])
				if dataSet.source == 'z':
					uvwCoords = uvUtils.computeUVW(antennaStands, HA=0.0, dec=equ.dec, site=obs, freq=self.refVal)
				else:
					HA = obs.sidereal_time() - dataSet.source.ra
					dec = dataSet.source.dec * 180/numpy.pi
					
					uvwCoords = uvUtils.computeUVW(antennaStands, HA=HA, dec=dec, site=obs, freq=self.refVal)
				uvwCoords *= 1.0 / self.refVal
				
				tempMList = {}
				for stokes in self.stokes:
					tempMList[stokes] = {}
					
			def _baselineSortMapped(x, y, mapper=mapper, shift=8):
				x1, x2 = splitBaseline(x, shift=16)
				y1, y2 = splitBaseline(y, shift=16)
				if mapper is not None:
					x1 = mapper[x1]
					x2 = mapper[x2]
					y1 = mapper[y1]
					y2 = mapper[y2]
				# Reconstruct the baseline in FITS IDI format
				xPrime = mergeBaseline(x1, x2, shift=shift)
				yPrime = mergeBaseline(y1, y2, shift=shift)
				
				if xPrime > yPrime:
					return 1
				elif xPrime < yPrime:
					return -1
				else:
					return 0
				
			# Loop over the data store in the dataDict and extract each baseline
			baselines = list(dataSet.dataDict.keys())
			baselines.sort(cmp=_baselineSortMapped)
			for baseline in baselines: 
				# validate baseline antenna ID's
				stand1, stand2 = splitBaseline(baseline, shift=16)
				if mapper is not None:
					stand1 = mapper[stand1]
					stand2 = mapper[stand2]
				# Reconstruct the baseline in FITS IDI format
				baselineMapped = mergeBaseline(stand1, stand2, shift=8)
				
				if (stand1 not in ids) or (stand2 not in ids):
					raise ValueError("baseline 0x%04x (%i to %i) contains unknown antenna IDs" % (baselineMapped, stand1, stand2))
					
				# validate data matrix shape
				visData = dataSet.dataDict[baseline]
				if (len(visData) != self.nChan):
					raise ValueError("data for baseline 0x%04x is the wrong size %s (!= %s)" % (baselineMapped, len(visData), self.nChan))
					
				# Calculate the UVW coordinates for antenna pair.  Catch auto-
				# correlations and set their UVW's to zeros.
				try:
					which = (numpy.where( uvwBaselines == baseline ))[0][0]
					uvw = numpy.squeeze(uvwCoords[which,:,0])
				except IndexError:
					uvw = numpy.zeros((3,))
					
				# Load the data into a matrix that splits the real and imaginary parts out
				matrix = numpy.zeros((2*self.nChan,), dtype=numpy.float32)
				matrix[0::2] = visData.real
				matrix[1::2] = visData.imag
				tempMList[dataSet.pol][baseline] = matrix.ravel()
				
				if dataSet.pol == self.stokes[0]:
					blineList.append(baselineMapped)
					rawList.append(baseline)
					uList.append(uvw[0])
					vList.append(uvw[1])
					wList.append(uvw[2])
					dateList.append(utc0)
					timeList.append(utc - utc0) 
					intTimeList.append(dataSet.intTime)
					sourceList.append(sourceID)
					nameList.append(name)
					
			if dataSet.pol == self.stokes[-1]:
				for bl in rawList:
					matrix = numpy.zeros((2*self.nStokes*self.nChan,), dtype=numpy.float32)
					for p in xrange(self.nStokes):
						try:
							matrix[(0+2*p)::(2*self.nStokes)] = tempMList[self.stokes[p]][bl][0::2]
							matrix[(1+2*p)::(2*self.nStokes)] = tempMList[self.stokes[p]][bl][1::2]
						except KeyError:
							stand1, stand2 = splitBaseline(bl, shift=16)
							newBL = mergeBaseline(stand2, stand1, shift=16)
							print 'WARNING: Keyerror', bl, bl in tempMList[self.stokes[p]], stand1, stand2, newBL, newBL in tempMList[self.stokes[p]]
							
					mList.append(matrix.ravel())
					
				rawList = []
		nBaseline = len(blineList)
		nSource = len(nameList)
		
		# Visibility Data
		c1 = pyfits.Column(name='FLUX', format='%iE' % (2*self.nStokes*self.nChan), unit='UNCALIB', 
						array=numpy.array(mList, dtype=numpy.float32))
		# Baseline number (first*256+second)
		c2 = pyfits.Column(name='BASELINE', format='1J', 
						array=numpy.array(blineList))
		# Julian date at 0h
		c3 = pyfits.Column(name='DATE', format='1D', unit='DAYS',
						array = numpy.array(dateList))
		# Time elapsed since 0h
		c4 = pyfits.Column(name='TIME', format='1D', unit = 'DAYS', 
						array = numpy.array(timeList))
		# Integration time (seconds)
		c5 = pyfits.Column(name='INTTIM', format='1D', unit='SECONDS', 
						array=numpy.array(intTimeList, dtype=numpy.float32))
		# U coordinate (light seconds)
		c6 = pyfits.Column(name='UU', format='1E', unit='SECONDS', 
						array=numpy.array(uList, dtype=numpy.float32))
		# V coordinate (light seconds)
		c7 = pyfits.Column(name='VV', format='1E', 
						array=numpy.array(vList, dtype=numpy.float32))
		# W coordinate (light seconds)
		c8 = pyfits.Column(name='WW', format='1E', 
						array=numpy.array(wList, dtype=numpy.float32))
		# Source ID number
		c9 = pyfits.Column(name='SOURCE', format='1J', 
						array=numpy.array(sourceList))
		# Frequency setup number
		c10 = pyfits.Column(name='FREQID', format='1J', 
						array=(numpy.zeros((nBaseline,), dtype=numpy.int32) + self.freq[0].id))
		# Filter number
		c11 = pyfits.Column(name='FILTER', format='1J', 
						array=numpy.zeros((nBaseline,), dtype=numpy.int32))
		# Gate ID number
		c12 = pyfits.Column(name='GATEID', format='1J', 
						array=numpy.zeros((nBaseline,), dtype=numpy.int32))
		# Weights
		c13 = pyfits.Column(name='WEIGHT', format='%iE' % self.nChan, 
						array=numpy.ones((nBaseline, self.nChan), dtype=numpy.float32))
						
		colDefs = pyfits.ColDefs([c6, c7, c8, c3, c4, c2, c11, c9, c10, c5, 
						c13, c12, c1])
						
		# Create the UV Data table and update its header
		uv = pyfits.new_table(colDefs)
		self._addCommonKeywords(uv.header, 'UV_DATA', 1)
		
		uv.header['NMATRIX'] = (1, 'number of UV data matricies')
		uv.header['MAXIS'] = (6, 'number of UV data matrix axes')
		uv.header['TMATX13'] = (True, 'axis 13 contains UV matrix')
		
		uv.header['MAXIS1'] = (2, 'number of pixels in COMPLEX axis')
		uv.header['CTYPE1'] = ('COMPLEX', 'axis 1 is COMPLEX axis')
		uv.header['CDELT1'] = 1.0
		uv.header['CRPIX1'] = 1.0
		uv.header['CRVAL1'] = 1.0
		
		uv.header['MAXIS2'] = (self.nStokes, 'number of pixels in STOKES axis')
		uv.header['CTYPE2'] = ('STOKES', 'axis 2 is STOKES axis (polarization)')
		if self.stokes[0] < 0:
			uv.header['CDELT2'] = -1.0
		else:
			uv.header['CDELT2'] = 1.0
		uv.header['CRPIX2'] = 1.0
		uv.header['CRVAL2'] = float(self.stokes[0])
		
		uv.header['MAXIS3'] = (self.nChan, 'number of pixels in FREQ axis')
		uv.header['CTYPE3'] = ('FREQ', 'axis 3 is FREQ axis (frequency)')
		uv.header['CDELT3'] = self.freq[0].chWidth
		uv.header['CRPIX3'] = self.refPix
		uv.header['CRVAL3'] = self.refVal
		
		uv.header['MAXIS4'] = (1, 'number of pixels in BAND (IF) axis')
		uv.header['CTYPE4'] = ('BAND', 'axis 4 is BAND axis')
		uv.header['CDELT4'] = 1.0
		uv.header['CRPIX4'] = 1.0
		uv.header['CRVAL4'] = 1.0
		
		uv.header['MAXIS5'] = (1, 'number of pixels in RA axis')
		uv.header['CTYPE5'] = ('RA', 'axis 5 is RA axis (position of phase center)')
		uv.header['CDELT5'] = 0.0
		uv.header['CRPIX5'] = 1.0
		uv.header['CRVAL5'] = 0.0
		
		uv.header['MAXIS6'] = (1, 'number of pixels in DEC axis')
		uv.header['CTYPE6'] = ('DEC', 'axis 6 is DEC axis (position of phase center)')
		uv.header['CDELT6'] = 0.0
		uv.header['CRPIX6'] = 1.0
		uv.header['CRVAL6'] = 0.0
		
		uv.header['TELESCOP'] = self.siteName
		uv.header['OBSERVER'] = 'ZASKY'
		uv.header['SORT'] = ('TB', 'data is sorted in [time,baseline] order')
		
		uv.header['VISSCALE'] = (1.0, 'UV data scale factor')
		
		uv.name = 'UV_DATA'
		self.FITS.append(uv)
		self.FITS.flush()
		
	def _writeMapper(self):
		"""
		Write a fits table that contains information about mapping stations 
		numbers to actual antenna numbers.  This information can be backed out of
		the names, but this makes the extraction more programmatic.
		"""
		
		c1 = pyfits.Column(name='ANNAME', format='A8', 
						array=numpy.array([ant.getName() for ant in self.array[0]['ants']]))
		c2 = pyfits.Column(name='NOSTA', format='1J', 
						array=numpy.array([self.array[0]['mapper'][ant.id] for ant in self.array[0]['ants']]))
		c3 = pyfits.Column(name='NOACT', format='1J', 
						array=numpy.array([ant.id for ant in self.array[0]['ants']]))
						
		colDefs = pyfits.ColDefs([c1, c2, c3])
		
		# Create the ID mapping table and update its header
		nsm = pyfits.new_table(colDefs)
		self._addCommonKeywords(nsm.header, 'NOSTA_MAPPER', 1)
		
		nsm.name = 'NOSTA_MAPPER'
		self.FITS.append(nsm)
		self.FITS.flush()
		
	def readArrayGeometry(self):
		"""
		Return a tuple with the array geodetic position and the local 
		positions for all antennas defined in the ARRAY_GEOMETRY table.
		"""
		
		try:
			ag = self.FITS['ARRAY_GEOMETRY']
		except IndexError:
			raise RuntimeError("File does not have an 'ARRAY_GEOMTRY' table.")
			
		# Array position
		arrayGeo = astro.rect_posn(ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ'])
		arrayGeo = astro.get_geo_from_rect(arrayGeo)
		
		# Antenna positions
		antennaGeo = {}
		antenna = ag.data
		antennaID = antenna.field('NOSTA')
		antennaPos = iter(antenna.field('STABXYZ'))
		for id in antennaID:
			antennaGeo[id] = antennaPos.next()
			
		# Return
		return (arrayGeo, antennaGeo)
		
	def readArrayMapper(self):
		"""
		Return a tuple with the array NOSTA mapper and inverse mapper (both
		dictionaries.  If the stand IDs have not been mapped, return None for
		both.
		"""
		
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


class AIPS(IDI):
	"""
	Sub-class of the FITS IDI writer for making files that *should* work 
	with AIPS nicely.  AIPS imposes a limit on antenna number of two digits
	(1-99).  This sub-class overwrite the setGeometry() function with one that
	enforces the two digit limit and maps accordingly.  It also sets the FITS
	`LWATYPE` keyword in the primary HDU to a value of `IDI-AIPS-ZA` to 
	distinguish files written by this writer from the standard IDI writer.
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
			return "L%03i" % self.id
			
	def setGeometry(self, site, antennas, bits=8):
		"""
		Given a station and an array of stands, set the relevant common observation
		parameters and add entries to the self.array list.
		
		.. versionchanged:: 0.4.0
			Switched over to passing in Antenna instances generated by the
			:mod:`lsl.common.stations` module instead of a list of stand ID
			numbers.
		"""
		
		# Make sure that we have been passed 99 or fewer stands
		if len(antennas) > 99:
			raise RuntimeError("FITS IDI for AIPS supports up to 99 antennas only, given %i" % len(antennas))
			
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
			
		# Create the stand mapper to deal with the fact that stands range from 
		# 1 to 258, not 1 to 99
		mapper = {}
		if stands.max() > 99:
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
			print "FITS IDI: stand ID mapping enabled"
			for key, value in mapper.iteritems():
				print "FITS IDI:  stand #%i -> mapped #%i" % (key, value)
				
		self.nAnt = len(ants)
		self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'enableMapper': enableMapper, 'inputAnts': antennas} )
		
	def _writePrimary(self):
		"""
		Write the primary HDU to file.
		"""
		
		primary = pyfits.PrimaryHDU()
		
		primary.header['NAXIS'] = (0, 'indicates IDI file')
		primary.header['EXTEND'] = (True, 'indicates IDI file')
		primary.header['GROUPS'] = (True, 'indicates IDI file')
		primary.header['GCOUNT'] = 0
		primary.header['PCOUNT'] = 0
		primary.header['OBJECT'] = 'BINARYTB'
		primary.header['TELESCOP'] = self.siteName
		primary.header['INSTRUME'] = self.siteName
		primary.header['OBSERVER'] = ('ZASKY', 'zenith all-sky image')
		primary.header['ORIGIN'] = 'LSL'
		primary.header['CORRELAT'] = ('LWASWC', 'Correlator used')
		primary.header['FXCORVER'] = ('1', 'Correlator version')
		primary.header['LWATYPE'] = ('IDI-AIPS-ZA', 'LWA FITS file type')
		primary.header['LWAMAJV'] = (IDIVersion[0], 'LWA FITS file format major version')
		primary.header['LWAMINV'] = (IDIVersion[1], 'LWA FITS file format minor version')
		primary.header['DATE-OBS'] = (self.refTime, 'IDI file data collection date')
		ts = str(astro.get_date_from_sys())
		primary.header['DATE-MAP'] = (ts.split()[0], 'IDI file creation date')
		
		self.FITS.append(primary)
		self.FITS.flush()


class ExtendedIDI(IDI):
	"""
	Sub-class of the FITS IDI writer for making files that support up to
	65,535 antennas.  This is done by changing the packing of baselines 
	stored in the UVDATA table.  The new packing for baseline (i,j) is:
	  (i << 16) & (j)
	It also sets the FITS `LWATYPE` keyword in the primary HDU to a value 
	of `IDI-EXTENDED-ZA` to distinguish files written by this writer 
	from the standard IDI writer.
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
			return "LWA%05i" % self.id
			
	def setGeometry(self, site, antennas, bits=8):
		"""
		Given a station and an array of stands, set the relevant common observation
		parameters and add entries to the self.array list.
		
		.. versionchanged:: 0.4.0
			Switched over to passing in Antenna instances generated by the
			:mod:`lsl.common.stations` module instead of a list of stand ID
			numbers.
		"""
		
		# Make sure that we have been passed 99 or fewer stands
		if len(antennas) > 65535:
			raise RuntimeError("Extended FITS IDI supports up to 65535 antennas only, given %i" % len(antennas))
			
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
			
		# 65,535 antennas should be enough for anybody
		mapper = {}
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
			print "FITS IDI: stand ID mapping enabled"
			for key, value in mapper.iteritems():
				print "FITS IDI:  stand #%i -> mapped #%i" % (key, value)
				
		self.nAnt = len(ants)
		self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'enableMapper': enableMapper, 'inputAnts': antennas} )
		
	def _writePrimary(self):
		"""
		Write the primary HDU to file.
		"""
		
		primary = pyfits.PrimaryHDU()
		
		primary.header['NAXIS'] = (0, 'indicates IDI file')
		primary.header['EXTEND'] = (True, 'indicates IDI file')
		primary.header['GROUPS'] = (True, 'indicates IDI file')
		primary.header['GCOUNT'] = 0
		primary.header['PCOUNT'] = 0
		primary.header['OBJECT'] = 'BINARYTB'
		primary.header['TELESCOP'] = self.siteName
		primary.header['INSTRUME'] = self.siteName
		primary.header['OBSERVER'] = ('ZASKY', 'zenith all-sky image')
		primary.header['ORIGIN'] = 'LSL'
		primary.header['CORRELAT'] = ('LWASWC', 'Correlator used')
		primary.header['FXCORVER'] = ('1', 'Correlator version')
		primary.header['LWATYPE'] = ('IDI-EXTENDED-ZA', 'LWA FITS file type')
		primary.header['LWAMAJV'] = (IDIVersion[0], 'LWA FITS file format major version')
		primary.header['LWAMINV'] = (IDIVersion[1], 'LWA FITS file format minor version')
		primary.header['DATE-OBS'] = (self.refTime, 'IDI file data collection date')
		ts = str(astro.get_date_from_sys())
		primary.header['DATE-MAP'] = (ts.split()[0], 'IDI file creation date')
		
		self.FITS.append(primary)
		self.FITS.flush()
		
	def _writeData(self):
		"""
		Define the UV_Data table (group 3, table 1).
		"""
		
		(arrPos, ag) = self.readArrayGeometry()
		(mapper, inverseMapper) = self.readArrayMapper()
		ids = ag.keys()
		
		obs = ephem.Observer()
		obs.lat = arrPos.lat * numpy.pi/180
		obs.lon = arrPos.lng * numpy.pi/180
		obs.elev = arrPos.elv * numpy.pi/180
		obs.pressure = 0
		
		# Retrieve the original list of Antenna objects and convert them to
		# a dictionary index by the stand ID number
		inputAnts = {}
		for ant in self.array[0]['inputAnts']:
			inputAnts[ant.stand.id] = ant
			
		mList = []
		uList = []
		vList = []
		wList = []
		timeList = []
		dateList = []
		intTimeList = []
		blineList = []
		rawList = []
		nameList = []
		sourceList = []
		for dataSet in self.data:
			if dataSet.pol == self.stokes[0]:
				utc = astro.taimjd_to_utcjd(dataSet.obsTime)
				date = astro.get_date(utc)
				date.hours = 0
				date.minutes = 0
				date.seconds = 0
				utc0 = date.to_jd()
				
				obs.date = utc - astro.DJD_OFFSET
				
				if dataSet.source == 'z':
					## Zenith pointings
					equ = astro.equ_posn( obs.sidereal_time()*180/numpy.pi, obs.lat*180/numpy.pi )
					
					# format 'source' name based on local sidereal time
					raHms = astro.deg_to_hms(equ.ra)
					(tsecs, secs) = math.modf(raHms.seconds)
					name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
				else:
					## Real-live sources (ephem.Body instances)
					name = dataSet.source.name
					
				sourceID = self._sourceTable.index(name) + 1
				
				# Compute the uvw coordinates of all baselines
				if inverseMapper is not None:
					standIDs = []
					for stand in self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA'):
						standIDs.append(inverseMapper[stand])
					standIDs = numpy.array(standIDs)
				else:
					standIDs = self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA')
				antennaStands = []
				for standID in standIDs:
					antennaStands.append( inputAnts[standID] )
					
				uvwBaselines = numpy.array([mergeBaseline(a1.stand.id, a2.stand.id, shift=16) for a1,a2 in uvUtils.getBaselines(antennaStands)])
				if dataSet.source == 'z':
					uvwCoords = uvUtils.computeUVW(antennaStands, HA=0.0, dec=equ.dec, site=obs, freq=self.refVal)
				else:
					HA = obs.sidereal_time() - dataSet.source.ra
					dec = dataSet.source.dec * 180/numpy.pi
					
					uvwCoords = uvUtils.computeUVW(antennaStands, HA=HA, dec=dec, site=obs, freq=self.refVal)
				uvwCoords *= 1.0 / self.refVal
				
				tempMList = {}
				for stokes in self.stokes:
					tempMList[stokes] = {}
					
			def _baselineSortMapped(x, y, mapper=mapper, shift=16):
				x1, x2 = splitBaseline(x, shift=16)
				y1, y2 = splitBaseline(y, shift=16)
				if mapper is not None:
					x1 = mapper[x1]
					x2 = mapper[x2]
					y1 = mapper[y1]
					y2 = mapper[y2]
				# Reconstruct the baseline in FITS IDI format
				xPrime = mergeBaseline(x1, x2, shift=shift)
				yPrime = mergeBaseline(y1, y2, shift=shift)
				
				if xPrime > yPrime:
					return 1
				elif xPrime < yPrime:
					return -1
				else:
					return 0
				
			# Loop over the data store in the dataDict and extract each baseline
			baselines = list(dataSet.dataDict.keys())
			baselines.sort(cmp=_baselineSortMapped)
			for baseline in baselines: 
				# validate baseline antenna ID's
				stand1, stand2 = splitBaseline(baseline, shift=16)
				if mapper is not None:
					stand1 = mapper[stand1]
					stand2 = mapper[stand2]
				# Reconstruct the baseline in FITS IDI Extended format
				baselineMapped = mergeBaseline(stand1, stand2, shift=16)
				
				if (stand1 not in ids) or (stand2 not in ids):
					raise ValueError("baseline 0x%04x (%i to %i) contains unknown antenna IDs" % (baselineMapped, stand1, stand2))
					
				# validate data matrix shape
				visData = dataSet.dataDict[baseline]
				if (len(visData) != self.nChan):
					raise ValueError("data for baseline 0x%04x is the wrong size %s (!= %s)" % (baselineMapped, len(visData), self.nChan))
					
				# Calculate the UVW coordinates for antenna pair.  Catch auto-
				# correlations and set their UVW's to zeros.
				try:
					which = (numpy.where( uvwBaselines == baseline ))[0][0]
					uvw = numpy.squeeze(uvwCoords[which,:,0])
				except IndexError:
					uvw = numpy.zeros((3,))
					
				# Load the data into a matrix that splits the real and imaginary parts out
				matrix = numpy.zeros((2*self.nChan,), dtype=numpy.float32)
				matrix[0::2] = visData.real
				matrix[1::2] = visData.imag
				tempMList[dataSet.pol][baseline] = matrix.ravel()
				
				if dataSet.pol == self.stokes[0]:
					blineList.append(baselineMapped)
					rawList.append(baseline)
					uList.append(uvw[0])
					vList.append(uvw[1])
					wList.append(uvw[2])
					dateList.append(utc0)
					timeList.append(utc - utc0) 
					intTimeList.append(dataSet.intTime)
					sourceList.append(sourceID)
					nameList.append(name)
					
			if dataSet.pol == self.stokes[-1]:
				for bl in rawList:
					matrix = numpy.zeros((2*self.nStokes*self.nChan,), dtype=numpy.float32)
					for p in xrange(self.nStokes):
						try:
							matrix[(0+2*p)::(2*self.nStokes)] = tempMList[self.stokes[p]][bl][0::2]
							matrix[(1+2*p)::(2*self.nStokes)] = tempMList[self.stokes[p]][bl][1::2]
						except KeyError:
							stand1, stand2 = splitBaseline(bl, shift=16)
							newBL = mergeBaseline(stand2, stand1, shift=16)
							print 'WARNING: Keyerror', bl, bl in tempMList[self.stokes[p]], stand1, stand2, newBL, newBL in tempMList[self.stokes[p]]
							
					mList.append(matrix.ravel())
					
				rawList = []
		nBaseline = len(blineList)
		nSource = len(nameList)
		
		# Visibility Data
		c1 = pyfits.Column(name='FLUX', format='%iE' % (2*self.nStokes*self.nChan), unit='UNCALIB', 
						array=numpy.array(mList, dtype=numpy.float32))
		# Baseline number (first*256+second)
		c2 = pyfits.Column(name='BASELINE', format='1J', 
						array=numpy.array(blineList))
		# Julian date at 0h
		c3 = pyfits.Column(name='DATE', format='1D', unit='DAYS',
						array = numpy.array(dateList))
		# Time elapsed since 0h
		c4 = pyfits.Column(name='TIME', format='1D', unit = 'DAYS', 
						array = numpy.array(timeList))
		# Integration time (seconds)
		c5 = pyfits.Column(name='INTTIM', format='1D', unit='SECONDS', 
						array=numpy.array(intTimeList, dtype=numpy.float32))
		# U coordinate (light seconds)
		c6 = pyfits.Column(name='UU', format='1E', unit='SECONDS', 
						array=numpy.array(uList, dtype=numpy.float32))
		# V coordinate (light seconds)
		c7 = pyfits.Column(name='VV', format='1E', 
						array=numpy.array(vList, dtype=numpy.float32))
		# W coordinate (light seconds)
		c8 = pyfits.Column(name='WW', format='1E', 
						array=numpy.array(wList, dtype=numpy.float32))
		# Source ID number
		c9 = pyfits.Column(name='SOURCE', format='1J', 
						array=numpy.array(sourceList))
		# Frequency setup number
		c10 = pyfits.Column(name='FREQID', format='1J', 
						array=(numpy.zeros((nBaseline,), dtype=numpy.int32) + self.freq[0].id))
		# Filter number
		c11 = pyfits.Column(name='FILTER', format='1J', 
						array=numpy.zeros((nBaseline,), dtype=numpy.int32))
		# Gate ID number
		c12 = pyfits.Column(name='GATEID', format='1J', 
						array=numpy.zeros((nBaseline,), dtype=numpy.int32))
		# Weights
		c13 = pyfits.Column(name='WEIGHT', format='%iE' % self.nChan, 
						array=numpy.ones((nBaseline, self.nChan), dtype=numpy.float32))
						
		colDefs = pyfits.ColDefs([c6, c7, c8, c3, c4, c2, c11, c9, c10, c5, 
						c13, c12, c1])
						
		# Create the UV Data table and update its header
		uv = pyfits.new_table(colDefs)
		self._addCommonKeywords(uv.header, 'UV_DATA', 1)
		
		uv.header['NMATRIX'] = (1, 'number of UV data matricies')
		uv.header['MAXIS'] = (6, 'number of UV data matrix axes')
		uv.header['TMATX13'] = (True, 'axis 13 contains UV matrix')
		
		uv.header['MAXIS1'] = (2, 'number of pixels in COMPLEX axis')
		uv.header['CTYPE1'] = ('COMPLEX', 'axis 1 is COMPLEX axis')
		uv.header['CDELT1'] = 1.0
		uv.header['CRPIX1'] = 1.0
		uv.header['CRVAL1'] = 1.0
		
		uv.header['MAXIS2'] = (self.nStokes, 'number of pixels in STOKES axis')
		uv.header['CTYPE2'] = ('STOKES', 'axis 2 is STOKES axis (polarization)')
		if self.stokes[0] < 0:
			uv.header['CDELT2'] = -1.0
		else:
			uv.header['CDELT2'] = 1.0
		uv.header['CRPIX2'] = 1.0
		uv.header['CRVAL2'] = float(self.stokes[0])
		
		uv.header['MAXIS3'] = (self.nChan, 'number of pixels in FREQ axis')
		uv.header['CTYPE3'] = ('FREQ', 'axis 3 is FREQ axis (frequency)')
		uv.header['CDELT3'] = self.freq[0].chWidth
		uv.header['CRPIX3'] = self.refPix
		uv.header['CRVAL3'] = self.refVal
		
		uv.header['MAXIS4'] = (1, 'number of pixels in BAND (IF) axis')
		uv.header['CTYPE4'] = ('BAND', 'axis 4 is BAND axis')
		uv.header['CDELT4'] = 1.0
		uv.header['CRPIX4'] = 1.0
		uv.header['CRVAL4'] = 1.0
		
		uv.header['MAXIS5'] = (1, 'number of pixels in RA axis')
		uv.header['CTYPE5'] = ('RA', 'axis 5 is RA axis (position of phase center)')
		uv.header['CDELT5'] = 0.0
		uv.header['CRPIX5'] = 1.0
		uv.header['CRVAL5'] = 0.0
		
		uv.header['MAXIS6'] = (1, 'number of pixels in DEC axis')
		uv.header['CTYPE6'] = ('DEC', 'axis 6 is DEC axis (position of phase center)')
		uv.header['CDELT6'] = 0.0
		uv.header['CRPIX6'] = 1.0
		uv.header['CRVAL6'] = 0.0
		
		uv.header['TELESCOP'] = self.siteName
		uv.header['OBSERVER'] = 'ZASKY'
		uv.header['SORT'] = ('TB', 'data is sorted in [time,baseline] order')
		
		uv.header['VISSCALE'] = (1.0, 'UV data scale factor')
		
		uv.name = 'UV_DATA'
		self.FITS.append(uv)
		self.FITS.flush()
