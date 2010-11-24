# -*- coding: utf-8 -*-

"""Module for writting correlator output to a FITS IDI file.  The classes and 
functions defined in this module are based heavily off the lwda_fits library.  
This module is still under active development."""

import os
import re
import sys
import math
import numpy
import pyfits
from datetime import datetime

from lsl import astro
from lsl.common import constants
from lsl.common import dp as dp_common
from lsl.common.stations import geo2ecef
from lsl.correlator import uvUtils
from lsl.misc import mathutil
from lsl.common.warns import warnExperimental

__version__ = '0.1'
__revision__ = '$ Revision: 8 $'
__all__ = ['IDI', 'StokesCodes', '__version__', '__revision__', '__all__']


IDIVersion = (2, 0)

StokesCodes = {'I': 1, 'Q': 2, 'U': 3, 'V': 4, 
			'RR': -1, 'LL': -2, 'RL': -3, 'LR': -4, 
			'XX': -5, 'YY': -6, 'XY': -7, 'YX': -8}


class IDI(object):
	class _Antenna(object):
		"""Holds information describing the location and properties of an antenna."""

		def __init__(self, id, x, y, z):
			self.id = id
			self.x = x
			self.y = y
			self.z = z
			self.levels = 1
			self.polA = {'Type': 'X', 'Angle': 0.0, 'Cal': 1.0}
			self.polB = {'Type': 'Y', 'Angle': 0.0, 'Cal': 1.0}

		def getName(self):
			return "LWA%03i" % self.id

	class _Frequency:
		"""Holds information about the frequency setup used in the file."""

		def __init__(self, offset, channelWidth, bandwidth):
			self.id = 1
			self.bandFreq = offset
			self.chWidth = channelWidth
			self.totalBW = bandwidth
			self.sideBand = 1
			self.baseBand = 0

	class _UVData(object):
		"""Represents one UV visibility data set for a given observation time."""
    
		def __init__(self, obsTime, intTime, dataDict):
			self.obsTime = obsTime
			self.intTime = intTime
			self.dataDict = dataDict
		
		def time(self):
			return self.obsTime

	def parseRefTime(self, refTime):
		"""Given a time as either a integer, float, string, or datetime object, 
		convert it to a string in the formation 'YYYY-MM-DDTHH:MM:SS'."""

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
		"""Convert a reference time string to an astro.date object."""

		dateStr = self.refTime.replace('T', '-').replace(':', '-').split('-')
		return astro.date(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), int(dateStr[3]), int(dateStr[4]), float(dateStr[5]))


	def __init__(self, filename, refTime=0.0):
		# File-specific information
		self.filename = filename
		self.hdulist = None

		# Observation-specifc information
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
		"""Given a list of Stokes parameters, update the object's parameters"""

		for pol in polList:
			if type(pol).__name__ == 'str':
				self.stokes.append(StokesCodes[pol.upper()])
			else:
				self.stokes.append(pol)

		self.nStokes = len(self.stokes)

	def setFrequency(self, freq):
		"""Given a numpy array of frequencies, set the relevant common observation
		parameters and add an entry to the self.freq list."""

		self.nChan = len(freq)
		self.refVal = freq[0]
		self.refPix = 1
		self.channelWidth = numpy.abs(freq[1] - freq[0])
		totalWidth = numpy.abs(freq[-1] - freq[0])

		freqSetup = self._Frequency(0.0, self.channelWidth, totalWidth)
		self.freq.append(freqSetup)

	def setGeometry(self, site, stands):
		"""Given a station and an array of stands, set the relevant common observation
		parameters and add entries to the self.array list."""

		# Make sure that we have been passed 255 or fewer stands
		if len(stands) > 255:
			raise RuntimeError("FITS IDI support up to 255 antennas only, given %i" % len(stands))

		arrayX, arrayY, arrayZ = site.getGeocentricLocation()
		xyz = uvUtils.getXYZ(stands)

		# Create the stand mapper to deal with the fact that stands range from 
		# 1 to 258, not 1 to 255
		mapper = {}
		if stands.max() > 255:
			enableMapper = True
		else:
			enableMapper = False

		ants = []
		topo2eci = site.getECEFTransform()
		for i in range(len(stands)):
			eci = numpy.dot(topo2eci, xyz[i,:])
			ants.append( self._Antenna(stands[i], eci[0], eci[1], eci[2]) )
			if enableMapper:
				mapper[stands[i]] = i+1
			else:
				mapper[stands[i]] = stands[i]

		# If the mapper has been enabled, tell the user about it
		if enableMapper:
			print "FITS IDI: stand ID mapping enabled"
			for key, value in mapper.iteritems():
				print "FITS IDI:  stand #%i -> mapped #%i" % (key, value)

		self.nAnt = len(ants)
		self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'enableMapper': enableMapper} )

	def addDataSet(self, obsTime, intTime, baselines, visibilities):
			"""Create a UVData object to store a collection of visibilites."""

			dataDict = {}
			for (stand1,stand2), visData in zip(baselines, visibilities):
				baseline = (stand1 << 16) | stand2
				dataDict[baseline] = visData
				
			self.data.append( self._UVData(obsTime, intTime, dataDict) )

	def write(self):
		"""Fill in the FITS-IDI file will all of the tables in the 
		correct order."""
		
		self.__writePrimary()
		self.__writeGeometry()
		self.__writeFrequency()
		self.__writeAntenna()
		self.__writeBandpass()
		self.__writeSource()
		self.__writeStars()
		self.__writeData()

	def __addCommonKeywords(self, hdr, name, revision):
		"""Added keywords common to all table headers."""

		hdr.update('EXTNAME', name, 'FITS-IDI table name')
		hdr.update('EXTVER', 1, 'table instance number') 
		hdr.update('TABREV', revision, 'table format revision number')
		hdr.update('NO_STKD', 1, 'number of Stokes parameters')
		hdr.update('STK_1', self.stokes[0], 'first Stokes parameter')
		hdr.update('NO_BAND', 1, 'number of frequency bands')
		hdr.update('NO_CHAN', self.nChan, 'number of frequency channels')
		hdr.update('REF_FREQ', self.refVal, 'reference frequency (Hz)')
		hdr.update('CHAN_BW', self.channelWidth, 'channel bandwidth (Hz)')
		hdr.update('REF_PIXL', float(self.refPix), 'reference frequency bin')
	
		date = self.refTime.split('-')
		name = "ZA%s%s%s" % (date[0][2:], date[1], date[2])
		hdr.update('OBSCODE', name, 'zenith all-sky image')
	
		hdr.update('ARRNAM', 'LWA-1')      
		hdr.update('RDATE', self.refTime, 'file data reference date')

	def __makeAppendTable(self, extension, AddRows=1):
		"""Private function to make a temporary table for appending data."""

		nrows = self.hdulist[extension].data.shape[0]
		tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+AddRows)
		for key in list(self.hdulist[extension].header.keys()):
			tempHDU.header.update(key, self.hdulist[extension].header[key])
	
		return tempHDU

	def __applyAppendTable(self, extension, tempHDU):
		"""Private function to replace the given extension with the temporary
		table."""

		self.hdulist[extension] = tempHDU
		self.flush()

	def __writePrimary(self):
		"""Write the primary HDU to file."""

		primary = pyfits.PrimaryHDU()

		primary.header.update('NAXIS', 0, 'indicates IDI file')
		primary.header.update('EXTEND', True, 'indicates IDI file')
		primary.header.update('GROUPS', True, 'indicates IDI file')
		primary.header.update('OBJECT', 'BINARYTB')
		primary.header.update('TELESCOP', 'LWA-1')
		primary.header.update('INSTRUME', 'LWA-1')
		primary.header.update('OBSERVER', 'ZASKY', 'zenith all-sky image')
		primary.header.update('ORIGIN', 'LSL')
		primary.header.update('LWATYPE', 'IDI-ZA', 'LWA FITS file type')
		primary.header.update('LWAMAJV', IDIVersion[0], 'LWA FITS file format major version')
		primary.header.update('LWAMINV', IDIVersion[1], 'LWA FITS file format minor version')
		primary.header.update('DATE-OBS', self.refTime, 'IDI file data collection date')
		ts = str(astro.get_date_from_sys())
		primary.header.update('DATE-MAP', ts.split()[0], 'IDI file creation date')
		
		self.FITS.append(primary)
		self.FITS.flush()

	def __writeGeometry(self):
		"""Define the Array_Geometry table (group 1, table 1)."""

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
		self.__addCommonKeywords(ag.header, 'ARRAY_GEOMETRY', 1)

		ag.header.update('EXTVER', 1, 'array ID')
		ag.header.update('ARRNAM', 'LWA-1')
		ag.header.update('FRAME', 'GEOCENTRIC', 'coordinate system')
		ag.header.update('NUMORB', 0, 'number of orbital parameters')
		ag.header.update('FREQ', self.refVal, 'reference frequency (Hz)')
		ag.header.update('TIMSYS', 'UTC', 'time coordinate system')

		date = self.refTime2AstroDate()
		utc0 = date.to_jd()
		gst0 = astro.get_apparent_sidereal_time(utc0)
		ag.header.update('GSTIA0', gst0 * 15, 'GAST (deg) at RDATE 0 hours')
		
		utc1 = utc0 + 1
		gst1 = astro.get_apparent_sidereal_time(utc1)
		if gst1 < gst0:
			gst1 += 24.0
		ds = gst1 - gst0
		deg = ds * 15.0      
		ag.header.update('DEGPDY', 360.0 + deg, 'rotation rate of the earth (deg/day)')
		
		ag.header.update('UT1UTC', 0.0, 'difference UT1 - UTC for reference date')
		ag.header.update('IATUTC', astro.leap_secs(utc0), 'TAI - UTC for reference date')
		ag.header.update('POLARX', 0.0)
		ag.header.update('POLARY', 0.0)

		ag.header.update('ARRAYX', self.array[0]['center'][0], 'array ECI X coordinate (m)')
		ag.header.update('ARRAYY', self.array[0]['center'][1], 'array ECI Y coordinate (m)')
		ag.header.update('ARRAYZ', self.array[0]['center'][2], 'array ECI Z coordinate (m)')

		ag.header.update('NOSTAMAP', int(self.array[0]['enableMapper']), 'Mapping enabled for stand numbers')

		ag.name = 'ARRAY_GEOMETRY'
		self.FITS.append(ag)
		self.FITS.flush()

		if self.array[0]['enableMapper']:
			self.__writeMapper()

	def __writeFrequency(self):
		"""Define the Frequency table (group 1, table 3)."""

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
		self.__addCommonKeywords(fq.header, 'FREQUENCY', 2)
		
		# Add the table to the file
		fq.name = 'FREQUENCY'
		self.FITS.append(fq)
		self.FITS.flush()

	def __writeAntenna(self):
		"""Define the Antenna table (group 2, table 1)."""

		# Central time of period covered by record in days
		c1 = pyfits.Column(name='TIME', unit='DAYS', format='1D', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.float64))
		# Durration of period covered by record in days
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
		c6 = pyfits.Column(name='FREQID', format='IJ', 
						array=(numpy.zeros((self.nAnt,), dtype=numpy.int32) + self.freq[0].id))
		# Number of digitizer levels
		c7 = pyfits.Column(name='NO_LEVELS', format='1J', 
						array=numpy.array([ant.levels for ant in self.array[0]['ants']]))
		# Feed A polarization label
		c8 = pyfits.Column(name='POLTYA', format='A1', 
						array=numpy.array([ant.polA['Type'] for ant in self.array[0]['ants']]))
		# Feed A orientation in degrees
		c9 = pyfits.Column(name='POLAA', format='1E', 
						array=numpy.array([ant.polA['Angle'] for ant in self.array[0]['ants']], dtype=numpy.int32))
		# Feed A polarization parameters
		c10 = pyfits.Column(name='POLCALA', format='1E', 
						array=numpy.array([ant.polA['Cal'] for ant in self.array[0]['ants']], dtype=numpy.int32))
		# Feed B polarization label
		c11 = pyfits.Column(name='POLTYB', format='A1', 
						array=numpy.array([ant.polB['Type'] for ant in self.array[0]['ants']]))
		# Feed B orientation in degrees
		c12 = pyfits.Column(name='POLAB', format='1E', 
						array=numpy.array([ant.polB['Angle'] for ant in self.array[0]['ants']], dtype=numpy.int32))
		# Feed B polarization parameters
		c13 = pyfits.Column(name='POLCALB', format='1E', 
						array=numpy.array([ant.polB['Cal'] for ant in self.array[0]['ants']], dtype=numpy.int32))

		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13])

		# Create the Antenna table and update it's header
		an = pyfits.new_table(colDefs)
		self.__addCommonKeywords(an.header, 'ANTENNA', 1)

		an.header.update('NOPCAL', 0, 'number of polarization parameters')
		an.header.update('POLTYPE', 'APPROX', 'polarization parameterization')
		
		an.name = 'ANTENNA'
		self.FITS.append(an)
		self.FITS.flush()

	def __writeBandpass(self):
		"""Define the Bandpass table (group 2, table 3)."""

		# Central time of period covered by record in days
		c1 = pyfits.Column(name='TIME', unit='DAYS', format='1D', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.float64))
		# Durration of period covered by record in days
		c2 = pyfits.Column(name='TIME_INTERVAL', unit='DAYS', format='1E',
						array=(2*numpy.ones((self.nAnt,), dtype=numpy.float32)))
		# Source ID
		c3 = pyfits.Column(name='SOURCE_ID', format='1J', 
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Antenna number
		c4 = pyfits.Column(name='ANTENNA_NO', format='1J', 
						array=self.FITS['ANTENNA'].data.field('ANTENNA_NO'))
		# Array number
		c5 = pyfits.Column(name='ARRAY', format='1J', 
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Frequency setup number
		c6 = pyfits.Column(name='FREQID', format='IJ',
						array=(numpy.zeros((self.nAnt,), dtype=numpy.int32) + self.freq[0].id))
		# Bandwidth in Hz
		c7 = pyfits.Column(name='BANDWIDTH', unit='HZ', format='1E',
						array=(numpy.zeros((self.nAnt,), dtype=numpy.float32)+self.freq[0].totalBW))
		# Band frequency in Hz
		c8 = pyfits.Column(name='BAND_FREQ', unit='HZ', format='1D',
						array=(numpy.zeros((self.nAnt,), dtype=numpy.float32)+self.freq[0].bandFreq))
		# Referance antenna number
		c9 = pyfits.Column(name='REFANT_1', format='1J',
						array=numpy.ones((self.nAnt,), dtype=numpy.int32))
		# Real part of the bandpass
		c10 = pyfits.Column(name='BREAL_1', format='%dE' % self.nChan,
						array=numpy.ones((self.nAnt,self.nChan), dtype=numpy.float32))
		# Imagniary part of the bandpass
		c11 = pyfits.Column(name='BIMAG_1', format='%dE' % self.nChan,
						array=numpy.zeros((self.nAnt,self.nChan), dtype=numpy.float32))

		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11])

		# Create the Bandpass table and update its header
		bp = pyfits.new_table(colDefs)
		self.__addCommonKeywords(bp.header, 'BANDPASS', 1)

		# Figure out how many polarization are present.  First, reverse the 
		# StokesCodes dictionary so we can convert codes to names.  Second, go 
		# through all of the codes in self.stokes and look at them.  The list 
		# of polarizations used goes into pols.
		codes = {}
		for key in list(StokesCodes.keys()):
			codes[StokesCodes[key]] = key
		pols = []
		for sp in self.stokes:
			name = codes[sp]
			if len(name) == 2:
				part1 = name[0]
				part2 = name[1]
			else:
				part1 = name
				part2 = name
			if part1 not in pols:
				pols.append(part1)
			if part2 not in pols:
				pol.append(part2)

		bp.header.update('NO_ANT', self.nAnt)
		bp.header.update('NO_POL', len(pols))
		bp.header.update('NO_BACH', self.nChan)
		bp.header.update('STRT_CHN', self.refPix)

		bp.name = 'BANDPASS'
		self.FITS.append(bp)
		self.FITS.flush()

	def __writeSource(self):
		"""Define the Source table (group 1, table 2)."""

		hrz = astro.hrz_posn(0, 90)
		(arrPos, ag) = self.readArrayGeometry()
		ids = ag.keys()

		nameList = []
		raList = []
		decList = []
		raPoList = []
		decPoList = []
		sourceID = 1
		for dataSet in self.data:
			utc = astro.taimjd_to_utcjd(dataSet.obsTime)
			date = astro.get_date(utc)
			date.hours = 0
			date.minutes = 0
			date.seconds = 0
			utc0 = date.to_jd()
			equ = hrz.to_equ(arrPos, utc)
			
			# current apparent zenith equatorial coordinates
			raList.append(equ.ra)
			decList.append(equ.dec)
			
			# J2000 zenith equatorial coordinates
			equPo = astro.get_equ_prec2(equ, utc, astro.J2000_UTC_JD)
			raPoList.append(equPo.ra)
			decPoList.append(equPo.dec)
			
			# format 'source' name based on local sidereal time
			raHms = astro.deg_to_hms(equ.ra)
			(tsecs, secs) = math.modf(raHms.seconds)
			name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), int(tsecs * 10.0))
			nameList.append(name)
			sourceID += 1
		nSource = len(nameList)

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
		# Mean equinox
		c12 = pyfits.Column(name='EPOCH', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64) + 2000.0)
		# Appearrent right ascension in degrees
		c13 = pyfits.Column(name='RAAPP', format='1D', 
						array=numpy.array(raList))
		# Apparent declination in degrees
		c14 = pyfits.Column(name='DECAPP', format='1D', 
						array=numpy.array(decList))
		# Right ascension at mean equinox in degrees
		c15 = pyfits.Column(name='RAEPO', format='1D', 
						array=numpy.array(raPoList))
		# Declination at mean equinox in degrees
		c16 = pyfits.Column(name='DECEPO', format='1D', 
						array=numpy.array(decPoList))
		# Systemic velocity in m/s
		c17 = pyfits.Column(name='SYSVEL', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64))
		# Velocity type
		c18 = pyfits.Column(name='VELTYP', format='A8', 
						array=numpy.array(('GEOCENTR',)).repeat(nSource))
		# Velocity definition
		c19 = pyfits.Column(name='VELDEF', format='A8', 
						array=numpy.array(('OPTICAL',)).repeat(nSource))
		# Line rest frequency in Hz
		c20 = pyfits.Column(name='RESTFREQ', format='1D', 
						array=(numpy.zeros((nSource,), dtype=numpy.float64) + self.refVal))
		# Proper motion in RA in degrees/day
		c21 = pyfits.Column(name='PMRA', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64))
		# Proper motion in Dec in degrees/day
		c22 = pyfits.Column(name='PMDEC', format='1D', 
						array=numpy.zeros((nSource,), dtype=numpy.float64))
		# Parallax of source in arc sec.
		c23 = pyfits.Column(name='PARALLAX', format='1E', 
						array=numpy.zeros((nSource,), dtype=numpy.float32))

		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, 
							c21, c22, c23])

		# Create the Source table and update its header
		sr = pyfits.new_table(colDefs)
		self.__addCommonKeywords(sr.header, 'SOURCE', 1)
		
		sr.name = 'SOURCE'
		self.FITS.append(sr)
		self.FITS.flush()

	def __writeStars(self):
		"""Define the STARS table (group ?, table ?)."""

		nSource = len(self.FITS['SOURCE'].data.field('SOURCE_ID'))

		sunLngList = []
		sunLatList = []
		for dataSet in self.data:
			utc = astro.taimjd_to_utcjd(dataSet.obsTime)
			
			 #get the galactic coordinates of Sun's current position
			sunGal = astro.get_solar_equ_coords(utc).to_gal(utc)
			sunLngList.append(sunGal.l)
			sunLatList.append(sunGal.b)
        
		# Source ID number
		sourceId = pyfits.Column(name = 'SOURCE', format = '1J', 
							array = numpy.arange(1, nSource + 1, dtype = numpy.int32))
		# Solar galactic latitude in degrees
		sunLng = pyfits.Column(name = 'SUNGALL', format = '1D', unit = 'DEGREES',
							array = numpy.array(sunLngList))
		# Solar galactic longitude in degrees
		sunLat = pyfits.Column(name = 'SUNGALB', format = '1D', unit = 'DEGREES',
							array = numpy.array(sunLatList))
		
		# Create the Stars table and update its header
		st = pyfits.new_table([sourceId, sunLng, sunLat])
		self.__addCommonKeywords(st.header, 'STARS', 1) 

		st.name = 'STARS'
		self.FITS.append(st)
		self.FITS.flush()

	def __writeData(self):
		"""Define the UV_Data table (group 3, table 1)."""

		hrz = astro.hrz_posn(0, 90)
		(arrPos, ag) = self.readArrayGeometry()
		(mapper, inverseMapper) = self.readArrayMapper()
		ids = ag.keys()

		mList = []
		uList = []
		vList = []
		wList = []
		timeList = []
		dateList = []
		intTimeList = []
		blineList = []
		nameList = []
		sourceList = []
		sourceID = 1
		for dataSet in self.data:
			utc = astro.taimjd_to_utcjd(dataSet.obsTime)
			date = astro.get_date(utc)
			date.hours = 0
			date.minutes = 0
			date.seconds = 0
			utc0 = date.to_jd()
			equ = hrz.to_equ(arrPos, utc)
			
			# format 'source' name based on local sidereal time
			raHms = astro.deg_to_hms(equ.ra)
			(tsecs, secs) = math.modf(raHms.seconds)
			name = "ZA%02d%02d%02d%01d" % (raHms.hours, raHms.minutes, int(secs), 
				int(tsecs * 10.0))
			nameList.append(name)
			
			# Compute the uvw coordinates of all baselines
			if inverseMapper is not None:
				standIDs = []
				for stand in self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA'):
					standIDs.append(inverseMapper[stand])
				standIDs = numpy.array(standIDs)
			else:
				standIDs = self.FITS['ARRAY_GEOMETRY'].data.field('NOSTA')
			uvwBaselines = numpy.array([(i << 16) | j for i,j in uvUtils.getBaselines(standIDs)])
			uvwCoords = uvUtils.computeUVW(standIDs, HA=0.0, dec=equ.dec, freq=self.refVal)
			uvwCoords *= 1.0 / self.refVal

			# Loop over the data store in the dataDict and extract each baseline
			baselines = list(dataSet.dataDict.keys())
			baselines.sort()
			for baseline in baselines: 
		
				# validate baseline antenna ID's
				stand1 = (baseline >> 16) & 65535
				stand2 = baseline & 65535
				if mapper is not None:
					stand1 = mapper[stand1]
					stand2 = mapper[stand2]
				# Recontruct the baseline in FITS IDI format
				baselineMapped = (stand1 << 8) | stand2

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
				##matrix = numpy.empty((self.nChan, 2), dtype=numpy.float32)
				##matrix[:,0] = visData.real
				##matrix[:,1] = visData.imag
				matrix = numpy.empty((2*self.nChan,), dtype=numpy.float32)
				matrix[0::2] = visData.real
				matrix[1::2] = visData.imag
			
				blineList.append(baselineMapped)
				mList.append(matrix.ravel())
				uList.append(uvw[0])
				vList.append(uvw[1])
				wList.append(uvw[2])
				dateList.append(utc0)
				timeList.append(utc - utc0) 
				intTimeList.append(dataSet.intTime)
				sourceList.append(sourceID)
			sourceID += 1
		nBaseline = len(blineList)
		nSource = len(nameList)

		# Visibility Data
		c1 = pyfits.Column(name='FLUX', format='%iD' % (2*self.nChan), unit='UNCALIB', 
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
		# Inegration time (seconds)
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
		
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13])

		# Create the UV Data table and update its header
		uv = pyfits.new_table(colDefs)
		self.__addCommonKeywords(uv.header, 'UV_DATA', 1)

		uv.header.update('NMATRIX', 1, 'number of UV data matricies')
		uv.header.update('MAXIS', 6, 'number of UV data matrix axes')
		uv.header.update('TMATX1', True, 'axis 1 contains UV matrix')
		
		uv.header.update('MAXIS1', 2, 'number of pixels in COMPLEX axis')
		uv.header.update('CTYPE1', 'COMPLEX', 'axis 1 is COMPLEX axis')
		uv.header.update('CDELT1', 1.0)
		uv.header.update('CRPIX1', 1.0)
		uv.header.update('CRVAL1', 1.0)
		
		uv.header.update('MAXIS2', 1, 'number of pixels in STOKES axis')
		uv.header.update('CTYPE2', 'STOKES', 'axis 2 is STOKES axis (polarization)')
		uv.header.update('CDELT2', -1.0)
		uv.header.update('CRPIX2', 1.0)
		uv.header.update('CRVAL2', float(self.stokes[0]))
		
		uv.header.update('MAXIS3', self.nChan, 'number of pixels in FREQ axis')
		uv.header.update('CTYPE3', 'FREQ', 'axis 3 is FREQ axis (frequency)')
		uv.header.update('CDELT3', self.freq[0].chWidth)
		uv.header.update('CRPIX3', self.refPix)
		uv.header.update('CRVAL3', self.refVal)
		
		uv.header.update('MAXIS4', 1, 'number of pixels in BAND (IF) axis')
		uv.header.update('CTYPE4', 'BAND', 'axis 4 is BAND axis')
		uv.header.update('CDELT4', 1.0)
		uv.header.update('CRPIX4', 1.0)
		uv.header.update('CRVAL4', 1.0)
		
		uv.header.update('MAXIS5', 1, 'number of pixels in RA axis')
		uv.header.update('CTYPE5', 'RA', 'axis 5 is RA axis (position of phase center)')
		uv.header.update('CDELT5', 0.0)
		uv.header.update('CRPIX5', 1.0)
		uv.header.update('CRVAL5', 0.0)
		
		uv.header.update('MAXIS6', 1, 'number of pixels in DEC axis')
		uv.header.update('CTYPE6', 'DEC', 'axis 6 is DEC axis (position of phase center)')
		uv.header.update('CDELT6', 0.0)
		uv.header.update('CRPIX6', 1.0)
		uv.header.update('CRVAL6', 0.0)
		
		uv.header.update('TELESCOP', 'LWA')
		uv.header.update('OBSERVER', 'ZASKY')
		uv.header.update('SORT', 'TB', 'data is sorted in [time,baseline] order')
		
		uv.header.update('VISSCALE', 1.0, 'UV data scale factor')
		
		uv.name = 'UV_DATA'
		self.FITS.append(uv)
		self.FITS.flush()

	def __writeMapper(self):
		"""Write a fits table that contains information about mapping stations 
		numbers to actual antenna numbers.  This information can be backed out of
		the names, but this makes the extraction more programatic."""

		c1 = pyfits.Column(name='ANNAME', format='A8', 
						array=numpy.array([ant.getName() for ant in self.array[0]['ants']]))
		c2 = pyfits.Column(name='NOSTA', format='1J', 
						array=numpy.array([self.array[0]['mapper'][ant.id] for ant in self.array[0]['ants']]))
		c3 = pyfits.Column(name='NOACT', format='1J', 
						array=numpy.array([ant.id for ant in self.array[0]['ants']]))

		colDefs = pyfits.ColDefs([c1, c2, c3])

		# Create the ID mapping table and update its header
		nsm = pyfits.new_table(colDefs)
		self.__addCommonKeywords(nsm.header, 'NOSTA_MAPPER', 1)

		nsm.name = 'NOSTA_MAPPER'
		self.FITS.append(nsm)
		self.FITS.flush()

	def readArrayGeometry(self):
		"""Return a tuple with the array geodetic position and the local 
		positions for all antennas defined in the ARRAY_GEOMETRY table."""

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
		"""Return a tuple with the array NOSTA mapper and inverse mapper (both
		dictionaries.  If the stand IDs have not been mapped, return None for
		both."""

		try:
			nsm = self.FITS['NOSTA_MAPPER']
		except IndexError:
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