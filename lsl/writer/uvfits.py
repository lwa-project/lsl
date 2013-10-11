# -*- coding: utf-8 -*-

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
from lsl.misc import geodesy

from lsl.writer.fitsidi import StokesCodes, NumericStokes

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['UV', 'StokesCodes', 'NumericStokes', '__version__', '__revision__', '__all__']


UVVersion = (1, 0)


def mergeBaseline(ant1, ant2):
	"""
	Merge two stand ID numbers into a single baseline.
	"""
	
	if ant1 > 255 or ant2 > 255:
		baseline = ant1*2048 + ant2 + 65536
	else:
		baseline = ant1*256 + ant2
		
	return baseline


def splitBaseline(baseline):
	"""
	Given a baseline, split it into it consistent stand ID numbers.
	"""
	
	if baseline >= 65536:
		ant1 = int((baseline - 65536) / 2048)
		ant2 = int((baseline - 65536) % 2048)
	else:
		ant1 = int(baseline / 256)
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
			
	class _UVData(object):
		"""
		Represents one UV visibility data set for a given observation time.
		"""
		
		def __init__(self, obsTime, intTime, dataDict, pol=StokesCodes['XX']):
			self.obsTime = obsTime
			self.intTime = intTime
			self.dataDict = dataDict
			self.pol = pol
			
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
		Initialize a new UVFITS object using a filename and a reference time 
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
			print "UVFITS: stand ID mapping enabled"
			for key, value in mapper.iteritems():
				print "UVFITS:  stand #%i -> mapped #%i" % (key, value)
				
		self.nAnt = len(ants)
		self.array.append( {'center': [arrayX, arrayY, arrayZ], 'ants': ants, 'mapper': mapper, 'enableMapper': enableMapper, 'inputAnts': antennas} )
		
	def addDataSet(self, obsTime, intTime, baselines, visibilities, pol='XX'):
		"""
		Create a UVData object to store a collection of visibilities.
		
		.. versionchanged:: 0.4.0
			Switched over to passing in Antenna instances generated by the
			:mod:`lsl.common.stations` module instead of a list of stand ID
			as part of the baselines.
		"""
		
		if type(pol).__name__ == 'str':
			numericPol = StokesCodes[pol.upper()]
		else:
			numericPol = pol
			
		dataDict = {}
		for bl, (ant1,ant2) in enumerate(baselines):
			baseline = mergeBaseline(ant1.stand.id, ant2.stand.id)
			dataDict[baseline] = visibilities[bl,:].squeeze()
			
		self.data.append( self._UVData(obsTime, intTime, dataDict, pol=numericPol) )
		
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
		
		hdr.update('EXTNAME', name, 'UVFITS table name')
		hdr.update('EXTVER', 1, 'table instance number') 
		hdr.update('TABREV', revision, 'table format revision number')
		
		date = self.refTime.split('-')
		name = "ZA%s%s%s" % (date[0][2:], date[1], date[2])
		hdr.update('OBSCODE', name, 'zenith all-sky image')
		
		hdr.update('ARRNAM', self.siteName)      
		hdr.update('RDATE', self.refTime, 'file data reference date')
		
	def _makeAppendTable(self, extension, AddRows=1):
		"""
		Private function to make a temporary table for appending data.
		"""
		
		nrows = self.hdulist[extension].data.shape[0]
		tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+AddRows)
		for key in list(self.hdulist[extension].header.keys()):
			tempHDU.header.update(key, self.hdulist[extension].header[key])
			
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
		
		self._writeGeometry(dummy=True)
		hrz = astro.hrz_posn(0, 90)
		(arrPos, ag) = self.readArrayGeometry(dummy=True)
		(mapper, inverseMapper) = self.readArrayMapper(dummy=True)
		ids = ag.keys()
		
		# Retrieve the original list of Antenna objects and convert them to
		# a dictionary index by the stand ID number
		inputAnts = {}
		for ant in self.array[0]['inputAnts']:
			inputAnts[ant.stand.id] = ant
			
		mList = []
		uList = []
		vList = []
		wList = []
		dateList = []
		blineList = []
		rawList = []
		for dataSet in self.data:
			if dataSet.pol == self.stokes[0]:
				utc = astro.taimjd_to_utcjd(dataSet.obsTime)
				date = astro.get_date(utc)
				date.hours = 0
				date.minutes = 0
				date.seconds = 0
				utc0 = date.to_jd()
				equ = hrz.to_equ(arrPos, utc)
				
				# format 'source' name based on local sidereal time
				raHms = astro.deg_to_hms(equ.ra)
				
				# Compute the uvw coordinates of all baselines
				if inverseMapper is not None:
					standIDs = []
					for stand in self.an.data.field('NOSTA'):
						standIDs.append(inverseMapper[stand])
					standIDs = numpy.array(standIDs)
				else:
					standIDs = self.an.data.field('NOSTA')
				antennaStands = []
				for standID in standIDs:
					antennaStands.append( inputAnts[standID] )
					
				uvwBaselines = numpy.array([mergeBaseline(a1.stand.id, a2.stand.id) for a1,a2 in uvUtils.getBaselines(antennaStands)])
				uvwCoords = uvUtils.computeUVW(antennaStands, HA=0.0, dec=equ.dec, freq=self.refVal)
				uvwCoords *= 1.0 / self.refVal
				
				tempMList = {}
				for stokes in self.stokes:
					tempMList[stokes] = {}
					
			# Loop over the data store in the dataDict and extract each baseline
			baselines = list(dataSet.dataDict.keys())
			baselines.sort()
			for baseline in baselines: 
				# validate baseline antenna ID's
				stand1, stand2 = splitBaseline(baseline)
				if mapper is not None:
					stand1 = mapper[stand1]
					stand2 = mapper[stand2]
				# Reconstruct the baseline in UVFITS/MIRIAD format
				baselineMapped = mergeBaseline(stand1, stand2)
				
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
				matrix = numpy.zeros((2,self.nChan,), dtype=numpy.float32)
				matrix[0,:] = visData.real
				matrix[1,:] = visData.imag
				tempMList[dataSet.pol][baseline] = matrix
				
				if dataSet.pol == self.stokes[0]:
					blineList.append(baselineMapped)
					rawList.append(baseline)
					uList.append(uvw[0])
					vList.append(uvw[1])
					wList.append(uvw[2])
					dateList.append(utc)
					
			if dataSet.pol == self.stokes[-1]:
				for bl in rawList:
					matrix = numpy.zeros((1,1,self.nChan,self.nStokes,2), dtype=numpy.float32)
					for p in xrange(self.nStokes):
						try:
							matrix[0,0,:,p,0] = tempMList[self.stokes[p]][bl][0,:]
							matrix[0,0,:,p,1] = tempMList[self.stokes[p]][bl][1,:]
						except KeyError:
							stand1, stand2 = splitBaseline(bl)
							newBL = mergeBaseline(stand2, stand1)
							print 'WARNING: Keyerror', bl, bl in tempMList[self.stokes[p]], stand1, stand2, newBL, newBL in tempMList[self.stokes[p]]
							
					mList.append(matrix)
					
				rawList = []
		nBaseline = len(blineList)
		
		# Create the UV Data table and update its header
		test = numpy.array(mList, dtype=numpy.float32)
		print test.shape, self.nStokes, self.nChan
		uv = pyfits.GroupData(numpy.array(mList, dtype=numpy.float32), parnames=['UU', 'VV', 'WW', 'BASELINE', 'DATE'], 
							pardata=[numpy.array(uList, dtype=numpy.float32), numpy.array(vList, dtype=numpy.float32), 
									numpy.array(wList, dtype=numpy.float32), numpy.array(blineList), 
									numpy.array(dateList)], bitpix=-32)
		primary = pyfits.GroupsHDU(uv)
		
		primary.header.update('EXTEND', True, 'indicates UVFITS file')
		primary.header.update('GROUPS', True, 'indicates UVFITS file')
		primary.header.update('OBJECT', 'BINARYTB')
		primary.header.update('TELESCOP', self.siteName)
		primary.header.update('INSTRUME', self.siteName)
		primary.header.update('OBSERVER', 'ZASKY', 'zenith all-sky image')
		primary.header.update('ORIGIN', 'LSL')
		primary.header.update('CORRELAT', 'LWASWC', 'Correlator used')
		primary.header.update('FXCORVER', '1', 'Correlator version')
		primary.header.update('LWATYPE', 'UV-ZA', 'LWA FITS file type')
		primary.header.update('LWAMAJV', UVVersion[0], 'LWA UVFITS file format major version')
		primary.header.update('LWAMINV', UVVersion[1], 'LWA UVFITS file format minor version')
		primary.header.update('DATE-OBS', self.refTime, 'UVFITS file data collection date')
		ts = str(astro.get_date_from_sys())
		primary.header.update('DATE-MAP', ts.split()[0], 'UVFITS file creation date')
		
		primary.header.update('CTYPE2', 'COMPLEX', 'axis 2 is COMPLEX axis')
		primary.header.update('CDELT2', 1.0)
		primary.header.update('CRPIX2', 1.0)
		primary.header.update('CRVAL2', 1.0)
		
		primary.header.update('CTYPE3', 'STOKES', 'axis 3 is STOKES axis (polarization)')
		if self.stokes[0] < 0:
			primary.header.update('CDELT3', -1.0)
		else:
			primary.header.update('CDELT3', 1.0)
		primary.header.update('CRPIX3', 1.0)
		primary.header.update('CRVAL3', float(self.stokes[0]))
		
		primary.header.update('CTYPE4', 'FREQ', 'axis 4 is FREQ axis (frequency)')
		primary.header.update('CDELT4', self.freq[0].chWidth)
		primary.header.update('CRPIX4', self.refPix)
		primary.header.update('CRVAL4', self.refVal)
		
		primary.header.update('CTYPE5', 'RA', 'axis 5 is RA axis (position of phase center)')
		primary.header.update('CDELT5', 0.0)
		primary.header.update('CRPIX5', 1.0)
		primary.header.update('CRVAL5', 0.0)
		
		primary.header.update('CTYPE6', 'DEC', 'axis 6 is DEC axis (position of phase center)')
		primary.header.update('CDELT6', 0.0)
		primary.header.update('CRPIX6', 1.0)
		primary.header.update('CRVAL6', 0.0)
		
		primary.header.update('TELESCOP', self.siteName)
		primary.header.update('OBSERVER', 'ZASKY')
		primary.header.update('SORT', 'TB', 'data is sorted in [time,baseline] order')
		
		primary.header.update('VISSCALE', 1.0, 'UV data scale factor')
		
		self.FITS.append(primary)
		self.FITS.flush()
		
	def _writeGeometry(self, dummy=False):
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
		c1 = pyfits.Column(name='ANNAME', format='A8', 
						array=numpy.array([ant.getName() for ant in self.array[0]['ants']]))
		# Station coordinates in meters
		c2 = pyfits.Column(name='STABXYZ', unit='METERS', format='3D', 
						array=xyz)
		# Station number
		c3 = pyfits.Column(name='NOSTA', format='1J', 
						array=numpy.array([self.array[0]['mapper'][ant.id] for ant in self.array[0]['ants']]))
		# Mount type (0 == alt-azimuth)
		c4 = pyfits.Column(name='MNTSTA', format='1J', 
						array=numpy.zeros((self.nAnt,), dtype=numpy.int32))
		# Axis offset in meters
		c5 = pyfits.Column(name='STAXOF', unit='METERS', format='3E', 
						array=numpy.zeros((self.nAnt,3), dtype=numpy.float32))
		# Feed A polarization label
		c6 = pyfits.Column(name='POLTYA', format='A1', 
						array=numpy.array([ant.polA['Type'] for ant in self.array[0]['ants']]))
		# Feed A orientation in degrees
		c7 = pyfits.Column(name='POLAA', format='1E', 
						array=numpy.array([ant.polA['Angle'] for ant in self.array[0]['ants']], dtype=numpy.float32))
		# Feed A polarization parameters
		c8 = pyfits.Column(name='POLCALA', format='2E', 
						array=numpy.array([ant.polA['Cal'] for ant in self.array[0]['ants']], dtype=numpy.float32))
		# Feed B polarization label
		c9 = pyfits.Column(name='POLTYB', format='A1', 
						array=numpy.array([ant.polB['Type'] for ant in self.array[0]['ants']]))
		# Feed B orientation in degrees
		c10 = pyfits.Column(name='POLAB', format='1E', 
						array=numpy.array([ant.polB['Angle'] for ant in self.array[0]['ants']], dtype=numpy.float32))
		# Feed B polarization parameters
		c11 = pyfits.Column(name='POLCALB', format='2E', 
						array=numpy.array([ant.polB['Cal'] for ant in self.array[0]['ants']], dtype=numpy.float32))
						
		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11])
		
		# Create the table and fill in the header
		ag = pyfits.new_table(colDefs)
		self._addCommonKeywords(ag.header, 'AIPS AN', 1)
		
		ag.header.update('EXTVER', 1, 'array ID')
		ag.header.update('ARRNAM', self.siteName)
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
		
		refDate = self.refTime2AstroDate()
		refMJD = refDate.to_jd() - astro.MJD_OFFSET
		eop = geodesy.getEOP(refMJD)
		if eop[0] is None:
			eop = [geodesy.EOP(mjd=refMJD),]
			
		ag.header.update('UT1UTC', eop[0].utDiff, 'difference UT1 - UTC for reference date')
		ag.header.update('IATUTC', astro.leap_secs(utc0), 'TAI - UTC for reference date')
		ag.header.update('POLARX', eop[0].x)
		ag.header.update('POLARY', eop[0].y)
		
		ag.header.update('ARRAYX', self.array[0]['center'][0], 'array ECI X coordinate (m)')
		ag.header.update('ARRAYY', self.array[0]['center'][1], 'array ECI Y coordinate (m)')
		ag.header.update('ARRAYZ', self.array[0]['center'][2], 'array ECI Z coordinate (m)')
		
		ag.header.update('NOSTAMAP', int(self.array[0]['enableMapper']), 'Mapping enabled for stand numbers')
		
		if dummy:
			self.an = ag
			if self.array[0]['enableMapper']:
				self._writeMapper(dummy=True)
				
		else:
			ag.name = 'AIPS AN'
			self.FITS.append(ag)
			self.FITS.flush()
			
			if self.array[0]['enableMapper']:
				self._writeMapper()
				
	def _writeMapper(self, dummy=False):
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
		
		if dummy:
			self.am = nsm
		else:
			nsm.name = 'NOSTA_MAPPER'
			self.FITS.append(nsm)
			self.FITS.flush()
			
	def readArrayGeometry(self, dummy=False):
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
			antennaGeo[id] = antennaPos.next()
			
		# Return
		return (arrayGeo, antennaGeo)
		
	def readArrayMapper(self, dummy=False):
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
