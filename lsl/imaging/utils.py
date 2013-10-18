# -*- coding: utf-8 -*-

"""
Module to support imaging correlated data.  This module provides utilities to read FITS IDI files
into data dictionaries (as described in :mod:`lsl.sim.vis`) and build AIPY ImgW instances from the 
data dictionaries.  Also included is a utility to sort data dictionaries by baselines.

.. versionadded:: 0.5.0

.. versionchanged:: 0.7.0
	Added support for UVFITS files and CASA measurement sets
"""

import os
import aipy
import pytz
import numpy
import pyfits
import string
from calendar import timegm
from datetime import datetime
from operator import itemgetter

from lsl import astro
from lsl.common import stations
from lsl.sim import vis as simVis
from lsl.writer.fitsidi import NumericStokes
from lsl.common.constants import c as vLight

__version__ = '0.4'
__revision__ = '$Rev$'
__all__ = ['baselineOrder', 'sortDataDict', 'pruneBaselineRange', 'rephaseData', 'CorrelatedData', 'CorrelatedDataIDI', 'CorrelatedDataUV', 'CorrelatedDataMS', 'buildGriddedImage', '__version__', '__revision__', '__all__']


def baselineOrder(bls):
	"""
	Like numpy.argsort(), but for a list of two-element tuples of baseline 
	pairs.  The resulting lists can then be used to sort a data dictionary
	a la sortDataDict().
	"""
	
	def __cmpBaseline(bl):
		return 1024*bl[0] + bl[1]
	
	return [i for (v, i) in sorted((v, i) for (i, v) in enumerate([__cmpBaseline(bl) for bl in bls]))]


def sortDataDict(dataDict, order=None):
	"""
	Sort a data dictionary by the specified order.  If no order is supplied, 
	the data dictionary is sorted by baseline using baselineOrder().
	"""
	
	if order is None:
		for pol in ['xx', 'yy', 'I']:
			try:
				if len(dataDict['bls'][pol]) == 0:
					continue
				order = baselineOrder(dataDict['bls'][pol])
				break
			except KeyError:
				pass
	
	for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
		for pol in dataDict[key].keys():
			try:
				newList = [dataDict[key][pol][i] for i in order]
				dataDict[key][pol] = newList
			except (KeyError, IndexError):
				pass
			
	return dataDict


def pruneBaselineRange(dataDict, uvMin=0, uvMax=numpy.inf):
	"""
	Prune baselines from a data dictionary that are less than uvMin or
	greater than or equal to uvMax.

	.. note::
		uvMin and uvMax should be specified in lambda
	"""

	# Force min to be less than max
	if uvMin > uvMax:
		temp = uvMin
		uvMin = uvMax
		uvMax = temp
		
	# Create the new output data dictionary
	newDict = {}
	for key in dataDict.keys():
		if key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
			newDict[key] = {}
			for pol in dataDict[key].keys():
				newDict[key][pol] = []
		else:
			newDict[key] = dataDict[key]

	# Find out the baseline lengths and create a list of good ones
	good = {}
	freq = dataDict['freq']
	for pol in dataDict['uvw'].keys():
		sizes = []
		for bl in dataDict['uvw'][pol]:
			uvw = bl[:,freq.size/2]
			sizes.append( numpy.sqrt((uvw**2).sum()) )
		sizes = numpy.array(sizes)
		good[pol] = list(numpy.where( (sizes >= uvMin) & (sizes < uvMax) )[0])
		
	# Prune
	for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
		for pol in dataDict[key].keys():
			lgp = len(good[pol])
			if lgp == 0:
				newDict[key][pol] = []
			elif lgp == 1:
				newDict[key][pol] = [dataDict[key][pol][good[pol][0]],]
			else:
				newDict[key][pol] = list(itemgetter(*good[pol])(dataDict[key][pol]))
				
	# Return
	return newDict


def rephaseData(aa, dataDict, currentPhaseCenter='z', newPhaseCenter='z'):
	"""
	Given an AntennaArray instance and a data dictionary, re-phase the data 
	to change the pointing center.
	"""
	
	# Load in basic information about the data
	freq = dataDict['freq']*1.0
	isMasked = dataDict['isMasked']
	pols = dataDict['bls'].keys()
	
	# Create the data dictionary that will hold the re-phased data
	dataDict2 = {}
	dataDict2['freq'] = freq 
	dataDict2['isMasked'] = isMasked
	for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
		dataDict2[key] = {}
		for p in pols:
			dataDict2[key][p] = []
			
	# Go!
	for p in pols:
		lastJD = None
		
		## Loop over baselines
		for k in xrange(len(dataDict['bls'][p])):
			### Load in the data
			i,j = dataDict['bls'][p][k]
			d   = dataDict['vis'][p][k]
			msk = dataDict['msk'][p][k]
			jd  = dataDict[ 'jd'][p][k]
			
			### Update the source positions if needed
			if jd != lastJD:
				# Update the time for the AntennaArray
				aa.set_jultime(jd)
				
				# Recompute
				if currentPhaseCenter is not 'z':
					currentPhaseCenter.compute(aa)
				if newPhaseCenter is not 'z':
					newPhaseCenter.compute(aa)
					
				lastJD = jd
				
			### Compute the uvw coordinates and the new phasing
			try:
				crd = aa.gen_uvw(j, i, src=newPhaseCenter)
				d = aa.unphs2src(d, currentPhaseCenter, j, i)
				d = aa.phs2src(d, newPhaseCenter, j, i)
			except aipy.phs.PointingError:
				raise RuntimeError("Rephasing center is below the horizon")
				
			### Save
			uvw = aa.gen_uvw(j, i, src=newPhaseCenter)
			uvw = numpy.squeeze(crd.compress(numpy.logical_not(msk), axis=2))
			vis = d.compress(numpy.logical_not(msk))
			wgt = numpy.ones_like(vis) * len(vis)

			dataDict2['bls'][p].append( (i,j) )
			dataDict2['uvw'][p].append( uvw )
			dataDict2['vis'][p].append( vis )
			dataDict2['wgt'][p].append( wgt )
			dataDict2['msk'][p].append( msk )
			dataDict2[ 'jd'][p].append( jd )
			
	# Done
	return dataDict2


def CorrelatedData(filename):
	"""
	Read in and work with FITS IDI and UVFITS files.  Returns either a 
	CorrelateDataIDI or CorrelatedDataUV instance.
	"""
	
	valid = False
	
	# FITS IDI
	try:
		return CorrelatedDataIDI(filename)
	except:
		pass
		
	# UVFITS
	try:
		return CorrelatedDataUV(filename)
	except:
		pass
		
	# Measurment Set
	try:
		return CorrelatedDataMS(filename)
	except:
		pass
		
	if not valid:
		raise RuntimeError("File '%s' does not appear to be either a FITS IDI file, UV FITS file, or MeasurmentSet" % filename)


class CorrelatedDataIDI(object):
	"""
	Class to make accessing information about a FITS IDI easy.  This wraps 
	all of the "messy" machinery needed to extract both the metadata and data 
	from the file and return them as common LSL objects.
	
	This class has three main attributes to interact with:
	  * getAntennaArray - Return a :class:`lsl.sim.vim.AntennaArray` instance
	    that represents the array where the data was obtained.  This is useful
	    for simulation proposes and computing source positions.
	  * getObserver - Return a ephem.Observer instance representing the array
	  * getDataSet - Return a data dictionary of all baselines for a given set
	    of observations
	    
	The class also includes a variety of useful metadata attributes:
	  * pols - Numpy array of polarization product codes
	  * freq - Numpy array of frequency channels in Hz
	  * station - LSL :class:`lsl.common.stations.LWAStation` instance for the
	    array
	  * dateObs - Datetime object for the reference date of the FIT IDI file
	  * antennas - List of :class:`lsl.common.stations.Antenna` instances
	  
	.. note::
		The CorrelatedData.antennas attribute should be used over 
		CorrelatedData.station.getAntennas() since the mapping in the FITS IDI
		file may not be the same as the digitizer order.
	"""
	
	def _createEmptyDataDict(self):
		"""
		Create an empty data dictionary that is appropriate for the current file.
		"""
		
		dataDict = {}
		dataDict['freq'] = self.freq 
		dataDict['isMasked'] = False
		for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
			dataDict[key] = {}
			
		for p in self.pols:
			name = NumericStokes[p]
			if len(name) == 2:
				name = name.lower()
			for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
				dataDict[key][name] = []
				
		return dataDict
	
	def __init__(self, filename):
		"""
		Initialize a new CorrelatedDataIDI instance from a FITS IDI file and 
		fill in the metadata.
		"""
		
		self.filename = filename
		
		# Open the file, check if it looks like FITS IDI, and pull out the UV_DATA table
		hdulist = pyfits.open(self.filename)
		tbls = [i.header['EXTNAME'] for i in hdulist[1:]]
		for tbl in ('ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'BANDPASS', 'SOURCE', 'UV_DATA'):
			if tbl not in tbls:
				raise RuntimeError("Cannot find table '%s' in '%s'" % (tbl, self.filename))
		
		self.extended = False
		if hdulist[0].header['LWATYPE'] == 'IDI-EXTENDED-ZA':
			self.extended = True
		ag = hdulist['ARRAY_GEOMETRY']
		uvData = hdulist['UV_DATA']
		
		# Antennas
		try:
			mapper = hdulist['NOSTA_MAPPER']
			
			nosta = mapper.data.field('NOSTA')
			noact = mapper.data.field('NOACT')
			stabxyz = ag.data.field('STABXYZ')
			anname = mapper.data.field('ANNAME')
		except KeyError:
			nosta = ag.data.field('NOSTA')
			noact = ag.data.field('NOSTA')
			stabxyz = ag.data.field('STABXYZ')
			anname = ag.data.field('ANNAME')
		
		# Station/telescope information
		self.telescope = hdulist[0].header['TELESCOP']
		self.dateObs = pytz.UTC.localize(datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S"))
		if self.telescope == 'LWA-1' or self.telescope == 'LWA1':
			self.station = stations.lwa1
		elif self.telescope == 'LWA-NA' or self.telescope == 'LWANA':
			self.station = stations.lwana
		else:
			geo = numpy.array([ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ']])
			site = stations.ecef2geo(*geo)
			
			lat  = site[0]
			ecii = numpy.array([[ 0.0,            1.0, 0.0           ],
							[-numpy.sin(lat), 0.0, numpy.cos(lat)],
							[ numpy.cos(lat), 0.0, numpy.sin(lat)]])
							
			antennas = []
			for line,act in zip(ag.data, noact):
				enz = numpy.dot(ecii, line['STABXYZ'])
				
				stand = stations.Stand(act, *enz)
				antennas.append(stations.Antenna(2*(stand.id-1)+1, stand=stand, pol=0))
				
			self.station = stations.LWAStation(ag.header['ARRNAM'], site[0]*180/numpy.pi, site[1]*180/numpy.pi, site[2], antennas=antennas)
		self.standMap = {}
		self.stands = []
		for sta, act in zip(nosta, noact):
			self.standMap[sta] = act
			self.stands.append(act)
			
		self.antennaMap = {}
		self.antennas = []
		for ant in self.station.getAntennas():
			if ant.stand.id in self.stands and ant.pol == 0:
				self.antennas.append(ant)
				self.antennaMap[ant.stand.id] = ant
				
		# Polarization and frequency
		self.pols  = numpy.arange(1, uvData.header['MAXIS2']+1) - uvData.header['CRPIX2']
		self.pols *= uvData.header['CDELT2'] 
		self.pols += uvData.header['CRVAL2']
		self.freq  = numpy.arange(1, uvData.header['NO_CHAN']+1) - uvData.header['REF_PIXL']
		self.freq *= uvData.header['CHAN_BW']
		self.freq += uvData.header['REF_FREQ']
		
		# Total baseline count
		self.totalBaselineCount = len(uvData.data['BASELINE'])
		
		# Close
		hdulist.close()
	
	def getAntennaArray(self):
		"""
		Return an AIPY AntennaArray instance for the array that made the 
		observations contained here.
		"""
		
		# Get the date of observations
		refJD = astro.unix_to_utcjd(timegm(self.dateObs.timetuple()))
		
		# Return
		return simVis.buildSimArray(self.station, self.antennas, self.freq/1e9, jd=refJD)
		
	def getObserver(self):
		"""
		Return a ephem.Observer instances for the array described in the file.
		"""
		
		return self.station.getObserver()
		
	def getDataSet(self, set, includeAuto=False, sort=True, uvMin=0, uvMax=numpy.inf):
		"""
		Return a baseline sorted data dictionary for the specified data set.  
		By default this excludes the autocorrelations.  To include 
		autocorrelations set the value of 'includeAuto' to True.  Setting the
		'sort' keyword to False will disable the baseline sorting.  Optionally,
		baselines with lengths between uvMin and uvMax can only be returned.

		.. note::
			uvMin and uvMax should be specified in lambda
		"""
		
		# Open the file
		hdulist = pyfits.open(self.filename)
		uvData = hdulist['UV_DATA']
		
		# We need this a lot...
		nPol = len(self.pols)
		
		# Define the dictionary to return
		dataDict = self._createEmptyDataDict()

		# Set the source ID to look for (this is LWA specific)
		sourceID = range(set,set+1)
		
		# Loop over data rows
		found = False
		for row in uvData.data:
			# If we are on the right set...
			if row['SOURCE'] in sourceID:
				found = True
				
				# Load it.
				bl = row['BASELINE']
				if not self.extended:
					i = self.standMap[(bl >> 8) & 255]
					j = self.standMap[bl & 255]
				else:
					i = self.standMap[(bl >> 16) & 65535]
					j = self.standMap[bl & 65535]
				if i == j and not includeAuto:
					## Skip auto-correlations
					continue
				ri = numpy.where(self.stands == i)[0][0]
				rj = numpy.where(self.stands == j)[0][0]

				uvw = numpy.array([row['UU'], row['VV'], row['WW']])
				
				jd = row['DATE'] + row['TIME']
				uvw = numpy.array([numpy.dot(uvw[0], self.freq), numpy.dot(uvw[1], self.freq), numpy.dot(uvw[2], self.freq)])
				flux = row['FLUX']
				
				for c,p in enumerate(self.pols):
					name = NumericStokes[p]
					if len(name) == 2:
						name = name.lower()
					
					vis = numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
					vis.real = flux[2*c+0::(2*nPol)]
					vis.imag = flux[2*c+1::(2*nPol)]
					wgt = numpy.ones(vis.size)
				
					dataDict['uvw'][name].append( uvw ) 
					dataDict['vis'][name].append( vis )
					dataDict['wgt'][name].append( wgt )
					dataDict['msk'][name].append( numpy.zeros(len(vis), dtype=numpy.int16) )
					dataDict['bls'][name].append( (ri,rj) )
					dataDict['jd' ][name].append( jd )
		
		# Close
		hdulist.close()
		
		# Make sure we found something
		if not found:
			raise RuntimeError("Cannot find baseline set %i in FITS IDI file", set)
		
		# Sort
		if sort:
			sortDataDict(dataDict)
			
		# Prune
		if uvMin != 0 or uvMax != numpy.inf:
			dataDict = pruneBaselineRange(dataDict, uvMin=uvMin, uvMax=uvMax)
			
		# Return
		return dataDict


class CorrelatedDataUV(object):
	"""
	Class to make accessing information about a UVFITS file easy.  This wraps 
	all of the "messy" machinery needed to extract both the metadata and data 
	from the file and return them as common LSL objects.
	
	This class has three main attributes to interact with:
	  * getAntennaArray - Return a :class:`lsl.sim.vim.AntennaArray` instance
	    that represents the array where the data was obtained.  This is useful
	    for simulation proposes and computing source positions.
	  * getObserver - Return a ephem.Observer instance representing the array
	  * getDataSet - Return a data dictionary of all baselines for a given set
	    of observations
	    
	The class also includes a variety of useful metadata attributes:
	  * pols - Numpy array of polarization product codes
	  * freq - Numpy array of frequency channels in Hz
	  * station - LSL :class:`lsl.common.stations.LWAStation` instance for the
	    array
	  * dateObs - Datetime object for the reference date of the FIT IDI file
	  * antennas - List of :class:`lsl.common.stations.Antenna` instances
	  
	.. note::
		The CorrelatedDataUV.antennas attribute should be used over 
		CorrelatedDataUV.station.getAntennas() since the mapping in the UVFITS
		file may not be the same as the digitizer order.
	"""
	
	def _createEmptyDataDict(self):
		"""
		Create an empty data dictionary that is appropriate for the current file.
		"""
		
		dataDict = {}
		dataDict['freq'] = self.freq 
		dataDict['isMasked'] = False
		for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
			dataDict[key] = {}
			
		for p in self.pols:
			name = NumericStokes[p]
			if len(name) == 2:
				name = name.lower()
			for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
				dataDict[key][name] = []
				
		return dataDict
	
	def __init__(self, filename):
		"""
		Initialize a new CorrelatedDataUV instance from a FITS IDI file and 
		fill in the metadata.
		"""
		
		self.filename = filename
		
		# Open the various tables that we need
		hdulist = pyfits.open(filename)
		
		uvData = hdulist[0]
		ag = hdulist['AIPS AN']
		
		# Station/telescope information
		self.telescope = hdulist[0].header['TELESCOP']
		dt = hdulist[0].header['DATE-OBS']
		dt = dt.rsplit('.', 1)[0]
		self.dateObs = pytz.UTC.localize(datetime.strptime(dt, "%Y-%m-%dT%H:%M:%S"))
		if self.telescope == 'LWA-1' or self.telescope == 'LWA1':
			self.station = stations.lwa1
		elif self.telescope == 'LWA-2' or self.telescope == 'LWA2':
			self.station = stations.lwa2
		else:
			geo = numpy.array([ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ']])
			site = stations.ecef2geo(*geo)
			
			lat  = site[0]
			ecii = numpy.array([[ 0.0,            1.0, 0.0           ],
							[-numpy.sin(lat), 0.0, numpy.cos(lat)],
							[ numpy.cos(lat), 0.0, numpy.sin(lat)]])
							
			antennas = []
			for line in ag.data:
				enz = numpy.dot(ecii, line['STABXYZ'])
				
				stand = stations.Stand(line['NOSTA'], *enz)
				antennas.append( stations.Antenna(2*(stand.id-1)-1, stand=stand) )
				
			self.station = stations.LWAStation(ag.header['ARRNAM'], site[0]*180/numpy.pi, site[1]*180/numpy.pi, site[2], antennas=antennas)
			
		# Antennas
		nosta = ag.data.field('NOSTA')
		noact = ag.data.field('NOSTA')
		anname = ag.data.field('ANNAME')
			
		self.standMap = {}
		self.stands = []
		for nosta, noact in zip(nosta, noact):
			self.standMap[nosta] = noact
			self.stands.append(noact)
			
		self.antennaMap = {}
		self.antennas = []
		if self.station is not None:
			for ant in self.station.getAntennas():
				if ant.stand.id in self.stands and ant.pol == 0:
					self.antennas.append(ant)
					self.antennaMap[ant.stand.id] = ant
		
		# Polarization and frequency
		self.pols  = numpy.arange(1, uvData.header['NAXIS3']+1) - uvData.header['CRPIX3']
		self.pols *= uvData.header['CDELT3'] 
		self.pols += uvData.header['CRVAL3']
		self.freq  = numpy.arange(1, uvData.header['NAXIS4']+1) - uvData.header['CRPIX4']
		self.freq *= uvData.header['CDELT4']
		self.freq += uvData.header['CRVAL4']
		
		# Total baseline count
		self.totalBaselineCount = len(hdulist[0].data['BASELINE'])
		
		# Close
		hdulist.close()
	
	def getAntennaArray(self):
		"""
		Return an AIPY AntennaArray instance for the array that made the 
		observations contained here.
		"""
		
		# Get the date of observations
		refJD = astro.unix_to_utcjd(timegm(self.dateObs.timetuple()))
		
		# Return
		return simVis.buildSimArray(self.station, self.antennas, self.freq/1e9, jd=refJD)
		
	def getObserver(self):
		"""
		Return a ephem.Observer instances for the array described in the file.
		"""
		
		return self.station.getObserver()
		
	def getDataSet(self, set, includeAuto=False, sort=True, uvMin=0, uvMax=numpy.inf):
		"""
		Return a baseline sorted data dictionary for the specified data set.  
		By default this excludes the autocorrelations.  To include 
		autocorrelations set the value of 'includeAuto' to True.  Setting the
		'sort' keyword to False will disable the baseline sorting.  Optionally,
		baselines with lengths between uvMin and uvMax can only be returned.

		.. note::
			uvMin and uvMax should be specified in lambda
		"""
		
		# Open the file
		hdulist = pyfits.open(self.filename)
		uvData = hdulist[0]
		
		# We need this a lot...
		nPol = len(self.pols)
		
		# Define the dictionary to return
		dataDict = self._createEmptyDataDict()

		# Set the source ID to look for (this is LWA specific)
		sourceID = range(set,set+1)
		
		# Loop over data rows
		found = False
		for row in uvData.data:
			# Load it.
			bl = int(row['BASELINE'])
			if bl >= 65536:
				a1 = int((bl - 65536) / 2048)
				a2 = int((bl - 65536) % 2048)
			else:
				a1 = int(bl / 256)
				a2 = int(bl % 256)
			i = self.standMap[a1]
			j = self.standMap[a2]
			
			if i == j and not includeAuto:
				## Skip auto-correlations
				continue
			ri = numpy.where(self.stands == i)[0][0]
			rj = numpy.where(self.stands == j)[0][0]
			
			uvw = numpy.array([row['UU'], row['VV'], row['WW']])
			
			jd = row['DATE']
			uvw = numpy.array([numpy.dot(uvw[0], self.freq), numpy.dot(uvw[1], self.freq), numpy.dot(uvw[2], self.freq)])
			flux = row['DATA'][0,0,:,:,:]
			
			for c,p in enumerate(self.pols):
				name = NumericStokes[p]
				if len(name) == 2:
					name = name.lower()
					
				vis = numpy.zeros(flux.shape[0], dtype=numpy.complex64)
				vis.real = flux[:,c,0]
				vis.imag = flux[:,c,1]
				wgt = numpy.ones(vis.size)
				msk = numpy.zeros(vis.size, dtype=numpy.int16)
			
				dataDict['uvw'][name].append( uvw ) 
				dataDict['vis'][name].append( vis )
				dataDict['wgt'][name].append( wgt )
				dataDict['msk'][name].append( msk )
				dataDict['bls'][name].append( (ri,rj) )
				dataDict['jd' ][name].append( jd )
			found = True
		
		# Close
		hdulist.close()
		
		# Make sure we found something
		if not found:
			raise RuntimeError("Cannot find baseline set %i in UVFITS file", set)
		
		# Sort
		if sort:
			sortDataDict(dataDict)
			
		# Prune
		if uvMin != 0 or uvMax != numpy.inf:
			dataDict = pruneBaselineRange(dataDict, uvMin=uvMin, uvMax=uvMax)
			
		# Return
		return dataDict


try:
	from pyrap.tables import *
	
	# Stokes codes for CASA Measurement Sets
	NumericStokesMS = {1:'I', 2:'Q', 3:'U', 4:'V', 
				    9:'XX', 10:'XY', 11:'YX', 12:'YY'}
	
	class CorrelatedDataMS(object):
		"""
		Class to make accessing information about a MS easy.  This wraps 
		all of the "messy" machinery needed to extract both the metadata and data 
		from the file and return them as common LSL objects.
		
		This class has three main attributes to interact with:
		* getAntennaArray - Return a :class:`lsl.sim.vim.AntennaArray` instance
		that represents the array where the data was obtained.  This is useful
		for simulation proposes and computing source positions.
		* getObserver - Return a ephem.Observer instance representing the array
		* getDataSet - Return a data dictionary of all baselines for a given set
		of observations
		
		The class also includes a variety of useful metadata attributes:
		* pols - Numpy array of polarization product codes
		* freq - Numpy array of frequency channels in Hz
		* station - LSL :class:`lsl.common.stations.LWAStation` instance for the
		array
		* dateObs - Datetime object for the reference date of the FIT IDI file
		* antennas - List of :class:`lsl.common.stations.Antenna` instances
		
		.. note::
			The CorrelatedDataMS.antennas attribute should be used over 
			CorrelatedDataMS.station.getAntennas() since the mapping in the MS
			may not be the same as the digitizer order.
		"""
		
		def _createEmptyDataDict(self):
			"""
			Create an empty data dictionary that is appropriate for the current file.
			"""
			
			dataDict = {}
			dataDict['freq'] = self.freq 
			dataDict['isMasked'] = False
			for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
				dataDict[key] = {}
				
			for p in self.pols:
				name = NumericStokesMS[p]
				if len(name) == 2:
					name = name.lower()
				for key in ('uvw', 'vis', 'wgt', 'msk', 'bls', 'jd'):
					dataDict[key][name] = []
					
			return dataDict
		
		def __init__(self, filename):
			"""
			Initialize a new CorrelatedData instance from a MS and fill 
			in the metadata.
			"""
			
			self.filename = filename
			
			# Open the various tables that we need
			data = table(self.filename, ack=False)
			try:
				ants = table(os.path.join(self.filename, 'ANTENNA'), ack=False)
			except:
				raise RuntimeError("Cannot find table 'ANTENNA' in '%s'" % self.filename)
			try:
				pols = table(os.path.join(self.filename, 'POLARIZATION'), ack=False)
			except:
				raise RuntimeError("Cannot find table 'POLARIZATION' in '%s'" % self.filename)
			try:
				obs = table(os.path.join(self.filename, 'OBSERVATION'), ack=False)
			except:
				raise RuntimeError("Cannot find table 'OBSERVATION' in '%s'" % self.filename)
			try:
				spw = table(os.path.join(self.filename, 'SPECTRAL_WINDOW'), ack=False)
			except:
				raise RuntimeError("Cannot find table 'SPECTRAL_WINDOW' in '%s'" % self.filename)
				
			# Station/telescope information
			self.telescope = obs.col('TELESCOPE_NAME')[0]
			if self.telescope == 'LWA-1' or self.telescope == 'LWA1':
				self.station = stations.lwa1
			elif self.telescope == 'LWA-2' or self.telescope == 'LWA2':
				self.station = stations.lwa2
			else:
				## Get latitude and longitude for all antennas
				lat = numpy.array([], dtype=numpy.float64)
				lng = numpy.array([], dtype=numpy.float64)
				elv = numpy.array([], dtype=numpy.float64)
				for row in ants.col('POSITION'):
					la,ln,el = stations.ecef2geo(*row)
					lat = numpy.append(lat, la*180/numpy.pi)
					lng = numpy.append(lng, ln*180/numpy.pi)
					elv = numpy.append(elv, el)
					
				## Estimate the center of the station
				sLat = robust.mean(lat)
				sLng = robust.mean(lng)
				sElv = robust.mean(elv)
				
				## Build a preliminayr represenation of the station
				self.station = stations.LWAStation(ants.col('STATION')[0], sLat, sLng, sElv)
				
				## Fill in the antennas instances
				antennas = []
				for i in xrange(lat.size):
					enz = self.station.getENZOffset((lat[i], lng[i], elv[i]))
					sid = int(ants.col('NAME')[i].translate(None, string.letters))
					
					stand = stations.Stand(sid, *enz)
					antennas.append( stations.Antenna(2*(stand.id-1)-1, stand=stand) )
				self.station.antennas = antennas
				
			# Antennas
			self.standMap = {}
			self.stands = []
			for i,noact in enumerate(ants.col('NAME')):
				noact = int(noact[3:])
				self.standMap[i] = noact
				self.stands.append(noact)
				
			self.antennaMap = {}
			self.antennas = []
			if self.station is not None:
				for ant in self.station.getAntennas():
					if ant.stand.id in self.stands and ant.pol == 0:
						self.antennas.append(ant)
						self.antennaMap[ant.stand.id] = ant
			
			# Polarization and frequency
			self.pols = pols.col('CORR_TYPE')[0]
			self.freq  = numpy.array( spw.col('CHAN_FREQ')[0] )
			
			# Total baseline count
			self.totalBaselineCount = data.nrows()
			
			# Data set times
			self._times = []
			for t in data.col('TIME'):
				if t not in self._times:
					self._times.append(t)
			jd = self._times[0] / 3600.0 / 24.0 + astro.MJD_OFFSET
			self.dateObs = pytz.UTC.localize(datetime.utcfromtimestamp(astro.utcjd_to_unix(jd)))
			
			# Close
			data.close()
			ants.close()
			pols.close()
			obs.close()
			spw.close()
		
		def getAntennaArray(self):
			"""
			Return an AIPY AntennaArray instance for the array that made the 
			observations contained here.
			"""
			
			# Get the date of observations
			refJD = astro.unix_to_utcjd(timegm(self.dateObs.timetuple()))
			
			# Return
			return simVis.buildSimArray(self.station, self.antennas, self.freq/1e9, jd=refJD)
			
		def getObserver(self):
			"""
			Return a ephem.Observer instances for the array described in the file.
			"""
			
			return self.station.getObserver()
			
		def getDataSet(self, set, includeAuto=False, sort=True, uvMin=0, uvMax=numpy.inf):
			"""
			Return a baseline sorted data dictionary for the specified data set.  
			By default this excludes the autocorrelations.  To include 
			autocorrelations set the value of 'includeAuto' to True.  Setting the
			'sort' keyword to False will disable the baseline sorting.  Optionally,
			baselines with lengths between uvMin and uvMax can only be returned.

			.. note::
				uvMin and uvMax should be specified in lambda
			"""
			
			# Open the data table
			data = table(self.filename, ack=False)
			
			# We need this a lot...
			nPol = len(self.pols)
			
			# Define the dictionary to return
			dataDict = self._createEmptyDataDict()

			# Load in something we can iterate over
			uvw  = data.col('UVW')
			ant1 = data.col('ANTENNA1')
			ant2 = data.col('ANTENNA2')
			vis  = data.col('DATA')
			time = data.col('TIME')
			
			# Set the time to look for
			targetTime = self._times[set-1]
			
			# Loop over data rows
			found = False
			for u,a1,a2,v,t in zip(uvw, ant1, ant2, vis, time):
				if t != targetTime:
					continue
				found = True

				if a1 == a2 and not includeAuto:
					## Skip auto-correlations
					continue
				r1 = numpy.where(self.stands == a1)
				r2 = numpy.where(self.stands == a2)
				
				jd = t / 3600.0 / 24.0 + astro.MJD_OFFSET
				u = numpy.array(u)
				v = numpy.array(v)
				w = numpy.ones(v.shape)
				m = numpy.zeros(v.shape, dtype=numpy.int16)
				
				u2 = numpy.zeros((u.size, v.shape[0]), dtype=u.dtype)
				for c in xrange(v.shape[0]):
					u2[:,c] = u * (self.freq[c] / vLight)
					
				for c,p in enumerate(self.pols):
					name = NumericStokesMS[p]
					if len(name) == 2:
						name = name.lower()
						
					dataDict['uvw'][name].append( u2 ) 
					dataDict['vis'][name].append( v[:,c] )
					dataDict['wgt'][name].append( w[:,c] )
					dataDict['msk'][name].append( m[:,c] )
					dataDict['bls'][name].append( (r1,r2) )
					dataDict['jd' ][name].append( jd )
			data.close()
			
			# Make sure we found something
			if not found:
				raise RuntimeError("Cannot find baseline set %i in MS", set)
			
			# Sort
			if sort:
				sortDataDict(dataDict)
				
			# Prune
			if uvMin != 0 or uvMax != numpy.inf:
				dataDict = pruneBaselineRange(dataDict, uvMin=uvMin, uvMax=uvMax)
				
			# Return
			return dataDict
			
except ImportError:
	import warnings
	warnings.warn('Cannot import pyrap.tables, MS support disabled')
	
	class CorrelatedMS(object):
		"""
		Class to make accessing information about a MS easy.  This wraps 
		all of the "messy" machinery needed to extract both the metadata and data 
		from the file and return them as common LSL objects.
		"""
		
		def __init__(self, filename):
			raise RuntimeError("Cannot import pyrap.tables, MS support disabled")


def buildGriddedImage(dataDict, MapSize=80, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None):
	"""
	Given a data dictionary, build an aipy.img.ImgW object of gridded uv data 
	which can be used for imaging.  The ImgW object itself is returned by this 
	function to make it more versatile.
	"""
	
	im = aipy.img.ImgW(size=MapSize, res=MapRes, wres=MapWRes)
	
	# Make sure we have the right polarization
	if pol not in dataDict['bls'].keys() and pol.lower() not in dataDict['bls'].keys():
		raise RuntimeError("Data dictionary does not have data for polarization '%s'" % pol)
		
	if chan is not None:
		# Make sure that `chan' is an array by trying to find its length
		try:
			junk = len(chan)
		except TypeError:
			chan = [chan]
			
		# Build up the data using only the specified channels
		uvw = []
		vis = []
		wgt = []
		for (i,j),d,u,w,m in zip(dataDict['bls'][pol], dataDict['vis'][pol], dataDict['uvw'][pol], dataDict['wgt'][pol], dataDict['msk'][pol]):
			u = u[:,chan]
			u.shape = (3, len(chan))
			
			uvw.append(u)
			vis.append(numpy.array([d[chan]]))
			wgt.append(numpy.array([w[chan]]))
	else:
		uvw = dataDict['uvw'][pol]
		vis = dataDict['vis'][pol]
		wgt = dataDict['wgt'][pol]
		
	uvw = numpy.concatenate(uvw, axis=1)
	vis = numpy.concatenate(vis)
	wgt = numpy.concatenate(wgt)
	
	uvw, vis, wgt = im.append_hermitian(uvw, vis, wgts=wgt)
	im.put(uvw, vis, wgts=wgt)
	
	return im
