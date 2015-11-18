# -*- coding: utf-8 -*-

"""
Module to support imaging correlated data.  This module provides utilities to 
read FITS IDI files into data dictionaries (as described in 
:mod:`lsl.sim.vis`) and build AIPY ImgW instances from the data dictionaries.
Also included is a utility to sort data dictionaries by baselines.

.. versionadded:: 0.5.0

.. versionchanged:: 1.0.0
	Added support for UVFITS files and CASA measurement sets

.. versionchanged:: 1.0.1
	Added the plotGriddedImage() function
	
.. versionchanged:: 1.1.0
	Added the getImageRADec() and getImageAzEl() functions to complement
	plotGriddedImage() and make it easier to work with phase centers 
	that are not at zenith.  Added in the ImgWPlus class to add support
	for imaging weighting and tapering.
"""

import os
import re
import sys
import aipy
import ephem
import numpy
import pyfits
import string
try:
	import cStringIO as StringIO
except ImportError:
	import StringIO
from calendar import timegm
from datetime import datetime
from operator import itemgetter

from lsl import astro
from lsl.common import stations
from lsl.sim import vis as simVis
from lsl.writer.fitsidi import NumericStokes
from lsl.common.constants import c as vLight

try:
	import pyfftw
	
	# Enable the PyFFTW cache
	if not pyfftw.interfaces.cache.is_enabled():
		pyfftw.interfaces.cache.enable()
		pyfftw.interfaces.cache.set_keepalive_time(60)
		
	fft2Function = lambda x: pyfftw.interfaces.numpy_fft.fft2(x)
	ifft2Function = lambda x: pyfftw.interfaces.numpy_fft.ifft2(x)

except ImportError:
	fft2Function = numpy.fft.fft2
	ifft2Function = numpy.fft.ifft2

__version__ = '0.7'
__revision__ = '$Rev$'
__all__ = ['baselineOrder', 'sortDataDict', 'pruneBaselineRange', 'rephaseData', 'CorrelatedData', 
		 'CorrelatedDataIDI', 'CorrelatedDataUV', 'CorrelatedDataMS', 'ImgWPlus', 'buildGriddedImage', 
		 'plotGriddedImage', 'getImageRADec', 'getImageAzEl', '__version__', '__revision__', '__all__']


# Regular expression for trying to get the stand number out of an antenna
# name
_annameRE = re.compile('^.*?(?P<id>\d{1,3})$')


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
		for pol in ['xx', 'yy', 'rr', 'll', 'I']:
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
				crd = aa.gen_uvw(j, i, src=newPhaseCenter)[:,0,:]
				d = aa.unphs2src(d, currentPhaseCenter, j, i)
				d = aa.phs2src(d, newPhaseCenter, j, i)
			except aipy.phs.PointingError:
				raise RuntimeError("Rephasing center is below the horizon")
				
			### Save
			if isMasked:
				crd = crd.compress(numpy.logical_not(msk), axis=2)
			vis = d.compress(numpy.logical_not(msk))
			wgt = numpy.ones_like(vis) * len(vis)

			dataDict2['bls'][p].append( (i,j) )
			dataDict2['uvw'][p].append( crd )
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
	except IOError as e:
		raise e
	except:
		pass
		
	# UVFITS
	try:
		return CorrelatedDataUV(filename)
	except IOError as e:
		raise e
	except:
		pass
		
	# Measurment Set
	try:
		return CorrelatedDataMS(filename)
	except IOError as e:
		raise e
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
		for tbl in ('ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'SOURCE', 'UV_DATA'):
			if tbl not in tbls:
				raise RuntimeError("Cannot find table '%s' in '%s'" % (tbl, self.filename))
		
		self.extended = False
		try:
			if hdulist[0].header['LWATYPE'] == 'IDI-EXTENDED-ZA':
				self.extended = True
		except KeyError:
			## Catch for LEDA64-NM data
			pass
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
		try:
			self.telescope = hdulist[0].header['TELESCOP']
			self.dateObs = datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
		except ValueError:
			## Catch for DiFX FITS-IDI files
			self.dateObs = datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%d")
		except KeyError:
			## Catch for LEDA64-NM data
			self.telescope = uvData.header['TELESCOP']
			self.dateObs = datetime.strptime(uvData.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
			
		## Extract the site position
		geo = numpy.array([ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ']])
		site = stations.ecef2geo(*geo)
		
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
			noact = numpy.array(noact2)
			
		## Create the ECI -> topocentric transform
		lat  = site[0]
		ecii = numpy.array([[ 0.0,            1.0, 0.0           ],
						[-numpy.sin(lat), 0.0, numpy.cos(lat)],
						[ numpy.cos(lat), 0.0, numpy.sin(lat)]])
						
		## Build up the list of antennas
		antennas = []
		for line,act in zip(ag.data, noact):
			enz = numpy.dot(ecii, line['STABXYZ'])
			
			stand = stations.Stand(act, *enz)
			antennas.append(stations.Antenna(2*(stand.id-1)+1, stand=stand, pol=0))
			
		## Build up the station
		self.station = stations.LWAStation(ag.header['ARRNAM'], site[0]*180/numpy.pi, site[1]*180/numpy.pi, site[2], antennas=antennas)
		
		self.standMap = {}
		self.stands = []
		for sta, act in zip(nosta, noact):
			self.standMap[sta] = act
			self.stands.append(act)
			
		self.antennaMap = {}
		self.antennas = []
		for stand in self.stands:
			for ant in self.station.getAntennas():
				if ant.stand.id == stand and ant.pol == 0:
					self.antennas.append(ant)
					self.antennaMap[ant.stand.id] = ant
					break
					
		# Polarization and frequency
		self.pols  = numpy.arange(1, uvData.header['MAXIS2']+1) - uvData.header['CRPIX2']
		self.pols *= uvData.header['CDELT2'] 
		self.pols += uvData.header['CRVAL2']
		self.freq  = numpy.arange(1, uvData.header['NO_CHAN']+1, dtype=numpy.float64) - uvData.header['REF_PIXL']
		self.freq *= uvData.header['CHAN_BW']
		self.freq += uvData.header['REF_FREQ']
		
		# Total baseline count
		self.totalBaselineCount = len(uvData.data['BASELINE'])
		
		# Integration count
		jd = uvData.data['DATE'] + uvData.data['TIME']
		self.integrationCount = len(numpy.unique(jd))
		
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
		
		# Get the date of observations
		refJD = astro.unix_to_utcjd(timegm(self.dateObs.timetuple()))
		
		obs = self.station.getObserver()
		obs.date = refJD - astro.DJD_OFFSET
		return obs
		
	def getDataSet(self, set, includeAuto=False, sort=True, uvMin=0, uvMax=numpy.inf):
		"""
		Return a baseline sorted data dictionary for the specified data set.  
		By default this excludes the autocorrelations.  To include 
		autocorrelations set the value of 'includeAuto' to True.  Setting the
		'sort' keyword to False will disable the baseline sorting.  Optionally,
		baselines with lengths between uvMin and uvMax can only be returned.

		.. note::
			uvMin and uvMax should be specified in lambda
			
		.. versionchanged:: 1.1.0
			'set' can now be either an integer or a list to pull back multiple 
			integrations.
		"""
		
		# Open the file
		hdulist = pyfits.open(self.filename)
		uvData = hdulist['UV_DATA']
		
		# We need this a lot...
		nPol = len(self.pols)
		
		# Define the dictionary to return
		dataDict = self._createEmptyDataDict()

		# Set the source ID to look for (this is LWA specific)
		if type(set) == list:
			sourceID = set
		else:
			sourceID = range(set,set+1)
			
		# Loop over data rows
		found = False
		setCounter = 0
		for row in uvData.data:
			if setCounter == 0:
				blCheck = row['baseline']
			if row['BASELINE'] == blCheck:
				setCounter += 1
				
			# If we've passed the correct set
			if setCounter > max(sourceID):
				break
				
			# If we are on the right set...
			if setCounter in sourceID:
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
				
				try:
					uvw = numpy.array([row['UU'], row['VV'], row['WW']])
				except KeyError:
					## Catch for DiFX FITS-IDI files that call them **---SIN
					uvw = numpy.array([row['UU---SIN'], row['VV---SIN'], row['WW---SIN']])
					
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
		Initialize a new CorrelatedDataUV instance from a UVFITS file and 
		fill in the metadata.
		"""
		
		self.filename = filename
		
		# Open the various tables that we need
		hdulist = pyfits.open(filename)
		
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
			self.dateObs = datetime.strptime(dt, "%Y-%m-%dT%H:%M:%S")
		except ValueError:
			## Catch for AIPS UVFITS files which only have a date set
			self.dateObs = datetime.strptime(dt, "%Y-%m-%d")
			
		## Extract the site position
		geo = numpy.array([ag.header['ARRAYX'], ag.header['ARRAYY'], ag.header['ARRAYZ']])
		site = stations.ecef2geo(*geo)
		
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
			noact = numpy.array(noact2)
			
		## Create the ECI -> topocentric transform
		lat  = site[0]
		ecii = numpy.array([[ 0.0,            1.0, 0.0           ],
						[-numpy.sin(lat), 0.0, numpy.cos(lat)],
						[ numpy.cos(lat), 0.0, numpy.sin(lat)]])
						
		## Build up the list of antennas
		antennas = []
		for line,act in zip(ag.data, noact):
			enz = numpy.dot(ecii, line['STABXYZ'])
			
			stand = stations.Stand(act, *enz)
			antennas.append(stations.Antenna(2*(stand.id-1)+1, stand=stand, pol=0))
			
		## Build up the station
		self.station = stations.LWAStation(ag.header['ARRNAM'], site[0]*180/numpy.pi, site[1]*180/numpy.pi, site[2], antennas=antennas)
		
		self.standMap = {}
		self.stands = []
		for nosta, noact in zip(nosta, noact):
			self.standMap[nosta] = noact
			self.stands.append(noact)
			
		self.antennaMap = {}
		self.antennas = []
		for stand in self.stands:
			for ant in self.station.getAntennas():
				if ant.stand.id == stand and ant.pol == 0:
					self.antennas.append(ant)
					self.antennaMap[ant.stand.id] = ant
					break
					
		# Polarization and frequency
		self.pols  = numpy.arange(1, uvData.header['NAXIS3']+1) - uvData.header['CRPIX3']
		self.pols *= uvData.header['CDELT3'] 
		self.pols += uvData.header['CRVAL3']
		self.freq  = numpy.arange(1, uvData.header['NAXIS4']+1, dtype=numpy.float64) - uvData.header['CRPIX4']
		self.freq *= uvData.header['CDELT4']
		self.freq += uvData.header['CRVAL4']
		
		# Total baseline count
		self.totalBaselineCount = len(hdulist[0].data['BASELINE'])
		
		# Integration count
		jd = hdulist[0].data['DATE']
		self.integrationCount = len(numpy.unique(jd))
		
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
		
		# Get the date of observations
		refJD = astro.unix_to_utcjd(timegm(self.dateObs.timetuple()))
		
		obs = self.station.getObserver()
		obs.date = refJD - astro.DJD_OFFSET
		return obs
		
	def getDataSet(self, set, includeAuto=False, sort=True, uvMin=0, uvMax=numpy.inf):
		"""
		Return a baseline sorted data dictionary for the specified data set.  
		By default this excludes the autocorrelations.  To include 
		autocorrelations set the value of 'includeAuto' to True.  Setting the
		'sort' keyword to False will disable the baseline sorting.  Optionally,
		baselines with lengths between uvMin and uvMax can only be returned.

		.. note::
			uvMin and uvMax should be specified in lambda
			
		.. versionchanged:: 1.1.0
			'set' can now be either an integer or a list to pull back multiple 
			integrations.
		"""
		
		# Open the file
		hdulist = pyfits.open(self.filename)
		uvData = hdulist[0]
		
		# We need this a lot...
		nPol = len(self.pols)
		
		# Define the dictionary to return
		dataDict = self._createEmptyDataDict()

		# Set the source ID to look for (this is LWA specific)
		if type(set) == list:
			sourceID = set
		else:
			sourceID = range(set,set+1)
			
		# Loop over data rows
		found = False
		setCounter = 0
		for row in uvData.data:
			if setCounter == 0:
				blCheck = row['baseline']
			if row['BASELINE'] == blCheck:
				setCounter += 1
				
			# If we've passed the correct set
			if setCounter > max(sourceID):
				break
				
			# If we are on the right set...
			if setCounter in sourceID:
				found = True
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
				
				try:
					uvw = numpy.array([row['UU'], row['VV'], row['WW']])
				except KeyError:
					### Catch for AIPS UVFITS data which calls them **---SIN
					uvw = numpy.array([row['UU---SIN'], row['VV---SIN'], row['WW---SIN']])
					
				jd = row['DATE']
				uvw = numpy.array([numpy.dot(uvw[0], self.freq), numpy.dot(uvw[1], self.freq), numpy.dot(uvw[2], self.freq)])
				if len(row['DATA'].shape) == 6:
					## Fix for AIPS UVFITS data which includes an 'IF' column
					flux = row['DATA'][0,0,0,:,:,:]
				else:
					flux = row['DATA'][0,0,:,:,:]
				if flux.shape[-1] == 3:
					## Catch for AIPS UVFITS which has a third entry in COMPLEX which includes the weight
					wgt = flux[:,:,2]
					flux = flux[:,:,:2]
				else:
					wgt = None
					
				for c,p in enumerate(self.pols):
					name = NumericStokes[p]
					if len(name) == 2:
						name = name.lower()
					
					vis = numpy.zeros(flux.shape[0], dtype=numpy.complex64)
					vis.real = flux[:,c,0]
					vis.imag = flux[:,c,1]
					if wgt is None:
						wgt = numpy.ones(vis.size)
					msk = numpy.zeros(vis.size, dtype=numpy.int16)
					
					dataDict['uvw'][name].append( uvw ) 
					dataDict['vis'][name].append( vis )
					dataDict['wgt'][name].append( wgt )
					dataDict['msk'][name].append( msk )
					dataDict['bls'][name].append( (ri,rj) )
					dataDict['jd' ][name].append( jd )
					
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
			for stand in self.stands:
				for ant in self.station.getAntennas():
					if ant.stand.id == stand and ant.pol == 0:
						self.antennas.append(ant)
						self.antennaMap[ant.stand.id] = ant
						break
						
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
			self.dateObs = datetime.utcfromtimestamp(astro.utcjd_to_unix(jd))
			
			# Integration count
			jd = numpy.array(self._times) / 3600.0 / 24.0 + astro.MJD_OFFSET
			self.integrationCount = len(numpy.unique(jd))
			
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
			
			# Get the date of observations
			refJD = astro.unix_to_utcjd(timegm(self.dateObs.timetuple()))
			
			obs = self.station.getObserver()
			obs.date = refJD - astro.DJD_OFFSET
			return obs
			
		def getDataSet(self, set, includeAuto=False, sort=True, uvMin=0, uvMax=numpy.inf):
			"""
			Return a baseline sorted data dictionary for the specified data set.  
			By default this excludes the autocorrelations.  To include 
			autocorrelations set the value of 'includeAuto' to True.  Setting the
			'sort' keyword to False will disable the baseline sorting.  Optionally,
			baselines with lengths between uvMin and uvMax can only be returned.

			.. note::
				uvMin and uvMax should be specified in lambda
				
			.. versionchanged:: 1.1.0
				'set' can now be either an integer or a list to pull back multiple 
				integrations.
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
	warnings.warn('Cannot import pyrap.tables, MS support disabled', ImportWarning)
	
	class CorrelatedMS(object):
		"""
		Class to make accessing information about a MS easy.  This wraps 
		all of the "messy" machinery needed to extract both the metadata and data 
		from the file and return them as common LSL objects.
		"""
		
		def __init__(self, filename):
			raise RuntimeError("Cannot import pyrap.tables, MS support disabled")


class ImgWPlus(aipy.img.ImgW):
	"""
	Sub-class of the aipy.img.ImgW class that adds support for different 
	visibility weighting scheme and uv plane tapering.  This class also
	adds in a couple of additional methods that help determine the size of
	the field of view and the pixels near the phase center.
	"""
	
	def getFieldOfView(self):
		"""
		Return the approximate size of the field of view in radians.  The 
		field of view calculate is based off the maximum and minimum values
		of L found for the inverted uv matrix.
		"""
		
		# Get the L and M coordinates
		l,m = self.get_LM()
		
		# Find the maximum and minimum values of L
		lMax = numpy.where( l == l.max() )
		lMin = numpy.where( l == l.min() )
		#print lMax, lMin, l[lMax], l[lMin], m[lMax], m[lMin]
		
		# Convert these locations into topocentric
		xMax, xMin = l.data[lMax], l.data[lMin]
		yMax, yMin = m.data[lMax], m.data[lMin]
		zMax, zMin = numpy.sqrt(1 - xMax**2 - yMax**2), numpy.sqrt(1 - xMin**2 - yMin**2)
		azAltMax = aipy.coord.top2azalt((xMax,yMax,zMax))
		azAltMin = aipy.coord.top2azalt((xMin,yMin,zMin))
		
		# Get the separation between the two
		d = 2*numpy.arcsin( numpy.sqrt( numpy.sin((azAltMax[1]-azAltMin[1])/2)**2+numpy.cos(azAltMax[1])*numpy.cos(azAltMin[1])*numpy.sin((azAltMax[0]-azAltMin[0])/2)**2 ) )
		
		return d.max()
		
	def getPixelSize(self):
		"""
		Return the approximate size of pixels at the phase center in radians.
		The pixel size is averaged over the four pixels that neighboor the 
		phase center.
		"""
		
		# Get the L and M coordinates
		l,m = self.get_LM()
		
		sizes = []
		x0, y0 = l[0,0], m[0,0]
		z0 = numpy.sqrt(1 - x0**2 - y0**2)
		for offX,offY in ((0,1), (1,0), (0,-1), (-1,0)):
			x1, y1 = l[offX,offY], m[offX,offY]
			z1 = numpy.sqrt(1 - x1**2 - y1**2)
			
			# Convert these locations into topocentric
			azAlt0 = aipy.coord.top2azalt((x0,y0,z0))
			azAlt1 = aipy.coord.top2azalt((x1,y1,z1))
			
			# Get the separation between the two
			d = 2*numpy.arcsin( numpy.sqrt( numpy.sin((azAlt1[1]-azAlt0[1])/2)**2+numpy.cos(azAlt0[1])*numpy.cos(azAlt1[1])*numpy.sin((azAlt1[0]-azAlt0[0])/2)**2 ) )
			
			# Save
			sizes.append(d)
		sizes = numpy.array(sizes)
		
		return sizes.mean()
		
	def _gen_img(self, data, center=(0,0), weighting='natural', localFraction=0.5, robust=0.0, taper=(0.0, 0.0)):
		"""
		Return the inverse FFT of the provided data, with the 0,0 point 
		moved to 'center'.  In the images return north is up and east is 
		to the left.
		
		There are a few keywords that control how the image is formed.  
		There are:
		 * weighting - The weighting scheme ('natural', 'uniform', or 
		               'briggs') used on the data;
		 * localFraction - The fraction of the uv grid that is consider 
		                   "local" for the 'uniform' and 'briggs' methods;
		 * robust - The value for the weighting robustness under the 
		            'briggs' method; and
		 * taper - The size of u and v Gaussian tapers at the 30% level.
		"""
		
		# Make sure that we have a valid weighting scheme to use
		if weighting not in ('natural', 'uniform', 'briggs'):
			raise ValueError("Unknown weighting scheme '%s'" % weighting)
			
		# Make sure that we have a valid localFraction value
		if localFraction <= 0 or localFraction > 1:
			raise ValueError("Invalid localFraction value")
			
		# Apply the weighting
		if weighting == 'natural':
			## Natural weighting - we already have it
			pass
		
		elif weighting == 'uniform':
			## Uniform weighting - we need to calculate it
			dens = numpy.abs(self.bm[0])
			size = dens.shape[0]
			
			from scipy.ndimage import uniform_filter
			dens = uniform_filter(dens, size=size*localFraction)
			dens /= dens.max()
			dens[numpy.where( dens < 1e-8 )] = 0
			
			data = data/dens
			data[numpy.where(dens == 0)] = 0.0
			
		elif weighting == 'briggs':
			## Robust weighting - we need to calculate it
			dens = numpy.abs(self.bm[0])
			size = dens.shape[0]
			
			from scipy.ndimage import uniform_filter
			dens = uniform_filter(dens, size=size*localFraction)
			dens /= dens.max()
			dens[numpy.where( dens < 1e-8 )] = 0
			
			f2 = (5*10**-robust)**2 / (dens**2).mean()
			dens = 1.0 / (1.0 + f2/dens)
			data = data/dens*dens.max()
			data[numpy.where(dens == 0)] = 0.0
			
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
				cu = numpy.log(0.3) / taper[0]**2
				taper1 = numpy.exp(cu*u**2)
				
			taper2 = 1.0
			if taper[1] > 0.0:
				cv = numpy.log(0.3) / taper[1]**2
				taper2 = numpy.exp(cv*v**2)
				
			data = data*taper1*taper2
			
		return aipy.img.recenter(ifft2Function(data).real.astype(numpy.float32), center)
		
	def image(self, center=(0,0), weighting='natural', localFraction=0.5, robust=0.0, taper=(0.0, 0.0)):
		"""Return the inverse FFT of the UV matrix, with the 0,0 point moved
		to 'center'.  In the images return north is up and east is 
		to the left.
		
		There are a few keywords that control how the image is formed.  
		There are:
		 * weighting - The weighting scheme ('natural', 'uniform', or 
		               'briggs') used on the data;
		 * localFraction - The fraction of the uv grid that is consider 
		                   "local" for the 'uniform' and 'briggs' methods;
		 * robust - The value for the weighting robustness under the 
		            'briggs' method; and
		 * taper - The size of u and v Gaussian tapers at the 30% level.
		"""
		
		return self._gen_img(self.uv, center=center, weighting=weighting, localFraction=localFraction, robust=robust, taper=taper)
		
	def bm_image(self, center=(0,0), term=None, weighting='natural', localFraction=0.5, robust=0.0, taper=(0.0, 0.0)):
		"""Return the inverse FFT of the sample weightings (for all mf_order
		terms, or the specified term if supplied), with the 0,0 point
		moved to 'center'.  In the images return north is up and east is 
		to the left.
		
		There are a few keywords that control how the image is formed.  
		There are:
		 * weighting - The weighting scheme ('natural', 'uniform', or 
		               'briggs') used on the data;
		 * localFraction - The fraction of the uv grid that is consider 
		                   "local" for the 'uniform' and 'briggs' methods;
		 * robust - The value for the weighting robustness under the 
		            'briggs' method; and
		 * taper - The size of u and v Gaussian tapers at the 30% level.
		"""
		
		if not term is None:
			return self._gen_img(self.bm[term], center=center, weighting=weighting, localFraction=localFraction, robust=robust, taper=taper)
		else:
			return [self._gen_img(b, center=center, weighting=weighting, localFraction=localFraction, robust=robust, taper=taper) for b in self.bm]


def buildGriddedImage(dataDict, MapSize=80, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None, verbose=True):
	"""
	Given a data dictionary, build an aipy.img.ImgW object of gridded uv data 
	which can be used for imaging.  The ImgW object itself is returned by this 
	function to make it more versatile.
	"""
	
	im = ImgWPlus(size=MapSize, res=MapRes, wres=MapWRes)
	
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
	
	if not verbose:
		sys.stdout = StringIO.StringIO()
		
	uvw, vis, wgt = im.append_hermitian(uvw, vis, wgts=wgt)
	im.put(uvw, vis, wgts=wgt)
	
	if not verbose:
		sys.stdout.close()
		sys.stdout = sys.__stdout__
	
	return im


def plotGriddedImage(ax, gimg, shifted=True, interpolation='nearest', **kwargs):
	"""
	Given a blank matplotlib axes instance and a gridded image generated by 
	the buildGriddedImage() function, plot the image on the axes and setup
	the basic coordinate system.  This function returns the NumPy array of
	the image that is plotted.
	
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
		img = numpy.roll(img, imgSize/2, axis=0)
		img = numpy.roll(img, imgSize/2, axis=1)
		
	# Plot
	ax.imshow(img, extent=(1,-1,-1,1), origin='lower', interpolation=interpolation, **kwargs)
	
	return img


def _radec_of(aa, az, alt):
	# az/el -> HA/dec
	HA = numpy.arctan2(numpy.sin(az-numpy.pi), (numpy.cos(az-numpy.pi)*numpy.sin(aa.lat) + numpy.tan(alt)*numpy.cos(aa.lat)))
	dec = numpy.arcsin(numpy.sin(aa.lat)*numpy.sin(alt) - numpy.cos(aa.lat)*numpy.cos(alt)*numpy.cos(az-numpy.pi))
	
	# HA -> RA
	RA = aa.sidereal_time() - HA
	
	# radians -> degrees
	RA = RA * 180.0/numpy.pi
	dec = dec * 180.0/numpy.pi
	
	# RA/dec -> astro.eqn_posn()
	pos = astro.equ_posn(RA, dec)
	
	# Correct for aberration
	pos2 = astro.get_equ_aber(pos, site.date+astro.DJD_OFFSET)
	dRA, dDec = pos2.ra - pos.ra, pos2.dec - pos.dec
	pos.ra = (pos.ra - dRA) % 360.0
	pos.ra %= 360.0
	pos.dec = pos.dec - dDec
	
	# Correct for nutation
	pos2 = astro.get_equ_nut(pos, site.date+astro.DJD_OFFSET)
	dRA, dDec = pos2.ra - pos.ra, pos2.dec - pos.dec
	pos.ra = (pos.ra - dRA) % 360.0
	pos.ra %= 360.0
	pos.dec = pos.dec - dDec
	
	# Precess back to J2000
	pos = astro.get_precession(aa.date+astro.DJD_OFFSET, pos, ephem.J2000+astro.DJD_OFFSET)
	RA, dec = pos.ra, pos.dec
	
	# degrees -> radians
	RA = RA * numpy.pi/180.0
	dec = dec * numpy.pi/180.0
	
	return RA, dec 


def getImageRADec(gimg, aa, phaseCenter='z', shifted=True):
	"""
	Given a gridded image generated by the buildGriddedImage() function
	and an AntennaArray instance, return a two-element tuple containing
	the RA and dec. values (in radians) for each pixel in the image.  
	
	The 'phaseCenter' keyword controls what the phase center of the image 
	is and defaults to zenith.
	
	.. versionadded: 1.1.0
	"""
	
	# Get the phase center
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)
		pcRA, pcDec = phaseCenter.ra, phaseCenter.dec
	else:
		pcRA, pcDec = aa.sidereal_time(), aa.lat
	rotInv = aipy.coord.top2eq_m(0, pcDec)
	
	# Extract the raw topocentric coordinates and convert to equatorial
	top = gimg.get_top()
	oldShape = top[0].shape
	top = (top[0].ravel(), top[1].ravel(), top[2].ravel())
	eq = numpy.dot(rotInv, top)
	eq = (eq[0].reshape(oldShape), eq[1].reshape(oldShape), eq[2].reshape(oldShape))
	
	# Over to RA/Dec
	ra, dec = aipy.coord.eq2radec(eq)
	
	# Correct for the phaseCenter
	ra += pcRA
	ra %= 2*numpy.pi
	
	# Shift, if needed
	if shifted:
		raSize = ra.shape[0]	# should be square
		ra = numpy.roll(ra, raSize/2, axis=0)
		ra = numpy.roll(ra, raSize/2, axis=1)
		
		decSize = dec.shape[0]	# should be square
		dec = numpy.roll(dec, decSize/2, axis=0)
		dec = numpy.roll(dec, decSize/2, axis=1)
		
	# Done
	return ra, dec


def getImageAzEl(gimg, aa, phaseCenter='z', shifted=True):
	"""
	Given a gridded image generated by the buildGriddedImage() function
	and an AntennaArray instance, return a two-element tuple containing
	the azimuth and elevation (altitude), both in radians, for each pixel
	in the image.
	
	The 'phaseCenter' keyword controls what the phase center of the image 
	is and defaults to zenith.
	
	.. versionadded: 1.1.0
	"""
	
	# Get the phase center
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)
		pcRA, pcDec = phaseCenter.ra, phaseCenter.dec
	else:
		pcRA, pcDec = aa.sidereal_time(), aa.lat
	rot = aipy.coord.eq2top_m(0, pcDec)
	
	# Get the RA and dec. coordinates for each pixel
	ra, dec = getImageRADec(gimg, aa, phaseCenter=phaseCenter, shifted=shifted)
	
	# Convert to azimuth and elevation using PyEphem
	bdy = aipy.amp.RadioFixedBody(0, 0)
	
	az, el = ra*0.0, dec*0.0
	for i in xrange(az.shape[0]):
		for j in xrange(az.shape[1]):
			bdy._ra = ra[i,j]
			bdy._dec = dec[i,j]
			bdy.compute(aa)
			az[i,j], el[i,j] = bdy.az, bdy.alt
	
	# Done
	return az, el
