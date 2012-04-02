# -*- coding: utf-8 -*-

"""
Module to support imaging correlated data.  This module provides utilities to read FITS IDI files
into data dictionaries (as described in :mod:`lsl.sim.vis`) and build AIPY ImgW instances from the 
data dictionaries.  Also included is a utility to sort data dictionaries by baselines.

.. versionadded:: 0.5.0
"""

import aipy
import pytz
import numpy
import pyfits
from calendar import timegm
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.sim import vis as simVis
from lsl.writer.fitsidi import NumericStokes

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['baselineOrder', 'sortDataDict', 'CorrelatedData', 'buildGriddedImage', '__version__', '__revision__', '__all__']


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


class CorrelatedData(object):
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
		Initialize a new CorrelatedData instance from a FITS IDI file and 
		fill in the metadata.
		"""
		
		self.filename = filename
		
		# Open the file, check if it looks like FITS IDI, and pull out the UV_DATA table
		hdulist = pyfits.open(self.filename)
		tbls = [i.header['EXTNAME'] for i in hdulist[1:]]
		for tbl in ('ARRAY_GEOMETRY', 'FREQUENCY', 'ANTENNA', 'BANDPASS', 'SOURCE', 'UV_DATA'):
			if tbl not in tbls:
				raise RuntimeError("Cannot find table '%s' in '%s'" % (tbl, self.filename))
		
		uvData = hdulist['UV_DATA']
		
		# Station/telescope information
		self.telescope = hdulist[0].header['TELESCOP']
		self.dateObs = pytz.UTC.localize(datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S"))
		if self.telescope == 'LWA-1':
			self.station = stations.lwa1
		elif self.telescope == 'LWA-2':
			self.station = stations.lwa2
		else:
			self.station = None
		
		# Antennas
		try:
			mapper = hdulist['NOSTA_MAPPER']
			
			nosta = mapper.data.field('NOSTA')
			noact = mapper.data.field('NOACT')
			anname = mapper.data.field('ANNAME')
		except KeyError:
			ag = hdulist['ARRAY_GEOMETRY']

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
		
	def getDataSet(self, set, includeAuto=False, sort=True):
		"""
		Return a baseline sorted data dictionary for the specified data set.  
		By default this excludes the autocorrelations.  To include 
		autocorrelations set the value of 'includeAuto' to True.  Setting the
		'sort' keyword to False will disable the baseline sorting.
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
				i = self.standMap[(bl >> 8) & 255]
				j = self.standMap[bl & 255]
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
		
		# Sort and return
		if sort:
			return sortDataDict(dataDict)
		else:
			return dataDict


def buildGriddedImage(dataDict, MapSize=80, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None):
	"""
	Given a data dictionary, build an aipy.img.ImgW object of gridded uv data 
	which can be used for imaging.  The ImgW object itself is returned by this 
	function to make it more versatile.
	"""

	im = aipy.img.ImgW(size=MapSize, res=MapRes, wres=MapWRes)

	if chan is not None:
		# Make sure that `chan' is an array by trying to find its length
		try:
			junk = len(chan)
		except TypeError:
			chan = [chan]
			
		# Make sure we have the right polarization
		if pol not in dataDict['bls'].keys() and pol.lower() not in dataDict['bls'].keys():
			raise RuntimeError("Data dictionary does not have data for polarization '%s'" % pol)

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