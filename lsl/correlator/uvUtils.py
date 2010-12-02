# -*- coding: utf-8 -*

"""This module stores various functions that are needed for computing UV 
coverage and time delays.  The functions in the module:
  * return the x, y, and z coordinates of a stand or array of stands
  * return the relative x, y, and z offsets between two stands
  * return the cable delays as a function of frequency for a stand
  * compute the u, v, and w coordinates of all baselines defined by an array 
    of stands
  * compute the track through the uv-plane of a collection of baselines as 
    the Earth rotates.
"""

import os
import csv
import math
import numpy

from lsl.common.paths import data as dataPath
from lsl.common.constants import *

__version__ = '0.3'
__revision__ = '$ Revision: 15 $'
__all__ = ['validateStand', 'getXYZ', 'getRelativeXYZ', 'PositionCache', 'cableDelay', 'CableCache', 'signalDelay', 'SignalCache', 'cableAttenuation', 'getBaselines', 'baseline2antenna', 'antenna2baseline', 'computeUVW', 'computeUVTrack', 'uvUtilsError', '__version__', '__revision__', '__all__']


class uvUtilsError(Exception):
	"""Base exception class for this module"""

	def __init__(self, strerror):
		self.strerror = strerror
	
	def __str__(self):
		return "%s" % self.strerror


def validateStand(stand, max=258):
	"""Make sure that the provided stand number is valid.  If it is not between
	1 and max (default of 258 for LWA-1), raise an error."""

	if stand < 1 or stand > max:
			raise uvUtilsError('Stand #%i is out of range (1-%i)' % (stand, max))

	return True


def _loadPositionData(filename='lwa1-positions.csv'):
	"""Private function to load in the stand location data CSV file.  Someday 
	this should be absorbed into this script so that there is not need to keep 
	yet another file floating around."""

	# lwa1-asbuilt.csv contains the stand (x,y,z) coordinates for all 257 stands +
	# the outlier (#258).  The data come from the "LWA-1 Antenna Position and 
	# Cable Data" memo, version 1.
	try:
		csvFH = open(os.path.join(dataPath, filename), 'r')
	except IOError, (errno, strerror):
		raise uvUtilsError("%s: stand position file '%s'" % (strerror, filename))
	csvData = csv.reader(csvFH, delimiter=',', quotechar='"')
	data = numpy.empty((258,3))
	
	i = 0
	for row in csvData:
		if row[0].find('ID') == -1:
			data[i,0] = float(row[1])
			data[i,1] = float(row[2])
			data[i,2] = float(row[3])
			i = i + 1

	csvFH.close()

	return data


def getXYZ(stands):
	"""Function to return a numpy array of the stand's x, y, and z coordinates
	(in meters) from the center pole.  If the coordinates of more than one 
	stand are needed, getXYZ can be called with a numpy array of stand 
	numbers.

	.. versionchanged:: 0.3
		Changed the numbering of the outlier from -1 to a more consistent 258
	"""

	try:
		junk = len(stands)
	except TypeError:
		stands = numpy.array([stands])

	standsXYZ = _loadPositionData()

	out = numpy.zeros((stands.shape[0],3))
	for i in range(stands.shape[0]):
		validateStand(stands[i])
		out[i,:] = standsXYZ[stands[i]-1,:]

	return out


def getRelativeXYZ(stand1, stand2):
	"""Function to return the relative coordinate difference between two 
	stands.  The output is a three-element numpy array."""

	xyzs = getXYZ(numpy.array([stand1, stand2]))
	xyz1 = numpy.squeeze(xyzs[0,:])
	xyz2 = numpy.squeeze(xyzs[1,:])

	deltaXYZ = xyz2 - xyz1

	return deltaXYZ


class PositionCache(object):
	"""PositionCache is a safe alternative to calling uvUtils.getXYZ or
	uvUtils.getRealtiveXYZ over and over again.  PositionCache loads the 
	stand positions from the included CSV file once and then uses that 
	stored information for all subsequent calls.

	.. warning::
		This assumes LWA-1 as the stations.  This will need to be changed some day.
	"""

	def __init__(self):
		"""Initialize the cache by loading in all stand positions."""

		self.standsXYZ = _loadPositionData()
		
	def __isValid(self, stand):
		if stand < 1 or stand > 258:
			raise uvUtilsError('Stand #%i is out of range (1-258)' % stand)
		
	def getXYZ(self, stands):
		"""Function to return a numpy array of the stand's x, y, and z coordinates
		(in meters) from the center pole.  If the coordinates of more than one 
		stand are needed, getXYZ can be called with a numpy array of stand 
		numbers."""

		try:
			junk = len(stands)
		except TypeError:
			stands = numpy.array([stands])
	
		out = numpy.zeros((stands.shape[0],3))
		for i in range(stands.shape[0]):
			self.__isValid(stands[i])
			out[i,:] = self.standsXYZ[stands[i]-1,:]
	
		return out
			
	def getRelativeXYZ(self, stand1, stand2):
		"""Function to return the relative coordinate difference between two 
		stands.  The output is a three-element numpy array."""

		xyzs = self.getXYZ(numpy.array([stand1, stand2]))
		xyz1 = numpy.squeeze(xyzs[0,:])
		xyz2 = numpy.squeeze(xyzs[1,:])
	
		deltaXYZ = xyz2 - xyz1
	
		return deltaXYZ
	
	def getZenithDelay(self, stand1, stand2):
		"""Return the geometrical delay in seconds for zenith as viewed by
		a baseline between stand1 and stand2."""

		import aipy
		
		b = self.getRelativeXYZ(stand1, stand2)
		top = aipy.coord.azalt2top(numpy.array([[numpy.pi/4],[numpy.pi/2]]))
		geoDelay = numpy.dot(b, top[:,0]) / c
		
		return geoDelay


def _loadDelayData(filename='lwa1-cables.csv'):
	"""Private function to load in the cable data CSV file.  Someday that file 
	should be absorbed into this script so that there is not need to keep another 
	file floating around."""

	# lwa1-cables.csv contains the stand cable lengths in meters for 256 of the
	# LWA-1 stands.  The data come from the "LWA-1 Antenna Position and Cable 
	# Data" memo, version 1.
	try:
		csvFH = open(os.path.join(dataPath, filename), 'r')	
	except IOError, (errno, strerror):
		raise uvUtilsError("%s: cable length file '%s'" % (strerror, filename))
	csvData = csv.reader(csvFH, delimiter=',', quotechar='"')
	data = numpy.empty(258)
	
	i = 0
	for row in csvData:
		if row[0].find('ID') == -1:
			data[i] = float(row[1])
			i = i + 1

	csvFH.close()

	return data


def cableDelay(stand, freq):
	"""For a given stands, return a numpy array of the cable delay in seconds 
	for a specific frequency (in Hz).  If delays for more than one frequency
	are needed, the frequencies can be passed in as a numpy array.
	"""

	# Stands start at 1, indices do not
	validateStand(stand)

	# Try this so that freq can be either a scalar or a list
	try:
		junk = len(freq)
	except TypeError:
		freq = numpy.array([freq])

	# Catch the outlier and load cable and delay information into holder variables.
	# This makes it easier to do one common set of operations with different cables.
	cableLengths = _loadDelayData()
	cableLength = cableLengths[stand-1]
	velocityFactor = 0.83
	tauD0 = 2.4e-9 # dispersion delay in seconds
	if stand == 258:
		velocityFactor = 0.85
		tauD0 = 4.6e-9	# dispersion delay in seconds

	# First, calculate the delay from cable length assuming a velocity factor.
	delay = cableLength / (velocityFactor * c)
	# Next, calculate the frequency-dependent dispersion delay and add it to 
	# the delay calculated above.  The form of the dispersion delay comes from
	# Section 4, Equation 25 of Steve Ellingson's "Dispersion in Coaxial Cables" 
	# memo.
	dispDelay = tauD0*(cableLength/100.0)/numpy.sqrt(freq/10.0e6)
	delays = delay + dispDelay

	return delays


class CableCache(object):
	"""CableCache is a safe alternative to calling uvUtils.cableDelay over 
	and over again.  CableCache loads the cable length file once and then
	uses that stored information for all subsequent calls.

	.. warning::
		This assumes LWA-1 as the stations.  This will need to be changed some day.
	"""

	def __init__(self, freq, applyDispersion=True):
		"""Initialize the cache by loading in all stand cable lengths and begin
		filling in the attributes."""

		self.cableLengths = _loadDelayData()
		try:
			junk = len(freq)
		except TypeError:
			freq = numpy.array([freq])
		self.freq = freq
		self.disp = 1.0/numpy.sqrt(numpy.abs(self.freq)/10.0e6)
		self.applyDispersion = applyDispersion
	
	def updateFreq(self, freq):
		"""Update the freq attribute of the cache to change which frequencies the
		delays are calculated for."""

		self.freq = freq
		self.disp = 1.0/numpy.sqrt(numpy.abs(self.freq)/10.0e6)

	def updateApplyDispersion(self, applyDispersion):
		"""Update the applyDispersion attribute of the cache to turn the cable
		dispersion on and off."""

		self.applyDispersion = applyDispersion

	def __isValid(self, stand):
		if stand < 1 or stand > 258:
			raise uvUtilsError('Stand #%i is out of range (1-258)' % stand)

	def cableDelay(self, stand, freq=None):
		"""For a given stands, return a numpy array of the cable delay in seconds 
		for a specific frequency (in Hz).  If delays for more than one frequency
		are needed, the frequencies can be passed in as a numpy array."""

		self.__isValid(stand)
	
		# Catch the outlier and load cable and delay information into holder variables.
		# This makes it easier to do one common set of operations with different cables.
		cableLength = self.cableLengths[stand-1]
		velocityFactor = 0.83
		tauD0 = 2.4e-9 # dispersion delay in seconds
		if stand == 258:
			velocityFactor = 0.85
			tauD0 = 4.6e-9	# dispersion delay in seconds
	
		# First, calculate the delay from cable length assuming a velocity factor.
		delay = cableLength / (velocityFactor * c)
		# Next, calculate the frequency-dependent dispersion delay and add it to 
		# the delay calculated above.  The form of the dispersion delay comes from
		# Section 4, Equation 25 of Steve Ellingson's "Dispersion in Coaxial Cables" 
		# memo.
		if freq is None:
			dispDelay = tauD0*(cableLength/100.0)*self.disp
		else:
			try:
				junk = len(freq)
			except TypeError:
				freq = numpy.array([freq])
			dispDelay = tauD0*(cableLength/100.0)/numpy.sqrt(numpy.abs(freq)/10.0e6)

		if not self.applyDispersion:
			dispDelay *= 0.0
				
		delays = delay + dispDelay

		return delays


def signalDelay(stand, freq, cache=None):
	"""Experimental function that wraps the cable delay for a stand along with 
	any other delays determined by phase fitting.  The list of stands that can be 
	corrected is incomplete.  Similar to cableDelay, a numpy array is returned.

	.. note::
		Currently all additional delays are set to 0.

	.. warning::
		This assumes LWA-1 as the stations.  This will need to be changed some day.
	"""

	# Additional delays in ns found from phase fitting.  There are currently no
	# other delays added in and signalDelay returns the same values as cableDelay.
	addDelay = {}
	#addDelay = {4: 0.0, -1: -208.6, 205: 96.5, 158: 197.2, 9: 249.8, 246: -153.4, 168: 372.4, 69: 124.4, 14: 314.9, 80: 339.7, 118: 401.3, 254: -81.1, 34: 239.7, 38: 292.5, 181: -123.6, 67: 79.2, 183: -108.6, 206: -118.2, 174: 467.5, 153: -93.2}

	# Stands start at 1, indices do not
	validateStand(stand)

	# Get the delays from the cable model and add in any additional delays found 
	# in the addDelay dictionary.
	if cache is None:
		delay = cableDelay(stand, freq)
	else:
		cache.updateFreq(freq)
		delay = cache.cableDelay(stand)
	if stand in addDelay.keys():
		delay += (addDelay[stand]*1e-9)

	return delay
	

class SignalCache(CableCache):
	"""Subclass of CableCache that does for signalDelay what CableCache did for
	cableDelay."""

	def signalDelay(self, stand):
		"""Experimental function that wraps the cable delay for a stand along with 
		any other delays determined by phase fitting.  The list of stands that can be 
		corrected is incomplete.  Similar to cableDelay, a numpy array is returned."""

		return signalDelay(stand, self.freq, cache=self)


def cableAttenuation(stand):
	"""For a given stands, return the factor needed to correct for cable 
	losses."""

	# Stands start at 1, indicies do not
	validateStand(stand)

	# Catch the outlier and load cable and delay information into holder variables.
	# This makes it easier to do one common set of operations with different cables.
	if stand == -1:
		# For the outlier, the attenuation comes from the Times Microwave LMR-400 
		# data sheet for 50 Mhz.
		cableLength = 315.344 # m
		dBperM = 0.029 # dB/m
	else:
		cableLengths = _loadDelayData()
		cableLength = cableLengths[stand-1]
		dBperM = 0.047 # dB / m

	dBLoss = dBperM*cableLength
	multFactor = 10.0**(dBLoss/10.0)

	return multFactor


def getBaselines(stands, IncludeAuto=False, Indicies=False):
	"""Generate a list of two-element tuples that describe which antennae
	compose each of the output uvw triplets from computeUVW/computeUVTrack.
	If the Indicies keyword is set to True, the two-element tuples 
	contain the indicies of the stands array used, rather than the actual
	stand numbers."""

	if IncludeAuto:
		offset = 0
	else:
		offset = 1

	try:
		N = stands.shape[0]
	except AttributeError:
		stands = numpy.array(stands)
		N = stands.shape[0]
	out = []
	for i in range(0, N-offset):
		for j in range(i+offset, N):
			if Indicies:
				out.append( (i, j) )
			else:
				out.append( (stands[i], stands[j]) )

	return out
	

def baseline2antenna(baseline, stands, BaselineList=None, IncludeAuto=False, Indicies=False):
	"""Given a baseline number, a list of stands, and options of how the base-
	line listed was generated, convert the baseline number to antenna numbers.
	Alternatively, use a list of baselines instead of generating a new list.  
	This utility is useful for figuring out what antennae comprise a baseline."""

	# Build up the list of baselines using the options provided
	if BaselineList is None:
		BaselineList = getBaselines(stands, IncludeAuto=IncludeAuto, Indicies=Indices)
	
	# Select the correct one and return based on the value of Indicies
	i,j = BaselineList[baseline]
	if Indicies:
		return i,j
	else:
		return stands[i],stands[j]


def antenna2baseline(ant1, ant2, stands, BaselineList=None, IncludeAuto=False, Indicies=False):
	"""Given two antenna numbers, a list of stands, and options to how the base-
	line listed was generated, convert the antenna pair to  a baseline number. 
	This utility is useful for picking out a particular pair from a list of
	baselines."""

	# Build up the list of baselines using the options provided
	if BaselineList is None:
		BaselineList = getBaselines(stands, IncludeAuto=IncludeAuto, Indicies=Indices)

	# If we don't want indicies, convert stand numbers to indicies
	if not Indicies:
		ant1 = (numpy.where( stands == ant1 ))[0][0]
		ant2 = (numpy.where( stands == ant2 ))[0][0]
	
	# Loop over the baselines until we find one that matches.  If we don't find 
	# one, return -1
	i = 0
	for baseline in BaselineList:
		if ant1 in baseline and ant2 in baseline:
			return i
		else:
			i = i + 1
	else:
			return -1


def computeUVW(stands, HA=0.0, dec=34.070, freq=49.0e6, IncludeAuto=False):
	"""Compute the uvw converate of a baselines formed by a collection of 
	stands.  The coverage is computed at a given HA (in hours) and 
	declination (in degrees) for LWA-1.  The frequency provided (in Hz) can 
	either as a scalar or as a numpy array.  If more than one frequency is 
	given, the output is a three dimensional with dimensions of baselines, 
	uvw, and frequencies."""

	# Try this so that freq can be either a scalar or a list
	try:
		junk = len(freq)
	except TypeError:
		freq = numpy.array([freq])

	N = stands.shape[0]
	baselines = getBaselines(stands, IncludeAuto=IncludeAuto, Indicies=True)
	Nbase = len(baselines)
	Nfreq = freq.shape[0]
	uvw = numpy.zeros((Nbase,3,Nfreq))

	# Phase center coordinates
	# Convert numbers to radians and, for HA, hours to degrees
	HA2 = HA * 15.0 * deg2rad
	dec2 = dec * deg2rad
	lat2 = 34.070 * deg2rad

	# Coordinate transformation matrices
	trans1 = numpy.matrix([[0, -math.sin(lat2), math.cos(lat2)],
					   [1, 0,               0],
					   [0, math.cos(lat2),  math.sin(lat2)]])
	trans2 = numpy.matrix([[math.sin(HA2),                 math.cos(HA2),                 0],
					   [-math.sin(dec2)*math.cos(HA2), math.sin(dec2)*math.sin(HA2),  math.cos(dec2)],
					   [math.cos(dec2)*math.cos(HA2),  -math.cos(dec2)*math.sin(HA2), math.sin(dec2)]])

	# Position cache
	posCache = PositionCache()

	count = 0
	for i,j in baselines:
		# Go from a east, north, up coordinate system to a celestial equation, 
		# east, north celestial pole system
		xyzPrime = posCache.getRelativeXYZ(stands[i], stands[j])
		xyz = trans1*numpy.matrix([[xyzPrime[0]],[xyzPrime[1]],[xyzPrime[2]]])
		
		# Go from CE, east, NCP to u, v, w
		temp = trans2*xyz
		for k in range(Nfreq):
			uvw[count,:,k] = numpy.squeeze(temp) * freq[k] / c
		count = count + 1

	return uvw


def computeUVTrack(stands, dec=34.070, freq=49.0e6):
	"""Whereas computeUVW provides the uvw coverage at a particular time, 
	computeUVTrack provides the complete uv plane track for a long 
	integration.  The output is a three dimensional numpy array with 
	dimensions baselines, uv, and 512 points along the track ellipses.  
	Unlike computeUVW, however, only a single frequency (in Hz) can be 
	specified."""

	# Try this so that freq can be either a scalar or a list
	try:
		junk = len(freq)
	except TypeError:
		freq = numpy.array([freq])

	N = stands.shape[0]
	Nbase = N*(N-1)/2
	uvTrack = numpy.zeros((Nbase,2,512))

	# Phase center coordinates
	# Convert numbers to radians and, for HA, hours to degrees
	dec2 = dec * math.pi/180.0
	lat2 = 34.070 * math.pi/180.0

	# Coordinate transformation matrices
	trans1 = numpy.matrix([[0, -math.sin(lat2), math.cos(lat2)],
					   [1, 0,               0],
					   [0, math.cos(lat2),  math.sin(lat2)]])

	# Position cache
	posCache = PositionCache()

	count = 0
	for i,j in getBaselines(stands, Indicies=True):
		# Go from a east, north, up coordinate system to a celestial equation, 
		# east, north celestial pole system
		xyzPrime = posCache.getRelativeXYZ(stands[i], stands[j])
		xyz = trans1*numpy.matrix([[xyzPrime[0]],[xyzPrime[1]],[xyzPrime[2]]])

		uRange = numpy.linspace(-math.sqrt(xyz[0]**2 + xyz[1]**2), math.sqrt(xyz[0]**2 + xyz[1]**2), num=256)
		vRange1 = numpy.sqrt(xyz[0]**2 + xyz[1]**2 - uRange**2)*sin(dec2) + xyz[2]*cos(dec2)
		vRange2 = -numpy.sqrt(xyz[0]**2 + xyz[1]**2 - uRange**2)*sin(dec2) + xyz[2]*cos(dec2)

		uvTrack[count,0,0:256] = uRange * freq / c
		uvTrack[count,1,0:256] = vRange1 * freq / c
		uvTrack[count,0,256:512] = uRange[::-1] * freq / c
		uvTrack[count,1,256:512] = vRange2[::-1] * freq / c
		count = count + 1
			
	return uvTrack
