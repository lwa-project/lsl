# -*- coding: utf-8 -*-

"""
Classes and methods to model sky brightness and visibility.
"""

### This module handles the skymaps from the 74 MHz skymap.
### David Munton, ARL:UT Jan 2007

from numpy import *
import os,sys

import pyfits

from lsl import astro

__version__   = '0.2'
__revision__ = '$Rev$'
__author__    = 'J. York'
__maintainer__ = 'Jayce Dowell'


### This code is the base class for the sky map. It takes as input a skymap file name and frequency to
### which the skymap corresponds.  It has the following methods:
###     _init_ - takes the array coordinate filename as an input argument.
###     NormalizePower - Converts the skymap powers (in Kelvin radiated into 4 pi ster) into a power seen at the antenna.
###     ComputeTotalPowerFromSky - Sums the power for all sources in the sky
###     ScaleSourcePowerstoFrequency - Scales the skymap from the base 73.8 MHz to the desired frequency.

class SkyMapError(Exception):
	"""
	Base class for SkyMap exceptions.
	"""
	
	pass


class SkyMapLoadError(SkyMapError):
	"""
	Exception raised when a SkyMap fails to load.
	"""
	
	pass


class SkyMapPowerError(SkyMapError):
	"""
	Exception raised when an operation on a normalized SkyMap fails to 
	compute (such as summing or multiplying over an empty data set).
	"""
	
	pass


class SkyMap(object):
	"""
	The class for handling the model sky brightness maps.  This code is the base 
	class for the sky map. It takes as input a skymap file name and frequency to
	which the skymap corresponds.  It has the following methods:
	  1. _init_ - takes the array coordinate filename as an input argument.
	  2. NormalizePower - Converts the skymap powers (in Kelvin radiated into 4 pi ster)
	  into a power seen at the antenna.
	  3. ComputeTotalPowerFromSky - Sums the power for all sources in the sky
	  4. ScaleSourcePowerstoFrequency - Scales the skymap from the base 73.8 MHz to the
	  desired frequency.
	"""
	
	# Class data 
	degToRad = (pi/180.)             # Usual conversion factor
	    
	def __init__(self, skyMapFileName=None, freqMHz=73.9):
		"""
		Initialize the SkyMap object with an optional full file path to 
		the skymap file.
		"""
		
		if skyMapFileName is None:
			from lsl.common.paths import data as dataPath
			skyMapFileName = os.path.join(dataPath, 'skymap', 'LFmap_73.9.fits')
			
		fits = pyfits.open(skyMapFileName)
		
		hdr = fits[0].header
		
		self.numPixelsX = hdr['NAXIS1']
		self.numPixelsY = hdr['NAXIS2']
		
		self.freq = float(hdr['FINALFRQ'].split()[0])
	  
		ra = (arange(1, self.numPixelsX + 1, dtype = float32) - hdr['CRPIX1']) * hdr['CDELT1']
		self.ra = ra.reshape(1, len(ra)).repeat(self.numPixelsY, 0).ravel()
		
		dec = (arange(1, self.numPixelsY + 1, dtype = float32) - hdr['CRPIX2']) * hdr['CDELT2']
		self.dec = dec.repeat(self.numPixelsX, -1)
		
		self._power = empty((self.numPixelsX * self.numPixelsY,), float32) 
		data = fits[0].data
		
		for i in range(self.numPixelsY):
			ii = i * self.numPixelsX
			for j in range(self.numPixelsX):
				self._power[ii + j] = data[i, j]
		
		fits.close()
		
	def NormalizePower(self):
		"""
		Compute the skymap power (total power radiated into 4 pi steradians) into 
		a power at antenna, based on pixel count.
		"""
		
		# The cosine term is the projection of the receiving area onto the direction 
		# of the source
		return self._power*cos(self.dec * self.degToRad) 
		
	def ComputeTotalPowerFromSky(self):
		"""
		Compute and return the the total power from the sky.
		"""
		
		if len(self._power) == 0:
			raise SkyMapPowerError("self._power contains 0 elements")
		return self.NormalizePower().sum()


class SkyMapGSM(SkyMap):
	"""
	Extension of the SkyMap class to use a Global Sky Model-based set of maps
	of T0 @ f0, alpha, and beta, where:
	  .. math:: T(f) = T_0 \\times \\left(\\frac{f}{f_0}\\right)^{\\alpha + \\beta*\\log{f/f_0}}
	
	For more information on the Global Sky Model, see: http://space.mit.edu/~angelica/gsm/index.html
	"""
	
	def __init__(self, skyMapFileName=None, freqMHz=73.9):
		"""
		Initialize the SkyMapGSM object with an optional full file path to 
		the skymap file.
		"""
		
		if skyMapFileName is None:
			from lsl.common.paths import data as dataPath
			skyMapFileName = os.path.join(dataPath, 'skymap', 'gsm.npz')
			
		# Since we are using a pre-computed GSM which is a NPZ file, read it
		# in with numpy.load.
		dataDict = load(skyMapFileName)
		
		# RA and dec. are stored in the 'eCoords' keyword as radians
		self.ra = dataDict['eCoords'][:,0].ravel() / self.degToRad
		self.dec = dataDict['eCoords'][:,1].ravel() / self.degToRad
		
		# Compute the temperature for the current frequency using the stroed 
		# values of T0, alpha (the spectral index), and beta (the spectral
		# curvature)
		fRatio = freqMHz / dataDict['f0']
		self._power = dataDict['T0']*fRatio**(dataDict['alpha']+dataDict['beta']*log10(fRatio))
		
	def NormalizePower(self):
		"""
		Compute the skymap power (total power radiated into 4 pi steradians) into 
		a power at antenna, based on pixel count.
		"""
		
		return self._power


class ProjectedSkyMap(object):
	"""
	The class for handling the model sky brightness maps over a particular site.
	This code is the base class for the sky map visible at a specific location. It 
	takes as input a skymap file name and frequency to which the skymap corresponds.
	It inherits from class SkyMap. It has the following methods:
	  1. _init_ - takes the array coordinate filename as an input argument.
	  2. ComputeDirectionCosines - Computes the direction cosines 
	  3. ComputeTotalPowerFromVisibleSky - Sums the power for all visible sources in 
	  the sky.
	"""
	
	def __init__(self, skyMapObject, lat, lon, utc_jd):
		"""
		Initialize the skymap at input lat,lon (decimal degrees) and time (in 
		UTC julian day).
		"""
		
		self.skyMapObject=skyMapObject
		
		assert -90<=lat<=90, ValueError('lat = %g not in [-90,90]' % lat)
		assert -360<=lon<=360, ValueError('lon = %g not in [-360,360]' % lon)
		self.lat=lat
		self.lon=lon
		self.time=utc_jd
		
		# Compute the ra and dec locations to a altitude and azimuth. This requires a lat/lon and a time (UTC).
		# Alt and Az are compressed to show only the visible source.
		## Extract the RA/dec values
		ra = self.skyMapObject.ra
		dec = self.skyMapObject.dec
		
		## Compute the LST at the lat/lon (hours)
		lst = astro.get_local_sidereal_time(self.lon, self.time)
		
		## RA -> HA (deg)
		ha = lst*15.0 - ra
		
		## HA, dec, and lat in degrees to alt and az in radians
		sinAlt = sin(dec*pi/180)*sin(self.lat*pi/180) + cos(dec*pi/180)*cos(lat*pi/180)*cos(ha*pi/180)
		alt = arcsin(sinAlt)
		cosAz = (sin(dec*pi/180)-sinAlt*sin(self.lat*pi/180))/(cos(alt)*cos(lat*pi/180))
		az = arccos(cosAz)
		swap = where(sin(ha*pi/180) < 0)
		az[swap] = 2*pi-az[swap]
		
		## Convert to alt and az to degrees
		self.alt = alt*180/pi
		self.az  = az*180/pi
		
		# Compress the az/alt so that only the visible sources are available. "
		visibleMask = self.alt > 0.0
		
		#Compress arrays to hide non-visible sources
		self.visibleAlt = compress(visibleMask, self.alt)
		self.visibleAz = compress(visibleMask, self.az)
		self.visiblePower = compress(visibleMask, self.skyMapObject._power)
		try:
			pixelSolidAngle = (2.0 * pi)**2/(self.skyMapObject.numPixelsX*self.skyMapObject.numPixelsY)
		except AttributeError:
			pixelSolidAngle = 2.5566346e-4
		fractionSolidAngle = (1.0/4.0*pi) * pixelSolidAngle
		
		# The cosine term is the projection of the receiving area onto the direction 
		# of the source
		normalizedPower = self.skyMapObject.NormalizePower() * fractionSolidAngle  
		self.visibleNormalizedPower = compress(visibleMask, normalizedPower)
		#self.visibleNormalizedPower = self.visiblePower*cos(self.visibleAlt * self.skyMapObject.degToRad)
		self.visibleRa = compress(visibleMask, self.skyMapObject.ra)
		self.visibleDec = compress(visibleMask, self.skyMapObject.dec)
		
	def ComputeDirectionCosines(self):
		"""
		Compute the direction cosines and return the tuple of arrays (l,m,n).
		"""
		
		altRad = self.visibleAlt*self.skyMapObject.degToRad
		azRad = self.visibleAz*self.skyMapObject.degToRad
		l = cos(altRad)*sin(azRad)
		m = cos(altRad)*cos(azRad)
		n = sin(altRad)
		return (l,m,n)
		
	def ComputeTotalPowerFromVisibleSky(self):
		"""
		Compute and return the the total power from visible portion of the sky.
		"""
		
		if len(self.visibleNormalizedPower) == 0:
			raise SkyMapPowerError("visibleNormalizedPower contains 0 elements")
		totalVisiblePower = sum(self.visibleNormalizedPower)
		return totalVisiblePower
