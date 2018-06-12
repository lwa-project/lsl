# -*- coding: utf-8 -*-

# Python3 compatibility
from __future__ import print_function
import sys
if sys.version_info > (3,):
	xrange = range
	long = int

"""
Module for creating object oriented representations of the LWA stations.
"""

import os
import re
import copy
import numpy
import ephem
import struct

from lsl.astro import DJD_OFFSET
from lsl.common.paths import data as dataPath
from lsl.common import mcs, mcsADP
from lsl.common.constants import c as speedOfLight
from lsl.misc.mathutil import to_dB, from_dB
from lsl.misc.total_sorting import cmp_to_total

__version__ = '2.2'
__revision__ = '$Rev$'
__all__ = ['geo2ecef', 'ecef2geo', 'LWAStation', 'Antenna', 'Stand', 'FEE', 'Cable', 'ARX', 'LSLInterface', 
		 'parseSSMIF', 'lwa1', 'lwavl', 'lwana', 'lwasv',  'getFullStations', 
		 'PrototypeStation', 'prototypeSystem', 
		 '__version__', '__revision__', '__all__']


_id2name = {'VL': 'LWA1', 'NA': 'LWANA', 'SV': 'LWASV'}


def geo2ecef(lat, lon, elev):
	"""
	Convert latitude (rad), longitude (rad), elevation (m) to earth-
	centered, earth-fixed coordinates.
	"""
	
	WGS84_a = 6378137.000000
	WGS84_b = 6356752.314245
	N = WGS84_a**2 / numpy.sqrt(WGS84_a**2*numpy.cos(lat)**2 + WGS84_b**2*numpy.sin(lat)**2)
	
	x = (N+elev)*numpy.cos(lat)*numpy.cos(lon)
	y = (N+elev)*numpy.cos(lat)*numpy.sin(lon)
	z = ((WGS84_b**2/WGS84_a**2)*N+elev)*numpy.sin(lat)
	
	return (x, y, z)


def ecef2geo(x, y, z):
	"""
	Convert earth-centered, earth-fixed coordinates to (rad), longitude 
	(rad), elevation (m) using Bowring's method.
	"""
	
	WGS84_a = 6378137.000000
	WGS84_b = 6356752.314245
	e2 = (WGS84_a**2 - WGS84_b**2) / WGS84_a**2
	ep2 = (WGS84_a**2 - WGS84_b**2) / WGS84_b**2
	
	# Distance from rotation axis
	p = numpy.sqrt(x**2 + y**2)
	
	# Longitude
	lon = numpy.arctan2(y, x)
	p = numpy.sqrt(x**2 + y**2)
	
	# Latitude (first approximation)
	lat = numpy.arctan2(z, p)
	
	# Latitude (refined using Bowring's method)
	psi = numpy.arctan2(WGS84_a*z, WGS84_b*p)
	num = z + WGS84_b*ep2*numpy.sin(psi)**3
	den = p - WGS84_a*e2*numpy.cos(psi)**3
	lat = numpy.arctan2(num, den)
	
	# Elevation
	N = WGS84_a**2 / numpy.sqrt(WGS84_a**2*numpy.cos(lat)**2 + WGS84_b**2*numpy.sin(lat)**2)
	elev = p / numpy.cos(lat) - N
	
	return lat, lon, elev


class LWAStationBase(object):
	"""
	Base object to hold information about the a LWA station.  Stores station:
	  * Name (name)
	  * ID code (id)
	  * List of Antenna instances (antennas)
	  
	.. versionchanged:: 1.2.0
		Added a new 'interface' keyword to reference a LSLInterface instance.
	
	.. versionadded:: 1.0.0
	"""
	
	def __init__(self, name, id='', antennas=None, interface=None):
		self.name = name
		self.id = id
		
		if antennas is None:
			self.antennas = []
		else:
			self.antennas = antennas
			
		if interface is None:
			self.interface = LSLInterface()
		else:
			self.interface = interface
			
	def __str__(self):
		return "%s (%s) with %i antennas" % (self.name, self.id, len(self.antennas))
		
	def __reduce__(self):
		return (LWAStationBase, (self.name, self.id, self.antennas, self.interface))
		
	def _sortAntennas(self, attr='digitizer'):
		"""
		Sort the antennas list by the specified attribute.  The default
		attribute is the digitizer number.
		"""
		
		# Build the sorting function
		def sortFcn(x, y):
			if getattr(x, attr) > getattr(y, attr):
				return 1
			elif getattr(x, attr) < getattr(y, attr):
				return -1
			else:
				return 0
				
		self.antennas.sort(cmp=sortFcn)


class LWAStation(ephem.Observer, LWAStationBase):
	"""
	Object to hold information about the a LWA station.  This object can
	create a ephem.Observer representation of itself and identify which stands
	were in use at a given time.  Stores station:
	  * Name (name)
	  * ID code (id)
	  * Latitiude in radians [but initialized as degrees] (N is positive, lat)
	  * Longitude in radians [but initialized as degrees] (W is negative, long)
	  * Elevation in meters (elev)
	  * List of Antenna instances (antennas)
	  
	LWAStation provides several functional attributes for dealing with the 
	station's location on Earth.  These include:
	  * getObserver: Return an ephem.Observer instance representing the station
	  * getAIPYLocation: Return a tuple for setting the location of an AIPY
	    AntennaArray instance
	  * getGeocentricLocation: Return a tuple of the EC-EF coordinates of the 
	    station
	  * getECITransform: Return a 3x3 transformation matrix to convert antenna
	    positions to ECI coordinates
	  * getECIInverseTransform: Return a 3x3 transformation matrix to convert antenna
	    positions from ECI coordinates
	  * getENZOffset: Return the east, north, and vertical offsets to a point on
	    the surface of the Earth
	  * getPointingAndDistance: Return the pointing direction and distance to 
	    another location on the surface of the Earth
	    
	LWAStation also provides several functional attributes for dealing with
	the station's antennas.  These include:
	  * getAntennas:  Return a list of antennas
	  * getStands: Return a list of stands
	  * getPols:  Return a list of polarizations
	  * getCables: Return a list of cables
	  
	.. versionchanged:: 1.2.0
		Added a new 'interface' attribute which provides referenves to various modules
		to help interface with the station.
	
	.. versionchanged:: 1.0.0
		Converted LWAStation to be an instance of LWAStationBase and ephem.Observer
		to make it easier to work with ephem.Body objects.
		
		Added additional functions for dealing with other locations.
		
		Changed getECEFTransform() to getECITransform() to make the function name
		consistent with its purpose.
	"""
	
	def __init__(self, name, lat, long, elev, id='', antennas=None, interface=None):
		LWAStationBase.__init__(self, name, id=id, antennas=antennas, interface=interface)
		ephem.Observer.__init__(self)
		
		self.lat = lat * numpy.pi/180.0
		self.long = long * numpy.pi/180.0
		self.elev = elev
		self.pressure = 0.0
		
	def __str__(self):
		return "%s (%s) at lat: %.3f, lng: %.3f, elev: %.1f m with %i antennas" % (self.name, self.id, self.lat*180.0/numpy.pi, self.long*180.0/numpy.pi, self.elev, len(self.antennas))
		
	def __reduce__(self):
		return (LWAStation, (self.name, self.lat*180/numpy.pi, self.long*180/numpy.pi, self.elev, self.id, self.antennas, self.interface))
		
	def compute(self, body):
		"""
		Update the provided ephem.Body instance with the current location as 
		viewed from the site.
		
		.. versionadded:: 1.0.0
		"""
		
		body.compute(self)
		
	def getObserver(self, date=None, JD=False):
		"""
		Return a ephem.Observer object for this site.
		"""
		
		oo = ephem.Observer()
		oo.lat = 1.0*self.lat
		oo.long = 1.0*self.long
		oo.elev = 1.0*self.elev
		oo.pressure = 0.0
		if date is not None:
			if JD:
				# If the date is Julian, convert to Dublin Julian Date 
				# which is used by ephem
				date -= DJD_OFFSET
			oo.date = date
			
		return oo
		
	def getAIPYLocation(self):
		"""
		Return a tuple that can be used by AIPY for specifying a array
		location.
		"""
		
		return (self.lat, self.long, self.elev)
		
	def getGeocentricLocation(self):
		"""
		Return a tuple with earth-centered, earth-fixed coordinates for the station.
		"""
		
		return geo2ecef(self.lat, self.long, self.elev)
		
	def getECITransform(self):
		"""
		Return a 3x3 transformation matrix that converts a baseline in 
		[east, north, elevation] to earth-centered inertial coordinates
		for that baseline [x, y, z].  Based off the 'local_to_eci' 
		function in the lwda_fits-dev library.
		"""
		
		return numpy.array([[0.0, -numpy.sin(self.lat), numpy.cos(self.lat)], 
						[1.0, 0.0,                  0.0], 
						[0.0, numpy.cos(self.lat),  numpy.sin(self.lat)]])
						
	def getECIInverseTransform(self):
		"""
		Return a 3x3 transformation matrix that converts a baseline in 
		earth-centered inertial coordinates [x, y, z] to [east, north, 
		elevation] for that baseline.
		"""
		
		return numpy.array([[ 0.0,                 1.0, 0.0                ],
						[-numpy.sin(self.lat), 0.0, numpy.cos(self.lat)],
						[ numpy.cos(self.lat), 0.0, numpy.sin(self.lat)]])
						
	def getENZOffset(self, locTo):
		"""
		Given another location on the surface of the Earth, either as a 
		LWAStation instance or a three-element tuple of latitude (deg.), 
		longitude (deg.), and elevation (m), return the topocentric offset
		in meter along the east, north, and vertical directions.
		"""
		
		ecefFrom = self.getGeocentricLocation()
		try:
			ecefTo = locTo.getGeocentricLocation()
		except AttributeError:
			ecefTo = geo2ecef(float(locTo[0])*numpy.pi/180, float(locTo[1])*numpy.pi/180, locTo[2])
			
		ecefFrom = numpy.array(ecefFrom)
		ecefTo = numpy.array(ecefTo)
		
		rho = ecefTo - ecefFrom
		rot = numpy.array([[ numpy.sin(self.lat)*numpy.cos(self.long), numpy.sin(self.lat)*numpy.sin(self.long), -numpy.cos(self.lat)], 
					    [-numpy.sin(self.long),                     numpy.cos(self.long),                      0                  ],
					    [ numpy.cos(self.lat)*numpy.cos(self.long), numpy.cos(self.lat)*numpy.sin(self.long),  numpy.sin(self.lat)]])
		sez = numpy.dot(rot, rho)
		
		# Convert from south, east, zenith to east, north, zenith
		enz = 1.0*sez[[1,0,2]]
		enz[1] *= -1.0
		
		return enz
		
	def getPointingAndDistance(self, locTo):
		"""
		Given another location on the surface of the Earth, either as a 
		LWAStation instance or a three-element tuple of latitude (deg.), 
		longitude (deg.), and elevation (m), return the bearing azimuth/
		elevation in radians and distance in meters to the location.
		
		.. versionchanged:: 1.0.1
			Renamed from getPointingAndDirection to getPointingAndDistance
		"""
		
		ecefFrom = self.getGeocentricLocation()
		try:
			ecefTo = locTo.getGeocentricLocation()
		except AttributeError:
			ecefTo = geo2ecef(float(locTo[0])*numpy.pi/180, float(locTo[1])*numpy.pi/180, locTo[2])
			
		ecefFrom = numpy.array(ecefFrom)
		ecefTo = numpy.array(ecefTo)
		
		rho = ecefTo - ecefFrom
		rot = numpy.array([[ numpy.sin(self.lat)*numpy.cos(self.long), numpy.sin(self.lat)*numpy.sin(self.long), -numpy.cos(self.lat)], 
					    [-numpy.sin(self.long),                     numpy.cos(self.long),                      0                  ],
					    [ numpy.cos(self.lat)*numpy.cos(self.long), numpy.cos(self.lat)*numpy.sin(self.long),  numpy.sin(self.lat)]])
		sez = numpy.dot(rot, rho)
		
		d = numpy.sqrt( (rho**2).sum() )
		el = numpy.arcsin(sez[2] / d)
		az = numpy.arctan2(sez[1], -sez[0])
		
		return az, el, d
		
	def getAntennas(self):
		"""
		Return the list of Antenna instances for the station, sorted by 
		digitizer number.
		"""
		
		# Sort and return
		self._sortAntennas()
		return self.antennas
		
	def getStands(self):
		"""
		Return a list of Stand instances for each antenna, sorted by 
		digitizer number.
		"""
		
		# Sort and return
		self._sortAntennas()
		return [ant.stand for ant in self.antennas]
	
	def getPols(self):
		"""
		Return a list of polarization (0 == N-S; 1 == E-W) for each antenna, 
		sorted by digitizer number.
		"""
		
		# Sort and return
		self._sortAntennas()
		return [ant.pol for ant in self.antennas]
		
	def getCables(self):
		"""
		Return a list of Cable instances for each antenna, sorted by
		digitizer number.
		"""
		
		# Sort and return
		self._sortAntennas()
		return [ant.cable for ant in self.antennas]


@cmp_to_total
class Antenna(object):
	"""
	Object to store the information about an antenna.  Stores antenna:
	  * ID number (id)
	  * ARX instance the antenna is attached to (arx)
	  * DP1/ROACH board number (board)
	  * DP1/ROACH  digitizer number (digitizer)
	  * DP/ADP rack input connector (input)
	  * Stand instance the antenna is part of (stand)
	  * Polarization (0 == N-S; pol)
	  * Antenna vertical mis-alignment in degrees (theta)
	  * Antenna rotation mis-alignment in degrees (phi)
	  * Fee instance the antenna is attached to (fee)
	  * Port of the FEE used for the antenna (feePort)
	  * Cable instance used to connect the antenna (cable)
	  * Status of the antenna (status)
	  
	Status codes are:
	  * 0 == Not installed
	  * 1 == Bad
	  * 2 == Suspect, possibly bad
	  * 3 == OK
	
	.. versionchanged:: 1.0.0
		Added an attribute to hold the DP rack input connector label.
	"""
	
	def __init__(self, id, arx=None, board=0, digitizer=0, input='', stand=None, pol=0, theta=0.0, phi=0.0, fee=None, feePort=1, cable=None, status=0):
		self.id = int(id)
		if arx is None:
			self.arx = ARX(0, 0, 0)
		else:
			self.arx = arx
		self.board = int(board)
		self.digitizer = int(digitizer)
		self.input = input
		
		if stand is None:
			self.stand = Stand(0, 0, 0, 0)
		else:
			self.stand = stand
			
		self.pol = int(pol)
		self.theta = float(theta)
		self.phi = float(phi)
		
		if fee is None:
			self.fee = FEE('', 0)
		else:
			self.fee = fee
		self.feePort = feePort
		
		if cable is None:
			self.cable = Cable('', 0)
		else:
			self.cable = cable
			
		self.status = int(status)
		
	def __str__(self):
		return "Antenna %i: stand=%i, polarization=%i; digitizer %i; status is %i" % (self.id, self.stand.id, self.pol, self.digitizer, self.status)
		
	def __reduce__(self):
		return (Antenna, (self.id, self.arx, self.board, self.digitizer, self.input, self.stand, self.pol, self.theta, self.phi, self.fee, self.feePort, self.cable, self.status))
		
	def __cmp__(self, y):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0
			
	def response(self, dB=False):
		"""
		Return a two-element tuple (freq in Hz, mis-match efficiency) for a model LWA1 
		antenna from Hicks et al. (2012, PASP, 124, 1090).
		
		.. versionadded:: 1.0.0
		"""
		
		# Find the filename to use
		filename = os.path.join(dataPath, 'BurnsZ.txt')
		
		# Read in the data
		data = numpy.loadtxt(filename)
		freq = data[:,0]*1e6
		ime = data[:,3]
		if dB:
			ime = to_dB(ime)
			
		# Return
		return (freq, ime)
		
	def getStatus(self):
		"""
		Return the combined antenna + FEE status as a two digit number 
		with the first digit representing the antenna status and the 
		second the FEE status.
		"""
		
		return 10*self.status + self.fee.status


@cmp_to_total
class Stand(object):
	"""
	Object to store the information (location and ID) about a stand.  
	Stores stand:
	  * ID number (id)
	  * Position relative to the center stake in meters (x,y,z)
	  
	The x, y, and z positions can also be accessed through subscripts:
	  Stand[0] = x
	  Stand[1] = y
	  Stand[2] = z
	  
	.. versionchanged:: 1.0.0
		Added the option to get the positions via subscripts.
	"""
	
	def __init__(self, id, x, y, z):
		self.id = int(id)
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		
	def __cmp__(self, y):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0
			
	def __str__(self):
		return "Stand %i:  x=%+.2f m, y=%+.2f m, z=%+.2f m" % (self.id, self.x, self.y, self.z)
		
	def __reduce__(self):
		return (Stand, (self.id, self.x, self.y, self.z))
		
	def __getitem__(self, key):
		if key == 0:
			return self.x
		elif key == 1:
			return self.y
		elif key == 2:
			return self.z
		else:
			raise ValueError("Subscript %i out of range" % key)
			
	def __setitem__(self, key, value):
		if key == 0:
			self.x = float(value)
		elif key == 1:
			self.y = float(value)
		elif key == 2:
			self.z = float(value)
		else:
			raise ValueError("Subscript %i out of range" % key)
			
	def __add__(self, std):
		try:
			# If its a Stand instance, do this
			out = (self.x+std.x, self.y+std.y, self.z+std.z)
		except AttributeError:
			try:
				# Maybe it is a list/tuple, so do this
				out = (self.x+std[0], self.y+std[1], self.z+std[2])
			except TypeError:
				out = (self.x+std, self.y+std, self.z+std)
		
		return out
		
	def __sub__(self, std):
		try:
			# If its a Stand instance, do this
			out = (self.x-std.x, self.y-std.y, self.z-std.z)
		except AttributeError:
			try:
				# Maybe it is a list/tuple, so do this
				out = (self.x-std[0], self.y-std[1], self.z-std[2])
			except TypeError:
				out = (self.x-std, self.y-std, self.z-std)
				
		return out


@cmp_to_total
class FEE(object):
	"""
	Object to store the information about a FEE.  Stores FEE:
	  * ID name (id)
	  * ID number (idNumber)
	  * Gain of port 1 (gain1)
	  * Gain of part 2 (gain2)
	  * Status (status)
	  
	Status codes are:
	  * 0 == Not installed
	  * 1 == Bad
	  * 2 == Suspect, possibly bad
	  * 3 == OK
	"""
	
	def __init__(self, id, idNumber, gain1=0, gain2=0, status=0):
		self.id = str(id)
		self.idNumber = int(idNumber)
		self.gain1 = float(gain1)
		self.gain2 = float(gain2)
		self.status = int(status)
		
	def __str__(self):
		return "FEE '%s': gain1=%.2f, gain2=%.2f; status is %i" % (self.id, self.gain1, self.gain2, self.status)
		
	def __reduce__(self):
		return (FEE, (self.id, self.idNumber, self.gain1, self.gain2, self.status))
		
	def __cmp__(self, y):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0
			
	def response(self, dB=False):
		"""
		Return a two-element tuple (freq in Hz, gain) for the frequency-
		dependent gain for a v1.7 FEE from LWA Memo #190, FEE0010, 
		Figure 3.
		
		.. versionadded:: 1.0.1
		"""
		
		# Find the filename to use
		filename = os.path.join(dataPath, 'fee.txt')
		
		# Read in the data
		data = numpy.loadtxt(filename)
		freq = data[:,0]*1e6
		gai = data[:,3]
		if not dB:
			gai = from_dB(gai)
			
		# Return
		return (freq, gai)


@cmp_to_total
class Cable(object):
	"""
	Object to store information about a cable.  Stores cable:
	  * ID name (id)
	  * Length in meters (length)
	  * Velocity factor (fractional, vf)
	  * Dispersive delay (seconds, dd)
	  * Gain term that goes as the square root of frequency (a0)
	  * Gain term that goes as frequency (a1)
	  * Gain term reference frequency (Hz, aFreq)
	  * Cable length stretch factor (stretch)
	  * Clock offset (seconds, clockOffset)
	
	The object also as a functional attribute named 'delay' that computes the
	cable delay for a particular frequency or collection of frequencies in 
	Hz.
	"""
	
	def __init__(self, id, length, vf=0, dd=0, a0=0.00428, a1=0.00000, aFreq=10e6, stretch=1.0):
		self.id = str(id)
		self.length = float(length)
		self.stretch = float(stretch)
		self.vf = float(vf)
		self.dd = float(dd)
		self.a0 = float(a0)
		self.a1 = float(a1)
		self.aFreq = float(aFreq)
		self.clockOffset = 0.0
		
	def __str__(self):
		return "Cable '%s' with length %.2f m (stretched to %.2f m)" % (self.id, self.length, self.length*self.stretch)
		
	def __reduce__(self):
		return (Cable, (self.id, self.length, self.vf, self.dd, self.a0, self.a1, self.aFreq, self.stretch))
		
	def __cmp__(self, y):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0
			
	def setClockOffset(self, offset):
		"""
		Add a clock offset (in seconds) to the cable model.
		"""
		
		self.clockOffset = float(offset)
		
	def clearClockOffset(self):
		"""
		Clear the clock offset of the cable model.
		"""
		
		self.clockOffset = 0.0
		
	def delay(self, frequency=49e6, ns=False):
		"""Get the delay associated with the cable in second (or nanoseconds 
		if the 'ns' keyword is set to True) for a particular frequency or 
		collection of frequencies in Hz."""
		
		# Bulk delay for the cable
		bulkDelay = self.length*self.stretch / (self.vf * speedOfLight)
		
		# Dispersion delay
		dispDelay = self.dd * (self.length*self.stretch / 100.0) / numpy.sqrt(frequency / numpy.array(self.aFreq))
		
		totlDelay = bulkDelay + dispDelay + self.clockOffset
		
		if ns:
			return totlDelay*1e9
		else:
			return totlDelay
			
	def attenuation(self, frequency=49e6, dB=False):
		"""Get the multiplicative factor needed to correct for the cable 
		loss for a specific frequency (in Hz).  If attenuations for more 
		than one frequency are needed, the frequencies can be passed in as 
		a numpy array.
		
		.. versionchanged:: 1.0.0
			Added the `dB' keyword to allow dB to be returned.
		"""
	
		atten = 2 * self.a0 * self.length*self.stretch * numpy.sqrt(frequency / numpy.array(self.aFreq)) 
		atten += self.a1 * self.length*self.stretch * (frequency / numpy.array(self.aFreq))
		atten = numpy.exp(atten)
		
		if dB:
			atten = to_dB(atten)
			
		return atten
		
	def gain(self, frequency=49e6, dB=False):
		"""Get the cable gain ("inverse loss") for a specific frequency (in 
		Hz).  If gains for more than one frequency are needed, the 
		frequencies can be passed in as a numpy array.
		
		.. versionchanged:: 1.0.0
			Added the `dB' keyword to allow dB to be returned.
		"""
		
		gai = 1.0 / self.attenuation(frequency=frequency, dB=False)
		
		if dB:
			gai = to_dB(gai)
			
		return gai
		
	def response(self, dB=False):
		"""
		Return a two-element tuple (freq in Hz, attenuation) for the cable
		using the model from LWA Memo #170.
		
		.. versionadded:: 1.0.1
		"""
		
		# Generate the frequencies to use
		freq = numpy.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
		freq = freq*1e6
		
		# Compute the attenuation
		attn = self.attenuation(frequency=freq, dB=dB)
		
		# Return
		return (freq, attn)


class ARX(object):
	"""
	Object to store information about a ARX board/channel combination.  Stores ARX:
	  * ID name (id)
	  * Channel number (channel; 1-16)
	  * ASP channel number (aspChannel; 1-520)
	  * ASP rack input connector label (input)
	  * ASP rack output connector label (output)
	
	The object also as a functional attribute named 'delay' that computes the
	cable delay for a particular frequency or collection of frequencies in 
	Hz.
	
	.. versionchanged:: 1.0.0
		Added attributes to hold the ASP rack input and output connector 
		labels.
	"""
	
	def __init__(self, id, channel=0, aspChannel=0, input='', output=''):
		self.id = id
		self.channel = int(channel)
		self.aspChannel = int(aspChannel)
		self.input = input
		self.output = output
		
	def __str__(self):
		return "ARX Board %s, channel %i (ASP Channel %i)" % (self.id, self.channel, self.aspChannel)
		
	def __reduce__(self):
		return (ARX, (self.id, self.channel, self.aspChannel, self.input, self.output))
		
	def response(self, filter='full', dB=True):
		"""
		Return a two-element tuple (freq in Hz, S21 magnitude in dB) for 
		the ARX response for the current board/channel from the "ARX0026" 
		memo on the "LWA Engineering Documents" wiki.
		
		Filter options are:
		  * 1 or 'full'
		  * 2 or 'reduced'
		  * 3 or 'split'
		  
		.. versionchanged:: 1.0.0
			Add an option to specify whether the magnitude should be 
			returned in dB or not.
		"""
		
		# Find the filename to use
		filename = 'ARX_board_%4s_filters_ch%i.npz' % (self.id, self.channel)
		filename = os.path.join(dataPath, 'arx', filename)
		
		# Read in the file and convert it to a numpy array
		try:
			dataDict = numpy.load(filename)
		except IOError:
			raise RuntimeError("Could not find the response data for ARX board #%s, channel %i" % (self.id, self.channel))

		freq = dataDict['freq']
		data = dataDict['data']
		try:
			dataDict.close()
		except AttributeError:
			pass
			
		if not dB:
			data = from_dB(data)
			
		# Return or raise an error
		if filter == 1 or filter == 'full':
			return (freq, data[:,0])
		elif filter == 2 or filter == 'reduced':
			return (freq, data[:,1])
		elif filter == 3 or filter == 'split':
			return (freq, data[:,2])
		else:
			raise ValueError("Unknown ARX filter '%s'" % filter)


class LSLInterface(object):
	"""
	Object to store information about how to work with the station in LSL.
	This includes names for the:
	  * Backend module to use (backend)
	  * MCS module to use (mcs)
	  * SDF module to use (sdf)
	  * Metadata module to use (metabundle)
	  * SDM module to use (sdm)
	
	.. versionadded:: 1.2.0
	"""
	
	def __init__(self, backend=None, mcs=None, sdf=None, metabundle=None, sdm=None):
		self.backend = backend
		self.mcs = mcs
		self.sdf = sdf
		self.metabundle = metabundle
		self.sdm = sdm
		
	def __str__(self):
		return "LSL Interfaces:\n Backend: %s\n MCS: %s\n SDF: %s\n Metadata: %s\n SDM: %s" % \
				(self.backend, self.mcs, self.sdf, self.metabundle, self.sdm)
		
	def __reduce__(self):
		return (LSLInterface, (self.backend, self.mcs, self.sdf, self.metabundle, self.sdm))
		
	def getModule(self, which):
		value = getattr(self, which)
		if value is None:
			raise RuntimeError("Unknown module for interface type '%s'" % which)
		exec('import %s as tempModule' % value, None, locals())
		return tempModule


def __parseTextSSMIF(filename):
	"""
	Given a human-readable (text) SSMIF file and return a collection of
	variables via locals() containing the files data.
	"""
	
	fh = open(filename, 'r')
	
	kwdRE = re.compile(r'(?P<keyword>[A-Z_0-9]+)(\[(?P<id1>[0-9]+?)\])?(\[(?P<id2>[0-9]+?)\])?(\[(?P<id3>[0-9]+?)\])?')
	
	# Loop over the lines in the file
	for line in fh:
		line = line.replace('\n', '')
		line = line.replace('\r', '')
		if len(line) == 0 or line.isspace():
			continue
		if line[0] == '#':
			continue
		
		keywordSection, value = line.split(None, 1)
		value = value.split('#', 1)[0]
		
		mtch = kwdRE.match(keywordSection)
		keyword = mtch.group('keyword')
		
		ids = [-1, -1, -1]
		for i in xrange(3):
			try:
				ids[i] = int(mtch.group('id%i' % (i+1)))
			except TypeError:
				pass
			
		#
		# Station Data
		#
			
		if keyword == 'STATION_ID':
			idn = str(value)
			continue
		if keyword == 'GEO_N':
			lat = float(value)
			continue
		if keyword == 'GEO_E':
			lon = float(value)
			continue
		if keyword == 'GEO_EL':
			elv = float(value)
			continue
		
		#
		# Stand & Antenna Data
		#
		
		if keyword == 'N_STD':
			nStand = int(value)
			
			stdPos = [[0.0, 0.0, 0.0] for n in xrange(nStand)]
			stdAnt = [n/2+1 for n in xrange(2*nStand)]
			stdOrie = [n % 2 for n in xrange(2*nStand)]
			stdStat = [3 for n in xrange(2*nStand)]
			stdTheta = [0.0 for n in xrange(2*nStand)]
			stdPhi = [0.0 for n in xrange(2*nStand)]
			
			stdDesi = [1 for x in xrange(2*nStand)]
			
			continue
			
		if keyword == 'STD_LX':
			stdPos[ids[0]-1][0] = float(value)
			continue
		if keyword == 'STD_LY':
			stdPos[ids[0]-1][1] = float(value)
			continue
		if keyword == 'STD_LZ':
			stdPos[ids[0]-1][2] = float(value)
			continue
		
		if keyword == 'ANT_STD':
			stdAnt[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ANT_ORIE':
			stdOrie[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ANT_STAT':
			stdStat[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ANT_THETA':
			stdTheta[ids[0]-1] = float(value)
			continue
		
		if keyword == 'ANT_PHI':
			stdPhi[ids[0]-1] = float(value)
			continue
		
		#
		# FEE, Cable, & SEP Data
		#
		
		if keyword == 'N_FEE':
			nFee = int(value)
			
			feeID = ["UNK" for n in xrange(nFee)]
			feeStat = [3 for n in xrange(nFee)]
			feeDesi = [1 for n in xrange(nFee)]
			feeGai1 = [35.7 for n in xrange(nFee)]
			feeGai2 = [35.7 for n in xrange(nFee)]
			feeAnt1 = [2*n+1 for n in xrange(nFee)]
			feeAnt2 = [2*n+2 for n in xrange(nFee)]
			
			continue
			
		if keyword == 'FEE_ID':
			feeID[ids[0]-1] = value
			continue
		
		if keyword == 'FEE_STAT':
			feeStat[ids[0]-1] = int(value)
			continue
			
		if keyword == 'FEE_DESI':
			feeDesi[ids[0]-1] = int(value)
			continue
			
		if keyword == 'FEE_GAI1':
			feeGai1[ids[0]-1] = float(value)
			continue
		if keyword == 'FEE_GAI2':
			feeGai2[ids[0]-1] = float(value)
			continue
			
		if keyword == 'FEE_ANT1':
			feeAnt1[ids[0]-1] = int(value)
			continue
		if keyword == 'FEE_ANT2':
			feeAnt2[ids[0]-1] = int(value)
			continue
		
		
		if keyword == 'N_RPD':
			nRPD = int(value)
			
			rpdID = ['UNK' for n in xrange(nRPD)]
			rpdStat = [3 for n in xrange(nRPD)]
			rpdLeng = [0.0 for n in xrange(nRPD)]
			rpdVF = [83.0 for n in xrange(nRPD)]
			rpdDD = [2.4 for n in xrange(nRPD)]
			rpdA0 = [0.00428 for n in xrange(nRPD)]
			rpdA1 = [0.00000 for n in xrange(nRPD)]
			rpdFre = [10e6 for n in xrange(nRPD)]
			rpdStr = [1.0 for n in xrange(nRPD)]
			rpdDesi = [1 for n in xrange(nRPD)]
			rpdAnt = [n+1 for n in xrange(nRPD)]
			
			continue
		
		if keyword == 'RPD_ID':
			rpdID[ids[0]-1] = value
			continue
		
		if keyword == 'RPD_STAT':
			rpdStat[ids[0]-1] = int(value)
			continue
		
		if keyword == 'RPD_LENG':
			rpdLeng[ids[0]-1] = float(value)
			continue
		
		if keyword == 'RPD_VF':
			if ids[0] == -1:
				rpdVF = [float(value) for n in xrange(nRPD)]
			else:
				rpdVF[ids[0]-1] = float(value)
			continue
		if keyword == 'RPD_DD':
			if ids[0] == -1:
				rpdDD = [float(value) for n in xrange(nRPD)]
			else:
				rpdDD[ids[0]-1] = float(value)
			continue
		
		if keyword == 'RPD_A0':
			if ids[0] == -1:
				rpdA0 = [float(value) for n in xrange(nRPD)]
			else:
				rpdA0[ids[0]-1] = float(value)
			continue
		
		if keyword == 'RPD_A1':
			if ids[0] == -1:
				rpdA1 = [float(value) for n in xrange(nRPD)]
			else:
				rpdA1[ids[0]-1] = float(value)
			continue
		
		if keyword == 'RPD_FREF':
			if ids[0] == -1:
				rpdFre = [float(value) for n in xrange(nRPD)]
			else:
				rpdFre[ids[0]-1] = float(value)
			continue
		
		if keyword == 'RPD_STR':
			if ids[0] == -1:
				rpdStr = [float(value) for n in xrange(nRPD)]
			else:
				rpdStr[ids[0]-1] = float(value)
			continue
		
		if keyword == 'RPD_DESI':
			rpdDesi[ids[0]-1] = value
			continue
		
		if keyword == 'RPD_ANT':
			rpdAnt[ids[0]-1] = int(value)
			continue
		
		
		if keyword == 'N_SEP':
			nSEP = int(value)
			
			sepCbl = ['UNK' for n in xrange(nSEP)]
			sepLeng = [0.0 for n in xrange(nSEP)]
			sepDesi = [1 for n in xrange(nSEP)]
			sepGain = [0.0 for n in xrange(nSEP)]
			sepAnt = [n+1 for n in xrange(nSEP)]
			
			continue
		
		if keyword == 'SEP_CABL':
			sepCbl[ids[0]-1] = value
			continue
		
		if keyword == 'SEP_LENG':
			sepLeng[ids[0]-1] = float(value)
			continue
		
		if keyword == 'SEP_DESI':
			sepDesi[ids[0]-1] = int(value)
			continue
		
		if keyword == 'SEP_GAIN':
			sepGain[ids[0]-1] = float(value)
			continue
		
		if keyword == 'SEP_ANT':
			sepAnt[ids[0]-1] = int(value)
			continue
		
		#
		# ARX (ARB) Data
		#
		
		if keyword == 'N_ARB':
			nARX = int(value)
			
			arxID = ["UNK" for n in xrange(nARX)]
			arxSlot = [0 for n in xrange(nARX)]
			arxDesi = [0 for n in xrange(nARX)]
			arxRack = [0 for n in xrange(nARX)]
			arxPort = [0 for n in xrange(nARX)]
			
			continue
		
		if keyword == 'N_ARBCH':
			nChanARX = int(value)
			
			arxStat = [[3 for c in xrange(nChanARX)] for n in xrange(nARX)]
			arxAnt = [[n*nChanARX+c+1 for c in xrange(nChanARX)] for n in xrange(nARX)]
			arxIn = [["UNK" for c in xrange(nChanARX)] for n in xrange(nARX)]
			arxOut = [["UNK" for c in xrange(nChanARX)] for n in xrange(nARX)]
			
			continue
		
		if keyword == 'ARB_ID':
			arxID[ids[0]-1] = value
			continue
		
		if keyword == 'ARB_SLOT':
			try:
				arxSlot[ids[0]-1] = int(value)
			except ValueError:
				arxSlot[ids[0]-1] = value
		
		if keyword == 'ARB_DESI':
			arxDesi[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ARB_RACK':
			arxRack[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ARB_PORT':
			arxRack[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ARB_STAT':
			arxStat[ids[0]-1][ids[1]-1] = int(value)
			continue
		
		if keyword == 'ARB_ANT':
			arxAnt[ids[0]-1][ids[1]-1] = int(value)
			continue
		
		if keyword == 'ARB_IN':
			arxIn[ids[0]-1][ids[1]-1] = value
			continue
		
		if keyword == 'ARB_OUT':
			arxOut[ids[0]-1][ids[1]-1] = value
			continue
		
		#
		# DP 1 & 2 Data - LWA1
		#
		
		if keyword == 'N_DP1':
			nDP1 = int(value)
			
			dp1ID = ["UNK" for n in xrange(nDP1)]
			dp1Slot = [0 for n in xrange(nDP1)]
			dp1Desi = [1 for n in xrange(nDP1)]
			
			continue
		
		if keyword == 'N_DP1CH':
			nChanDP1 = int(value)
			
			dp1Stat = [[3 for c in xrange(nChanDP1)] for n in xrange(nDP1)]
			dp1InR = [["UNK" for c in xrange(nChanDP1)] for n in xrange(nDP1)]
			dp1InC = [["UNK" for c in xrange(nChanDP1)] for n in xrange(nDP1)]
			dp1Ant = [[n*nChanDP1+c+1 for c in xrange(nChanDP1)] for n in xrange(nDP1)]
			
			continue
		
		if keyword == 'DP1_ID':
			dp1ID[ids[0]-1] = value
			continue
		
		if keyword == 'DP1_SLOT':
			dp1Slot[ids[0]-1] = value
			continue
		
		if keyword == 'DP1_DESI':
			dp1Desi[ids[0]-1] = int(value)
			continue
		
		if keyword == 'DP1_STAT':
			dp1Stat[ids[0]-1][ids[1]-1] = int(value)
			continue
		
		if keyword == 'DP1_INR':
			dp1InR[ids[0]-1][ids[1]-1] = value
			continue
		
		if keyword == 'DP1_INC':
			dp1InC[ids[0]-1][ids[1]-1] = value
			continue
		
		if keyword == 'DP1_ANT':
			dp1Ant[ids[0]-1][ids[1]-1] = int(value)
		
		
		if keyword == 'N_DP2':
			nDP2 = int(value)
			
			dp2ID = ["UNK" for n in xrange(nDP2)]
			dp2Slot = ["UNK" for n in xrange(nDP2)]
			dp2Stat = [3 for n in xrange(nDP2)]
			dp2Desi = [1 for n in xrange(nDP2)]
			
			continue
		
		if keyword == 'DP2_ID':
			dp2ID[ids[0]-1] = value
			continue
		
		if keyword == 'DP2_SLOT':
			dp2Slot[ids[0]-1] = value
			continue
		
		if keyword == 'DP2_STAT':
			dp2Stat[ids[0]-1] = int(value)
			continue
		
		if keyword == 'DP2_DESI':
			dp2Desi[ids[0]-1] = int(value)
			continue
		
		#
		# ROACH & Server Data - LWA-SV
		#
		
		if keyword == 'N_ROACH':
			nRoach = int(value)
			
			roachID = ["UNK" for n in xrange(nRoach)]
			roachSlot = [0 for n in xrange(nRoach)]
			roachDesi = [1 for n in xrange(nRoach)]
			
			continue
		
		if keyword == 'N_ROACHCH':
			nChanRoach = int(value)
			
			roachStat = [[3 for c in xrange(nChanRoach)] for n in xrange(nRoach)]
			roachInR = [["UNK" for c in xrange(nChanRoach)] for n in xrange(nRoach)]
			roachInC = [["UNK" for c in xrange(nChanRoach)] for n in xrange(nRoach)]
			roachAnt = [[n*nChanRoach+c+1 for c in xrange(nChanRoach)] for n in xrange(nRoach)]
			
			continue
		
		if keyword == 'ROACH_ID':
			roachID[ids[0]-1] = value
			continue
		
		if keyword == 'ROACH_SLOT':
			roachSlot[ids[0]-1] = value
			continue
		
		if keyword == 'ROACH_DESI':
			roachDesi[ids[0]-1] = int(value)
			continue
		
		if keyword == 'ROACH_STAT':
			roachStat[ids[0]-1][ids[1]-1] = int(value)
			continue
		
		if keyword == 'ROACH_INR':
			roachInR[ids[0]-1][ids[1]-1] = value
			continue
		
		if keyword == 'ROACH_INC':
			roachInC[ids[0]-1][ids[1]-1] = value
			continue
		
		if keyword == 'ROACH_ANT':
			roachAnt[ids[0]-1][ids[1]-1] = int(value)
		
		
		if keyword == 'N_SERVER':
			nServer = int(value)
			
			serverID = ["UNK" for n in xrange(nServer)]
			serverSlot = ["UNK" for n in xrange(nServer)]
			serverStat = [3 for n in xrange(nServer)]
			serverDesi = [1 for n in xrange(nServer)]
			
			continue
		
		if keyword == 'SERVER_ID':
			serverID[ids[0]-1] = value
			continue
		
		if keyword == 'SERVER_SLOT':
			serverSlot[ids[0]-1] = value
			continue
		
		if keyword == 'SERVER_STAT':
			serverStat[ids[0]-1] = int(value)
			continue
		
		if keyword == 'SERVER_DESI':
			serverDesi[ids[0]-1] = int(value)
			continue
		
		#
		# DR Data
		#
		
		if keyword == 'N_DR':
			nDR = int(value)
			
			drStat = [0 for n in xrange(nDR)]
			drID = ["UNK" for n in xrange(nDR)]
			drShlf = [0 for n in xrange(nDR)]
			drPC = ["UNK" for n in xrange(nDR)]
			drDP = [0 for n in xrange(nDR)]
			
			continue
		
		if keyword == 'DR_STAT':
			drStat[ids[0]-1] = int(value)
			continue
		
		if keyword == 'DR_ID':
			drID[ids[0]-1] = value
			continue
		
		if keyword == 'DR_SHLF':
			drShlf[ids[0]-1] = int(value)
			continue
		
		if keyword == 'DR_PC':
			drPC[ids[0]-1] = value
			continue
		
		if keyword == 'DR_DP':
			drDP[ids[0]-1] = int(value)
			continue
		
	fh.close()
	
	return locals()


def __parseBinarySSMIF(filename):
	"""
	Given a binary packed SSMIF file and return a collection of
	variables via locals() containing the files data.
	"""
	
	fh = open(filename, 'rb')
	
	# Read in the first four bytes to get the version code and go from there
	version = fh.read(4)
	version = struct.unpack('<i', version)[0]
	fh.seek(0)
	
	if version in (8,):
		## ADP
		mode = mcsADP
	else:
		## DP
		mode = mcs
	bssmif = mode.parseCStruct(mode.SSMIF_STRUCT, charMode='int', endianness='little')
	bsettings = mode.parseCStruct(mode.STATION_SETTINGS_STRUCT, endianness='little')
	
	fh.readinto(bssmif)
	
	#
	# Station Data
	#
	idn = [chr(i) for i in bssmif.sStationID]
	idn = ''.join([i for i in idn if i != '\x00'])
	lat = bssmif.fGeoN
	lon = bssmif.fGeoE
	elv = bssmif.fGeoEl
	
	#
	# Stand & Antenna Data
	#
	stdPos   = [list(i) for i in zip(bssmif.fStdLx, bssmif.fStdLy, bssmif.fStdLz)]
	stdAnt   = list(bssmif.iAntStd)
	stdOrie  = list(bssmif.iAntOrie)
	stdStat  = list(bssmif.iAntStat)
	stdTheta = list(bssmif.fAntTheta)
	stdPhi   = list(bssmif.fAntPhi)
	stdDesi  = list(bssmif.eAntDesi)
	
	#
	# FEE, Cable, & SEP Data
	#
	feeID   = mcs.single2multi([chr(i) for i in bssmif.sFEEID], *bssmif.dims['sFEEID'])
	feeID   = [''.join([k for k in i if k != '\x00']) for i in feeID]
	feeStat = list(bssmif.iFEEStat)
	feeDesi = list(bssmif.eFEEDesi)
	feeGai1 = list(bssmif.fFEEGai1)
	feeGai2 = list(bssmif.fFEEGai2)
	feeAnt1 = list(bssmif.iFEEAnt1)
	feeAnt2 = list(bssmif.iFEEAnt2)
	
	rpdID   = mcs.single2multi([chr(i) for i in bssmif.sRPDID], *bssmif.dims['sRPDID'])
	rpdID   = [''.join([k for k in i if k != '\x00']) for i in rpdID]
	rpdStat = list(bssmif.iRPDStat)
	rpdDesi = list(bssmif.eRPDDesi)
	rpdLeng = list(bssmif.fRPDLeng)
	rpdVF   = list(bssmif.fRPDVF)
	rpdDD   = list(bssmif.fRPDDD)
	rpdA0   = list(bssmif.fRPDA0)
	rpdA1   = list(bssmif.fRPDA1)
	rpdFre  = list(bssmif.fRPDFref)
	rpdStr  = list(bssmif.fRPDStr)
	rpdAnt  = list(bssmif.iRPDAnt)
	
	sepCbl  = mcs.single2multi([chr(i) for i in bssmif.sSEPCabl], *bssmif.dims['sSEPCabl'])
	sepCbl  = [''.join([k for k in i if k != '\x00']) for i in sepCbl]
	sepLeng = list(bssmif.fSEPLeng)
	sepDesi = list(bssmif.eSEPDesi)
	sepGain = list(bssmif.fSEPGain)
	sepAnt  = list(bssmif.iSEPAnt)
	
	#
	# ARX (ARB) Data
	#
	nChanARX = bssmif.nARBCH
	arxID    = mcs.single2multi([chr(i) for i in bssmif.sARBID], *bssmif.dims['sARBID'])
	arxID    = [''.join([k for k in i if k != '\x00']) for i in arxID]
	arxSlot  = list(bssmif.iARBSlot)
	arxDesi  = list(bssmif.eARBDesi)
	arxRack  = list(bssmif.iARBRack)
	arxPort  = list(bssmif.iARBPort)
	arxStat  = mcs.single2multi(bssmif.eARBStat, *bssmif.dims['eARBStat'])
	arxAnt   = mcs.single2multi(bssmif.iARBAnt, *bssmif.dims['iARBAnt'])
	arxIn    = mcs.single2multi([chr(i) for i in bssmif.sARBIN], *bssmif.dims['sARBIN'])
	arxIn    = [[''.join(i) for i in j] for j in arxIn]
	arxOut   = mcs.single2multi([chr(i) for i in bssmif.sARBOUT], *bssmif.dims['sARBOUT'])
	arxOUt   = [[''.join(i) for i in j] for j in arxOut]
	
	try:
		#
		# DP 1 & 2 Data
		#
		dp1ID   = mcs.single2multi([chr(i) for i in bssmif.sDP1ID], *bssmif.dims['sDP1ID'])
		dp1ID   = [''.join([k for k in i if k != '\x00']) for i in dp1ID]
		dp1Slot = mcs.single2multi([chr(i) for i in bssmif.sDP1Slot], *bssmif.dims['sDP1Slot'])
		dp1Slot = [''.join([k for k in i if k != '\x00']) for i in dp1Slot]
		dp1Desi = list(bssmif.eDP1Desi)
		dp1Stat = list(bssmif.eDP1Stat)
		dp1InR  = mcs.single2multi([chr(i) for i in bssmif.sDP1INR], *bssmif.dims['sDP1INR'])
		dp1InR  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in dp1InR]
		dp1InC  = mcs.single2multi([chr(i) for i in bssmif.sDP1INC], *bssmif.dims['sDP1INC'])
		dp1InC  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in dp1InC]
		dp1Ant  = mcs.single2multi(bssmif.iDP1Ant, *bssmif.dims['iDP1Ant'])
		
		dp2ID   = mcs.single2multi([chr(i) for i in bssmif.sDP2ID], *bssmif.dims['sDP2ID'])
		dp2ID   = [''.join([k for k in i if k != '\x00']) for i in dp2ID]
		dp2Slot = mcs.single2multi([chr(i) for i in bssmif.sDP2Slot], *bssmif.dims['sDP2Slot'])
		dp2Slot = [''.join([k for k in i if k != '\x00']) for i in dp2Slot]
		dp2Stat = list(bssmif.eDP2Stat)
		dp2Desi = list(bssmif.eDP2Desi)
	except AttributeError:
		#
		# ROACH & Server Data
		#
		roachID   = mcs.single2multi([chr(i) for i in bssmif.sRoachID], *bssmif.dims['sRoachID'])
		roachID   = [''.join([k for k in i if k != '\x00']) for i in roachID]
		roachSlot = mcs.single2multi([chr(i) for i in bssmif.sRoachSlot], *bssmif.dims['sRoachSlot'])
		roachSlot = [''.join([k for k in i if k != '\x00']) for i in roachSlot]
		roachDesi = list(bssmif.eRoachDesi)
		roachStat = list(bssmif.eRoachStat)
		roachInR  = mcs.single2multi([chr(i) for i in bssmif.sRoachINR], *bssmif.dims['sRoachINR'])
		roachInR  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in roachInR]
		roachInC  = mcs.single2multi([chr(i) for i in bssmif.sRoachINC], *bssmif.dims['sRoachINC'])
		roachInC  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in roachInC]
		roachAnt  = mcs.single2multi(bssmif.iRoachAnt, *bssmif.dims['iRoachAnt'])
		
		serverID   = mcs.single2multi([chr(i) for i in bssmif.sServerID], *bssmif.dims['sServerID'])
		serverID   = [''.join([k for k in i if k != '\x00']) for i in serverID]
		serverSlot = mcs.single2multi([chr(i) for i in bssmif.sServerSlot], *bssmif.dims['sServerSlot'])
		serverSlot = [''.join([k for k in i if k != '\x00']) for i in serverSlot]
		serverStat = list(bssmif.eServerStat)
		serverDesi = list(bssmif.eServerDesi)
		
	#
	# DR Data
	#
	drStat = list(bssmif.eDRStat)
	drID   = mcs.single2multi([chr(i) for i in bssmif.sDRID], *bssmif.dims['sDRID'])
	drID   = [''.join([k for k in i if k != '\x00']) for i in drID]
	drShlf = [0 for i in xrange(bssmif.nDR)]
	drPC   = mcs.single2multi([chr(i) for i in bssmif.sDRPC], *bssmif.dims['sDRPC'])
	drPC   = [''.join([k for k in i if k != '\x00']) for i in drPC]
	drDP   = list(bssmif.iDRDP)
	
	fh.readinto(bsettings)
	
	fh.close()
	
	return locals()


def parseSSMIF(filename):
	"""
	Given a SSMIF file, return a fully-filled LWAStation instance.  This function
	supports both human-readable files (filenames with '.txt' extensions) or 
	binary packed files (filenames with '.dat' extensions).
	"""
	
	# Find out if we have a .txt or .dat file and process accordingly
	base, ext = os.path.splitext(filename)
	
	# Read in the ssmif to a dictionary of variables
	if ext == '.dat':
		ssmifDataDict = __parseBinarySSMIF(filename)
	elif ext == '.txt':
		ssmifDataDict = __parseTextSSMIF(filename)
	else:
		raise ValueError("Unknown file extension '%s', cannot tell if it is text or binary" % ext)
	
	# Unpack the dictionary into the current variable scope
	## Site
	idn = ssmifDataDict['idn']
	lat = ssmifDataDict['lat']
	lon = ssmifDataDict['lon']
	elv = ssmifDataDict['elv']
	## Stands
	stdPos = ssmifDataDict['stdPos']
	## FEEs
	feeID   = ssmifDataDict['feeID']
	feeGai1 = ssmifDataDict['feeGai1']
	feeGai2 = ssmifDataDict['feeGai2']
	feeStat = ssmifDataDict['feeStat']
	feeAnt1 = ssmifDataDict['feeAnt1']
	feeAnt2 = ssmifDataDict['feeAnt2']
	## Cables
	rpdID   = ssmifDataDict['rpdID']
	rpdLeng = ssmifDataDict['rpdLeng']
	rpdVF   = ssmifDataDict['rpdVF']
	rpdDD   = ssmifDataDict['rpdDD']
	rpdA0   = ssmifDataDict['rpdA0']
	rpdA1   = ssmifDataDict['rpdA1']
	rpdStr  = ssmifDataDict['rpdStr']
	rpdFre  = ssmifDataDict['rpdFre']
	rpdAnt  = ssmifDataDict['rpdAnt']
	## Antennas
	stdAnt   = ssmifDataDict['stdAnt']
	stdOrie  = ssmifDataDict['stdOrie']
	stdTheta = ssmifDataDict['stdTheta']
	stdPhi   = ssmifDataDict['stdPhi']
	stdStat  = ssmifDataDict['stdStat']
	## ARX
	nChanARX = ssmifDataDict['nChanARX']
	arxID    = ssmifDataDict['arxID']
	arxAnt   = ssmifDataDict['arxAnt']
	arxIn    = ssmifDataDict['arxIn']
	arxOut   = ssmifDataDict['arxOut']
	### DP/ADP
	try:
		isDP = True
		dp1Ant = ssmifDataDict['dp1Ant']
		dp1InR = ssmifDataDict['dp1InR']
	except KeyError:
		isDP = False
		roachAnt = ssmifDataDict['roachAnt']
		roachInR = ssmifDataDict['roachInR']
		
	# Build up a list of Stand instances and load them with data
	i = 1
	stands = []
	for pos in stdPos:
		stands.append(Stand(i, *pos))
		i += 1
	
	# Build up a list of FEE instances and load them with data
	i = 1
	fees = []
	for id,gain1,gain2,stat in zip(feeID, feeGai1, feeGai2, feeStat):
		fees.append(FEE(id, i, gain1=gain1, gain2=gain2, status=stat))
		i += 1
	
	# Build up a list of Cable instances and load them with data
	i = 1
	cables = []
	for id,length,vf,dd,a0,a1,stretch,aFreq in zip(rpdID, rpdLeng, rpdVF, rpdDD, rpdA0, rpdA1, rpdStr, rpdFre):
		cables.append(Cable(id, length, vf=vf/100.0, dd=float(dd)*1e-9, a0=a0, a1=a1, aFreq=aFreq, stretch=stretch))
		i += 1
	
	# Build up a list of Antenna instances and load them with antenna-level
	# data
	i = 1
	antennas = []
	for ant,pol,theta,phi,stat in zip(stdAnt, stdOrie, stdTheta, stdPhi, stdStat):
		antennas.append(Antenna(i, stand=stands[ant-1], pol=pol, theta=theta, phi=phi, status=stat))
		i += 1
	
	# Associate FEEs with Antennas and set the FEE port numbers
	i = 1
	for fee,ant in zip(fees, feeAnt1):
		antennas[ant-1].fee = fee
		antennas[ant-1].feePort = 1
		i += 1
	i = 1
	for fee,ant in zip(fees, feeAnt2):
		antennas[ant-1].fee = fee
		antennas[ant-1].feePort = 2
		i += 1
	
	# Associate Cables with Antennas
	i = 1
	for cbl,ant in zip(cables, rpdAnt):
		antennas[ant-1].cable = cbl
		i += 1
	
	# Associate ARX boards/channels with Antennas
	for i in xrange(len(arxAnt)):
		for j in xrange(len(arxAnt[i])):
			ant = arxAnt[i][j]
			if ant == 0 or ant > 520:
				continue
			
			boardID = arxID[i]
			channel = j + 1
			antennas[ant-1].arx = ARX(boardID, channel=channel, aspChannel=i*nChanARX + j + 1, input=arxIn[i][j], output=arxOut[i][j])
			
	if isDP:
		# Associate DP 1 board and digitizer numbers with Antennas - DP1 boards are 2-14 and 16-28 
		# with DP2 boards at 1 and 15.
		i = 1
		j = 1
		for brd,inp in zip(dp1Ant,dp1InR):
			for ant,con in zip(brd,inp):
				antennas[ant-1].board = i + 1 + (i/14)
				antennas[ant-1].digitizer = j
				antennas[ant-1].input = con
				j += 1
			i += 1
	else:
		# Associate ROACH board and digitizer numbers with Antennas.
		i = 1
		j = 1
		for brd,inp in zip(roachAnt,roachInR):
			for ant,con in zip(brd,inp):
				antennas[ant-1].board = i
				antennas[ant-1].digitizer = j
				antennas[ant-1].input = con
				j += 1
			i += 1
			
	# Build a Station
	try:
		## LSL interface support
		if idn in ('VL',):
			interface = LSLInterface(backend='lsl.common.dp', 
								mcs='lsl.common.mcs', 
								sdf='lsl.common.sdf', 
								metabundle='lsl.common.metabundle', 
								sdm='lsl.common.sdm')
		elif idn in ('SV',):
			interface = LSLInterface(backend='lsl.common.adp', 
								mcs='lsl.common.mcsADP', 
								sdf='lsl.common.sdfADP', 
								metabundle='lsl.common.metabundleADP', 
								sdm='lsl.common.sdmADP')
		else:
			interface = None
			
		station = LWAStation(_id2name[idn], lat, lon, elv, id=idn, antennas=antennas, interface=interface)
	except KeyError:
		station = LWAStation('New LWA Station', lat, lon, elv, id=idn, antennas=antennas)
		
	# And return it
	return station


# LWAVL
_ssmifvl = os.path.join(dataPath, 'lwa1-ssmif.txt')
lwavl = parseSSMIF(_ssmifvl)

# LWAVL is also known as LWA1
lwa1 = lwavl

# LWANA
_ssmifna = os.path.join(dataPath, 'lwana-ssmif.txt')
lwana = parseSSMIF(_ssmifna)

# LWASV
_ssmifsv = os.path.join(dataPath, 'lwasv-ssmif.txt')
lwasv = parseSSMIF(_ssmifsv)


def getFullStations():
	"""
	Function to return a list of full stations.
	
	.. versionadded:: 1.2.0
	"""
	
	return [lwa1, lwasv]


class PrototypeStation(LWAStation):
	"""
	LWAStation class to provide backward compatiability to the 20-stand prototype 
	system.
	"""
	
	def __init__(self, base):
		pLat = base.lat * 180 / numpy.pi
		pLng = base.long * 180 / numpy.pi
		
		super(PrototypeStation, self).__init__('%s prototype' % base.name, pLat, pLng, base.elev, id='PS', antennas=copy.deepcopy(base.antennas))
	
	def __standsList(self, date):
		"""
		For a given datetime object, return a list of two-element tuples (stand ID 
		and polarization) for what antennas are connected to what digitizers.
		"""
		
		from datetime import datetime
		
		if date < datetime(2010, 11, 23):
			return [(258,0), (4,0), (158,0), (205,0), (246,0), (9,0), (69,0), (168,0), (80,0), (14,0), 
					(254,0), (118,0), (38,0), (34,0), (67,0), (181,0), (206,0), (183,0), (153,0), (174,0)]
		elif date < datetime(2011, 1, 12):
			return [(206,0), (183,0), (153,0), (174,0), (38,0), (34,0), (67,0), (181,0), (80,0), (14,0), 
					(254,0), (118,0), (246,0), (9,0), (69,0), (168,0), (258,0), (4,0), (158,0), (205,0)]
		elif date < datetime(2011, 2, 11):
			return [(206,0), (183,0), (153,0), (174,0), (38,0), (34,0), (67,0), (187,0), (80,0), (14,0), 
					(254,0), (118,0), (246,0), (9,0), (69,0), (168,0), (258,0), (4,0), (158,0), (205,0)]
		elif date < datetime(2011, 2, 22):
			return [(157,0), (157,1), (155,0), (155,1), (185,0), (185,1), (159,0), (159,1), (63,0), (63,1), 
					(67,0), (67,1), (65,0), (65,1), (127,0), (127,1), (123,0), (123,1), (125,0), (125,1)]
		elif date < datetime(2011, 3, 12):
			return [(183,0), (183,1), (200,0), (200,1), (212,0), (212,1), (14,0), (14,1), (210,0), (210,1), 
					(80,0), (80,1), (228,0), (228,1), (15,0), (15,1)]
		else:
			return [(183,0), (183,1), (200,0), (200,1), (35,0), (35,1), (108,0), (108,1), (99,0), (99,1), 
					(97,0), (97,1), (8,0), (8,1), (168,0), (168,1), (129,0), (129,1), (155,0), (155,1)]
		
	def __filterAntennas(self, date):
		"""
		For a given datetime object, return a list of Antenna instances with the 
		correct digitizer mappings.
		"""
		
		# Build the sorting function
		def sortFcn(x, y):
			if x.digitizer > y.digitizer:
				return 1
			elif x.digitizer < y.digitizer:
				return -1
			else:
				return 0
		
		# Get the list of stand,polarization tuples that we are looking for
		stdpolList = self.__standsList(date)
		
		# Build up the list of actual antennas connected to the prototype system 
		# and fix the digitizer numbers
		protoAnts = []
		for ant in self.antennas:
			stdpolPair = (ant.stand.id,ant.pol)
			if stdpolPair in stdpolList:
				protoAnts.append(ant)
				protoAnts[-1].digitizer = stdpolList.index(stdpolPair) + 1
				protoAnts[-1].board = 1
				
		# Sort the list of prototype antennas and return
		protoAnts.sort(cmp=sortFcn)
		return protoAnts
	
	def getAntennas(self, date):
		"""Return the list of Antenna instances for the station, sorted by 
		digitizer number, for a given DateTime instance.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return protoAnts
		
	def getStands(self, date):
		"""
		Return a list of Stand instances for each antenna, sorted by 
		digitizer number, for a given DateTime instance.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return [ant.stand for ant in protoAnts]
	
	def getPols(self, date):
		"""
		Return a list of polarization (0 == N-S; 1 == E-W) for each antenna, 
		sorted by digitizer number, for a given DateTime instance.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return [ant.pol for ant in protoAnts]
		
	def getCables(self, date):
		"""
		Return a list of Cable instances for each antenna, sorted by
		digitizer number, for a given DateTime instance.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return [ant.cable for ant in protoAnts]


prototypeSystem = PrototypeStation(lwa1)
