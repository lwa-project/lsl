#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import numpy
import ephem

from lsl.common.paths import data as dataPath
from lsl.common.constants import *

__version__ = "0.5"
__revision__ = "$ Revision: 27 $"
__all__ = ['status2string', 'geo2ecef', 'LWAStation', 'Antenna', 'Stand', 'FEE', 'Cable', 'parseSSMIF', 'lwa1', 'PrototypeStation', 'prototypeSystem', '__version__', '__revision__', '__all__']


_id2name = {'VL': 'LWA-1'}

def status2string(code):
	"""
	Convert a numerical SSMIF status code to a string.
	"""
	
	# Make sure we have an integer
	code = int(code)
	
	# Loop through the options
	if code == 0:
		return "Not installed"
	elif code == 1:
		return "Bad"
	elif code == 2:
		return "Suspect, possibly bad"
	elif code == 3:
		return "OK"
	else:
		return "Unknown status code '%i'" % code


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


class LWAStation(object):
	"""
	Object to hold information about the a LWA station.  This object can
	create a ephem.Observer representation of itself and identify which stands
	were in use at a given time.  Stores station:
	  * Name (name)
	  * ID code (id)
	  * Latitiude in radians [but initialized as degrees] (N is positive, lat)
	  * Longitude in radians [but initialized as degrees] (W is negative, lon)
	  * Elevation in meters (elev)
	  * List of Antenna instances (antennas)
	  
	LWAStation also provides several functional attributes for dealing with
	the station's location on Earth.  These include:
	  * getObserver: Return an ephem.Observer instance representing the station
	  * getAIPYLocation: Return a tuple for setting the location of an AIPY
	    AntennaArray instance
	  * getGeocentricLocation: Return a tuple of the EC-EF coordinates of the 
	    station
	  * getECEFTransform: Return a 3x3 tranformation matrix to convert antenna
	    positions to EC-EF coordinates
	"""

	def __init__(self, name, lat, long, elev, id='', antennas=None):
		self.name = name
		self.id = id
		self.lat = lat*numpy.pi/180.0
		self.long = long*numpy.pi/180.0
		self.elev = elev
		
		if antennas is None:
			self.antennas = []
		else:
			self.antennas = antennas

	def getObserver(self, date=None, JD=False):
		"""
		Return a ephem.Observer object for this site.
		"""

		oo = ephem.Observer()
		oo.lat = self.lat
		oo.long = self.long
		oo.elev = self.elev
		oo.pressure = 0.0
		if date is not None:
			if JD:
				# If the date is Julian, convert to Dublin Julian Date 
				# which is used by ephem
				date -= 2415020.0
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

	def getECEFTransform(self):
		"""
		Return a 3x3 tranformation matrix that converts a baseline in 
		[east, north, elevation] to earh-centered, earth-fixed coordinates
		for that baseline [x, y, z].  Based off the 'local_to_eci' function
		in the lwda_fits-dev library.
		"""

		return numpy.array([[0.0, -numpy.sin(self.lat), numpy.cos(self.lat)], 
						[1.0, 0.0,                  0.0], 
						[0.0, numpy.cos(self.lat),  numpy.sin(self.lat)]])
						
	def __sortAntennas(self, attr='digitizer'):
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
		
	def getAntennas(self):
		"""
		Return the list of Antenna instances for the station, sorted by 
		digitizer number.
		"""
		
		# Sort and return
		self.__sortAntennas()
		return self.antennas
		
	def getStands(self):
		"""
		Return a list of Stand instances for each antenna, sorted by 
		digitizer number.
		"""
		
		# Sort and return
		self.__sortAntennas()
		return [ant.stand for ant in self.antennas]
	
	def getPols(self):
		"""
		Return a list of polarization (0 == N-S; 1 == E-W) for each antenna, 
		sorted by digitizer number.
		"""
		
		# Sort and return
		self.__sortAntennas()
		return [ant.pol for ant in self.antennas]
		
	def getCables(self):
		"""
		Return a list of Cable instances for each antenna, sorted by
		digitizer number.
		"""
		
		# Sort and return
		self.__sortAntennas()
		return [ant.cable for ant in self.antennas]


class Antenna(object):
	"""
	Object to store the information about an antenna.  Stores antenna:
	  * ID number (id)
	  * DP1 board number (board)
	  * DP1 digitizer number (digiziter)
	  * Stand instance the antenna is part of (stand)
	  * Polarization (0 == N-S; pol)
	  * Antenna mis-alignment in degrees (phi)
	  * Fee instance the antenna is attached to (fee)
	  * Port of the FEE used for the antenna (feePort)
	  * Cable instance used to connect the antenna (cable)
	  * Status of the antenna (status)
	  
	Status codes are:
	  * 0 == Not installed
	  * 1 == Bad
	  * 2 == Suspect, possibly bad
	  * 3 == OK
	"""
	
	def __init__(self, id, board=0, digitizer=0, stand=None, pol=0, theta=0.0, phi=0.0, fee=None, feePort=1, cable=None, status=0):
		self.id = int(id)
		self.board = int(board)
		self.digitizer = int(digitizer)
		
		if stand is None:
			self.stand = Stand(0, 0, 0, 0)
		else:
			self.stand = stand
			
		self.pol = int(pol)
		self.theta = float(theta)
		self.phi = float(phi)
		
		if fee is None:
			self.fee = FEE('')
		else:
			self.fee = fee
		self.feePort = feePort
		
		if cable is None:
			self.cable = Cable('', 0)
		else:
			self.cable = cable
			
		self.status = int(status)
		
	def __cmp__(self, y):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0
		
	def __str__(self):
		return "Antenna %i: stand=%i, polarization=%i; digitizer %i; status is %i" % (self.id, self.stand.id, self.pol, self.digitizer, self.status)
		
	def getStatus(self):
		"""
		Return the combined antenna + FEE status as a two digit number 
		with the first digit representing the antenna status and the 
		second the FEE status.
		"""
		
		return 10*self.status + self.fee.status


class Stand(object):
	"""
	Object to store the information (location and ID) about a stand.  
	Stores stand:
	  * ID number (id)
	  * Position relative to the center stake in meters (x,y,z)
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


class FEE(object):
	"""
	Object to store the information about a FEE.  Stores FEE:
	  * ID name (id)
	  * Gain of port 1 (gain1)
	  * Gain of part 2 (gain2)
	  * Status (status)
	  
	Status codes are:
	  * 0 == Not installed
	  * 1 == Bad
	  * 2 == Suspect, possibly bad
	  * 3 == OK
	"""
	
	def __init__(self, id, gain1=0, gain2=0, status=0):
		self.id = str(id)
		self.gain1 = float(gain1)
		self.gain2 = float(gain2)
		self.status = int(status)
	
	def __cmp__(self, y):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0
	
	def __str__(self):
		return "Fee '%s': gain1=%.2f, gain2=%.2f; status is %i" % (self.id, self.gain1, self.gain2, self.status)


class Cable(object):
	"""
	Object to store information about a cable.  Stores cable:
	  * ID name (id)
	  * Length in meters (length)
	  * Velocity factor (fractional, vf)
	  * Dispersive delay (seconds, dd)
	
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

	def __cmp__(self):
		if self.id > y.id:
			return 1
		elif self.id < y.id:
			return -1
		else:
			return 0

	def __str__(self):
		return "Cable '%s' with length %.2f m (stretched to %.2f m)" % (self.id, self.length, self.length*self.stretch)
		
	def delay(self, frequency=49e6, ns=False):
		"""Get the delay associated with the cable in second (or nanoseconds 
		if the 'ns' keyword is set to True) for a particular frequency or 
		collection of frequencies in Hz."""
		
		# Bulk delay for the cable
		bulkDelay = self.length*self.stretch / (self.vf * c)
		
		# Dispersize delay in 
		dispDelay = self.dd * (self.length*self.stretch / 100.0) * numpy.sqrt(frequency / numpy.array(10e6))
		
		totlDelay = bulkDelay + dispDelay
		
		if ns:
			return totlDelay*1e9
		else:
			return totlDelay

	def attenuation(self, frequency=49e6):
		"""Get the multiplicative factor needed to correct for the cable 
		loss for a specific frequency (in Hz).  If attenuations for more 
		than one frequency are needed, the frequencies can be passed in as 
		a numpy array.

		.. note::
			This function assumes a cable of type KSR200DB from SLC0014, version 3
		"""
	
		atten = 2 * self.a0 * self.length*self.stretch * numpy.sqrt(frequency / numpy.array(self.aFreq)) 
		atten += self.a1 * self.length*self.stretch * (frequency / numpy.array(self.aFreq))
		atten = numpy.exp(atten)
		return atten
		
	def gain(self, frequency=49e6):
		"""Get the cable gain ("inverse loss") for a specific frequency (in 
		Hz).  If gains for more than one frequency are needed, the 
		frequencies can be passed in as a numpy array.

		.. note::
			This function assumes a cable of type KSR200DB from SLC0014, version 3
		"""
		
		return 1.0 / self.attenuation(frequency=frequency)


def parseSSMIF(filename):
	"""
	Given a SSMIF file, return a fully-filled LWAStation instance.
	"""
	
	fh = open(filename, 'r')
	
	kwdRE = re.compile(r'(?P<keyword>[A-Z_0-9]+)(\[(?P<id1>[0-9]+?)\])?(\[(?P<id2>[0-9]+?)\])?(\[(?P<id3>[0-9]+?)\])?')
	
	# Loop over the lines in the file
	for line in fh:
		line = line.replace('\n', '')
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
		if keyword == 'RDP_DD':
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
		
		if keyword == 'RDP_DESI':
			rpdDesi[ids[0]-1] = value
			continue
		
		if keyword == 'RDP_ANT':
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
		# ARX Data
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
		# DP 1 & 2 Data
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
			dp1Slot[ids[0]-1] = int(value)
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
		fees.append(FEE(id, gain1=gain1, gain2=gain2, status=stat))
	
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
	
	
	# Associate DP 1 board and digitizer numbers with Antennas
	i = 1
	j = 1
	for brd in dp1Ant:
		for ant in brd:
			antennas[ant-1].board = i
			antennas[ant-1].digitizer = j
			j += 1
		i += 1
		
	# Build a Station
	try:
		station = LWAStation(_id2name[idn], lat, lon, elv, id=idn, antennas=antennas)
	except KeyError:
		station = LWAStation('New LWA Station', lat, lon, elv, id=idn, antennas=antennas)
	
	# And return it
	return station


_ssmif = os.path.join(dataPath, 'lwa1-ssmif.txt')
lwa1 = parseSSMIF(_ssmif)


class PrototypeStation(LWAStation):
	"""
	LWAStation class to provide backward compatiability to the 20-stand prototype 
	system.
	"""
	
	def __init__(self, base):
		pLat = base.lat * 180 / numpy.pi
		pLng = base.long * 180 / numpy.pi
	
		super(PrototypeStation, self).__init__('%s prototype' % base.name, pLat, pLng, base.elev, id='', antennas=base.antennas)
	
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
		digitizer number.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return protoAnts
		
	def getStands(self):
		"""
		Return a list of Stand instances for each antenna, sorted by 
		digitizer number.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return [ant.stand for ant in protoAnts]
	
	def getPols(self):
		"""
		Return a list of polarization (0 == N-S; 1 == E-W) for each antenna, 
		sorted by digitizer number.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return [ant.pol for ant in protoAnts]
		
	def getCables(self):
		"""
		Return a list of Cable instances for each antenna, sorted by
		digitizer number.
		"""
		
		# Sort and return
		protoAnts = self.__filterAntennas(date)
		return [ant.cable for ant in protoAnts]


prototypeSystem = PrototypeStation(lwa1)
