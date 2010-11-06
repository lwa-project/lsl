# -*- coding: utf-8 -*-

"""Module that contains common values about the LWA station sites."""

import math
import ephem
import numpy


__version__ = '0.2'
__revision__ = '$ Revision: 8 $'
__all__ = ['geo2ecef', 'LWAStation', 'lwa1']


def geo2ecef(lat, lon, elev):
	"""Convert latitude (rad), longitude (rad), elevation (m) to earth-
	centered, earth-fixed coordinates."""

	WGS84_a = 6378137.000000
	WGS84_b = 6356752.314245
	N = WGS84_a**2 / math.sqrt(WGS84_a**2*math.cos(lat)**2 + WGS84_b**2*math.sin(lon)**2)

	x = (N+elev)*math.cos(lat)*math.cos(lon)
	y = (N+elev)*math.cos(lat)*math.sin(lon)
	z = ((WGS84_b**2/WGS84_a**2)*N+elev)*math.sin(lat)

	return (x, y, z)


class LWAStation(object):
	"""Object to hold information about the a LWA station.  This object can
	create a ephem.Observer representation of itself and identify which stands
	were in use at a given time."""

	def __init__(self, name, lat, long, elev):
		self.name = name
		self.lat = lat*math.pi/180.0
		self.long = long*math.pi/180.0
		self.elev = elev

	def getObserver(self, date=None, JD=False):
		"""Return a ephem.Observer object for this site."""

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
		"""Return a tuple that can be used by AIPY for specifying a array
		location."""

		return (self.lat, self.long, self.elev)

	def getGeocentricLocation(self):
		"""Return a tuple with earth-centered, earth-fixed coordinates for the station."""

		return geo2ecef(self.lat, self.long, self.elev)

	def getECEFTransform(self):
		"""Return a 3x3 tranformation matrix that converts a baseline in 
		[east, north, elevation] to earh-centered, earth-fixed coordinates
		for that baseline [x, y, z].  Based off the 'local_to_eci' function
		in the lwda_fits-dev library."""

		return numpy.array([[0.0, -math.sin(self.lat), math.cos(self.lat)], 
						[1.0, 0.0,                 0.0], 
						[0.0, math.cos(self.lat), math.sin(self.lat)]])


class lwa1(LWAStation):
	"""Object to hold information about the first LWA station.  This object can
	create a ephem.Observer representation of itself and identify which stands
	were in use at a given time."""

	def __init__(self):
		super(lwa1, self).__init__('LWA-1', 34.070, -107.628, 2133.6)

	def getStands(self, date=None, JD=False):
		"""Return a numpy array that contains the stands used at a given 
		point in time."""

		oo = self.getObserver(date=date, JD=JD)
		fDate = float(oo.date)

		if fDate >= 40475.0:
			# Current as of 10/29/2010 for both TBN and TBW
			stands = numpy.array([-1, 4, 158, 205, 246, 9, 69, 168, 80, 14, 254, 118, 38, 34, 67, 181, 206, 183, 153, 174])
		elif fDate >= 40435.0:
			# Current as of 9/15/2010
			# Note:  This is for the input reversed TBW data.  Any TBN data should have 
			# the stand numbers exchanged in pairs.  For example, TBN data would be:
			# -1, 4, 158, 205, etc.
			stands = numpy.array([4, -1, 205, 158, 9, 246, 168, 69, 14, 80, 118, 254, 34, 38, 181, 67, 183, 206, 174, 153])
		elif fDate >= 40414.0:
			# Between 8/25/2010 and 9/15/2010
			stands = numpy.array([214, 212, 228, 206, 127, 157, 187, 208, 123, 125])
		else:
			# Before 8/25/2010
			stands = numpy.array([214, 212, 228, 204, 225, 210, 187, 174, 123, 125])
		
		return stands