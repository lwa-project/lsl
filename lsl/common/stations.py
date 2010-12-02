# -*- coding: utf-8 -*-

"""Module that contains information about the LWA station sites.  This 
information includes::
  * latitude, 
  * longitude, 
  * elevation, and 
  * stand list.

The information for each site is stored as a LWAStation object.  This object 
has a variety of attributes for dealing with station-related parameters, such
as creating an ephem.Observer object for the station, getting the station's
geocentric location, and creating a geocentric coordinate transformation 
matrix for the station's baselines.

Currently, data exists for the following LWA stations::
  * lwa1 - LWA 1 near the VLA center.
"""

import math
import ephem
import numpy


__version__ = '0.2'
__revision__ = '$ Revision: 10 $'
__all__ = ['geo2ecef', 'LWAStation', 'lwa1', '__version__', '__revision__', '__all__']


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

		if fDate >= ephem.date('2010/11/23'):
			# Current as of 11/30/2010
			stands = numpy.array([206, 183, 174, 153, 38, 34, 181, 67, 80, 14, 118, 254, 246, 9, 168, 69, 258, 4, 205, 158])
		elif fDate >= ephem.Date('2010/10/29'):
			# Current as of 10/29/2010 for both TBN and TBW
			stands = numpy.array([258, 4, 158, 205, 246, 9, 69, 168, 80, 14, 254, 118, 38, 34, 67, 181, 206, 183, 153, 174])
		elif fDate >= ephem.Date('2010/09/15'):
			# Current as of 9/15/2010
			# Note:  This is for the input reversed TBW data.  Any TBN data should have 
			# the stand numbers exchanged in pairs.  For example, TBN data would be:
			# 258, 4, 158, 205, etc.
			stands = numpy.array([4, 258, 205, 158, 9, 246, 168, 69, 14, 80, 118, 254, 34, 38, 181, 67, 183, 206, 174, 153])
		elif fDate >= ephem.Date('2010/08/25'):
			# Between 8/25/2010 and 9/15/2010
			stands = numpy.array([214, 212, 228, 206, 127, 157, 187, 208, 123, 125])
		else:
			# Before 8/25/2010
			stands = numpy.array([214, 212, 228, 204, 225, 210, 187, 174, 123, 125])
		
		return stands