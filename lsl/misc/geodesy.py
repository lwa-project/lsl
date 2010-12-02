# -*- coding: utf-8 -*-

"""Module for querying the earth orientation parameters for a given date/list
of dates."""

import os
import gzip
import ephem
import numpy
import urllib

import lsl.astro as astro
import lsl.common.paths as paths

__version__ = '0.1'
__revision__ = '$ Revision: 3 $'
__all__ = ['EOP', 'getEOP', 'getEOPRange', '__version__', '__revision__', '__all__']

# Last valid MJD in the historic EOP data included with LSL
# MJD 55518 == November 18, 2010.
Historic1973Limit = 55518.0


class EOP(object):
	"""Object for storing the geodetic parameters relevant for DiFX input 
	files:
	  * mjd - modified Julian Date of the measurement/prediction
	  * x - difference between the celestial ephemeric pole (CEP) and the 
	    international reference pole (IRP) in the direction of the IERS
	    meridian [arc seconds]
	  * y - difference between the CEP and IRP in the direction of 90 degrees 
	    west longitude [arc seconds]
	  * UT1-UTC - difference between rotation angle about the pole and UTC
	    [seconds]
	  * type - whether the values for the given MJD are observed (final/IERS) or
	    predicted.
	"""

	def __init__(self, mjd=0.0, x=0.0, y=0.0, utDiff=0.0, type='final'):
		self.mjd = mjd
		self.x = x
		self.y = y
		self.utDiff = utDiff
		self.type = type

		self.__setDate()

	def fromMAIA(self, line):
		"""Given a line from a MAIA standard rapiod EOP data (IAU2000) file, fill
		in the object's values with the needed information."""

		self.mjd = float(line[7:15])
		self.x = float(line[18:27])
		self.y = float(line[38:46])
		self.utDiff = float(line[58:68])
		if line[57] == 'I':
			self.type = 'final'
		else:
			self.type = 'prediction'

		self.__setDate()

	def __setDate(self):
		"""Use the ephem.Data object to get an easy-to-use date into the structure."""

		# Catch statement for working is version of LSL older than 0.3
		try:
			self.date = ephem.Date(self.mjd + astro.MJD_OFFSET - astro.DJD_OFFSET)
		except AttributeError:
			self.date = ephem.Date(self.mjd + astro.MJD_OFFSET - 2415020.0)

	def __str__(self):
		"""Create a string representation of the EOP object that shows the MJD, x, 
		y, and UT1-UTC values."""

		return "%.1f (%s): x=%.6f y=%.6f UT1-UTC=%.6f (%s)" % (self.mjd, str(self.date), self.x, self.y, self.utDiff, self.type)

	def __eq__(self, y):
		"""Determine if MJDs of two EOP objects are equal, or if the MJD of a EOP 
		object equal data of a numeric MJD."""

		tX = self.mjd
		try:
			tY = y.mjd
		except:
			tY = float(y)

		if tX == tY:
			return True
		else:
			return False

	def __cmp__(self, y):
		"""Method for soring EOP objects based on their MJDs."""

		tX = float(self.date)
		try:
			tY = float(y.date)
		except AttributeError:
			tY = float(y)
		if tY > tX:
			return -1
		elif tX > tY:
			return 1
		else:
			return 0


def __loadHistoric1973():
	"""Load in historical values.  The file included with LSL contains values 
	from January 2, 1973 to November 18, 2010."""

	# Open the file and read in via the gzip.GzipFile object
	heopFH = gzip.GzipFile(os.path.join(paths.data, 'astro', 'eop-historical.dat.gz'), 'rb')
	lines = heopFH.readlines()
	heopFH.close()

	heops = []
	for line in lines:
		try:
			newEOP = EOP()
			newEOP.fromMAIA(line)
			# Only include "final" values, not predictions
			if newEOP.type == 'final':
				heops.append(newEOP)
		except ValueError:
			pass
	
	return heops


def __loadHistoric1992():
	"""Load in historical values from the web.  The downloaded file includes 
	values from January 1, 1992 until today (usually)."""

	eopFH = urllib.urlopen('http://maia.usno.navy.mil/ser7/finals2000A.daily')
	lines = eopFH.readlines()
	eopFH.close()

	eops = []
	for line in lines:
		newEOP = EOP()
		newEOP.fromMAIA(line) 
		# Only include "final" values, not predictions
		if newEOP.type == 'final':
			eops.append(newEOP)
	
	return eops


def __loadCurrent90():
	"""Load data for the current 90-day period from MAIA via the web."""

	eopFH = urllib.urlopen('http://maia.usno.navy.mil/ser7/finals2000A.daily')
	lines = eopFH.readlines()
	eopFH.close()

	eops = []
	for line in lines:
		newEOP = EOP()
		newEOP.fromMAIA(line) 
		eops.append(newEOP)
	
	return eops


def getEOP(mjd=None):
	"""Return a list of earth orientation parameter objects for the specified 
	MJDs.  A MJD of 'None' returns the values for today's date."""
	
	try:
		junk = len(mjd)
	except TypeError:
		if mjd is None:
			try:
				mjd = [int(float(ephem.now()) + astro.DJD_OFFSET - astro.MJD_OFFSET)]
			except AttributeError:
				mjd = [int(float(ephem.now()) + 2415020.0 - astro.MJD_OFFSET)]
		else:
			mjd = [int(mjd)]
	mjd = numpy.array(mjd)

	oldEOPs = __loadHistoric1973()
	if mjd.max() > Historic1973Limit:
		newEOPs = __loadCurrent90()
	else:
		newEOPs = []
	if mjd.min() > Historic1973Limit and mjd.min() < newEOPs[0]:
		oldEOPs.extend(__loadHistoric1992())
		oldEOPs.sort()

	outEOPs = []
	for day in mjd:
		if day in oldEOPs:
			outEOPs.append(oldEOPs[oldEOPs.index(day)])
		elif day in newEOPs:
			outEOPs.append(newEOPs[newEOPs.index(day)])
		else:
			outEOPs.append(None)
	
	return outEOPs
		

def getEOPRange(start=None, stop=None):
	"""Return a list of orientation parameter objects that span the start and 
	stop (inclusive) MJDs provided.  Values of 'None' indicate today's date."""

	if start is None:
		try:
			start = int(float(ephem.now()) + astro.DJD_OFFSET - astro.MJD_OFFSET)
		except AttributeError:
			start = int(float(ephem.now()) + 2415020.0 - astro.MJD_OFFSET)
	if stop is None:
		try:
			stop = int(float(ephem.now()) + astro.DJD_OFFSET - astro.MJD_OFFSET)
		except AttributeError:
			stop = int(float(ephem.now()) + 2415020.0 - astro.MJD_OFFSET)

	mjdList = numpy.arange(start, stop+1)
	return getEOP(mjdList)
