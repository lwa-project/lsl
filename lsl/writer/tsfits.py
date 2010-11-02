# -*- coding: utf-8 -*-

"""Modules to take TBW/TBN time series data and write it to a FITS file composed 
of binary tables.  This module is still under active development."""

import os
import sys
import time
import numpy
import pyfits
from datetime import datetime, timedelta, tzinfo

from ..common import dp as dp_common
from ..reader.tbw import filterCodes as tbnCodes
from ..reader.warnings import warnDeprecated

__version__ = '0.5'
__revision__ = '$ Revision: 15 $'
__all__ = ['UTC', 'TSFITS', 'TBW', 'TBN', '__version__', '__revision__', '__all__']


class UTC(tzinfo):
    """tzinfo object for UTC time."""

    def utcoffset(self, dt):
        return timedelta(0)

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return timedelta(0)


class TSFITS(object):
	def __init__(self, filename, mode, Overwrite=False, UseQueue=True):
		assert(mode in ['TBW', 'TBN'])

		self.filename = filename
		self.mode = mode
		self.hdulist = None
		self.UseQueue = UseQueue
		self.queue = {}
		self.queueLimit = 10000
		self.site = 'LWA'
		self.firstSamples = {}

		if os.path.exists(filename) and not Overwrite:
			self.hdulist = pyfits.open(self.filename, mode="update", memmap=0)
			self.standCount = self.hdulist[0].header['nstand']
		else:
			if os.path.exists(filename):
				os.path.delete(filename)

			self.standCount = 0
			primary = pyfits.PrimaryHDU()
			primary.header.update('OBJECT', 'zenith')
			primary.header.update('TELESCOP', self.site)
			primary.header.update('OBSMODE', self.mode)
			primary.header.update('NSTAND', self.standCount)

			hdulist = pyfits.HDUList([primary])
			hdulist.writeto(filename)

			self.hdulist = pyfits.open(self.filename, mode="update", memmap=0)

	def info(self):
		self.hdulist.info()

	def flush(self):
		self.hdulist.flush()

	def setSite(self, site):
		self.site = site.name

		self.hdulist[0].header.update('TELESCOP', self.site)

		self.flush()

	def __findExtension(self, stand):
		extension = None

		i = 1
		for hdu in self.hdulist[1:]:
			if int(hdu.header['stand']) == int(stand):
				extension = i
				break
			i = i + 1

		return extension

	def __makeAppendTable(self, extension, AddRows=1):
		"""Private function to make a temporary table for appending data."""

		nrows = self.hdulist[extension].data.shape[0]
		tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+AddRows)
		tempHDU.header = self.hdulist[extension].header.copy()
	
		return tempHDU

	def __applyAppendTable(self, extension, tempHDU):
		"""Private function to replace the given extension with the temporary
		table."""

		self.hdulist[extension] = tempHDU
		self.flush()

	def __addDataSingle(self, frame):
		if self.mode == 'TBW':
			stand = frame.parseID()
		else:
			stand, pol = frame.parseID()
		extension = self.__findExtension(stand)

		if extension is None:
			print "Stand '%i' not found, creating new binary table extension" % stand
			self.standCount = self.standCount + 1

			# Data
			if self.mode == 'TBW':
				c1 = pyfits.Column(name='data', format='%iI' % frame.data.xy.shape[1], array=frame.data.xy.astype(numpy.int16))
			else:
				c1 = pyfits.Column(name='data', format='512C', array=frame.data.iq.astype(numpy.csingle))
			# Polarization
			c2 = pyfits.Column(name='pol', format='1I')
			# Time
			c3 = pyfits.Column(name='time', format='1K')

			# Define the collection of columns
			colDefs = pyfits.ColDefs([c1, c2, c3])

			# Get the time of the first sample and convert it to a datetime
			self.firstSamples[stand] = long(frame.data.timeTag / dp_common.fS)
			firstSample = datetime.utcfromtimestamp(self.firstSamples[stand])

			tsfits = pyfits.new_table(colDefs)
			tsfits.header.update('EXTNAME', 'TIME SERIES', after='tfields')
			tsfits.header.update('EXTVER', self.standCount, after='EXTNAME')
			tsfits.header.update('STAND', stand, after='EXTVER')
			tsfits.header.update('DATE-OBS', firstSample.isoformat())

			self.hdulist.append(tsfits)
			self.hdulist[0].header.update('NSTAND', self.standCount)
			self.flush()
			
			self.hdulist[-1].data.field('pol')[0] = 0
			self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]
			if self.mode == 'TBW':
				self.hdulist[-1].data.field('pol')[1] = 1
				self.hdulist[-1].data.field('time')[1] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]
		else:
			# Make sure that we have a first sample time to reference to if 
			# we are adding on to the end of a file
			if stand not in self.firstSamples.keys():
				firstSample = datetime.strptime(self.hdulist[extension].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
				self.firstSamples[stand] = long(time.mktime(firstSample.timetuple()))

			nrows = self.hdulist[extension].data.shape[0]
			if self.mode == 'TBW':
				tempHDU = self.__makeAppendTable(extension, AddRows=2)
				tempHDU.data.field('data')[nrows:] = frame.data.xy.astype(numpy.int16)
				tempHDU.data.field('pol')[nrows:] = numpy.array([0,1])
				tempHDU.data.field('time')[nrows:] = numpy.array([frame.data.timeTag, frame.data.timeTag]) / dp_common.fS - self.firstSamples[stand]
			else:
				tempHDU = self.__makeAppendTable(extension, AddRows=1)
				tempHDU.data.field('data')[nrows:] = frame.data.iq.astype(numpy.csingle)
				tempHDU.data.field('pol')[nrows:] = numpy.array([pol])
				tempHDU.data.field('time')[nrows:] = numpy.array([frame.data.timeTag]) / dp_common.fS - self.firstSamples[stand]

			self.__applyAppendTable(extension, tempHDU)

		self.flush()

	def __addDataQueue(self, frame):
		if self.mode == 'TBW':
			stand = frame.parseID()
		else:
			stand, pol = frame.parseID()
		if stand not in self.queue.keys():
			self.queue[stand] = []
		self.queue[stand].append(frame)

		if len(self.queue[stand]) >= self.queueLimit:
			start = 0
			extension = self.__findExtension(stand)

			if extension is None:
				start = 1
				frame = self.queue[stand][0]
				
				print "Stand '%i' not found, creating new binary table extension" % stand
				self.standCount = self.standCount + 1

				# Data
				if self.mode == 'TBW':
					c1 = pyfits.Column(name='data', format='%iI' % frame.data.xy.shape[1], array=frame.data.xy.astype(numpy.int16))
				else:
					c1 = pyfits.Column(name='data', format='512C', array=frame.data.iq.astype(numpy.csingle))
				# Polarization
				c2 = pyfits.Column(name='pol', format='1I')
				# Time
				c3 = pyfits.Column(name='time', format='1K')

				# Define the collection of columns
				colDefs = pyfits.ColDefs([c1, c2, c3])

				# Get the time of the first sample and convert it to a datetime
				self.firstSamples[stand] = long(frame.data.timeTag / dp_common.fS)
				firstSample = datetime.utcfromtimestamp(self.firstSamples[stand])

				tsfits = pyfits.new_table(colDefs)
				tsfits.header.update('EXTNAME', 'TIME SERIES', after='tfields')
				tsfits.header.update('EXTVER', self.standCount, after='EXTNAME')
				tsfits.header.update('STAND', stand, after='EXTVER')
				tsfits.header.update('DATE-OBS', firstSample.isoformat('T'))

				self.hdulist.append(tsfits)
				self.hdulist[0].header.update('NSTAND', self.standCount)
				self.flush()
				
				self.hdulist[-1].data.field('pol')[0] = 0
				self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]
				if self.mode == 'TBW':
					self.hdulist[-1].data.field('time')[1] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]
					self.hdulist[-1].data.field('pol')[1] = 1

				self.flush()
				extension = self.__findExtension(stand)

			# Make sure that we have a first sample time to reference to if 
			# we are adding on to the end of a file
			if stand not in self.firstSamples.keys():
				firstSample = datetime.strptime(self.hdulist[extension].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
				self.firstSamples[stand] = long(time.mktime(firstSample.timetuple()))

			nrows = self.hdulist[extension].data.shape[0]
			if self.mode == 'TBW':
				tempHDU = self.__makeAppendTable(extension, AddRows=2*(self.queueLimit-start))
				for count,frame in zip(range(len(self.queue[stand][start:])), self.queue[stand][start:]):
					tempHDU.data.field('data')[nrows+2*count] = numpy.squeeze(frame.data.xy.astype(numpy.int16)[0,:])
					tempHDU.data.field('pol')[nrows+2*count] = 0
					tempHDU.data.field('time')[nrows+2*count] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]
					tempHDU.data.field('data')[nrows+2*count+1] = numpy.squeeze(frame.data.xy.astype(numpy.int16)[1,:])
					tempHDU.data.field('pol')[nrows+2*count+1] = 1
					tempHDU.data.field('time')[nrows+2*count+1] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]
			else:
				tempHDU = self.__makeAppendTable(extension, AddRows=(self.queueLimit-start))
				for count,frame in zip(range(len(self.queue[stand][start:])), self.queue[stand][start:]):
					stand,pol = frame.parseID()
					tempHDU.data.field('data')[nrows+count] = frame.data.iq.astype(numpy.csingle)
					tempHDU.data.field('pol')[nrows+count] = pol
					tempHDU.data.field('time')[nrows+count] = frame.data.timeTag / dp_common.fS - self.firstSamples[stand]

			self.__applyAppendTable(extension, tempHDU)
			del(self.queue[stand])

	def addStandData(self, frame):
		if self.mode == 'TBW':
			try:
				self.hdulist[0].header['TBWBITS']
			except:
				self.hdulist[0].header.update('TBWBITS', frame.getDataBits())
				self.flush()
		else:
			try:
				self.hdulist[0].header['FILTER']
			except:
				self.hdulist[0].header.update('FILTER', frame.getFilterCode())
				if frame.data.gain is not None:
					self.hdulist[0].header.update('GAIN', frame.data.gain)
				if frame.data.centralFreq is not None:
					self.hdulist[0].header.update('FREQ', frame.data.centralFreq)
				self.flush()

		if self.UseQueue:
			self.__addDataQueue(frame)
		else:
			self.__addDataSingle(frame)

	def getStandData(stand):
		extension = self.__findExtension(stand)
		
		if extension is None:
			output = {}
			for key in ['data', 'pol', 'time']:
				output[key] = None
		else:
			output = {}
			for key in ['data', 'pol', 'time']:
				output[key] = self.hdulist[extension].field(key)
		return output


class TBW(TSFITS):
	def __init__(self, filename, Overwrite=False, UseQueue=True):
		super(TBW, self).__init__(filename, 'TBW', Overwrite=Overwrite, UseQueue=UseQueue)
		

class TBN(TSFITS):
	def __init__(self, filename, Overwrite=False, UseQueue=True):
		super(TBN, self).__init__(filename, 'TBN', Overwrite=Overwrite, UseQueue=UseQueue)
