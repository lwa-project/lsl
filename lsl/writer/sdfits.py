# -*- coding: utf-8 -*-

"""Modules to take DRX/TBW/TBN data and write it to a SDFITS file.  This module
is still under active development."""

import os
import sys
import numpy
import pyfits

from lsl.common import dp as dp_common
from lsl.common.warns import warnDeprecated
from lsl.correlator import fx as correlate

__version__ = '0.3'
__revision__ = '$ Revision: 15 $'
__all__ = ['SDFITS', 'TBW', 'TBN', 'DRX', '__version__', '__revision__', '__all__']


class SDFITS(object):
	"""Class that holds TSFITS data until it is ready to be writen to disk."""

	def __init__(self, filename, mode, LFFT=128, Overwrite=False, UseQueue=True):
		"""Initialize a SDFITS object using a filename, an observation mode (TBW,  
		TBN, or DRX), and a FFT length in channels.  Optionally, SDFITS can be 
		told to overwrite the file if it already exists using the 'Overwrite' keyword."""

		assert(mode in ['TBW', 'TBN', 'DRX'])

		self.filename = filename
		self.mode = mode
		self.hdulist = None
		self.fftLength = LFFT
		self.UseQueue = UseQueue
		self.queue = {}
		self.queueLimit = 10000
		self.site = 'LWA'
		self.sampleRate = dp_common.fS
		self.firstSamples = {}

		if os.path.exists(filename) and not Overwrite:
			self.hdulist = pyfits.open(self.filename, mode="update", memmap=0)
			self.standCount = self.hdulist[0].header['nstand']
		else:
			if os.path.exists(filename):
				os.path.delete(filename)

			self.standCount = 0
			primary = pyfits.PrimaryHDU()
			primary.header.update('TELESCOP', self.site)
			primary.header.update('NSTAND', self.standCount)

			hdulist = pyfits.HDUList([primary])
			hdulist.writeto(filename)

			self.hdulist = pyfits.open(self.filename, mode="update", memmap=0)

	def info(self):
		"""Short-cut to the pyfits.info() function on an opened FITS file."""

		self.hdulist.info()

	def flush(self):
		"""Short-cut to the pyfits.flush() function on an opened FITS file."""

		self.hdulist.flush()

	def close(self):
		"""Empty the data queue (if it exists) and write all changes to disk 
		using flush."""

		if self.UseQueue:
			self.__emptyQueue()
		self.flush()

	def setSite(self, site):
		"""Set the TELESCOP keyword in the primary HDU using an lsl.common.stations
		object."""

		self.site = site.name

		for hdu in self.hdulist:
			hdu.header.update('TELESCOP', self.site)

		self.flush()

	def setSampleRate(self, sampleRate):
		"""Set the sample rate in Hz of the data for TBN and DRX observations."""

		if self.mode != 'TBW':
			self.sampleRate = sampleRate

	def setCentralFrequency(self, centralFreq):
		"""Set the central frequency in Hz of the data for TBN and DRX 
		observations."""

		if self.mode != 'TBW':
			self.centralFreq = centralFreq

	def __findExtension(self, stand):
		"""Private function to find out which extension stores the stand in 
		question.  None is returned is that stand is not in the SDFITS file."""

		extension = None

		i = 1
		for hdu in self.hdulist[1:]:
			if int(hdu.header['stand']) == int(stand):
				extension = i
				break
			i = i + 1

		return extension

	def __addDataSingle(self, frame):
		"""Private function to add a single data entry to a SDFITS file.  This 
		method is not particular fast since the existing file needs to be copied
		to append data to the end of a particular extension."""

		stand = frame.parseID()
		extension = self.__findExtension(stand)

		freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength, SampleRate=self.sampleRate)

		if extension is None:
			print "Stand '%i' not found, creating new binary table extension" % stand
			self.standCount = self.standCount + 1

			# Data - power spectrum
			c1 = pyfits.Column(name='data', format='%iD' % self.fftLength, array=framePS.astype(numpy.float_))
			# Polarization
			c2 = pyfits.Column(name='pol', format='1I')
			# Time
			c3 = pyfits.Column(name='time', format='1D')

			# Define the collection of columns
			colDefs = pyfits.ColDefs([c1, c2, c3])

			# Get the time of the first sample and convert it to a datetime
			self.firstSamples[stand] = long(frame.data.timeTag / dp_common.fS)
			firstSample = datetime.utcfromtimestamp(self.firstSamples[stand])

			sdfits = pyfits.new_table(colDefs)
			sdfits.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
			sdfits.header.update('EXTVER', self.standCount, after='EXTNAME')
			sdfits.header.update('STAND', stand, after='EXTVER')

			# Define static fields - core
			sdfits.header.update('OBJECT', 'zenith')
			sdfits.header.update('TELESCOP', self.site)
			sdfits.header.update('BANDWID', 78.0e6)
			sdfits.header.update('DATE-OBS', firstSample.isoformat('T'))
			sdfits.header.update('EXPOSURE', 0.0)
			sdfits.header.update('TSYS', 1.0)
			# Define static fields - virtual columns
			sdfits.header.update('CTYPE1', 'FREQ-OBS')
			sdfits.header.update('CRVAL1', freq[0])
			sdfits.header.update('CRPIX1', 1)
			sdfits.header.update('CDELT1', (freq[1]-freq[0]))
			sdfits.header.update('CTYPE2', 'HA')
			sdfits.header.update('CRVAL2', 0.0, unit='degrees')
			sdfits.header.update('CTYPE3', 'DEC')
			sdfits.header.update('CRVAL3', 0.0, unit='degrees')
			sdfits.header.update('OBSMODE', self.mode)
			sdfits.header.update('NCHAN', len(freq))

			self.hdulist.append(tsfits)
			self.hdulist[0].header.update('NSTAND', self.standCount)
			self.flush()
			
			self.hdulist[-1].data.field('pol')[0] = 0
			self.hdulist[-1].data.field('pol')[1] = 1
			self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS
			self.hdulist[-1].data.field('time')[1] = frame.data.timeTag / dp_common.fS
		else:
			# Make sure that we have a first sample time to reference to if 
			# we are adding on to the end of a file
			if stand not in self.firstSamples.keys():
				firstSample = long(datetime.strptime(self.hdulist[extension].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S"))
				self.firstSamples[stand] = time.mktime(firstSample.timetuple())

			nrows = self.hdulist[extension].data.shape[0]
			tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+2)
			for key in ['EXTNAME', 'EXTVER', 'STAND']:
				tempHDU.header.update(key, self.hdulist[extension].header[key])
			tempHDU.data.field('data')[nrows:] = framePS.astype(numpy.float_)
			tempHDU.data.field('pol')[nrows:] = numpy.array([0,1])
			tempHDU.data.field('time')[nrows:] = numpy.array([1,1])*frame.data.timeTag / dp_common.fS

			self.hdulist[extension] = tempHDU

		self.flush()

	def __addDataQueue(self, frame):
		"""Private function similar to __addDataSignle, but it saves the data to 
		memory (self.queue) to be written once the queue fills up."""

		stand = frame.parseID()
		if stand not in self.queue.keys():
			self.queue[stand] = []
		self.queue[stand].append(frame)

		if len(self.queue[stand]) >= self.queueLimit:
			start = 0
			extension = self.__findExtension(stand)

			if extension is None:
				start = 1
				frame = self.queue[stand][0]
				freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength, SampleRate=self.sampleRate)
				
				print "Stand '%i' not found, creating new binary table extension" % stand
				self.standCount = self.standCount + 1

				# Data - power spectrum
				c1 = pyfits.Column(name='data', format='%iD' % self.fftLength, array=framePS.astype(numpy.float_))
				# Polarization
				c2 = pyfits.Column(name='pol', format='1I')
				# Time
				c3 = pyfits.Column(name='time', format='1D')

				# Define the collection of columns
				colDefs = pyfits.ColDefs([c1, c2, c3])

				# Get the time of the first sample and convert it to a datetime
				self.firstSamples[stand] = long(frame.data.timeTag / dp_common.fS)
				firstSample = datetime.utcfromtimestamp(self.firstSamples[stand])

				sdfits = pyfits.new_table(colDefs)
				sdfits.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
				sdfits.header.update('EXTVER', self.standCount, after='EXTNAME')
				sdfits.header.update('STAND', stand, after='EXTVER')

				# Define static fields - core
				sdfits.header.update('OBJECT', 'zenith')
				sdfits.header.update('TELESCOP', 'LWA-1')
				sdfits.header.update('BANDWID', 78.0e6)
				sdfits.header.update('DATE-OBS', firstSample.isoformat('T'))
				sdfits.header.update('EXPOSURE', 0.0)
				sdfits.header.update('TSYS', 1.0)
				# Define static fields - virtual columns
				sdfits.header.update('CTYPE1', 'FREQ-OBS')
				sdfits.header.update('CRVAL1', freq[0])
				sdfits.header.update('CRPIX1', 1)
				sdfits.header.update('CDELT1', (freq[1]-freq[0]))
				sdfits.header.update('CTYPE2', 'HA')
				sdfits.header.update('CRVAL2', 0.0, unit='degrees')
				sdfits.header.update('CTYPE3', 'DEC')
				sdfits.header.update('CRVAL3', 0.0, unit='degrees')
				sdfits.header.update('OBSMODE', self.mode)
				sdfits.header.update('NCHAN', len(freq))

				self.hdulist.append(tsfits)
				self.hdulist[0].header.update('NSTAND', self.standCount)
				self.flush()
				
				self.hdulist[-1].data.field('pol')[0] = 0
				self.hdulist[-1].data.field('pol')[1] = 1
				self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS
				self.hdulist[-1].data.field('time')[1] = frame.data.timeTag / dp_common.fS

				self.flush()
				extension = self.__findExtension(stand)

			# Make sure that we have a first sample time to reference to if 
			# we are adding on to the end of a file
			if stand not in self.firstSamples.keys():
				firstSample = long(datetime.strptime(self.hdulist[extension].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S"))
				self.firstSamples[stand] = time.mktime(firstSample.timetuple())

			nrows = self.hdulist[extension].data.shape[0]
			tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+2*(self.queueLimit-start))
			for key in self.hdulist[extension].header.keys():
				tempHDU.header.update(key, self.hdulist[extension].header[key])
			for count,frame in zip(range(len(self.queue[stand][start:])), self.queue[stand][start:]):
				freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength)

				tempHDU.data.field('data')[nrows+2*count] = numpy.squeeze(framePS.astype(numpy.float_)[0,:])
				tempHDU.data.field('pol')[nrows+2*count] = 0
				tempHDU.data.field('time')[nrows+2*count] = frame.data.timeTag / dp_common.fS
				tempHDU.data.field('data')[nrows+2*count+1] = numpy.squeeze(framePS.astype(numpy.float_)[1,:])
				tempHDU.data.field('pol')[nrows+2*count+1] = 1
				tempHDU.data.field('time')[nrows+2*count+1] = frame.data.timeTag / dp_common.fS

			self.hdulist[extension] = tempHDU
			self.flush()
			del(self.queue[stand])

	def __emptyQueue(self):
		"""Private function to empty a empty the self.queue dictionary on demand
		if it exists."""

		for stand in self.queue.keys():
			if len(self.queue[stand]) == 0:
				continue
			start = 0
			extension = self.__findExtension(stand)

			if extension is None:
				start = 1
				frame = self.queue[stand][0]
				freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength, SampleRate=self.sampleRate)
				
				print "Stand '%i' not found, creating new binary table extension" % stand
				self.standCount = self.standCount + 1

				# Data - power spectrum
				c1 = pyfits.Column(name='data', format='%iD' % self.fftLength, array=framePS.astype(numpy.float_))
				# Polarization
				c2 = pyfits.Column(name='pol', format='1I')
				# Time
				c3 = pyfits.Column(name='time', format='1D')

				# Define the collection of columns
				colDefs = pyfits.ColDefs([c1, c2, c3])

				# Get the time of the first sample and convert it to a datetime
				self.firstSamples[stand] = long(frame.data.timeTag / dp_common.fS)
				firstSample = datetime.utcfromtimestamp(self.firstSamples[stand])

				sdfits = pyfits.new_table(colDefs)
				sdfits.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
				sdfits.header.update('EXTVER', self.standCount, after='EXTNAME')
				sdfits.header.update('STAND', stand, after='EXTVER')

				# Define static fields - core
				sdfits.header.update('OBJECT', 'zenith')
				sdfits.header.update('TELESCOP', 'LWA-1')
				sdfits.header.update('BANDWID', 78.0e6)
				sdfits.header.update('DATE-OBS', firstSample.isoformat('T'))
				sdfits.header.update('EXPOSURE', 0.0)
				sdfits.header.update('TSYS', 1.0)
				# Define static fields - virtual columns
				sdfits.header.update('CTYPE1', 'FREQ-OBS')
				sdfits.header.update('CRVAL1', freq[0])
				sdfits.header.update('CRPIX1', 1)
				sdfits.header.update('CDELT1', (freq[1]-freq[0]))
				sdfits.header.update('CTYPE2', 'HA')
				sdfits.header.update('CRVAL2', 0.0, unit='degrees')
				sdfits.header.update('CTYPE3', 'DEC')
				sdfits.header.update('CRVAL3', 0.0, unit='degrees')
				sdfits.header.update('OBSMODE', self.mode)
				sdfits.header.update('NCHAN', len(freq))

				self.hdulist.append(tsfits)
				self.hdulist[0].header.update('NSTAND', self.standCount)
				self.flush()
				
				self.hdulist[-1].data.field('pol')[0] = 0
				self.hdulist[-1].data.field('pol')[1] = 1
				self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS
				self.hdulist[-1].data.field('time')[1] = frame.data.timeTag / dp_common.fS

				self.flush()
				extension = self.__findExtension(stand)

			# Make sure that we have a first sample time to reference to if 
			# we are adding on to the end of a file
			if stand not in self.firstSamples.keys():
				firstSample = long(datetime.strptime(self.hdulist[extension].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S"))
				self.firstSamples[stand] = time.mktime(firstSample.timetuple())

			nrows = self.hdulist[extension].data.shape[0]
			tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+2*(self.queueLimit-start))
			for key in self.hdulist[extension].header.keys():
				tempHDU.header.update(key, self.hdulist[extension].header[key])
			for count,frame in zip(range(len(self.queue[stand][start:])), self.queue[stand][start:]):
				freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength)

				tempHDU.data.field('data')[nrows+2*count] = numpy.squeeze(framePS.astype(numpy.float_)[0,:])
				tempHDU.data.field('pol')[nrows+2*count] = 0
				tempHDU.data.field('time')[nrows+2*count] = frame.data.timeTag / dp_common.fS
				tempHDU.data.field('data')[nrows+2*count+1] = numpy.squeeze(framePS.astype(numpy.float_)[1,:])
				tempHDU.data.field('pol')[nrows+2*count+1] = 1
				tempHDU.data.field('time')[nrows+2*count+1] = frame.data.timeTag / dp_common.fS

			self.hdulist[extension] = tempHDU
			self.flush()
			del(self.queue[stand])

	def addStandData(self, frame):
		"""Add a frame object to the SDFITS file.  This function takes care of 
		figuring out which extension the data goes to and how it should be formated.
		This function also updates the primary HDU with information about the data, 
		i.e., TBW data bits."""

		if self.mode == 'TBW':
			try:
				self.hdulist[0].header['TBWBITS']
			except:
				self.hdulist[0].header.update('TBWBITS', frame.getDataBits())
				self.flush()

		if self.UseQueue:
			self.__addDataQueue(frame)
		else:
			self.__addDataSingle(frame)

	def getStandData(stand):
		"""Retrieve a dictionary of all data stored for a particular stand.  The 
		dictionary keys are:
		  * *data* - numpy array of data
		  * *pol* - numpy array of polarizations
		  * *time* - numpy array of times in samples at f_S since DATE-OBS
		"""

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


class TBW(SDFITS):
	"""Sub-class of SDFITS for dealing with TBW data in particular."""

	def __init__(self, filename, LFFT=128, Overwrite=False, UseQueue=True):
		super(TBW, self).__init__(filename, 'TBW', LFFT=LFFT, Overwrite=Overwrite, UseQueue=UseQueue)


class TBN(SDFITS):
	"""Sub-class of SDFITS for dealing with TBN data in particular."""

	def __init__(self, filename, LFFT=128, Overwrite=False, UseQueue=True):
		super(TBN, self).__init__(filename, 'TBN', LFFT=LFFT, Overwrite=Overwrite, UseQueue=UseQueue)


class DRX(SDFITS):
	"""Sub-class of SDFITS for dealing with DRX data in particular."""

	def __init__(self, filename, Overwrite=False, UseQueue=True):
		super(DRX, self).__init__(filename, 'DRX', LFFT=4096, Overwrite=Overwrite, UseQueue=UseQueue)
