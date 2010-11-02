# -*- coding: utf-8 -*-

"""Modules to take DRX/TBW/TBN data and write it to a SDFITS file.  This module
is still under active development."""

import os
import sys
import numpy
import pyfits

from ..common import dp as dp_common
from ..reader.warngins import warnDeprecated
from ..correlator import fx as correlate

__version__ = '0.3'
__revision__ = '$ Revision: 15 $'
__all__ = ['SDFITS', 'TBW', 'TBN', 'writeSDFITS', '__version__', '__revision__', '__all__']


class SDFITS(object):
	def __init__(self, filename, mode, LFFT=128, Overwrite=False, UseQueue=True):
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
		self.hdulist.info()

	def flush(self):
		self.hdulist.flush()

	def setSite(self, site):
		self.site = site.name

		for hdu in self.hdulist:
			hdu.header.update('TELESCOP', self.site)

		self.flush()

	def setSampleRate(self, sampleRate):
		if self.mode != 'TBW':
			self.sampleRate = sampleRate

	def setCentralFrequency(self, centralFreq):
		if self.mode != 'TBW':
			self.centralFreq = centralFreq

	def __findExtension(self, stand):
		extension = None

		i = 1
		for hdu in self.hdulist[1:]:
			if int(hdu.header['stand']) == int(stand):
				extension = i
				break
			i = i + 1

		return extension

	def __addDataSingle(self, frame):
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

	def addStandData(self, frame):
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
	def __init__(self, filename, LFFT=128, Overwrite=False, UseQueue=True):
		super(TBW, self).__init__(filename, 'TBW', LFFT=LFFT, Overwrite=Overwrite, UseQueue=UseQueue)

class TBN(SDFITS):
	def __init__(self, filename, LFFT=128, Overwrite=False, UseQueue=True):
		super(TBN, self).__init__(filename, 'TBN', LFFT=LFFT, Overwrite=Overwrite, UseQueue=UseQueue)

class DRX(SDFITS):
	def __init__(self, filename, Overwrite=False, UseQueue=True):
		super(DRX, self).__init__(filename, 'DRX', LFFT=4096, Overwrite=Overwrite, UseQueue=UseQueue)



	def __init__(self, filename, LFFT=128, Overwrite=False, UseQueue=True):
		self.filename = filename
		self.hdulist = None
		self.fftLength = LFFT
		self.UseQueue = UseQueue
		self.queue = {}
		self.queueLimit = 10000

		if os.path.exists(filename) and not Overwrite:
			self.hdulist = pyfits.open(self.filename, mode="update", memmap=1)
			self.standCount = self.hdulist[0].header['nstand']
		else:
			if os.path.exists(filename):
				os.path.delete(filename)

			self.standCount = 0
			primary = pyfits.PrimaryHDU()
			primary.header.update('TELESCOP', 'LWA-1')
			primary.header.update('NSTAND', self.standCount)

			hdulist = pyfits.HDUList([primary])
			hdulist.writeto(filename)

			self.hdulist = pyfits.open(self.filename, mode="update", memmap=1)

	def info(self):
		self.hdulist.info()

	def flush(self):
		self.hdulist.flush()

	def __findExtension(self, stand):
		extension = None

		i = 1
		for hdu in self.hdulist[1:]:
			if int(hdu.header['stand']) == int(stand):
				extension = i
				break
			i = i + 1

		return extension

	def __addDataSingle(self, frame, SampleRate=100000, CentralFreq=40.0e6):
		stand,pol = frame.parseID()
		extension = self.__findExtension(stand)

		freq, framePS = correlate.calcSpectra(frame.data.iq, LFFT=self.fftLength, SampleRate=SampleRate)

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

			sdfits = pyfits.new_table(colDefs)
			sdfits.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
			sdfits.header.update('EXTVER', self.standCount, after='EXTNAME')
			sdfits.header.update('STAND', stand, after='EXTVER')

			# Define static fields - core
			sdfits.header.update('OBJECT', 'zenith')
			sdfits.header.update('TELESCOP', 'LWA-1')
			sdfits.header.update('BANDWID', (freq.max()-freq.min()))
			sdfits.header.update('DATE-OBS', '0000-00-00T00:00:00')
			sdfits.header.update('EXPOSURE', 0.061)
			sdfits.header.update('TSYS', 1.0)
			# Define static fields - virtual columns
			sdfits.header.update('CTYPE1', 'FREQ-OBS')
			sdfits.header.update('CRVAL1', CentralFreq)
			sdfits.header.update('CRPIX1', 1)
			sdfits.header.update('CDELT1', (freq[1]-freq[0]))
			sdfits.header.update('CTYPE2', 'HA')
			sdfits.header.update('CRVAL2', 0.0, unit='degrees')
			sdfits.header.update('CTYPE3', 'DEC')
			sdfits.header.update('CRVAL3', 0.0, unit='degrees')
			sdfits.header.update('OBSMODE', 'TBN')
			sdfits.header.update('NCHAN', len(freq))

			self.hdulist.append(tsfits)
			self.hdulist[0].header.update('NSTAND', self.standCount)
			self.flush()
			
			self.hdulist[-1].data.field('pol')[0] = 0
			self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS
		else:
			nrows = self.hdulist[extension].data.shape[0]
			tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+1)
			for key in ['EXTNAME', 'EXTVER', 'STAND']:
				tempHDU.header.update(key, self.hdulist[extension].header[key])
			tempHDU.data.field('data')[nrows:] = framePS.astype(numpy.float_)
			tempHDU.data.field('pol')[nrows:] = numpy.array([pol])
			tempHDU.data.field('time')[nrows:] = numpy.array([frame.data.timeTag / dp_common.fS])

			self.hdulist[extension] = tempHDU

		self.flush()

	def __addDataQueue(self, frame, SampleRate=100000, CentralFreq=40.0e6):
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
				freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength, SampleRate=SampleRate)
				
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

				sdfits = pyfits.new_table(colDefs)
				sdfits.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
				sdfits.header.update('EXTVER', self.standCount, after='EXTNAME')
				sdfits.header.update('STAND', stand, after='EXTVER')

				# Define static fields - core
				sdfits.header.update('OBJECT', 'zenith')
				sdfits.header.update('TELESCOP', 'LWA-1')
				sdfits.header.update('BANDWID', (freq.max()-freq.min()))
				sdfits.header.update('DATE-OBS', '0000-00-00T00:00:00')
				sdfits.header.update('EXPOSURE', 0.061)
				sdfits.header.update('TSYS', 1.0)
				# Define static fields - virtual columns
				sdfits.header.update('CTYPE1', 'FREQ-OBS')
				sdfits.header.update('CRVAL1', CentralFreq)
				sdfits.header.update('CRPIX1', 1)
				sdfits.header.update('CDELT1', (freq[1]-freq[0]))
				sdfits.header.update('CTYPE2', 'HA')
				sdfits.header.update('CRVAL2', 0.0, unit='degrees')
				sdfits.header.update('CTYPE3', 'DEC')
				sdfits.header.update('CRVAL3', 0.0, unit='degrees')
				sdfits.header.update('OBSMODE', 'TBN')
				sdfits.header.update('NCHAN', len(freq))

				self.hdulist.append(tsfits)
				self.hdulist[0].header.update('NSTAND', self.standCount)
				self.flush()
				
				self.hdulist[-1].data.field('pol')[0] = 0
				self.hdulist[-1].data.field('time')[0] = frame.data.timeTag / dp_common.fS

				self.flush()
				extension = self.__findExtension(stand)

			nrows = self.hdulist[extension].data.shape[0]
			tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+(self.queueLimit-start))
			for key in ['EXTNAME', 'EXTVER', 'STAND', 'OBJECT', 'TELESCOP', 'BANDWID', 'DATE-OBS', 'EXPOSURE', 'TSYS', 'CTYPE1', 'CRVAL1', 'CRPIX1', 'CDELT1', 'CTYPE2', 'CRVAL2', 'CTYPE3', 'CRVAL3', 'OBSMODE', 'NCHAN']:
				tempHDU.header.update(key, self.hdulist[extension].header[key])
			for count,frame in zip(range(len(self.queue[stand][start:])), self.queue[stand][start:]):
				freq, framePS = correlate.calcSpectra(frame.data.xy, LFFT=self.fftLength)

				tempHDU.data.field('data')[nrows+count] = framePS.astype(numpy.float_)
				tempHDU.data.field('pol')[nrows+count] = 0
				tempHDU.data.field('time')[nrows+count] = frame.data.timeTag / dp_common.fS

			self.hdulist[extension] = tempHDU
			self.flush()
			del(self.queue[stand])

	def addStandData(self, frame, SampleRate=100000, CentralFreq=40.0e6):
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
			self.__addDataQueue(frame, SampleRate=SampleRate, CentralFreq=CentralFreq)
		else:
			self.__addDataSingle(frame, SampleRate=SampleRate, CentralFreq=CentralFreq)

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


def __writeDRX(data, time, Beam=1, Polarization=1, Tunning=1):
	# Data
	nTime = data.shape[0]
	c1 = pyfits.Column(name='data', format='4096C', array=numpy.zeros((nTime, 4096), dtype=numpy.complex64))

	# Frequency
	c2 = pyfits.Column(name='ctype1', format='A8')
	c3 = pyfits.Column(name='crval1', format='1D')
	c4 = pyfits.Column(name='cdelt1', format='1D')
	c5 = pyfits.Column(name='crpix1', format='1D')

	# RA/Dec
	c6 = pyfits.Column(name='ctype2', format='A8')
	c7 = pyfits.Column(name='crval2', format='1D', unit='degrees')
	c8 = pyfits.Column(name='ctype3', format='A8')
	c9 = pyfits.Column(name='crval3', format='1D', unit='degrees')

	# Beam/Tunning/Polarization
	c10 = pyfits.Column(name='ctype4', format='A8')
	c11 = pyfits.Column(name='crval4', format='1D')
	c12 = pyfits.Column(name='ctype5', format='A8')
	c13 = pyfits.Column(name='crval5', format='1D')
	c14 = pyfits.Column(name='ctype6', format='A8')
	c15 = pyfits.Column(name='crval6', format='1D')

	# Time
	c16 = pyfits.Column(name='time', format='1D')

	# Obs. Mode/Channels
	c17 = pyfits.Column(name='obsmode', format='A3')
	c18 = pyfits.Column(name='nchan', format='1D')

	colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13, c14, c15, c16, c17, c18])

	prim = pyfits.PrimaryHDU()
	sdfits11 = pyfits.new_table(colDefs)
	sdfits11.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
	sdfits12 = pyfits.new_table(colDefs)
	sdfits12.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
	sdfits21 = pyfits.new_table(colDefs)
	sdfits21.header.update('EXTNAME', 'SINGLE DISH', after='tfields')
	sdfits22 = pyfits.new_table(colDefs)
	sdfits22.header.update('EXTNAME', 'SINGLE DISH', after='tfields')

	for i in range(data.shape[0]):
		print "%i of %s" % (i, data.shape[0])
		sdfits11.data.field('data')[i,:] = numpy.reshape(data[i,:], (1,4096))
		sdfits11.data.field('time')[i] = time[i]
	sdfits11.data.field('ctype1')[:] = 'FREQ-OBS'
	sdfits11.data.field('crval1')[:] = 0
	sdfits11.data.field('cdelt1')[:] = 1
	sdfits11.data.field('crpix1')[:] = 1
	sdfits11.data.field('ctype2')[:] = 'RA'
	sdfits11.data.field('crval2')[:] = 180.0
	sdfits11.data.field('ctype3')[:] = 'DEC'
	sdfits11.data.field('crval3')[:] = -35.0
	sdfits11.data.field('ctype4')[:] = 'BEAM'
	sdfits11.data.field('crval4')[:] = Beam
	sdfits11.data.field('ctype5')[:] = 'TUNE'
	sdfits11.data.field('crval5')[:] = Tunning
	sdfits11.data.field('ctype6')[:] = 'POL'
	sdfits11.data.field('crval6')[:] = Polarization
	sdfits11.data.field('obsmode')[:] = 'DRX'
	sdfits11.data.field('nchan')[:] = data.shape[1]

	print sdfits11.header.ascardlist()

	masterList = pyfits.HDUList([prim, sdfits11, sdfits12, sdfits21, sdfits22])

	if os.path.exists('test-sdfits.fits'):
		os.remove('test-sdfits.fits')

	masterList.writeto('test-sdfits.fits')


# Deprecated functions
def readfits(*args, **kwargs):
	warnDeprecated('readfits', memo='should be called via readSDFITS')
	readSDFITS(*args, **kwargs)


def writefits(*args, **kwargs):
	warnDeprecated('writefits', memo="should be called via writeSDFITS with mode='DRX'")
	writeSDFITS('DRX', *args, **kwargs)

