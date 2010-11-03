# -*- coding: utf-8 -*-

"""Module for writting correlator output to a FITS IDI file.  This module is 
still under active development."""

import os
import sys
import numpy
import pyfits
from datetime import datetime

from lsl.common import dp as dp_common
from lsl.common.stations import geo2ecef
from lsl.correlator import uvUtils
from lsl.reader.warnings import warnDeprecated

__version__ = '0.1'
__revision__ = '$ Revision: 1 $'
__all__ = ['IDI', 'StokesCodes', '__version__', '__revision__', '__all__']


StokesCodes = {'I': 1, 'Q': 2, 'U': 3, 'V': 4, 
			'RR': -1, 'LL': -2, 'RL': -3, 'LR': -4, 
			'XX': -5, 'YY': -6, 'XY': -7, 'YX': -8}


class IDI(object):
	def __init__(self, filename, Overwrite=False):
		# File-specific information
		self.filename = filename
		self.hdulist = None

		# Observation-specifc information
		self.nAnt = 0
		self.nPol = 0
		self.nChan = 0
		self.nStokes = 0

		if os.path.exists(filename) and not Overwrite:
			self.hdulist = pyfits.open(self.filename, mode="update", memmap=1)
		else:
			if os.path.exists(filename):
				os.path.delete(filename)

			now = datetime.now()
			primary = pyfits.PrimaryHDU()
			
			primary.header.update('NAXIS', 0)
			primary.header.update('EXTEND', True)
			primary.header.update('GROUPS', True)
			primary.header.update('GCOUNT', 0)
			primary.header.update('PCOUNT', 0)
			primary.header.update('OBJECT', 'BINARYTB')
			primary.header.update('TELESCOP', 'LWA-1')
			primary.header.update('INSTRUME', 'LWA-1')
			primary.header.update('OBSERVER', 'LWA-1')
			primary.header.update('DATE-OBS', '0000-00-00T00:00:00')
			primary.header.update('DATE-MAP', now.isoformat())

			hdulist = pyfits.HDUList([primary])
			hdulist.writeto(filename)

			self.hdulist = pyfits.open(self.filename, mode="update", memmap=1)
			#self.initGeometry()
			#self.initSource()
			#self.initFrequency()
			#self.initAntenna()
			#self.initData()

	def info(self):
		self.hdulist.info()

	def flush(self):
		self.hdulist.flush()

	def fillFromAA(self, aa):
		self.nAnt = len(aa.ants)
		self.nChan = len(aa.get_afreqs())
		self.nPol = 1

	def __makeAppendTable(self, extension, AddRows=1):
		"""Private function to make a temporary table for appending data."""

		nrows = self.hdulist[extension].data.shape[0]
		tempHDU = pyfits.new_table(self.hdulist[extension].columns, nrows=nrows+AddRows)
		for key in self.hdulist[extension].header.keys():
			tempHDU.header.update(key, self.hdulist[extension].header[key])
	
		return tempHDU

	def __applyAppendTable(self, extension, tempHDU):
		"""Private function to replace the given extension with the temporary
		table."""

		self.hdulist[extension] = tempHDU
		self.flush()

	def initGeometry(self):
		"""Define the Array_Geometry table (group 1, table 1)."""

		# Antenna name
		c1 = pyfits.Column(name='anname', format='A8')
		# Station coordinates in meters
		c2 = pyfits.Column(name='stabxyz', unit='METERS', format='3D')
		# First order derivative of station coordinates in m/s
		c3 = pyfits.Column(name='derxyz', unit='METERS/S', format='3E')
		# Orbital elements
		c4 = pyfits.Column(name='orbparm', format='1D')
		# Station number
		c5 = pyfits.Column(name='nosta', format='1J')
		# Mount type (0 == alt-azimuth)
		c6 = pyfits.Column(name='mntsta', format='1J')
		# Axis offset in meters
		c7 = pyfits.Column(name='staxof', unit='METERS', format='3E')

		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7])

		ag = pyfits.new_table(colDefs)
		ag.header.update('EXTNAME', 'ARRAY_GEOMETRY', after='tfields')
		ag.header.update('TABREV', 1, after='EXTNAME')
		ag.header.update('EXTVER', 1, after='TABREV')
		ag.header.update('ARRNAM', 'TEST')
		ag.header.update('FRAME', 'GEOCENTRIC')
		ag.header.update('ARRAYX', 0.0)
		ag.header.update('ARRAYY', 0.0)
		ag.header.update('ARRAYZ', 0.0)
		ag.header.update('NUMORB', 0)
		ag.header.update('FREQ', 0.0)
		ag.header.update('TIMSYS', 'UTC')
		ag.header.update('RDATE', '0000-00-00')
		ag.header.update('GSTIA0', 0.0)
		ag.header.update('DEGPDY', 360.0)
		ag.header.update('UT1UTC', 0.0)
		ag.header.update('IATUTC', 0.0)
		ag.header.update('POLARX', 0.0)
		ag.header.update('POLARY', 0.0)
		ag.header.update('OBSCODE', 0.0)
		ag.header.update('NO_STKD', 0.0)
		ag.header.update('STK_1', 0.0)
		ag.header.update('NO_BAND', 0)
		ag.header.update('NO_CHAN', 0)
		ag.header.update('REF_FREQ', 0.0)
		ag.header.update('CHAN_BW', 0.0)
		ag.header.update('REF_PIXL', 0)

		self.hdulist.append(ag)
		self.flush()

	def initSource(self):
		"""Define the Source table (group 1, table 2)."""

		# Source ID number
		c1 = pyfits.Column(name='source_id', format='1J')
		# Source name
		c2 = pyfits.Column(name='source', format='A16')
		# Source qualifier
		c3 = pyfits.Column(name='qual', format='1J')
		# Calibrator code
		c4 = pyfits.Column(name='calcode', format='A4')
		# Frequency group ID
		c5 = pyfits.Column(name='freqid', format='1J')
		# Stokes I flux density in Jy
		c6 = pyfits.Column(name='iflux', format='1E')
		# Stokes I flux density in Jy
		c7 = pyfits.Column(name='qflux', format='1E')
		# Stokes I flux density in Jy
		c8 = pyfits.Column(name='uflux', format='1E')
		# Stokes I flux density in Jy
		c9 = pyfits.Column(name='vflux', format='1E')
		# Spectral index
		c10 = pyfits.Column(name='alpha', format='1E')
		# Frequency offset in Hz
		c11 = pyfits.Column(name='freqoff', format='1E')
		# Right ascension at mean equinox in degrees
		c12 = pyfits.Column(name='raepo', format='1D')
		# Declination at mean equinox in degrees
		c13 = pyfits.Column(name='decepo', format='1D')
		# Mean equinox
		c14 = pyfits.Column(name='equinox', format='A8')
		# Appearrent right ascension in degrees
		c15 = pyfits.Column(name='rapp', format='1D')
		# Apparent declination in degrees
		c16 = pyfits.Column(name='decapp', format='1D')
		# Systemic velocity in m/s
		c17 = pyfits.Column(name='sysvel', format='1D')
		# Velocity type
		c18 = pyfits.Column(name='veltyp', format='A8')
		# Velocity definition
		c19 = pyfits.Column(name='veldef', format='A8')
		# Line rest frequency in Hz
		c20 = pyfits.Column(name='restfreq', format='1D')
		# Proper motion in RA in degrees/day
		c21 = pyfits.Column(name='pmra', format='1D')
		# Proper motion in Dec in degrees/day
		c22 = pyfits.Column(name='pmdec', format='1D')
		# Parallax of source in arc sec.
		c23 = pyfits.Column(name='parallax', format='1E')

		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, 
							c21, c22, c23])

		sr = pyfits.new_table(colDefs)
		sr.header.update('EXTNAME', 'SOURCE', after='tfields')
		sr.header.update('TABREV', 1, after='EXTNAME')
		sr.header.update('OBSCODE', 0.0)
		sr.header.update('NO_STKD', 0.0)
		sr.header.update('STK_1', 0.0)
		sr.header.update('NO_BAND', 0)
		sr.header.update('NO_CHAN', 0)
		sr.header.update('REF_FREQ', 0.0)
		sr.header.update('CHAN_BW', 0.0)
		sr.header.update('REF_PIXL', 0)

		self.hdulist.append(sr)
		self.flush()

	def initFrequency(self):
		"""Define the Frequency table (group 1, table 3)."""

		# Frequency setup number
		c1 = pyfits.Column(name='freqid', format='1J')
		# Frequency offsets in Hz
		c2 = pyfits.Column(name='bandfreq', format='1D')
		# Channel width in Hz
		c3 = pyfits.Column(name='ch_width', format='1R')
		# Total bandwidths of bands
		c4 = pyfits.Column(name='total_bandwidth', format='1R')
		# Sideband flag
		c5 = pyfits.Column(name='sideband', format='1J')
		# Baseband channel
		c6 = pyfits.Column(name='bb_chan', format='1J')

		# Define the collection of columns
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6])

		fq = pyfits.new_table(colDefs)
		fq.header.update('EXTNAME', 'FREQUENCY', after='tfields')
		fq.header.update('TABREV', 2, after='EXTNAME')
		fq.header.update('OBSCODE', 0.0)
		fq.header.update('NO_STKD', 0.0)
		fq.header.update('STK_1', 0.0)
		fq.header.update('NO_BAND', 0)
		fq.header.update('NO_CHAN', 0)
		fq.header.update('REF_FREQ', 0.0)
		fq.header.update('CHAN_BW', 0.0)
		fq.header.update('REF_PIXL', 0)

		self.hdulist.append(fq)
		self.flush()

	def initAntenna(self):
		"""Define the Antenna table (group 2, table 1)."""

		# Central time of period covered by record in days
		c1 = pyfits.Column(name='time', units='DAYS', format='1D')
		# Durration of period covered by record in days
		c2 = pyfits.Column(name='time_interval', unit='DAYS', format='1E')
		# Antenna name
		c3 = pyfits.Column(name='anname', format='A8')
		# Antenna number
		c4 = pyfits.Column(name='antenna_no', format='1J')
		# Array number
		c5 = pyfits.Column(name='array', format='1J')
		# Frequency setup number
		c6 = pyfits.Column(name='freqid', format='IJ')
		# Number of digitizer levels
		c7 = pyfits.Column(name='no_levels', format='1J')
		# Feed A polarization label
		c8 = pyfits.Column(name='poltya', format='A1')
		# Feed A orientation in degrees
		c9 = pyfits.Column(name='polaa', format='1E')
		# Feed A polarization parameters
		c10 = pyfits.Column(name='polcala', format='1E')
		# Feed B polarization label
		c11 = pyfits.Column(name='poltyb', format='A1')
		# Feed B orientation in degrees
		c12 = pyfits.Column(name='polab', format='1E')
		# Feed B polarization parameters
		c13 = pyfits.Column(name='polcalb', format='1E')

		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13])

		an = pyfits.new_table(colDefs)
		an.header.update('EXTNAME', 'ANTENNA', after='tfields')
		an.header.update('TABREV', 1, after='EXTNAME')
		an.header.update('NOPCAL', 1)
		an.header.update('OBSCODE', 0.0)
		an.header.update('NO_STKD', 0.0)
		an.header.update('STK_1', 0.0)
		an.header.update('NO_BAND', 0)
		an.header.update('NO_CHAN', 0)
		an.header.update('REF_FREQ', 0.0)
		an.header.update('CHAN_BW', 0.0)
		an.header.update('REF_PIXL', 0)
		an.header.update('POLTYPE', 'X-Y LIN')
		
		self.hdulist.append(an)
		self.flush()

	def initBandpass(self):
		"""Define the Bandpass table (group 2, table 3)."""

		# Central time of period covered by record in days
		c1 = pyfits.Column(name='time', units='DAYS', format='1D')
		# Durration of period covered by record in days
		c2 = pyfits.Column(name='time_interval', unit='DAYS', format='1E')
		# Source ID
		c3 = pyfit.Column(name='source_id', format='1J')
		# Antenna number
		c4 = pyfits.Column(name='antenna_no', format='1J')
		# Array number
		c5 = pyfits.Column(name='array', format='1J')
		# Frequency setup number
		c6 = pyfits.Column(name='freqid', format='IJ')
		# Bandwidth in Hz
		c7 = pyfits.Column(name='bandwidth', unit='HZ', format='1E')
		# Band frequency in Hz
		c8 = pyfits.Column(name='band_freq', unit='HZ', format='1D')
		# Referance antenna number
		c9 = pyfits.Column(name='ref_ant1', format='1J')
		# Real part of the bandpass
		c10 = pyfits.Column(name='breal_1', format='%dE' % nChan)
		# Imagniary part of the bandpass
		c11 = pyfits.Column(name='bimag_1', format='%dE' % nChan)

		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11])

		bp = pyfits.new_table(colDefs)
		bp.header.update('EXTNAME', 'BANDPASS', after='tfields')
		bp.header.update('TABREV', 1, after='EXTNAME')
		bp.header.update('OBSCODE', 0.0)
		bp.header.update('NO_STKD', 0.0)
		bp.header.update('STK_1', 0.0)
		bp.header.update('NO_BAND', 0)
		bp.header.update('NO_CHAN', 0)
		bp.header.update('REF_FREQ', 0.0)
		bp.header.update('CHAN_BW', 0.0)
		bp.header.update('REF_PIXL', 0)
		bp.header.update('NO_ANT', 0)
		bp.header.update('NO_POL', 0)
		bp.header.update('NO_BACH', nChan)
		bp.header.update('STRT_CHN', 1)

		self.hdulist.append(bp)
		self.flush()

	def initData(self):
		"""Define the UV_Data table (group 3, table 1)."""

		# Complex visibility data (real, imag, weight)
		c1 = pyfits.Column(name='time', format='3D')
		# Stokes parameter
		c2 = pyfits.Column(name='stokes', format='1I')
		# Frequency channel number (needs corresponding CRPIX/CRDELT)
		c3 = pyfits.Column(name='freq', format='1J')
		# Band
		c4 = pyfits.Column(name='band', format='1I')
		# RA of phase center in degrees
		c5 = pyfits.Column(name='ra', format='1D')
		# Dec. of phase center in degrees
		c6 = pyfits.Column(name='dec', format='1D')
		# U coordinate (light seconds)
		c7 = pyfits.Column(name='uu--sin', format='1D')
		# V coordinate (light seconds)
		c8 = pyfits.Column(name='vv--sin', format='1D')
		# W coordinate (light seconds)
		c9 = pyfits.Column(name='ww--sin', format='1D')
		# Julian date at 0h
		c10 = pyfits.Column(name='date', format='1D')
		# Time elapsed since 0h
		c11 = pyfits.Column(name='time', format='1D')
		# Baseline number (first*256+second)
		c12 = pyfits.Column(name='baseline', format='1J')
		# Array number
		c13 = pyfits.Column(name='array', format='1J')
		# Source ID number
		c14 = pyfits.Column(name='source_id', format='1J')
		# Frequency setup number
		c15 = pyfits.Column(name='freqid', format='1J')
		# Inegration time (seconds)
		c16 = pyfits.Column(name='intim', format='1D')
		# Weights
		c17 = pyfits.Column(name='weight', format='1E')
		
		colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, 
							c11, c12, c13, c14, c15, c16, c17])

		uv = pyfits.new_table(colDefs)
		uv.header.update('EXTNAME', 'UV_DATA', after='tfields')
		uv.header.update('TABREV', 2, after='EXTNAME')
		uv.header.update('OBSCODE', 0.0)
		uv.header.update('NO_STKD', 0.0)
		uv.header.update('STK_1', 0.0)
		uv.header.update('NO_BAND', 0)
		uv.header.update('NO_CHAN', 0)
		uv.header.update('REF_FREQ', 0.0)
		uv.header.update('CHAN_BW', 0.0)
		uv.header.update('REF_PIXL', 0)
		
		self.hdulist.append(uv)
		self.flush()

	def setStokes(self, polList):
		"""Given a list of numerical Stokes parameters, update the headers 
		of all tables for NO_STKD and STK_1."""

		for hdu in self.hdulist:
			hdu.header.update('NO_STKD', len(polList))
			hdu.header.update('STK_1', polList[0])

		self.flush()

	def setFrequency(self, freq):
		"""Given a numpy array of frequencies, set the relavant header 
		keywords in all of the tables (NO_CHAN, REF_FREQ, etc.) and
		define a frequency setup.  The frequency setup ID is returned."""

		nBand = self.hdulist[1].header['NO_BAND'] + 1
		nChan = len(freq)
		refVal = freq[0]
		refPix = 0
		BW = numpy.abs(freq[1]-freq[0])

		for hdu in self.hdulist:
			hdu.header.update('NO_BAND', nBand)
			hdu.header.update('NO_CHAN', nChan)
			hdu.header.update('REF_FREQ', refVal)
			hdu.header.update('REF_PIXL', refPix)
			hdu.header.update('CHAN_BW', BW)

		self.flush()

		# Fixed frequency ID of 1
		freqID = 1
		bandFreq = 0.0
		BW = BW
		totalBW = BW*nChan
		sideband = 1

		# Append the data to the end of the table
		nrows = self.hdulist[3].data.shape[0]
		tempHDU = self.__makeAppendTable(3, AddRows=1)
		tempHDU.data.field('freqid')[nrows+1] = freqID
		tempHDU.data.field('bandfreq')[nrows+1] = bandFreq
		tempHDU.data.field('ch_width')[nrows+1] = BW
		tempHDU.data.field('total_bandwidth')[nrows+1] = totalBW
		tempHDU.data.field('sideband')[nrows+1] = 1

		self.__applyAppendTable(3, tempHDU)

		return freqID

	def setGeometry(self, site, stands):
		"""Given a numpy array of stands, update the ARRAY_GEOMETRY table."""

		arrayX, arrayY, arrayZ = site.getGeocentricLocation()
		self.hdulist[1].header.update('arrayx', arrayX)
		self.hdulist[1].header.update('arrayy', arrayY)
		self.hdulist[1].header.update('arrayz', arrayZ)

		for stand in stands:
			pass

	def baseline2baseline(stand1, stand2):
		return 256*stand1 + stand2
		
		
