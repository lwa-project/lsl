# -*- coding: utf-8 -*-

"""
Module that contains common values found in the MCS Joint Release 5 header file
src/exec/me.h and other functions useful for working with the MCS metadata.  
The header file values are:
  * ME_SSMIF_FORMAT_VERSION - SSMIF format version code
  * ME_MAX_NSTD - Maximum number of stands that can be described
  * ME_MAX_NFEE - Maximum number of FEEs that can be described
  * ME_MAX_FEEID_LENGTH - Number of characters in FEE ID name
  * ME_MAX_RACK - Maximum number of racks?
  * ME_MAX_PORT - Maximum number of ports?
  * ME_MAX_NRPD - Maxmimum number of RPD cables
  * ME_MAX_RPDID_LENGTH - Number of characters in the RPD ID name
  * ME_MAX_NSEP - Maximum number of SEP cable connections
  * ME_MAX_SEPID_LENGTH - Number of characters in the SEP ID name
  * ME_MAX_SEPCABL_LENGTH - Number of characters in the SEP cable ID name
  * ME_MAX_NARB - Maximum number of ARX boards
  * ME_MAX_NARBCH - Number of ARX channels per board
  * ME_MAX_ARBID_LENGTH - Number of characters in the ARX ID name
  * ME_MAX_NDP1 - Maximum number of DP1 boards
  * ME_MAX_NDP1CH - Number of channels per DP1 board
  * ME_MAX_DP1ID_LENGTH - Number of characters in the DP1 board ID name
  * ME_MAX_NDP2 - Maximum number of DP2 boards
  * ME_MAX_DP2ID_LENGTH - Number of characters in the DP2 board ID name
  * ME_MAX_NDR - Maximum number of data recorders
  * ME_MAX_DRID_LENGTH - Number of characters in the DR ID name
  * ME_MAX_NPWRPORT - Maximum number of power ports
  * ME_MAX_SSNAME_LENGTH - Number of characters in the power port ID names, for 
    codes used for PWR_NAME
  * LWA_MAX_NSTD - Maximum number of stands for the LWA
  * MIB_REC_TYPE_BRANCH - eType for MIB branch entries
  * MIB_REC_TYPE_VALUE - etype for MIB value entries
  * MIB_INDEX_FIELD_LENGTH - Number of characters in a MIB index field
  * MIB_LABEL_FIELD_LENGTH - Number of characters in a MIB label field
  * MIB_VAL_FIELD_LENGTH - Number of characters in a MIB value field
  
The other functions:
  * Convert MCS codes (status, summary, command ID, etc.) to human readable strings, 
  * Parse the binary packed metadata, 
  * Ready delays and gains for a custom beamforming SDF, and
  * Convert datetime instances to MJD and MPM values.
  * Apply a rotation axis-based pointing corretion to a azimuth/elevation pair
"""

import re
import math
import numpy
import struct
import anydbm
from ctypes import *
from aipy import coord
from datetime import datetime

from lsl.common.dp import gaintoDPG, DPGtogain, delaytoDPD, DPDtodelay


__version__ = '0.8'
__revision__ = '$Rev$'
__all__ = ['ME_SSMIF_FORMAT_VERSION', 'ME_MAX_NSTD', 'ME_MAX_NFEE', 'ME_MAX_FEEID_LENGTH', 'ME_MAX_RACK', 'ME_MAX_PORT', 
			'ME_MAX_NRPD', 'ME_MAX_RPDID_LENGTH', 'ME_MAX_NSEP', 'ME_MAX_SEPID_LENGTH', 'ME_MAX_SEPCABL_LENGTH', 
			'ME_MAX_NARB', 'ME_MAX_NARBCH', 'ME_MAX_ARBID_LENGTH', 'ME_MAX_NDP1', 'ME_MAX_NDP1CH', 'ME_MAX_DP1ID_LENGTH', 
			'ME_MAX_NDP2', 'ME_MAX_DP2ID_LENGTH', 'ME_MAX_NDR', 'ME_MAX_DRID_LENGTH', 'ME_MAX_NPWRPORT', 
			'ME_MAX_SSNAME_LENGTH', 'LWA_MAX_NSTD', 'MIB_REC_TYPE_BRANCH', 'MIB_REC_TYPE_VALUE', 'MIB_INDEX_FIELD_LENGTH', 
			'MIB_LABEL_FIELD_LENGTH', 'MIB_VAL_FIELD_LENGTH', 'IS_32BIT_PYTHON', 
			'delaytoMCSD', 'MCSDtodelay', 'gaintoMCSG', 'MCSGtogain',
			'mjdmpm2datetime', 'datetime2mjdmpm', 'status2string', 'summary2string', 'sid2string', 'cid2string', 
			'mode2string', 'parseCStruct', 'single2multi', 'applyPointingCorrection', 'MIB', 'MIBEntry', 
			'__version__', '__revision__', '__all__']


ME_SSMIF_FORMAT_VERSION = 5	# SSMIF format version code
ME_MAX_NSTD = 260			# Maximum number of stands that can be described
ME_MAX_NFEE = 260			# Maximum number of FEEs that can be described
ME_MAX_FEEID_LENGTH = 10		# Number of characters in FEE ID name
ME_MAX_RACK = 6			# Maximum number of racks?
ME_MAX_PORT = 50			# Maximum number of ports?
ME_MAX_NRPD = 520			# Maxmimum number of RPD cables
ME_MAX_RPDID_LENGTH = 25		# Number of characters in the RPD ID name
ME_MAX_NSEP = 520			# Maximum number of SEP cable connections
ME_MAX_SEPID_LENGTH = 25		# Number of characters in the SEP ID name
ME_MAX_SEPCABL_LENGTH = 25	# Number of characters in the SEP cable ID name
ME_MAX_NARB = 33			# Maximum number of ARX boards
ME_MAX_NARBCH = 16			# Number of ARX channels per board
ME_MAX_ARBID_LENGTH = 10		# Number of characters in the ARX ID name
ME_MAX_NDP1 = 26			# Maximum number of DP1 boards
ME_MAX_NDP1CH = 20			# Number of channels per DP1 board
ME_MAX_DP1ID_LENGTH = 10		# Number of characters in the DP1 board ID name
ME_MAX_NDP2 = 2			# Maximum number of DP2 boards
ME_MAX_DP2ID_LENGTH = 10		# Number of characters in the DP2 board ID name
ME_MAX_NDR = 5				# Maximum number of data recorders
ME_MAX_DRID_LENGTH = 10		# Number of characters in the DR ID name
ME_MAX_NPWRPORT = 50		# Maximum number of power ports
ME_MAX_SSNAME_LENGTH = 3		# Number of characters in the power port ID names, for codes used for PWR_NAME
LWA_MAX_NSTD = 260			# Maximum number of stands for the LWA
MIB_REC_TYPE_BRANCH = 0 		# eType for MIB branch entries
MIB_REC_TYPE_VALUE = 1 		# etype for MIB value entries
MIB_INDEX_FIELD_LENGTH = 12	# Number of characters in a MIB index field
MIB_LABEL_FIELD_LENGTH = 32	# Number of characters in a MIB label field
MIB_VAL_FIELD_LENGTH = 8192	# Number of characters in a MIB value field


IS_32BIT_PYTHON = True if struct.calcsize("P") == 4 else False


_cDecRE = re.compile(r'(?P<type>[a-z][a-z \t]+)[ \t]+(?P<name>[a-zA-Z_0-9]+)(\[(?P<d1>[\*\+A-Z_\d]+)\])?(\[(?P<d2>[\*\+A-Z_\d]+)\])?(\[(?P<d3>[\*\+A-Z_\d]+)\])?(\[(?P<d4>[\*\+A-Z_\d]+)\])?;')


def _twoByteSwap(value):
	return ((value & 0xFF) << 8) | ((value >> 8) & 0xFF)


def delaytoMCSD(delay):
	"""
	Given a delay in ns, convert it to a course and fine portion and into 
	the form expected by MCS in a custom beamforming SDF (little endian 
	16.12 unsigned integer).

	.. versionadded:: 0.6.3
	"""
	
	return _twoByteSwap( delaytoDPD(delay) )


def MCSDtodelay(delay):
	"""
	Given delay value from an OBS_BEAM_DELAY field in a custom beamforming 
	SDF, return the delay in ns.

	.. versionadded:: 0.6.3
	"""
	
	return DPDtodelay( _twoByteSwap(delay) )


def gaintoMCSG(gain):
	"""
	Given a gain (between 0 and 1), convert it to a gain in the form 
	expected by MCS in a custom beamforming SDF (little endian 16.1 
	signed integer).

	.. versionadded::0.6.3
	"""
	
	return _twoByteSwap( gaintoDPG(gain) )


def MCSGtogain(gain):
	"""
	Given a gain value from an OBS_BEAM_GAIN field in a custom beamforming
	SDF, return the decimal equivalent.

	.. versionadded:: 0.6.3
	"""
	
	return DPGtogain( _twoByteSwap(gain) )


def mjdmpm2datetime(mjd, mpm):
	"""
	Convert a MJD, MPM pair to a UTC-aware datetime instance.
	
	.. versionadded:: 0.5.2
	"""
	
	unix = mjd*86400.0 + mpm/1000.0 - 3506716800.0
	return datetime.utcfromtimestamp(unix)


def datetime2mjdmpm(dt):
	"""
	Convert a UTC datetime instance to a MJD, MPM pair (returned as a 
	two-element tuple).
	
	Based off: http://paste.lisp.org/display/73536
	
	.. versionadded:: 0.5.2
	"""
	
	year        = dt.year             
	month       = dt.month      
	day         = dt.day    

	hour        = dt.hour
	minute      = dt.minute
	second      = dt.second     
	millisecond = dt.microsecond / 1000

	# compute MJD         
	# adapted from http://paste.lisp.org/display/73536
	# can check result using http://www.csgnetwork.com/julianmodifdateconv.html
	a = (14 - month) // 12
	y = year + 4800 - a          
	m = month + (12 * a) - 3                    
	p = day + (((153 * m) + 2) // 5) + (365 * y)   
	q = (y // 4) - (y // 100) + (y // 400) - 32045
	mjdi = int(math.floor( (p+q) - 2400000.5))
	mjd = mjdi

	# compute MPM
	mpmi = int(math.floor( (hour*3600 + minute*60 + second)*1000 + millisecond ))
	mpm = mpmi
	return (mjd, mpm)


def status2string(code):
	"""
	Convert a numerical MCS status code to a string.
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


def summary2string(code):
	"""
	Convert a numerical MCS overall status code to an explination.
	"""
	
	if code < 0 or code > 6:
		raise ValueError("Invalid code %i" % code)
	
	if code == 0:
		return "Not normally used"
	elif code == 1:
		return "Normal"
	elif code == 2:
		return "Warning - issue(s) found, but still fully operational"
	elif code == 3:
		return "Error - problems found which limit or prevent proper function"
	elif code == 4:
		return "Booting - initializing; not yet fully operational"
	elif code == 5:
		return "Shutdown - shutting down; not ready for operation"
	else:
		return "Status is unknown"


def sid2string(sid):
	"""
	Convert a MCS subsystem ID code into a string.
	"""
	
	if sid < 1 or sid > 18:
		raise ValueError("Invalid sid code %i" % sid)
	
	if sid < 9:
		return "Null subsystem #%i" % sid
	elif sid == 10:
		return "MCS"
	elif sid == 11:
		return "SHL"
	elif sid == 12:
		return "ASP"
	elif sid == 13:
		return "DP"
	elif sid == 14:
		return "DR #1"
	elif sid == 15:
		return "DR #2"
	elif sid == 16:
		return "DR #3"
	elif sid == 17:
		return "DR #4"
	else:
		return "DR #5"


def cid2string(cid):
	"""
	Convert a MCS command code into a string.
	"""
	
	if cid < 0 or cid > 39:
		raise ValueError("Invalid cid code %i" % cid)
	
	if cid == 0:
		return "MCSSHT"
	elif cid == 1:
		return "PNG"
	elif cid == 2:
		return "RPT"
	elif cid == 3:
		return "SHT"
	elif cid == 4:
		return "INI"
	elif cid == 5:
		return "TMP"
	elif cid == 6:
		return "DIF"
	elif cid == 7:
		return "PWR"
	elif cid == 8:
		return "FIL"
	elif cid == 9:
		return "AT1"
	elif cid == 10:
		return "AT2"
	elif cid == 11:
		return "ATS"
	elif cid == 12:
		return "FPW"
	elif cid == 13:
		return "RXP"
	elif cid == 14:
		return "FEP"
	elif cid == 15:
		return "TBW"
	elif cid == 16:
		return "TBN"
	elif cid == 17:
		return "DRX"
	elif cid == 18:
		return "BAM"
	elif cid == 19:
		return "FST"
	elif cid == 20:
		return "CLK"
	elif cid == 21:
		return "REC"
	elif cid == 22:
		return "DEL"
	elif cid == 23:
		return "STP"
	elif cid == 24:
		return "GET"
	elif cid == 25:
		return "CPY"
	elif cid == 26:
		return "DMP"
	elif cid == 27:
		return "FMT"
	elif cid == 28:
		return "DWN"
	elif cid == 29:
		return "UP_"
	elif cid == 30:
		return "SEL"
	elif cid == 31:
		return "SYN"
	elif cid == 32:
		return "TST"
	elif cid == 33:
		return "BUF"
	elif cid == 34:
		return "NUL"
	elif cid == 35:
		return "ESN"
	elif cid == 36:
		return "ESF"
	elif cid == 37:
		return "OBS"
	elif cid == 38:
		return "OBE"
	else:
		return "SPC"


def mode2string(mode):
	"""
	Convert a MCS numeric observing mode into a string.
	"""
	
	if mode < 1 or mode > 6:
		raise ValueError("Invalid observing mode %i" % mode)
	
	if mode == 1:
		return "TRK_RADEC"
	elif mode == 2:
		return "TRK_SOL"
	elif mode == 3:
		return "TRK_JOV"
	elif mode == 4:
		return "STEPPED"
	elif mode == 5:
		return "TBW"
	else:
		return "TBN"


def parseCStruct(cStruct, charMode='str', endianness='native'):
	"""
	Function to take a C structure declaration and build a ctypes.Structure out 
	of it with the appropriate alignment, character interpretation*, and endianness
	(little, big, network, or native).
	
	*:  ctypes converts character arrays to Python strings until the first null is
	incountered.  This behavior causes problems for multi-dimension arrays of null
	filled strings.  By setting charMode to 'int', all char types are retuned as 
	bytes which can be converted to strings via chr().
	"""
	
	# Figure out how to deal with character arrays
	if charMode not in ('str', 'int'):
		raise RuntimeError("Unknown character mode: '%s'" % charMode)
	if charMode == 'str':
		baseCharType = c_char
	else:
		baseCharType = c_byte
	
	# Hold the basic fields and dimensions
	fields = []
	dims2 = {}
	
	# Split into lines and go!
	cStruct = cStruct.split('\n')
	for line in cStruct:
		## Skip structure declaration and lines too short to hold a declaration
		if '{' in line or '}' in line:
			continue
		if len(line) < 5:
			continue
		
		## RegEx the line to find the type, name, and dimensions (if needed) for
		## the next structure variable
		mtch = _cDecRE.search(line)
		if mtch is None:
			raise RuntimeError("Unparseable line: '%s'" % line)
		
		dec = mtch.group('type')
		dec = dec.rstrip()
		name = mtch.group('name')
		
		try:
			d1 = mtch.group('d1')
			if d1 is not None:
				d1 = eval(d1)
			d2 = mtch.group('d2')
			if d2 is not None:
				d2 = eval(d2)
			d3 = mtch.group('d3')
			if d3 is not None:
				d3 = eval(d3)
			d4 = mtch.group('d4')
			if d4 is not None:
				d4 = eval(d4)
		except NameError:
			raise RuntimeError("Unknown value in array index: '%s'" % line)
		
		## Basic data types
		if dec in ('signed int', 'int'):
			typ = c_int
		elif dec == 'unsigned int':
			typ = c_uint
		elif dec in ('signed short int', 'signed short', 'short int', 'short'):
			typ = c_short
		elif dec in ('unsigned short int', 'unsigned short'):
			typ = c_ushort
		elif dec in ('signed long int', 'signed long', 'long int', 'long'):
			if IS_32BIT_PYTHON:
				typ = c_longlong
			else:
				typ = c_long
		elif dec in ('unsigned long int', 'unsigned long'):
			if IS_32BIT_PYTHON:
				typ = c_uint64
			else:
				typ = c_ulong
		elif dec in ('signed long long', 'long long'):
			typ = c_longlong
		elif dec == 'unsigned long long':
			typ = c_uint64
		elif dec == 'float':
			typ = c_float
		elif dec == 'double':
			typ = c_double
		elif dec == 'char':
			typ = baseCharType
		elif dec == 'signed char':
			typ = c_byte
		elif dec == 'unsigned char':
			typ = c_ubyte
		else:
			raise RuntimeError("Unparseable line: '%s' -> type: %s, name: %s, dims: %s, %s, %s %s" % (line, dec, name, d1, d2, d3, d4))
		
		## Array identification and construction
		dims2[name] = []
		if d1 is not None:
			count = d1
			dims2[name].append(d1)
			
			if d2 is not None:
				count *= d2
				dims2[name].append(d2)
			if d3 is not None:
				count *= d3
				dims2[name].append(d3)
			if d4 is not None:
				count *= d4
				dims2[name].append(d4)
				
			typ *= count
		
		## Append
		fields.append( (name, typ) )
	
	# ctypes creation - endianess
	endianness = endianness.lower()
	if endianness not in ('little', 'big', 'network', 'native'):
		raise RuntimeError("Unknown endianness: '%s'" % endianness)
	
	if endianness == 'little':
		endianness = LittleEndianStructure
	elif endianness == 'big':
		endianness = BigEndianStructure
	elif endianness == 'network':
		endianness = BigEndianStructure
	else:
		endiannes = Structure
	
	# ctypes creation - actual
	class MyStruct(endianness):
		"""
		ctypes.Structure of the correct endianness for the provided
		C structure.  
		
		In addition to the standard attributes defined for a ctypes.Structure 
		instance there are a few additional attributes related to the parsed C
		structure.  They are:
		  * origC - String containing the original C structure
		  * dims  - Dictionary of the dimensionallity of the data, if needed, 
		            keyed by variable name
		"""
		
		origC = '\n'.join(cStruct)
		
		_fields_ = fields
		_pack_ = 8	# Pack it like we are 64-bit
		dims = dims2
		
		def __str__(self):
			"""
			Print out the structure in a nice easy-to-read formation that
			captures the various structure elements, their data types, and 
			their values.
			"""
			
			out = ''
			for f,d in self._fields_:
				out += '%s (%s): %s\n' % (f, d, eval("self.%s" % f))
			return out
			
		def returnDict(self):
			"""
			Return the structure as a simple Python dictionary keyed off the
			structure elements.
			"""
			
			output = {}
			for f,d in self._fields_:
				output[f] = eval("self.%s" % f)
			return output
	
	# Create and return
	return MyStruct()


def single2multi(inputList, dim1, dim2=None, dim3=None, dim4=None):
	if dim4 is not None:
		return _single2four(inputList, dim1, dim2, dim3, dim4)
		
	elif dim3 is not None:
		return _single2three(inputList, dim1, dim2, dim3)

	elif dim2 is not None:
		return _single2two(inputList, dim1, dim2)
		
	else:
		return list(inputList)


def _single2two(inputList, dim1, dim2):
	"""
	Convert a flatten list into a two-dimensional list.  This is useful
	for converting flat lists of ctypes into their two-dimensional forms.
	"""
	
	if dim1*dim2 < len(inputList):
		raise ValueError("Incompatiable dimensions: input=%i, output=%i by %i" % (len(inputList), dim1, dim2))
	
	outputList = []
	for i in xrange(dim1):
		outputList.append( [None for k in xrange(dim2)] )
		for j in xrange(dim2):
			try:
				outputList[i][j] = inputList[dim2*i+j]
			except IndexError:
				pass
	
	return outputList


def _single2three(inputList, dim1, dim2, dim3):
	"""
	Convert a flatten list into a three-dimensional list.  This is useful
	for converting flat lists of ctypes into their three-dimensional forms.
	"""
	
	if dim1*dim2*dim3 < len(inputList):
		raise ValueError("Incompatiable dimensions: input=%i, output=%i by %i by %i" % (len(inputList), dim1, dim2, dim3))
	
	outputList = []
	for i in xrange(dim1):
		outputList.append( [[None for l in xrange(dim3)] for k in xrange(dim2)] )
		for j in xrange(dim2):
			for k in xrange(dim3):
				try:
					outputList[i][j][k] = inputList[dim2*dim3*i+dim3*j+k]
				except IndexError:
					pass
	
	return outputList


def _single2four(inputList, dim1, dim2, dim3, dim4):
	"""
	Convert a flatten list into a four-dimensional list.  This is useful
	for converting flat lists of ctypes into their four-dimensional forms.
	"""
	
	if dim1*dim2*dim3*dim4 < len(inputList):
		raise ValueError("Incompatiable dimensions: input=%i, output=%i by %i by %i by %i" % (len(inputList), dim1, dim2, dim3, dim4))
	
	outputList = []
	for i in xrange(dim1):
		outputList.append( [[[None for m in xrange(dim4)] for l in xrange(dim3)] for k in xrange(dim2)] )
		for j in xrange(dim2):
			for k in xrange(dim3):
				for l in xrange(dim4):
					try:
						outputList[i][j][k][l] = inputList[dim2*dim3*dim4*i+dim3*dim4*j+dim4*k+l]
					except IndexError:
						pass
	
	return outputList


def _getRotationMatrix(theta, phi, psi, degrees=True):
	"""
	Generate the rotation matrix for a rotation about a specified axis.
	
	http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
	"""
	
	if degrees:
		theta = theta * numpy.pi/180.0
		phi   = phi * numpy.pi/180.0
		psi   = psi * numpy.pi/180.0
		
	# Axis
	u = numpy.array([numpy.cos(phi)*numpy.sin(theta), numpy.sin(phi)*numpy.sin(theta), numpy.cos(theta)])
	ux, uy, uz = u
	
	# Rotation matrix
	rot  = numpy.eye(3)*numpy.cos(psi) 
	rot += numpy.sin(psi)*numpy.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]]) 
	rot += (1-numpy.cos(psi))*numpy.tensordot(u, u, axes=0)
	
	return rot


def applyPointingCorrection(az, el, theta, phi, psi, lat=34.070, degrees=True):
	"""
	Given a azimuth and elevation pair, and an axis to rotate about, 
	perform the rotation.
	
	http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
	"""
	
	# Get the rotation matrix
	rot = _getRotationMatrix(theta, phi, psi, degrees=degrees)
	
	# Convert the az,alt coordinates to the unit vector
	if degrees:
		xyz = coord.azalt2top((az*numpy.pi/180.0, el*numpy.pi/180.0))
	else:
		xyz = coord.azalt2top((az, el))
		
	# Rotate
	xyzP = numpy.dot(rot, xyz)
	
	azP, elP = coord.top2azalt(xyzP)
	if degrees:
		azP *= 180/numpy.pi
		elP *= 180/numpy.pi
		
	return azP, elP


class MIB(object):
	"""
	Class to represent an entire MCS MIB.
	
	.. versionadded:: 0.6.5
	"""
	
	def __init__(self):
		self.entries = {}
		self.mapper = {}
		self.invMapper = {}
		
	def __str__(self):
		"""
		Describe the MIB as a string.
		"""
		
		nEntry = len(self.entries.keys())
		times = [self.entries[k].updateTime for k in self.entries.keys()]
		return "MIB containing %i entries from %s to %s" % (nEntry, min(times), max(times))
		
	def __getitem__(self, key):
		"""
		Allow access to the MIB through either the index of name.
		"""
		
		try:
			return self.entries[key]
		except KeyError:
			return self.entries[self.invMapper[key]]
			
	def keys(self, name=False):
		"""
		Return a list of entry indicies (or names if the 'name' keyword is set 
		to True) for the MIB.
		"""
		
		if name:
			# Index -> Name
			output = []
			for k in self.entries.keys():
				try:
					output.append( self.mapper[k] )
				except KeyError:
					output.append(k)
			return output
		else:
			# Index
			return self.entries.keys()
			
	def parseInitFile(self, filename):
		"""
		Given a MCS MIB initialization file, i.e., ASP_MIB_init.dat, create a 
		dictionary that maps indicies to human-readable names that can be used
		to clarify a MIBEntry.  Return this dictionary.
		"""
		
		# Open the file
		fh = open(filename)
		
		# Go!
		self.mapper = {}
		self.invMapper = {}
		for line in fh:
			## Parse out the line
			line = line.replace('\n', '')
			eType, index, name, default, dbmType, icdType = line.split(None, 5)
			
			## Skip over branches
			if eType == 'B':
				continue
				
			## Remember the mapping
			self.mapper[index] = name
			self.invMapper[name] = index
			
		# Done
		fh.close()
		return True
		
	def fromFile(self, filename, initFilename=None):
		"""
		Given the name of a GDBM database file, initialize the MIB.  
		Optionally, use the name of the MCS MIB initialization file to
		help convert indicies to names.
		"""
		
		# Parse the init. file (if we have one)
		if initFilename is not None:
			self.parseInitFile(initFilename)
			
		# Make sure we have the .pag file
		if filename[-4:] == '.dir':
			filename = filename.replace('.dir', '.pag')
			
		# Open the database
		db = anydbm.open(filename, 'ru')
		
		# Go!
		entry = db.firstkey()
		while entry is not None:
			value = db[entry]
			
			try:
				record = MIBEntry()
				record.fromEntry(value)
				self.entries[record.index] = record
				
			except ValueError:
				pass
				
			entry = db.nextkey(entry)
		db.close()
		
		# Done
		return True


class MIBEntry(object):
	"""
	Class for accessing and storing MCS MIB information contained in a GDBM 
	database.
	
	.. versionadded:: 0.6.5
	"""
	
	def __init__(self):
		"""
		Create the MIBEntry instance and fill with dummy values.
		"""
		
		# Basic information straight out of the mcs.h/record structure
		self.eType = 0
		self.index = ''
		self.value = ''
		self.dbmType = ''
		self.icdType = ''
		self._tv = (0, 0)
		
		# Additional information determined from the basic information
		self.updateTime = datetime.utcfromtimestamp(0)
		
	def __str__(self):
		"""
		Represent the class as a string for display purposes.
		"""
		
		return "Index: %s; Value: %s; Updated at %s" % (self.index, self.value, self.updateTime)
		
	def _parseValue(self, value, dataType):
		"""
		Convert an encoded value to something Pythonic (if possible).
		
		Understood codes:
		  * NUL:   No data stored (e.g., branch head entry)
		  * a####: printable (i.e., ASCII minus escape codes), #### = number of characters
		           e.g., "a3" means 3 printable ASCII-encoded characters
		  * r####: raw data (not printable), #### = number of bytes
		           e.g., "r1024" means 1024 bytes of raw data
		  * i1u:   integer, 1 byte,  unsigned, big-endian (=uint8)
		  * i2u:   integer, 2 bytes, unsigned, big-endian (=uint16)
		  * i4u:   integer, 4 bytes, unsigned, big-endian (=uint32)
		  * f4:    float, 4 bytes, big-endian (=float32)
		"""
		
		if dataType == 'NUL':
			return str(value)
		elif dataType[0] == 'a':
			return str(value)
		elif dataType[0] == 'r':
			return str(value)
		elif dataType[:3] == 'i1u':
			try:
				return struct.unpack('>1B', value)[0]
			except struct.error:
				return 0
		elif dataType[:3] == 'i2u':
			try:
				return struct.unpack('>1H', value)[0]
			except struct.error:
				return 0
		elif dataType[:3] == 'i4u':
			try:
				return struct.unpack('>1I', value)[0]
			except struct.error:
				return 0
		elif dataType[:2] == 'f4':
			try:
				return struct.unpack('>1f', value)[0]
			except struct.error:
				return 0.0
		else:
			raise ValueError("Unknown data type '%s'" % dataType)
			
	def fromEntry(self, value):
		"""
		Given an MIB entry straight out of a GDBM database, populate the 
		MIBEntry instance.
		
		.. note::
			The MIBEntry class does not currently support branch entries.  
			Entries of this type will raise a ValueError.
			
		.. note::
			The MIBEntry class currently only support entries with indicies
			that start with a number.  All other indicies will raise a 
			ValueError.
		"""
		
		# Initialize a new structure to parse the binary string
		record = parseCStruct("""
			int eType;
			char index[%i];
			char val[%i];
			char type_dbm[6];
			char type_icd[6];
			long tv[2];
			""" % (MIB_INDEX_FIELD_LENGTH, MIB_VAL_FIELD_LENGTH), endianness='little')
			
		# Squeeze the binary string into the structure using ctypes magic
		temp = create_string_buffer(value)
		memmove(addressof(record), addressof(temp), sizeof(record))
		
		# Validate
		if record.eType == MIB_REC_TYPE_BRANCH:
			raise ValueError("Cannot interpret MIB branch entries")
		if record.index[0] not in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
			raise ValueError("Entry index '%s' does not appear to be numeric" % record.index)
			
		# Basic information
		self.eType = int(record.eType)
		self.index = str(record.index)
		self.value = self._parseValue(record.val, record.type_dbm)
		self.dbmType = str(record.type_dbm)
		self.icdType = str(record.type_icd)
		self._tv = (int(record.tv[0]), int(record.tv[1]))
		
		# Time
		self.updateTime = datetime.utcfromtimestamp(record.tv[0] + record.tv[1]/1e9)
		
		return True