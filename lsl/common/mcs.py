# -*- coding: utf-8 -*-

"""
Module that contains common values found in the MCS Joint Release 5 header file
src/exec/me.h.  The values are:
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
"""

import re
import struct
from ctypes import *


__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['ME_SSMIF_FORMAT_VERSION', 'ME_MAX_NSTD', 'ME_MAX_NFEE', 'ME_MAX_FEEID_LENGTH', 'ME_MAX_RACK', 'ME_MAX_PORT', 
			'ME_MAX_NRPD', 'ME_MAX_RPDID_LENGTH', 'ME_MAX_NSEP', 'ME_MAX_SEPID_LENGTH', 'ME_MAX_SEPCABL_LENGTH', 
			'ME_MAX_NARB', 'ME_MAX_NARBCH', 'ME_MAX_ARBID_LENGTH', 'ME_MAX_NDP1', 'ME_MAX_NDP1CH', 'ME_MAX_DP1ID_LENGTH', 
			'ME_MAX_NDP2', 'ME_MAX_DP2ID_LENGTH', 'ME_MAX_NDR', 'ME_MAX_DRID_LENGTH', 'ME_MAX_NPWRPORT', 
			'ME_MAX_SSNAME_LENGTH', 'LWA_MAX_NSTD', 'status2string', 'summary2string', 
			'sid2string', 'cid2string', 'mode2string', 'parseCStruct', 'single2multi', 
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


cDeclRE = re.compile(r'(?P<type>[a-z][a-z \t]+)[ \t]+(?P<name>[a-zA-Z_0-9]+)(\[(?P<d1>[\*\+A-Z_\d]+)\])?(\[(?P<d2>[\*\+A-Z_\d]+)\])?(\[(?P<d3>[\*\+A-Z_\d]+)\])?(\[(?P<d4>[\*\+A-Z_\d]+)\])?;')


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
		mtch = cDeclRE.search(line)
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
			typ = c_long
		elif dec in ('unsigned long int', 'unsigned long'):
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
