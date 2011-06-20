# -*- coding: utf-8 -*-

"""
Module for reading in an interpreting binary-packed Station Dynamic MIB (SDM) 
files (as defined in MCS0031, v5).
"""

import copy
import struct

from lsl.common.mcs import *

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['SubSystemStatus', 'SubSubSystemStatus', 'StationsSettings', 'SDM', 'parseSDM', '__version__', '__revision__', '__all__']


def _guidedBinaryRead(fh, fmt):
	"""
	Function to wrap reading in packed binary data directrly from an open file
	handle.  This function calls struct.unpack() and struct.calcsize() to figure 
	out what to read and how.
	
	Return either a single item if a single item is requested or a list of items.
	
	Used by SDM.binaryFill()
	"""
	
	
	data = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))
	if len(data) == 1:
		return data[0]
	else:
		return list(data)


class SubSystemStatus(object):
	"""
	Python object that holds the status for a particular subsystem in a SDM 
	file.
	"""
	
	def __init__(self, name, summary=6, info='UNK', time=0):
		self.name = name
		self.summary = int(summary)
		self.info = str(info)
		self.time = float(time)
		
	def __str__(self):
		return "%s at %i: %s [%i]" % (self.name, self.time, self.info, self.summary)
		
	def _binaryRead(self, fh):
		"""
		Given an open file handle, interpret it in the context of a 
		subsystem_status_struct C structure and update the Python instance accordingly.
		"""
		
		self.summary = _guidedBinaryRead(fh, "<i")
		self.info = _guidedBinaryRead(fh, "<256s")
		ts, tu = _guidedBinaryRead(fh, "<2l")
		self.time = ts + tu/1.0e6


class SubSubSystemStatus(object):
	"""
	Python object that holds the status for the sub-subsystems in a SDM file.
	"""
	
	def __init__(self, feeStat=None, rpdStat=None, sepStat=None, arbStat=None, dp1Stat=None, dp2Stat=None, drStat=None):
		
		if feeStat is None:
			self.feeStat = [0 for n in xrange(ME_MAX_NFEE)]
		else:
			self.feeStat = feeStat
		
		if rpdStat is None:
			self.rpdStat = [0 for n in xrange(ME_MAX_NRPD)]
		else:
			self.rpdStat = rpdStat
			
		if sepStat is None:
			self.sepStat = [0 for n in xrange(ME_MAX_NSEP)]
		else:
			self.sepStat = sepStat
		
		if arbStat is None:
			self.arbStat = [0 for n in xrange(ME_MAX_NARB)]
		else:
			self.arbStat = arbStat
			
		if dp1Stat is None:
			self.dp1Stat = [0 for n in xrange(ME_MAX_NDP1)]
		else:
			self.dp1Stat = dp1Stat
		
		if dp2Stat is None:
			self.dp2Stat = [0 for n in xrange(ME_MAX_NDP2)]
		else:
			self.dp2Stat = dp2Stat
			
		if drStat is None:
			self.drStat = [0 for n in xrange(ME_MAX_NDR)]
		else:
			self.drStat = drStat
	
	def _binaryRead(self, fh):
		"""
		Given an open file handle, interpret it in the context of a 
		subsubsystem_status_struct C structure and update the Python instance accordingly.
		"""
		
		self.feeStat = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NFEE)
		self.rpdStat = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NRPD)
		self.sepStat = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NSEP)
		self.arbStat = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NARB)
		self.dp1Stat = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NDP1)
		self.dp2Stat = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NDP2)
		self.drStat  = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NDR)


class StationSettings(object):
	"""
	Python object that holds the status for the sub-subsystems in a SDM file.
	"""
	
	def __init__(self, report=None, update=None, fee=None, aspFlt=None, aspAT1=None, aspAT2=None, aspATS=None, 
				tbnGain=-1, drxGain=-1, nSTD=260):
		ME_MAX_NSTD = nSTD
					
		if report is None:
			self.report = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
		else:
			self.report = report
			
		if update is None:
			self.update = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
		else:
			self.update = update
			
		if fee is None:
			self.fee = [0 for n in xrange(ME_MAX_NSTD)]
		else:
			self.fee = fee
			
		if aspFlt is None:
			self.aspFlt = [0 for n in xrange(ME_MAX_NSTD)]
		else:
			self.aspFlt = aspFlt
			
		if aspAT1 is None:
			self.aspAT1 = [0 for n in xrange(ME_MAX_NSTD)]
		else:
			self.aspAT1 = aspAT1
			
		if aspAT2 is None:
			self.aspAT2 = [0 for n in xrange(ME_MAX_NSTD)]
		else:
			self.aspAT2 = aspAT2
			
		if aspATS is None:
			self.aspATS = [0 for n in xrange(ME_MAX_NSTD)]
		else:
			self.aspATS = aspATS
			
		self.tbnGain = tbnGain
		self.drxGain = drxGain
		
	def _binaryRead(self, fh):
		"""
		Given an open file handle, interpret it in the context of a 
		station_settings_struct C structure and update the Python instance accordingly.
		"""

		for ss in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
			self.report[ss] = _guidedBinaryRead(fh, "<h")
		for ss in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
			self.update[ss] = _guidedBinaryRead(fh, "<h")
			
		self.fee = _guidedBinaryRead(fh, "<%ih" % ME_MAX_NSTD)

		self.aspFlt = _guidedBinaryRead(fh, "<%ih" % ME_MAX_NSTD)
		self.aspAT1 = _guidedBinaryRead(fh, "<%ih" % ME_MAX_NSTD)
		self.aspAT2 = _guidedBinaryRead(fh, "<%ih" % ME_MAX_NSTD)
		self.aspATS = _guidedBinaryRead(fh, "<%ih" % ME_MAX_NSTD)
		
		self.tbnGain = _guidedBinaryRead(fh, "<h")
		self.drxgain = _guidedBinaryRead(fh, "<h")


class SDM(object):
	"""
	Python object that holds the various bits of information in a binary-packed Station
	Dynamic MIB file (SDM file).
	"""
	
	def __init__(self, station=None, shl=None, asp=None, dp=None, dr=None, status=None, antStatus=None, dpoStatus=None, settings=None):
		if station is None:
			self.station = SubSystemStatus('station')
		else:
			self.station = station
		if shl is None:
			self.shl = SubSystemStatus('shl')
		else:
			self.shl = shl
		if asp is None:
			self.asp = SubSystemStatus('asp')
		else:
			self.asp = asp
		if dp is None:
			self.dp  = SubSystemStatus('dp')
		else:
			self.dp = dp
		if dr is None:
			self.dr  = [SubSystemStatus('dr%i' % (n+1,)) for n in xrange(ME_MAX_NDR)]
		else:
			self.dr = dr
		
		if status is None:
			self.status = SubSubSystemStatus()
		else:
			self.status = status
		
		if antStatus is None:
			self.antStatus = [[0,0] for n in xrange(ME_MAX_NSTD)]
		else:
			self.antStatus = antStatus
		if dpoStatus is None:
			self.dpoStatus = [0 for n in xrange(ME_MAX_NDR)]
		else:
			self.dpoStatus = dpoStatus
			
		if settings is None:
			self.settings = StationSettings()
		else:
			self.settings = settings
		
	def updateAntennas(self, antennas):
		"""
		Given a list of :mod:`ls.common.stations.Antenna` instances, return a new list 
		of Antenna instances with updated antenna status codes.
		"""
		
		updatedAntennas = []
		for ant in antennas:
			updatedAntennas.apppend(ant)
			
			index = self.antStatus.index(ant.id)
			updatedAntennas[-1].status = self.antStatus[index]
			
		return updatedAntennas


def parseSDM(filename):
	"""
	Given a filename, read the file's contents into the SDM instance and return
	that instance.
	"""
		
	fh = open(filename, 'rb')
	
	# Create a new SDM instance
	dynamic = SDM()
	
	# Sub-system status sections
	dynamic.station._binaryRead(fh)
	dynamic.shl._binaryRead(fh)
	dynamic.asp._binaryRead(fh)
	dynamic.dp._binaryRead(fh)
	for n in xrange(ME_MAX_NDR):
			dynamic.dr[n]._binaryRead(fh)
	
	# Sub-sub-system status section
	dynamic.status._binaryRead(fh)
	
	# Antenna status and data path status
	dynamic.antStatus = _guidedBinaryRead(fh, "<%ii" % (2*ME_MAX_NSTD,))
	dynamic.dpoStatus = _guidedBinaryRead(fh, "<%ii" % ME_MAX_NDR)
	
	# Station settings section
	dynamic.settings._binaryRead(fh)
	
	fh.close()
	
	return dynamic

