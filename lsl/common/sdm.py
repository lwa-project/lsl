# -*- coding: utf-8 -*-

"""
Module for reading in an interpreting binary-packed Station Dynamic MIB (SDM) 
files (as defined in MCS0031, v5).
"""

import copy
import struct

__version__ = '0.1'
__revision__ = '$ Revision: 1 $'
__all__ = ['SubSystemStatus', 'SubSubSystemStatus', 'StationsSettings', 'SDM', '__version__', '__revision__', '__all__']


class SubSystemStatus(object):
	"""
	Python object that holds the status for a particular subsystem in a SDM 
	file.
	"""
	
	nBytes = 268
	
	def __init__(self, name, summary=6, info='UNK', time=0):
		self.name = name
		self.summary = int(summary)
		self.info = str(info)
		self.time = float(time)
		
	def __str__(self):
		return "%s at %i: %s [%i]" % (self.name, self.time, self.info, self.summary)
		
	def binaryFill(self, binaryData):
		"""
		Given a binary string, interpret it in the context of a 
		subsystem_status_struct C structure and update the Python instance accordingly.
		"""
		
		data = struct.unpack('<i256c2l', binaryData)
		
		self.summary = data[0]
		self.info = ''.join(data[1:-3])
		self.time = data[-2] + data[-1]/1.0e6


class SubSubSystemStatus(object):
	"""
	Python object that holds the status for the sub-subsystems in a SDM file.
	"""
	
	nBytes = 5464
	
	def __init__(self, feeStat=None, rpdStat=None, sepStat=None, arbStat=None, dp1Stat=None, dp2Stat=None, drStat=None, 
				nFEE=260, nRPD=520, nSEP=520, nARB=33, nDP1=26, nDP2=2, nDR=5):
		
		self.nFEE = nFEE
		self.nRPD = RPD
		self.nSEP = nSEP
		self.nARB = nARB
		self.nDP1 = nDP1
		self.nDP2 = nDP2
		self.nDR  = nDR
		
		if feeStat is None:
			self.feeStat = [0 for n in xrange(self.nFEE)]
		else:
			self.feeStat = feeStat
		
		if rpdStat is None:
			self.rpdStat = [0 for n in xrange(self.nRPD)]
		else:
			self.rpdStat = rpdStat
			
		if sepStat is None:
			self.sepStat = [0 for n in xrange(self.nSEP)]
		else:
			self.sepStat = sepStat
		
		if arbStat is None:
			self.arbStat = [0 for n in xrange(self.nARB)]
		else:
			self.arbStat = arbStat
			
		if dp1Stat is None:
			self.dp1Stat = [0 for n in xrange(self.nDP1)]
		else:
			self.dp1Stat = dp1Stat
		
		if dp2Stat is None:
			self.dp2Stat = [0 for n in xrange(self.nDP2)]
		else:
			self.dp2Stat = dp2Stat
			
		if drStat is None:
			self.drStat = [0 for n in xrange(self.nDR)]
		else:
			self.drStat = drStat
			
		self.updateBytes()
		
	def updateBytes(self):
		"""
		Update the size of the object based on the values of nFEE, nRPD, etc.
		"""
		
		self.nBytes = (self.nFEE + self.nRPD + self.nSEP + self.nARB + self.nDP1 + self.nDP2 + self.nDR)*4
	
	def binaryFill(self, binaryData):
		"""
		Given a binary string, interpret it in the context of a 
		subsubsystem_status_struct C structure and update the Python instance accordingly.
		"""
		
		nInts = self.nFEE + self.nRPD + self.nSEP + self.nARB + self.nDP1 + self.nDP2 + self.nDR
		data = struct.unpack('<%ii' % nInts, binaryData)
		
		self.feeStat = list(data[0:self.nFEE])
		del data[0:self.nFEE]
		
		self.rpdStat = list(data[0:self.nRPD])
		del data[0:self.nRPD]
		
		self.sepStat = list(data[0:self.nSEP])
		del data[0:self.nSEP]
		
		self.arbStat = list(data[0:self.nARB])
		del data[0:self.nARB]
		
		self.dp1Stat = list(data[0:self.nDP1])
		del data[0:self.nDP1]
		
		self.dp2Stat = list(data[0:self.nDP2])
		del data[0:self.nDP2]
		
		self.drStat  = list(data[0:self.nDR])
		del data[0:self.nDR]


class StationSettings(object):
	"""
	Python object that holds the status for the sub-subsystems in a SDM file.
	"""
	
	nBytes = 2640
	
	def __init__(self, report=None, update=None, fee=None, aspFlt=None, aspAT1=None, aspAT2=None, aspATS=None, 
				tbnGain=-1, drxGain=-1, nSTD=260):
		self.nSTD = nSTD
					
		if report is None:
			self.report = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
		else:
			self.report = report
			
		if update is None:
			self.update = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
		else:
			self.update = update
			
		if fee is None:
			self.fee = [0 for n in xrange(self.nSTD)]
		else:
			self.fee = fee
			
		if aspFlt is None:
			self.aspFlt = [0 for n in xrange(self.nSTD)]
		else:
			self.aspFlt = aspFlt
			
		if aspAT1 is None:
			self.aspAT1 = [0 for n in xrange(self.nSTD)]
		else:
			self.aspAT1 = aspAT1
			
		if aspAT2 is None:
			self.aspAT2 = [0 for n in xrange(self.nSTD)]
		else:
			self.aspAT2 = aspAT2
			
		if aspATS is None:
			self.aspATS = [0 for n in xrange(self.nSTD)]
		else:
			self.aspATS = aspATS
			
		self.tbnGain = tbnGain
		self.drxGain = drxGain
		
		self.updateBytes()
		
	def updateBytes(self):
		"""
		Update the size of the object based on the value of nSTD.
		"""
		
		self.nBytes = (18 + 5*self.nSTD + 2)*2
		
	def binaryFill(self, binaryData):
		"""
		Given a binary string, interpret it in the context of a 
		station_settings_struct C structure and update the Python instance accordingly.
		"""
		
		nShorts = 18 + 5*self.nSTD + 2
		data = struct.unpack('<%ih' % nShorts, binaryData)
		
		for ss in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
			self.report[ss] = data[0]
			del data[0]
		for ss in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
			self.update[ss] = data[0]
			del data[0]
			
		self.fee = list(data[0:self.nSTD])
		del data[0:self.nStand]
		
		self.aspFlt = list(data[0:self.nSTD])
		del data[0:self.nStand]
		
		self.aspAT1 = list(data[0:self.nSTD])
		del data[0:self.nStand]
		
		self.aspAT2 = list(data[0:self.nSTD])
		del data[0:self.nStand]
		
		self.aspATS = list(data[0:self.nSTD])
		del data[0:self.nStand]
		
		self.tbnGain = data[0]
		del data[0]
		self.drxgain = data[0]
		del data[0]

class SDM(object):
	"""
	Python object that holds the various bits of information in a binary-packed Station
	Dynamic MIB file (SDM file).
	"""
	
	nBytes = 268*9 + 5464 + 2080 + 20 + 2640
	
	def __init__(self, filename='', nSTD=260, nFEE=260, nRPD=520, nSEP=520, nARB=33, nDP1=26, nDP2=2, nDR=5):
		self.filename = filename
		
		self.nSTD = nSTD
		self.nFEE = nFEE
		self.nRPD = nRPD
		self.nSEP = nSEP
		self.nARB = nARB
		self.nDP1 = nDP1
		self.nDP2 = nDP2
		self.nDR  = nDR
		
		self.station = SubSystemStatus('station')
		self.shl = SubSystemStatus('shl')
		self.asp = SubSystemStatus('sdp')
		self.dp  = SubSystemStatus('dp')
		self.dr  = [SubSystemStatus('dr%i' % (n+1,)) for n in xrange(nDR)]
		self.status = SubSubSystemStatus(nFEE=self.nFEE, nRPD=self.nRPD, nSEP=self.nSEP, nARB=self.nARB, nDP1=self.nDP1, nDP2=self.nDP2, nDR=self.nDR)
		self.antStatus = [[0,0] for n in xrange(self.nSTD)]
		self.dpoStatus = [0 for n in xrange(self.nDR)]
		self.settings = StationSettings(nSTD=self.nSTD)
		
		self.updateBytes()
		
		if self.filename != '':
			self.load()
			
	def updateBytes(self):
		"""
		Update the size of the object based on the values of nSTD, nFEE, etc.
		
		.. note::
			updateBytes() makes sure the the relevant child instances of 
			SubSubSystemStatus and StationSettings have the right values and updates
			their sizes as well.
		"""
		
		# Update the values in self.status
		self.status.nFEE = self.nFEE
		self.status.nRPD = self.nRPD
		self.status.nSEP = self.nSEP
		self.status.nARB = self.nARB
		self.status.nDP1 = self.nDP1
		self.status.nDP2 = self.nDP2
		self.status.nDR  = self.nDR
		self.status.updateBytes()
		
		# Update the values in self.settings
		self.settings.nSTD = self.nSTD
		self.settings.updateBytes()
		
		# Build up the new size for the SDM instance
		self.nBytes  = self.station.nBytes + self.shl.nBytes + self.asp.nBytes + self.nDP.nBytes
		for n in xrange(self.nDR):
			self.nBytes += self.dr[n].nBytes
		self.nBytes += self.settings.status.nBytes
		self.nBytes += (self.nSTD*2*4)
		self.nBytes += (self.nDR*4)
		self.nBytes += self.settings.nBytes()
	
	def load(self, filename=''):
		"""
		Given a filename or using the filename the SDM instance was created with, 
		read the file's contents into the SDM instance.
		"""
		
		if filename is not '':
			self.filename = filename
		elif self.filename == '':
			raise RuntimeError('No filename provided')
		else:
			pass
		
		fh = open(self.filename, 'rb')
		data = fh.read(nBytes)
		fh.close()
		
		self.binaryFill(data)
	
	def binaryFill(self, binaryData):
		"""
		Given a binary string, interpret it in the context of a sdm_struct C 
		structure and update the Python instances accordingly.
		"""
		
		data = copy.deepcoy(binaryData)
		
		self.station.binaryFill(data[0:self.station.nBytes])
		del data[0:self.station.nBytes]
		
		self.shl.binaryFill(data[0:self.shl.nBytes])
		del data[0:self.shl.nBytes]
		
		self.asp.binaryFill(data[0:self.asp.nBytes])
		del data[0:self.asp.nBytes]
		
		self.dp.binaryFill(data[0:self.dp.nBytes])
		del data[0:self.dp.nBytes]
		
		for n in xrange(self.nDR):
			self.dr[n].binaryFill(data[0:self.dr[n].nBytes])
			del data[0:self.dr[n].nBytes]
		
		self.status.binaryFill(data[0:self.status.nBytes])
		del data[0:self.status.nBytes]
		
		#temp = list(struct.unpack('<%ii' % (2*self.nSTD,), data[0:2*self.nSTD*4]))
		#self.antStatus = [[temp[2*n], temp[2*n+1]] for n in xrange(self.nSTD)]
		self.antStatus = list(struct.unpack('<%ii' % (2*self.nSTD,), data[0:2*self.nSTD*4]))
		del data[0:(2*self.nSTD*4)]
		
		self.dpoStatus = list(struct.unpack('<%ii' % self.nDR, data[0:self.nDR*4]))
		del data[0:self.nDR*4]
		
		self.settings.binaryFill(data[0:self.settings.nBytes])
		del data[0:self.settings.nByts]
		
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
