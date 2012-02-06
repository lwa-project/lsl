# -*- coding: utf-8 -*-

"""
Module for reading in an interpreting binary-packed Station Dynamic MIB (SDM) 
files (as defined in MCS0031, v5).
"""

import copy

from lsl.common.mcs import *
from datetime import datetime

__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['SubSystemStatus', 'SubSubSystemStatus', 'StationsSettings', 'SDM', 'parseSDM', '__version__', '__revision__', '__all__']


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
		return "%s at %s: %s [%i = %s]" % (self.name, datetime.utcfromtimestamp(self.time), self.info, self.summary, summary2string(self.summary))
		
	def _binaryRead(self, fh):
		"""
		Given an open file handle, interpret it in the context of a 
		subsystem_status_struct C structure and update the Python instance accordingly.
		"""
		
		sssStruct = parseCStruct("""
		int summary;
		char info[256];
		long tv[2];
		""", endianness='little')
		
		fh.readinto(sssStruct)
		
		self.summary = sssStruct.summary
		self.info = sssStruct.info
		ts, tu = sssStruct.tv
		self.time = ts + tu/1.0e6


class SubSubSystemStatus(object):
	"""
	Python object that holds the status for the sub-subsystems in a SDM file.
	"""
	
	def __init__(self, fee=None, rpd=None, sep=None, arb=None, dp1=None, dp2=None, dr=None):
		
		if fee is None:
			self.fee = [0 for n in xrange(ME_MAX_NFEE)]
		else:
			self.fee = fee
		
		if rpd is None:
			self.rpd = [0 for n in xrange(ME_MAX_NRPD)]
		else:
			self.rpd = rpd
			
		if sep is None:
			self.sep = [0 for n in xrange(ME_MAX_NSEP)]
		else:
			self.sep = sep
		
		if arb is None:
			self.arb = [0 for n in xrange(ME_MAX_NARB)]
		else:
			self.arb = arb
			
		if dp1 is None:
			self.dp1 = [0 for n in xrange(ME_MAX_NDP1)]
		else:
			self.dp1 = dp1
		
		if dp2 is None:
			self.dp2 = [0 for n in xrange(ME_MAX_NDP2)]
		else:
			self.dp2 = dp2
			
		if dr is None:
			self.dr = [0 for n in xrange(ME_MAX_NDR)]
		else:
			self.dr = dr
	
	def _binaryRead(self, fh):
		"""
		Given an open file handle, interpret it in the context of a 
		subsubsystem_status_struct C structure and update the Python instance accordingly.
		"""
		
		ssssStruct = parseCStruct("""
		int    eFEEStat[ME_MAX_NFEE];                /* FEE_STAT[] */
		int    eRPDStat[ME_MAX_NRPD];                /* RPD_STAT[] */
		int    eSEPStat[ME_MAX_NSEP];                /* SEP_STAT[] */
		int    eARBStat[ME_MAX_NARB][ME_MAX_NARBCH]; /* ARB_STAT[][] */
		int    eDP1Stat[ME_MAX_NDP1][ME_MAX_NDP1CH]; /* DP1_STAT[][] */
		int    eDP2Stat[ME_MAX_NDP2];                /* DP2_STAT[] */
		int    eDRStat[ME_MAX_NDR];                  /* DR_STAT[] */
		""", endianness='little')
		
		fh.readinto(ssssStruct)
		
		self.fee = list(ssssStruct.eFEEStat)
		self.rpd = list(ssssStruct.eRPDStat)
		self.sep = list(ssssStruct.eSEPStat)
		self.arb = single2multi(ssssStruct.eARBStat, *ssssStruct.dims['eARBStat'])
		self.dp1 = single2multi(ssssStruct.eDP1Stat, *ssssStruct.dims['eDP1Stat'])
		self.dp2 = list(ssssStruct.eDP2Stat)
		self.dr  = list(ssssStruct.eDRStat)


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

		ssStruct = parseCStruct("""
		signed short int mrp_asp; // SESSION_MRP_ASP // MRP_ASP
		signed short int mrp_dp;  // SESSION_MRP_DP_ // MRP_DP_
		signed short int mrp_dr1; // SESSION_MRP_DR1 // MRP_DR1
		signed short int mrp_dr2; // SESSION_MRP_DR2 // MRP_DR2
		signed short int mrp_dr3; // SESSION_MRP_DR3 // MRP_DR3
		signed short int mrp_dr4; // SESSION_MRP_DR4 // MRP_DR4
		signed short int mrp_dr5; // SESSION_MRP_DR5 // MRP_DR5
		signed short int mrp_shl; // SESSION_MRP_SHL // MRP_SHL
		signed short int mrp_mcs; // SESSION_MRP_MCS // MRP_MCS
		signed short int mup_asp; // SESSION_MUP_ASP // MUP_ASP
		signed short int mup_dp;  // SESSION_MUP_DP_ // MUP_DP_
		signed short int mup_dr1; // SESSION_MUP_DR1 // MUP_DR1
		signed short int mup_dr2; // SESSION_MUP_DR2 // MUP_DR2
		signed short int mup_dr3; // SESSION_MUP_DR3 // MUP_DR3
		signed short int mup_dr4; // SESSION_MUP_DR4 // MUP_DR4
		signed short int mup_dr5; // SESSION_MUP_DR5 // MUP_DR5
		signed short int mup_shl; // SESSION_MUP_SHL // MUP_SHL
		signed short int mup_mcs; // SESSION_MUP_MCS // MUP_MCS
		signed short int fee[LWA_MAX_NSTD];     // OBS_FEE[LWA_MAX_NSTD][2]  // FEE[LWA_MAX_NSTD]
		signed short int asp_flt[LWA_MAX_NSTD]; // OBS_ASP_FLT[LWA_MAX_NSTD] // ASP_FLT[LWA_MAX_NSTD]
		signed short int asp_at1[LWA_MAX_NSTD]; // OBS_ASP_AT1[LWA_MAX_NSTD] // ASP_AT1[LWA_MAX_NSTD]
		signed short int asp_at2[LWA_MAX_NSTD]; // OBS_ASP_AT2[LWA_MAX_NSTD] // ASP_AT2[LWA_MAX_NSTD]
		signed short int asp_ats[LWA_MAX_NSTD]; // OBS_ASP_ATS[LWA_MAX_NSTD] // ASP_ATS[LWA_MAX_NSTD]
		signed short int tbn_gain; // OBS_TBN_GAIN // TBN_GAIN
		signed short int drx_gain; // OBS_DRX_GAIN // DRX_GAIN
		""", endianness='little')
		
		fh.readinto(ssStruct)

		self.report['ASP'] = ssStruct.mrp_asp
		self.report['DP_'] = ssStruct.mrp_dp
		self.report['DR1'] = ssStruct.mrp_dr1
		self.report['DR2'] = ssStruct.mrp_dr2
		self.report['DR3'] = ssStruct.mrp_dr3
		self.report['DR4'] = ssStruct.mrp_dr4
		self.report['DR5'] = ssStruct.mrp_dr5
		self.report['SHL'] = ssStruct.mrp_shl
		self.report['MCS'] = ssStruct.mrp_mcs
		
		self.update['ASP'] = ssStruct.mup_asp
		self.update['DP_'] = ssStruct.mup_dp
		self.update['DR1'] = ssStruct.mup_dr1
		self.update['DR2'] = ssStruct.mup_dr2
		self.update['DR3'] = ssStruct.mup_dr3
		self.update['DR4'] = ssStruct.mup_dr4
		self.update['DR5'] = ssStruct.mup_dr5
		self.update['SHL'] = ssStruct.mup_shl
		self.update['MCS'] = ssStruct.mup_mcs
		
		self.fee = single2multi(ssStruct.fee, *ssStruct.dims['fee'])

		self.aspFlt = list(ssStruct.asp_flt)
		self.aspAT1 = list(ssStruct.asp_at1)
		self.aspAT2 = list(ssStruct.asp_at2)
		self.aspATS = list(ssStruct.asp_ats)
		
		self.tbnGain = ssStruct.tbn_gain
		self.drxgain = ssStruct.drx_gain


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
	adpsStruct = parseCStruct("""
	int ant_stat[ME_MAX_NSTD][2]; /* corresponds to sc.Stand[i].Ant[k].iSS, but dynamically updated */
	int dpo_stat[ME_MAX_NDR];     /* corresponds to sc.DPO[i].iStat, but dynamically updated */
	""", endianness='little')
	
	fh.readinto(adpsStruct)
	
	dynamic.antStatus = single2multi(adpsStruct.ant_stat, *adpsStruct.dims['ant_stat'])
	dynamic.dpoStatus = single2multi(adpsStruct.dpo_stat, *adpsStruct.dims['dpo_stat'])
	
	# Station settings section
	dynamic.settings._binaryRead(fh)
	
	fh.close()
	
	return dynamic

