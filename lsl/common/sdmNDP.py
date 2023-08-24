"""
Module for reading in an interpreting binary-packed Station Dynamic MIB (SDM) 
files (as defined in MCS0031, v5).
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
from datetime import datetime

from lsl.common.mcsNDP import summary_to_string, parse_c_struct, flat_to_multi, \
                        STATION_SETTINGS_STRUCT, SUBSYSTEM_STATUS_STRUCT, SUBSUBSYSTEM_STATUS_STRUCT, \
                        ME_MAX_NSTD, ME_MAX_NFEE, ME_MAX_NRPD, ME_MAX_NSEP, ME_MAX_NARB, \
                        ME_MAX_NSNAP, ME_MAX_NSERVER, ME_MAX_NDR

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.1'
__all__ = ['SubSystemStatus', 'SubSubSystemStatus', 'StationSettings', 'SDM', 'parse_sdm']


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
        return "%s at %s: %s [%i = %s]" % (self.name, datetime.utcfromtimestamp(self.time), self.info, self.summary, summary_to_string(self.summary))
        
    def binary_read(self, fh):
        """
        Given an open file handle, interpret it in the context of a 
        subsystem_status_struct C structure and update the Python instance accordingly.
        """
        
        ## They are the same size so this really doesn't matter
        sssStruct = parse_c_struct(SUBSYSTEM_STATUS_STRUCT, endianness='little')
        
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
            self.fee = [0 for n in range(ME_MAX_NFEE)]
        else:
            self.fee = fee
        
        if rpd is None:
            self.rpd = [0 for n in range(ME_MAX_NRPD)]
        else:
            self.rpd = rpd
            
        if sep is None:
            self.sep = [0 for n in range(ME_MAX_NSEP)]
        else:
            self.sep = sep
        
        if arb is None:
            self.arb = [0 for n in range(ME_MAX_NARB)]
        else:
            self.arb = arb
            
        if dp1 is None:
            self.dp1 = [0 for n in range(ME_MAX_NSNAP)]
        else:
            self.dp1 = dp1
        
        if dp2 is None:
            self.dp2 = [0 for n in range(ME_MAX_NSERVER)]
        else:
            self.dp2 = dp2
            
        if dr is None:
            self.dr = [0 for n in range(ME_MAX_NDR)]
        else:
            self.dr = dr
            
    def binary_read(self, fh, overrides=None):
        """
        Given an open file handle, interpret it in the context of a 
        subsubsystem_status_struct C structure and update the Python instance accordingly.
        """
        
        # Figure out what to do
        ssssStruct = parse_c_struct(SUBSUBSYSTEM_STATUS_STRUCT, endianness='little', overrides=overrides)
        
        # Read
        fh.readinto(ssssStruct)
        
        # Parse and save
        self.fee = list(ssssStruct.eFEEStat)
        self.rpd = list(ssssStruct.eRPDStat)
        self.sep = list(ssssStruct.eSEPStat)
        self.arb = flat_to_multi(ssssStruct.eARBStat, *ssssStruct.dims['eARBStat'])
        self.dp1 = flat_to_multi(ssssStruct.eSnapStat, *ssssStruct.dims['eSnapStat'])
        self.dp2 = list(ssssStruct.eServerStat)
        self.dr  = list(ssssStruct.eDRStat)


class StationSettings(object):
    """
    Python object that holds the status for the sub-subsystems in a SDM file.
    """
    
    def __init__(self, report=None, update=None, fee_power=None, asp_filter=None, asp_atten_1=None, asp_atten_2=None, asp_atten_split=None, 
                drx_gain=-1, tbf_gain=-1):
        if report is None:
            self.report = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
        else:
            self.report = report
            
        if update is None:
            self.update = {'ASP': -1, 'DP_': -1, 'DR1': -1, 'DR2': -1, 'DR3': -1, 'DR4': -1, 'DR5': -1, 'SHL': -1, 'MCS': -1}
        else:
            self.update = update
            
        if fee_power is None:
            self.fee_power = [0 for n in range(ME_MAX_NSTD)]
        else:
            self.fee_power = fee_power
            
        if asp_filter is None:
            self.asp_filter = [0 for n in range(ME_MAX_NSTD)]
        else:
            self.asp_filter = asp_filter
            
        if asp_atten_1 is None:
            self.asp_atten_1 = [0 for n in range(ME_MAX_NSTD)]
        else:
            self.asp_atten_1 = asp_atten_1
            
        if asp_atten_2 is None:
            self.asp_atten_2 = [0 for n in range(ME_MAX_NSTD)]
        else:
            self.asp_atten_2 = asp_atten_2
            
        if asp_atten_split is None:
            self.asp_atten_split = [0 for n in range(ME_MAX_NSTD)]
        else:
            self.asp_atten_split = asp_atten_split
            
        self.drx_gain = drx_gain
        self.tbf_gain = tbf_gain
        
    def binary_read(self, fh):
        """
        Given an open file handle, interpret it in the context of a 
        station_settings_struct C structure and update the Python instance accordingly.
        """
        
        # Figure out what to do
        ssStruct = parse_c_struct(STATION_SETTINGS_STRUCT, endianness='little')
        
        # Read
        fh.readinto(ssStruct)
        
        # Parse and save
        ## Common
        self.report['ASP'] = ssStruct.mrp_asp
        self.report['DP_'] = ssStruct.mrp_dp
        self.report['DR1'] = ssStruct.mrp_dr1
        self.report['DR2'] = ssStruct.mrp_dr2
        self.report['DR3'] = ssStruct.mrp_dr3
        self.report['DR4'] = ssStruct.mrp_dr4
        self.report['SHL'] = ssStruct.mrp_shl
        self.report['MCS'] = ssStruct.mrp_mcs
        
        self.update['ASP'] = ssStruct.mup_asp
        self.update['DP_'] = ssStruct.mup_dp
        self.update['DR1'] = ssStruct.mup_dr1
        self.update['DR2'] = ssStruct.mup_dr2
        self.update['DR3'] = ssStruct.mup_dr3
        self.update['DR4'] = ssStruct.mup_dr4
        self.update['SHL'] = ssStruct.mup_shl
        self.update['MCS'] = ssStruct.mup_mcs
        
        self.fee_power = flat_to_multi(ssStruct.fee, *ssStruct.dims['fee'])
        
        self.asp_filter = list(ssStruct.asp_flt)
        self.asp_atten_1 = list(ssStruct.asp_at1)
        self.asp_atten_2 = list(ssStruct.asp_at2)
        self.asp_atten_split = list(ssStruct.asp_ats)
        
        self.drx_gain = ssStruct.drx_gain
        self.tbf_gain = ssStruct.tbf_gain


class SDM(object):
    """
    Python object that holds the various bits of information in a binary-packed Station
    Dynamic MIB file (SDM file).
    """
    
    def __init__(self, station=None, shl=None, asp=None, dp=None, dr=None, status=None, ant_status=None, dpo_status=None, settings=None):
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
            self.dr  = [SubSystemStatus('dr%i' % (n+1,)) for n in range(ME_MAX_NDR)]
        else:
            self.dr = dr
        
        if status is None:
            self.status = SubSubSystemStatus()
        else:
            self.status = status
        
        if ant_status is None:
            self.ant_status = [[0,0] for n in range(ME_MAX_NSTD)]
        else:
            self.ant_status = ant_status
        if dpo_status is None:
            self.dpo_status = [0 for n in range(ME_MAX_NDR)]
        else:
            self.dpo_status = dpo_status
            
        if settings is None:
            self.settings = StationSettings()
        else:
            self.settings = settings
            
    def update_antennas(self, antennas):
        """
        Given a list of :class:`lsl.common.stations.Antenna` instances, return a new list 
        of Antenna instances with updated antenna status codes.
        """
        
        updatedAntennas = []
        for ant in antennas:
            updatedAntennas.append(ant)
            
            index = self.ant_status.index(ant.id)
            updatedAntennas[-1].status = self.ant_status[index]
            
        return updatedAntennas


def parse_sdm(filename):
    """
    Given a filename, read the file's contents into the SDM instance and return
    that instance.
    """
    
    # Open the file
    with open(filename, 'rb') as fh:
        # Create a new SDM instance
        dynamic = SDM()
        
        # Sub-system status sections
        dynamic.station.binary_read(fh)
        dynamic.shl.binary_read(fh)
        dynamic.asp.binary_read(fh)
        dynamic.dp.binary_read(fh)
        mark = fh.tell()
        try:
            for n in range(ME_MAX_NDR):
                dynamic.dr[n].binary_read(fh)
                
            # Sub-sub-system status section
            dynamic.status.binary_read(fh)
            
            # Antenna status and data path status
            ndpsStruct = parse_c_struct("""
            int ant_stat[ME_MAX_NSTD][2]; /* corresponds to sc.Stand[i].Ant[k].iSS, but dynamically updated */
            int dpo_stat[ME_MAX_NDR];     /* corresponds to sc.DPO[i].iStat, but dynamically updated */
            """, endianness='little')
            
            fh.readinto(ndpsStruct)
            
        except IOError:
            fh.seek(mark)
            overrides = {'ME_MAX_NDR': 3}
            
            for n in range(3):
                dynamic.dr[n].binary_read(fh)
                
            # Sub-sub-system status section
            dynamic.status.binary_read(fh, overrides=overrides)
            
            # Antenna status and data path status
            ndpsStruct = parse_c_struct("""
            int ant_stat[ME_MAX_NSTD][2]; /* corresponds to sc.Stand[i].Ant[k].iSS, but dynamically updated */
            int dpo_stat[ME_MAX_NDR];     /* corresponds to sc.DPO[i].iStat, but dynamically updated */
            """, endianness='little', overrides=overrides)
            
            fh.readinto(ndpsStruct)
            
        dynamic.ant_status = flat_to_multi(ndpsStruct.ant_stat, *ndpsStruct.dims['ant_stat'])
        dynamic.dpo_status = flat_to_multi(ndpsStruct.dpo_stat, *ndpsStruct.dims['dpo_stat'])
        
        # Station settings section
        dynamic.settings.binary_read(fh)
        
    return dynamic
