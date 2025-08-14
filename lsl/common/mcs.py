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
 * ME_MAX_NSNAP - Maximum number of SNAP boards
 * ME_MAX_NSNAPCH - Number of channels per SNAP board
 * ME_MAX_SNAPID_LENGTH - Number of characters in the SNAP board ID name
 * ME_MAX_NSERVER - Maximum number of servers
 * ME_MAX_SERVERID_LENGTH - Number of characters in the server ID name
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
 * SSMIF_STRUCT - String representing the C structure of the binary SSMIF

The other functions:
 * Parse the binary packed metadata, 
"""

import re
import dbm
import enum
import math
import pytz
import numpy as np
import ctypes
import struct
from functools import reduce
from datetime import datetime

from astropy import units as astrounits
from astropy.coordinates import SphericalRepresentation, CartesianRepresentation

from lsl.common import ndp as ndpCommon

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.4'
__all__ = ['ME_SSMIF_FORMAT_VERSION', 'ME_MAX_NSTD', 'ME_MAX_NFEE', 'ME_MAX_FEEID_LENGTH', 'ME_MAX_RACK', 'ME_MAX_PORT', 
           'ME_MAX_NRPD', 'ME_MAX_RPDID_LENGTH', 'ME_MAX_NSEP', 'ME_MAX_SEPID_LENGTH', 'ME_MAX_SEPCABL_LENGTH', 
           'ME_MAX_NARB', 'ME_MAX_NARBCH', 'ME_MAX_ARBID_LENGTH', 'ME_MAX_NSNAP', 'ME_MAX_NSNAPCH', 'ME_MAX_SNAPID_LENGTH', 
           'ME_MAX_NSERVER', 'ME_MAX_SERVERID_LENGTH', 'ME_MAX_NDR', 'ME_MAX_DRID_LENGTH', 'ME_MAX_NPWRPORT', 
           'ME_MAX_SSNAME_LENGTH', 'LWA_MAX_NSTD', 'MIB_REC_TYPE_BRANCH', 'MIB_REC_TYPE_VALUE', 'MIB_INDEX_FIELD_LENGTH', 
           'MIB_LABEL_FIELD_LENGTH', 'MIB_VAL_FIELD_LENGTH', 
           'SSMIF_STRUCT', 'STATION_SETTINGS_STRUCT', 'SUBSYSTEM_STATUS_STRUCT', 'SUBSUBSYSTEM_STATUS_STRUCT', 
           'SSF_STRUCT', 'OSF_STRUCT', 'OSFS_STRUCT', 'BEAM_STRUCT', 'OSF2_STRUCT', 
           'delay_to_mcsd', 'mcsd_to_delay', 'gain_to_mcsg', 'mcsg_to_gain',
           'mjdmpm_to_datetime', 'datetime_to_mjdmpm', 'StatusCode', 'SummaryCode', 'SubsystemID', 'CommandID', 
           'ObservingMode', 'parse_c_struct', 'flat_to_multi', 'apply_pointing_correction', 'MIB', 'MIBEntry']


ME_SSMIF_FORMAT_VERSION = 10  # SSMIF format version code
ME_MAX_NSTD = 256             # Maximum number of stands that can be described
ME_MAX_NFEE = 256             # Maximum number of FEEs that can be described
ME_MAX_FEEID_LENGTH = 10      # Number of characters in FEE ID name
ME_MAX_RACK = 8               # Maximum number of racks?
ME_MAX_PORT = 50              # Maximum number of ports?
ME_MAX_NRPD = 512             # Maxmimum number of RPD cables
ME_MAX_RPDID_LENGTH = 25      # Number of characters in the RPD ID name
ME_MAX_NSEP = 512             # Maximum number of SEP cable connections
ME_MAX_SEPID_LENGTH = 25      # Number of characters in the SEP ID name
ME_MAX_SEPCABL_LENGTH = 25    # Number of characters in the SEP cable ID name
ME_MAX_NARB = 32              # Maximum number of ARX boards
ME_MAX_NARBCH = 16            # Number of ARX channels per board
ME_MAX_ARBID_LENGTH = 10      # Number of characters in the ARX ID name
ME_MAX_NSNAP = 16             # Maximum number of SNAP boards
ME_MAX_NSNAPCH = 64           # Number of channels per SNAP board
ME_MAX_SNAPID_LENGTH = 10     # Number of characters in the SNAP board ID name
ME_MAX_NSERVER = 9            # Maximum number of server
ME_MAX_SERVERID_LENGTH = 10   # Number of characters in the server ID name
ME_MAX_NDR = 4                # Maximum number of data recorders
ME_MAX_DRID_LENGTH = 10       # Number of characters in the DR ID name
ME_MAX_NPWRPORT = 50          # Maximum number of power ports
ME_MAX_SSNAME_LENGTH = 3      # Number of characters in the power port ID names, for codes used for PWR_NAME
LWA_MAX_NSTD = 256            # Maximum number of stands for the LWA
TPSS_FORMAT_VERSION = 6       # MCS0030 format version code
MIB_REC_TYPE_BRANCH = 0       # eType for MIB branch entries
MIB_REC_TYPE_VALUE = 1        # etype for MIB value entries
MIB_INDEX_FIELD_LENGTH = 12   # Number of characters in a MIB index field
MIB_LABEL_FIELD_LENGTH = 32   # Number of characters in a MIB label field
MIB_VAL_FIELD_LENGTH = 8192   # Number of characters in a MIB value field


SSMIF_STRUCT = """
    int    iFormatVersion;           /* FORMAT_VERSION */
    char   sStationID[3];            /* STATION_ID */
    double fGeoN;                    /* GEO_N */
    double fGeoE;                    /* GEO_E */
    double fGeoEl;                   /* GEO_EL */
    int    nStd;                     /* N_STD */
    double fStdLx[ME_MAX_NSTD];      /* STD_LX[] */
    double fStdLy[ME_MAX_NSTD];      /* STD_LY[] */
    double fStdLz[ME_MAX_NSTD];      /* STD_LZ[] */
    int    iAntStd[2*ME_MAX_NSTD];   /* ANT_STD[] */
    int    iAntOrie[2*ME_MAX_NSTD];  /* ANT_ORIE[] */
    int    iAntStat[2*ME_MAX_NSTD];  /* ANT_STAT[] */
    float  fAntTheta[2*ME_MAX_NSTD]; /* ANT_THETA[] */
    float  fAntPhi[2*ME_MAX_NSTD];   /* ANT_PHI[] */
    int    eAntDesi[2*ME_MAX_NSTD];  /* ANT_DESI[] */
    int    nFEE;                     /* N_FEE */
    char   sFEEID[ME_MAX_NFEE][ME_MAX_FEEID_LENGTH+1]; /* FEE_ID[] */
    int    iFEEStat[ME_MAX_NFEE];    /* FEE_STAT[] */
    int    eFEEDesi[ME_MAX_NFEE];    /* FEE_DESI[] */
    float  fFEEGai1[ME_MAX_NFEE];    /* FEE_GAI1[] */
    float  fFEEGai2[ME_MAX_NFEE];    /* FEE_GAI2[] */
    int    iFEEAnt1[ME_MAX_NFEE];    /* FEE_ANT1[] */
    int    iFEEAnt2[ME_MAX_NFEE];    /* FEE_ANT2[] */
    int    iFEERack[ME_MAX_NFEE];    /* FEE_RACK[] */
    int    iFEEPort[ME_MAX_NFEE];    /* FEE_PORT[] */
    int    nRPD;                     /* N_RPD */
    char   sRPDID[ME_MAX_NRPD][ME_MAX_RPDID_LENGTH+1]; /* RPD_ID[] */
    int    iRPDStat[ME_MAX_NRPD];    /* RPD_STAT[] */
    int    eRPDDesi[ME_MAX_NRPD];    /* RPD_DESI[] */
    float  fRPDLeng[ME_MAX_NRPD];    /* RPD_LENG[] */
    float  fRPDVF[ME_MAX_NRPD];      /* RPD_VF[] */
    float  fRPDDD[ME_MAX_NRPD];      /* RPD_DD[] */
    float  fRPDA0[ME_MAX_NRPD];      /* RPD_A0[] */
    float  fRPDA1[ME_MAX_NRPD];      /* RPD_A1[] */
    float  fRPDFref[ME_MAX_NRPD];    /* RPD_FREF[] */
    float  fRPDStr[ME_MAX_NRPD];     /* RPD_STR[] */
    int    iRPDAnt[ME_MAX_NRPD];     /* RPD_ANT[] */
    int    nSEP;                     /* N_SEP */
    char   sSEPID[ME_MAX_NSEP][ME_MAX_SEPID_LENGTH+1]; /* SEP_ID[] */
    int    iSEPStat[ME_MAX_NSEP];    /* SEP_STAT[] */
    char   sSEPCabl[ME_MAX_NSEP][ME_MAX_SEPCABL_LENGTH+1]; /* SEP_Cabl[] */
    float  fSEPLeng[ME_MAX_NSEP];    /* SEP_LENG[] */
    int    eSEPDesi[ME_MAX_NSEP];    /* SEP_DESI[] */
    float  fSEPGain[ME_MAX_NSEP];    /* SEP_GAIN[] */
    int    iSEPAnt[ME_MAX_NSEP];     /* SEP_ANT[] */
    int    nARB;                     /* N_ARB */
    int    nARBCH;                   /* N_ARBCH */
    char   sARBID[ME_MAX_NARB][ME_MAX_ARBID_LENGTH+1]; /* ARB_ID[] */
    int    iARBSlot[ME_MAX_NARB];    /* ARB_SLOT[] */
    int    eARBDesi[ME_MAX_NARB];    /* ARB_DESI[] */
    int    iARBRack[ME_MAX_NARB];    /* ARB_RACK[] */
    int    iARBPort[ME_MAX_NARB];    /* ARB_PORT[] */
    int    eARBStat[ME_MAX_NARB][ME_MAX_NARBCH];       /* ARB_STAT[][] */
    float  fARBGain[ME_MAX_NARB][ME_MAX_NARBCH];        /* ARB_GAIN[][] */
    int    iARBAnt[ME_MAX_NARB][ME_MAX_NARBCH];        /* ARB_ANT[][] */
    char   sARBIN[ME_MAX_NARB][ME_MAX_NARBCH][ME_MAX_ARBID_LENGTH+1]; /* ARB_IN[][] */
    char   sARBOUT[ME_MAX_NARB][ME_MAX_NARBCH][ME_MAX_ARBID_LENGTH+1]; /* ARB_OUT[][] */
    int    nSnapCh;                  /* N_SNAPCH */
    char   sSnapID[ME_MAX_NSNAP][ME_MAX_SNAPID_LENGTH+1]; /* SNAP_ID[] */
    char   sSnapSlot[ME_MAX_NSNAP][ME_MAX_SNAPID_LENGTH+1]; /* SNAP_SLOT[] */
    int    eSnapDesi[ME_MAX_NSNAP]; /* SNAP_DESI[] */
    int    eSnapStat[ME_MAX_NSNAP][ME_MAX_NSNAPCH];       /* SNAP_STAT[][] */
    char   sSnapINR[ME_MAX_NSNAP][ME_MAX_NSNAPCH][ME_MAX_SNAPID_LENGTH+1]; /* SNAP_INR[][] */
    char   sSnapINC[ME_MAX_NSNAP][ME_MAX_NSNAPCH][ME_MAX_SNAPID_LENGTH+1]; /* SNAP_INC[][] */
    int    iSnapAnt[ME_MAX_NSNAP][ME_MAX_NSNAPCH];        /* SNAP_ANT[][] */
    int    nServer;                  /* N_SERVER */
    char   sServerID[ME_MAX_NSERVER][ME_MAX_SERVERID_LENGTH+1]; /* SERVER_ID[] */
    char   sServerSlot[ME_MAX_NSERVER][ME_MAX_SERVERID_LENGTH+1]; /* SERVER_SLOT[] */
    int    eServerStat[ME_MAX_NSERVER];       /* SERVER_STAT[] */
    int    eServerDesi[ME_MAX_NSERVER];       /* SERVER_DESI[] */
    int    nDR;                     /* N_DR */
    int    eDRStat[ME_MAX_NDR];       /* DR_STAT[] */
    char   sDRID[ME_MAX_NDR][ME_MAX_DRID_LENGTH+1]; /* DR_ID[] */
    char   sDRPC[ME_MAX_NDR][ME_MAX_DRID_LENGTH+1]; /* DR_PC[] */
    int    iDRNDP[ME_MAX_NDR];       /* DR_NDP[] */
    int    nPwrRack;                /* N_PWR_RACK */
    int    nPwrPort[ME_MAX_RACK];   /* N_PWR_PORT[] */
    int    ePwrSS[ME_MAX_RACK][ME_MAX_NPWRPORT]; /* PWR_SS[][], converted to a LWA_SID_ value */
    char   sPwrName[ME_MAX_RACK][ME_MAX_NPWRPORT][ME_MAX_SSNAME_LENGTH+1]; /* PWR_NAME[][] */
    int    eCRA;                /* MCS_CRA */
    float  fPCAxisTh; /* PC_AXIS_TH */
    float  fPCAxisPh; /* PC_AXIS_PH */
    float  fPCRot;    /* PC_ROT */
"""


STATION_SETTINGS_STRUCT = """
    signed short int mrp_asp; // SESSION_MRP_ASP // MRP_ASP
    signed short int mrp_ndp; // SESSION_MRP_NDP // MRP_NDP
    signed short int mrp_dr1; // SESSION_MRP_DR1 // MRP_DR1
    signed short int mrp_dr2; // SESSION_MRP_DR2 // MRP_DR2
    signed short int mrp_dr3; // SESSION_MRP_DR3 // MRP_DR3
    signed short int mrp_dr4; // SESSION_MRP_DR4 // MRP_DR4
    signed short int mrp_dr5; // SESSION_MRP_DR5 // MRP_DR5
    signed short int mrp_shl; // SESSION_MRP_SHL // MRP_SHL
    signed short int mrp_mcs; // SESSION_MRP_MCS // MRP_MCS
    signed short int mup_asp; // SESSION_MUP_ASP // MUP_ASP
    signed short int mup_ndp; // SESSION_MUP_DP_ // MUP_NDP
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
    signed short int drx_gain; // OBS_DRX_GAIN // DRX_GAIN
"""


SUBSYSTEM_STATUS_STRUCT = """
    int summary;
    char info[256];
    long tv[2];
"""


SUBSUBSYSTEM_STATUS_STRUCT = """
    int    eFEEStat[ME_MAX_NFEE];                      /* FEE_STAT[] */
    int    eRPDStat[ME_MAX_NRPD];                      /* RPD_STAT[] */
    int    eSEPStat[ME_MAX_NSEP];                      /* SEP_STAT[] */
    int    eARBStat[ME_MAX_NARB][ME_MAX_NARBCH];       /* ARB_STAT[][] */
    int    eSnapStat[ME_MAX_NSNAP][ME_MAX_NSNAPCH];    /* SNAP_STAT[][] */
    int    eServerStat[ME_MAX_NSERVER];                /* SERVER_STAT[] */
    int    eDRStat[ME_MAX_NDR];                        /* DR_STAT[] */
"""

SSF_STRUCT = """
    unsigned short int FORMAT_VERSION;
    char PROJECT_ID[9];
    unsigned int SESSION_ID;
    unsigned short int SESSION_CRA;
    signed short int SESSION_DRX_BEAM;
    char SESSION_SPC[32];
    unsigned long int SESSION_START_MJD;
    unsigned long int SESSION_START_MPM;
    unsigned long int SESSION_DUR;
    unsigned int SESSION_NOBS;
    signed short int SESSION_MRP_ASP;
    signed short int SESSION_MRP_NDP;
    signed short int SESSION_MRP_DR1;
    signed short int SESSION_MRP_DR2;
    signed short int SESSION_MRP_DR3;
    signed short int SESSION_MRP_DR4;
    signed short int SESSION_MRP_DR5;
    signed short int SESSION_MRP_SHL;
    signed short int SESSION_MRP_MCS;
    signed short int SESSION_MUP_ASP;
    signed short int SESSION_MUP_NDP;
    signed short int SESSION_MUP_DR1;
    signed short int SESSION_MUP_DR2;
    signed short int SESSION_MUP_DR3;
    signed short int SESSION_MUP_DR4;
    signed short int SESSION_MUP_DR5;
    signed short int SESSION_MUP_SHL;
    signed short int SESSION_MUP_MCS;
    signed char SESSION_LOG_SCH;
    signed char SESSION_LOG_EXE;
    signed char SESSION_INC_SMIB;
    signed char SESSION_INC_DES;
"""


OSF_STRUCT = """
    unsigned short int FORMAT_VERSION;
    char               PROJECT_ID[9];
    unsigned int       SESSION_ID;
    signed short int   SESSION_DRX_BEAM;
    char               SESSION_SPC[32];
    unsigned int       OBS_ID; 
    unsigned long int  OBS_START_MJD;
    unsigned long int  OBS_START_MPM;
    unsigned long int  OBS_DUR;
    unsigned short int OBS_MODE;
    char               OBS_BDM[32];  /* added 140310 */
    float              OBS_RA;
    float              OBS_DEC;
    unsigned short int OBS_B;
    unsigned int       OBS_FREQ1;
    unsigned int       OBS_FREQ2;
    unsigned short int OBS_BW;
    unsigned int       OBS_STP_N;
    unsigned short int OBS_STP_RADEC;
"""


OSFS_STRUCT = """
    float              OBS_STP_C1;
    float              OBS_STP_C2;
    unsigned int       OBS_STP_T;
    unsigned int       OBS_STP_FREQ1;
    unsigned int       OBS_STP_FREQ2;
    unsigned short int OBS_STP_B;
"""


BEAM_STRUCT = """
    unsigned short int OBS_BEAM_DELAY[2*LWA_MAX_NSTD];
    signed short int   OBS_BEAM_GAIN[LWA_MAX_NSTD][2][2];
"""


OSF2_STRUCT = """
    signed short int   OBS_FEE[LWA_MAX_NSTD][2];
    signed short int   OBS_ASP_FLT[LWA_MAX_NSTD];
    signed short int   OBS_ASP_AT1[LWA_MAX_NSTD];
    signed short int   OBS_ASP_AT2[LWA_MAX_NSTD];
    signed short int   OBS_ASP_ATS[LWA_MAX_NSTD];
    unsigned int       OBS_TBT_SAMPLES;
    signed short int   OBS_DRX_GAIN;
    unsigned int       alignment;
"""


_cDecRE = re.compile(r'(?P<type>[a-z][a-z \t]+)[ \t]+(?P<name>[a-zA-Z_0-9]+)(\[(?P<d1>[\*\+A-Z_\d]+)\])?(\[(?P<d2>[\*\+A-Z_\d]+)\])?(\[(?P<d3>[\*\+A-Z_\d]+)\])?(\[(?P<d4>[\*\+A-Z_\d]+)\])?;')


def parse_c_struct(cStruct, char_mode='str', endianness='native', overrides=None):
    """
    Function to take a C structure declaration and build a ctypes.Structure out 
    of it with the appropriate alignment, character interpretation*, and endianness
    (little, big, network, or native).
    
    .. note::  ctypes converts character arrays to Python strings until the first null is
    incountered.  This behavior causes problems for multi-dimension arrays of null
    filled strings.  By setting char_mode to 'int', all char types are retuned as 
    bytes which can be converted to strings via chr().
    """
    
    # Process the macro overrides dictionary
    ndp_macros = {a: globals()[a] for a in __all__}
    if overrides is not None:
        ndp_macros.update(overrides)
    not_int = []
    for k in ndp_macros:
        if not isinstance(ndp_macros[k], (int, np.integer)):
            not_int.append(k)
    for k in not_int:
        del ndp_macros[k]
        
    # Figure out how to deal with character arrays
    if char_mode not in ('str', 'int'):
        raise RuntimeError(f"Unknown character mode: '{char_mode}'")
    if char_mode == 'str':
        baseCharType = ctypes.c_char
    else:
        baseCharType = ctypes.c_byte
    
    # Hold the basic fields and dimensions
    fields = []
    dims2 = {}
    
    # Split into lines and go!
    cStruct = cStruct.split('\n')
    for line in cStruct:
        ## Skip structure declaration, blank lines, comments, and lines too short to hold a 
        ## declaration
        line = line.strip().rstrip()
        if '{' in line or '}' in line:
            continue
        if line[:2] == '//':
            continue
        if len(line) < 5:
            continue
            
        ## RegEx the line to find the type, name, and dimensions (if needed) for
        ## the next structure variable
        mtch = _cDecRE.search(line)
        if mtch is None:
            raise RuntimeError(f"Unparseable line: '{line}'")
        
        dec = mtch.group('type')
        dec = dec.rstrip()
        name = mtch.group('name')
        
        try:
            d1 = mtch.group('d1')
            if d1 is not None:
                d1 = eval(d1, None, ndp_macros)
            d2 = mtch.group('d2')
            if d2 is not None:
                d2 = eval(d2, None, ndp_macros)
            d3 = mtch.group('d3')
            if d3 is not None:
                d3 = eval(d3, None, ndp_macros)
            d4 = mtch.group('d4')
            if d4 is not None:
                d4 = eval(d4, None, ndp_macros)
        except NameError:
            raise RuntimeError(f"Unknown value in array index: '{line}'")
        
        ## Basic data types
        if dec in ('signed int', 'int'):
            typ = ctypes.c_int
        elif dec == 'unsigned int':
            typ = ctypes.c_uint
        elif dec in ('signed short int', 'signed short', 'short int', 'short'):
            typ = ctypes.c_short
        elif dec in ('unsigned short int', 'unsigned short'):
            typ = ctypes.c_ushort
        elif dec in ('signed long int', 'signed long', 'long int', 'long'):
            typ = ctypes.c_long
        elif dec in ('unsigned long int', 'unsigned long'):
            typ = ctypes.c_ulong
        elif dec in ('signed long long', 'long long'):
            typ = ctypes.c_longlong
        elif dec == 'unsigned long long':
            typ = ctypes.c_uint64
        elif dec == 'float':
            typ = ctypes.c_float
        elif dec == 'double':
            typ = ctypes.c_double
        elif dec == 'char':
            typ = baseCharType
        elif dec == 'signed char':
            typ = ctypes.c_byte
        elif dec == 'unsigned char':
            typ = ctypes.c_ubyte
        else:
            raise RuntimeError(f"Unparseable line: '{line}' -> type: {dec}, name: {name}, dims: {d1}, {d2}, {d3}, {d4}")
        
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
    
    # ctypes creation - endianness
    endianness = endianness.lower()
    if endianness not in ('little', 'big', 'network', 'native'):
        raise RuntimeError(f"Unknown endianness: '{endianness}'")
    
    if endianness == 'little':
        endianness = ctypes.LittleEndianStructure
    elif endianness == 'big':
        endianness = ctypes.BigEndianStructure
    elif endianness == 'network':
        endianness = ctypes.BigEndianStructure
    else:
        endianness = ctypes.Structure
    
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
        _pack_ = 8    # Pack it like we are 64-bit
        dims = dims2
        
        def __str__(self):
            """
            Print out the structure in a nice easy-to-read formation that
            captures the various structure elements, their data types, and 
            their values.
            """
            
            out = ''
            for f,d in self._fields_:
                out += f"{f} ({d}): "+str(getattr(self, f))+'\n'
            return out
            
        def sizeof(self):
            """
            Return the size, in bytes, of the structure.
            """
            
            return ctypes.sizeof(self)
            
        def as_dict(self):
            """
            Return the structure as a simple Python dictionary keyed off the
            structure elements.
            """
            
            output = {}
            for f,d in self._fields_:
                output[f] = getattr(self, f)
            return output
    
    # Create and return
    return MyStruct()


def _two_byte_swap(value):
    return ((value & 0xFF) << 8) | ((value >> 8) & 0xFF)


def delay_to_mcsd(delay):
    """
    Given a delay in ns, convert it to a course and fine portion and into 
    the form expected by MCS in a custom beamforming SDF (little endian 
    16.12 unsigned integer).
    
    .. versionadded:: 0.6.3
    """
    
    return _two_byte_swap( ndpCommon.delay_to_dpd(delay) )


def mcsd_to_delay(delay):
    """
    Given delay value from an OBS_BEAM_DELAY field in a custom beamforming 
    SDF, return the delay in ns.
    
    .. versionadded:: 0.6.3
    """
    
    return ndpCommon.dpd_to_delay( _two_byte_swap(delay) )


def gain_to_mcsg(gain):
    """
    Given a gain (between 0 and 1), convert it to a gain in the form 
    expected by MCS in a custom beamforming SDF (little endian 16.1 
    signed integer).
    
    .. versionadded::0.6.3
    """
    
    return _two_byte_swap( ndpCommon.gain_to_dpg(gain) )


def mcsg_to_gain(gain):
    """
    Given a gain value from an OBS_BEAM_GAIN field in a custom beamforming
    SDF, return the decimal equivalent.

    .. versionadded:: 0.6.3
    """
    
    return ndpCommon.dpg_to_gain( _two_byte_swap(gain) )


def mjdmpm_to_datetime(mjd, mpm, tz=None):
    """
    Convert a MJD, MPM pair to a naive datetime instance.  If `tz` is not None
    the value is converted to a time zone-aware instance in the specified time
    zone.
    
    .. versionchanged:: 2.0.1
        Added the `tz` keyword and fixed the documentation.
        
    .. versionadded:: 0.5.2
    """
    
    unix = mjd*86400.0 + mpm/1000.0 - 3506716800.0
    dt = datetime.utcfromtimestamp(unix)
    if tz is not None:
        dt = pytz.utc.localize(dt)
        dt = dt.astimezone(tz)
    return dt


def datetime_to_mjdmpm(dt):
    """
    Convert a UTC datetime instance to a MJD, MPM pair (returned as a 
    two-element tuple).
    
    Based off: http://paste.lisp.org/display/73536
    
    .. versionchanged:: 2.0.1
        Better support for time zone-aware datetime instances
    
    .. versionadded:: 0.5.2
    """
    
    if dt.tzinfo is not None:
        dt = dt.astimezone(pytz.utc)
        
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


class StatusCode(enum.Enum):
    """
    MCS component status code.
    """
    
    NOTINSTALLED = 0
    BAD          = 1
    SUSPECT      = 2
    OK           = 3


class SummaryCode(enum.Enum):
    """
    MCS subsystem summary codes.
    """
    
    NULL    = 0
    NORMAL  = 1
    WARNING = 2
    ERROR   = 3
    BOOTING = 4
    SHUTDWN = 5
    UNK     = 6
    
    @property
    def description(self):
        if self.value == 0:
            return "Not normally used"
        elif self.value == 1:
            return "Normal"
        elif self.value == 2:
            return "Warning - issue(s) found, but still fully operational"
        elif self.value == 3:
            return "Error - problems found which limit or prevent proper function"
        elif self.value == 4:
            return "Booting - initializing; not yet fully operational"
        elif self.value == 5:
            return "Shutdown - shutting down; not ready for operation"
        elif self.value == 6:
            return "Unknown"
        else:
            raise ValueError(f"Unknown summary code '{code}'")


class SubsystemID(enum.Enum):
    """
    MCS subsystem IDs.
    """
    
    NU1 = 1
    NU2 = 2
    NU3 = 3
    NU4 = 4
    NU5 = 5
    NU6 = 6
    NU7 = 7
    NU8 = 8
    NU9 = 9
    MCS = 10
    SHL = 11
    ASP = 12
    DP  = 13
    DR1 = 14
    DR2 = 15
    DR3 = 16
    DR4 = 17
    DR5 = 18
    ADP = 19
    NDP = 20


class CommandID(enum.Enum):
    """
    MCS Command IDs.
    """
    
    MCSSHT = 0
    PNG    = 1
    RPT    = 2
    SHT    = 3
    INI    = 4
    TMP    = 5
    DIF    = 6
    PWR    = 7
    FIL    = 8
    AT1    = 9
    AT2    = 10
    ATS    = 11
    FPW    = 12
    RXP    = 13
    FEP    = 14
    TBW    = 15
    TBN    = 16
    DRX    = 17
    BAM    = 18
    FST    = 19
    CLK    = 20
    REC    = 21
    DEL    = 22
    STP    = 23
    GET    = 24
    CPY    = 25
    DMP    = 26
    FMT    = 27
    DWN    = 28
    UP_    = 29
    SEL    = 30
    SYN    = 31
    TST    = 32
    BUF    = 33
    NUL    = 34
    ESN    = 35
    ESF    = 36
    OBS    = 37
    OBE    = 38
    SPC    = 39
    TBF    = 40
    COR    = 41
    TBT    = 42
    TBS    = 43


class ObservingMode(enum.Enum):
    """
    MCS Observing Modes.
    """
    
    TRK_RADEC = 1
    TRK_SOL   = 2
    TRK_JOV   = 3
    STEPPED   = 4
    TBW       = 5
    TBN       = 6
    DIAG1     = 7
    TBF       = 8
    TRK_LUN   = 9
    TBT       = 10
    TBS       = 11


def flat_to_multi(inputList, *shape):
    """
    Convert a 1-D list into a multi-dimension list of shape 'shape'.
    
    .. versionchanged:: 3.0.0
        Switch from explicit 'dim1', 'dim2', etc. to a catch-all 'shape'.
    """
    
    if reduce(lambda x, y: x*y, shape, 1) != len(inputList):
        raise ValueError(f"Incompatiable dimensions: input={len(inputList)}, output={' by '.join([str(s) for s in shape])}")
        
    outputList = list(inputList)
    while len(shape) > 1:
        outputList = [outputList[i:i+shape[-1]] for i in range(0, len(outputList), shape[-1])]
        shape = shape[:-1]
        
    return outputList


def _get_rotation_matrix(theta, phi, psi, degrees=True):
    """
    Generate the rotation matrix for a rotation about a specified axis.
    
    http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    """
    
    if degrees:
        theta = theta * np.pi/180.0
        phi   = phi * np.pi/180.0
        psi   = psi * np.pi/180.0
        
    # Axis
    u = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
    
    # Rotation matrix
    rot  = np.eye(3)*np.cos(psi) 
    rot += np.sin(psi)*np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]]) 
    rot += (1-np.cos(psi))*np.tensordot(u, u, axes=0)
    
    return rot


def apply_pointing_correction(az, el, theta, phi, psi, lat=34.070, degrees=True):
    """
    Given a azimuth and elevation pair, and an axis to rotate about, 
    perform the rotation.
    
    http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    """
    
    # Get the rotation matrix
    rot = _get_rotation_matrix(theta, phi, psi, degrees=degrees)
    
    # Convert the az,alt coordinates to the unit vector
    if degrees:
        az = az*astrounits.deg
        el = el*astrounits.deg
    else:
        az = az*astrounits.rad
        el = el*astrounits.rad
    sph = SphericalRepresentation(az, el, 1)
    xyz = sph.to_cartesian()
    
    # Rotate
    yxz = xyz.xyz[[1,0,2]]
    xyzP = CartesianRepresentation((rot @ yxz)[[1,0,2]])
    
    sph = SphericalRepresentation.from_cartesian(xyzP)
    if degrees:
        azP = sph.lon.deg
        elP = sph.lat.deg
    else:
        azP = sph.lon.rad
        elP = sph.lat.rad
        
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
        out = f"MIB containing {nEntry} entries"
        if times:
            out += f" from {min(times)} to {max(times)}"
        return out
        
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
            
    def parse_init_file(self, filename):
        """
        Given a MCS MIB initialization file, i.e., ASP_MIB_init.dat, create a 
        dictionary that maps indicies to human-readable names that can be used
        to clarify a MIBEntry.  Return this dictionary.
        """
        
        # Open the file
        with open(filename, 'r') as fh:
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
        return True
        
    @classmethod
    def from_file(klass, filename, init_filename=None):
        """
        Given the name of a GDBM database file, initialize the MIB.  
        Optionally, use the name of the MCS MIB initialization file to
        help convert indicies to names.
        """
        
        # Parse the init. file (if we have one)
        mib = klass()
        if init_filename is not None:
            mib.parse_init_file(init_filename)
            
        # Make sure we have the .pag file
        if filename[-4:] == '.dir':
            filename = filename.replace('.dir', '.pag')
            
        # Open the database
        with dbm.open(filename, 'ru') as db:
            # Go!
            entry = db.firstkey()
            while entry is not None:
                value = db[entry]
                
                try:
                    record = MIBEntry.from_entry(value)
                    mib.entries[record.index] = record
                    
                except ValueError:
                    pass
                    
                entry = db.nextkey(entry)
                
        # Done
        return mib


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
        
        return f"Index: {self.index}; Value: {self.value}; Updated at {self.updateTime}"
        
    @staticmethod
    def _parse_value(value, dataType):
        """
        Convert an encoded value to something Pythonic (if possible).
        
        Understood codes:
         * NUL:   No data stored (e.g., branch head entry)
         * a####: printable (i.e., ASCII minus escape codes), #### = number of characters
                  e.g., "a3" means 3 printable ASCII-encoded characters
         * r####: raw data (not printable), #### = number of bytes
                  e.g., "r1024" means 1024 bytes of raw data
         * i1u:   integer, 1 byte,  unsigned, little-endian (=uint8)
         * i2u:   integer, 2 bytes, unsigned, litte-endian (=uint16)
         * i4u:   integer, 4 bytes, unsigned, little-endian (=uint32)
         * i4s:   integer, 4 bytes, signed, little-endian (=int32)
         * i8u:   integer, 8 bytes, unsigned, litte-endian (=uint64)
         * f4:    float, 4 bytes, little-endian (=float32)
         * f4r:   float, 4 bytes, big-ending (=float32)
         * f8:    float, 8 bytes, little-endian (=float64)
         * f8r:   float, 8 bytes, big-ending (=float64)
        """
        
        if dataType == 'NUL' or dataType == '':
            try:
                value = str(value, 'utf-8')
            except TypeError:
                value = str(value)
            return value
        elif dataType[0] == 'a':
            try:
                value = str(value, 'utf-8')
            except TypeError:
                value = str(value)
            return value
        elif dataType[0] == 'r':
            try:
                value = str(value, 'utf-8')
            except TypeError:
                value = str(value)
            return value

        elif dataType[:3] == 'i1u':
            try:
                return struct.unpack('<1B', value)[0]
            except struct.error:
                return 0
        elif dataType[:3] == 'i2u':
            try:
                return struct.unpack('<1H', value)[0]
            except struct.error:
                return 0
        elif dataType[:3] == 'i4u':
            try:
                return struct.unpack('<1I', value)[0]
            except struct.error:
                return 0
        elif dataType[:3] == 'i4s':
            try:
                return struct.unpack('<1i', value)[0]
            except struct.error:
                return 0
        elif dataType[:3] == 'i8u':
            try:
                return struct.unpack('<1Q', value)[0]
            except struct.error:
                return 0
        elif dataType[:3] == 'f4r':
            try:
                return struct.unpack('>1f', value)[0]
            except struct.error:
                return 0.0
        elif dataType[:2] == 'f4':
            try:
                return struct.unpack('<1f', value)[0]
            except struct.error:
                return 0.0
        elif dataType[:3] == 'f8r':
            try:
                return struct.unpack('>1d', value)[0]
            except struct.error:
                return 0.0
        elif dataType[:2] == 'f8':
            try:
                return struct.unpack('<1d', value)[0]
            except struct.error:
                return 0.0
        else:
            raise ValueError(f"Unknown data type '{dataType}'")
            
    @classmethod
    def from_entry(klass, value):
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
        record = parse_c_struct("""
            int eType;
            char index[%i];
            char val[%i];
            char type_dbm[6];
            char type_icd[6];
            long tv[2];
            """ % (MIB_INDEX_FIELD_LENGTH, MIB_VAL_FIELD_LENGTH), endianness='little')
            
        # Squeeze the binary string into the structure using ctypes magic
        temp = ctypes.create_string_buffer(value)
        ctypes.memmove(ctypes.addressof(record), ctypes.addressof(temp), ctypes.sizeof(record))
        
        # Bytes to Unicode for Python3
        try:
            index   = str(record.index, 'utf-8')
            dbmType = str(record.type_dbm, 'utf-8')
            icdType = str(record.type_icd, 'utf-8')
        except Exception as e:
            index   = str(record.index)
            dbmType = str(record.type_dbm)
            icdType = str(record.type_icd)
            
        # Validate
        if record.eType == MIB_REC_TYPE_BRANCH:
            raise ValueError("Cannot interpret MIB branch entries")
        try:
            if index[0] not in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
                raise ValueError(f"Entry index '{record.index}' does not appear to be numeric")
        except IndexError:
            raise ValueError(f"Entry index '{record.index}' does not appear to be numeric")
            
        # Basic information
        mibe = klass()
        mibe.eType = int(record.eType)
        mibe.index = index
        mibe.value = mibe._parse_value(record.val, dbmType)
        mibe.dbmType = dbmType
        mibe.icdType = icdType
        mibe._tv = (int(record.tv[0]), int(record.tv[1]))
        
        # Time
        mibe.updateTime = datetime.utcfromtimestamp(record.tv[0] + record.tv[1]/1e9)
        
        return mibe
