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
"""

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['ME_SSMIF_FORMAT_VERSION', 'ME_MAX_NSTD', 'ME_MAX_NFEE', 'ME_MAX_FEEID_LENGTH', 'ME_MAX_RACK', 'ME_MAX_PORT', 
			'ME_MAX_NRPD', 'ME_MAX_RPDID_LENGTH', 'ME_MAX_NSEP', 'ME_MAX_SEPID_LENGTH', 'ME_MAX_SEPCABL_LENGTH', 
			'ME_MAX_NARB', 'ME_MAX_NARBCH', 'ME_MAX_ARBID_LENGTH', 'ME_MAX_NDP1', 'ME_MAX_NDP1CH', 'ME_MAX_DP1ID_LENGTH', 
			'ME_MAX_NDP2', 'ME_MAX_DP2ID_LENGTH', 'ME_MAX_NDR', 'ME_MAX_DRID_LENGTH', 'ME_MAX_NPWRPORT', 
			'ME_MAX_SSNAME_LENGTH', '__version__', '__revision__', '__all__']

ME_SSMIF_FORMAT_VERSION = 4	# SSMIF format version code
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
