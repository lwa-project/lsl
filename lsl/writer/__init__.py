# -*- coding: utf-8 -*

"""
lsl.writer - Writers for exporting LWA data to various file formats.  
Formats include:
 * sdfits        - SDFITS filtes for storing single-dish data such as DRX, 
 * fitsidi       - FITS IDI writer for storing correlator output,
 * uvfits        - UVFITS writer also for storing correlator output,
 * measurmentset - CASA MS writer also for storing correlator output, and
 * vdif          - write data to the VLBI Data Interchange Format.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range



