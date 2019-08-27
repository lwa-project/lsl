# -*- coding: utf-8 -*
    
"""
lsl.common - Common information for the LSL package.  Including:
 * paths         - paths to the module and its data directories, 
 * stations      - information about LWA stations, 
 * dp            - information about the DP system, 
 * mcs           - information about DP-compatible stations at large, 
 * sdf           - read in and interpret SDF files, 
 * metabundle    - read in and work with MCS metadata,
 * adp           - information about the ADP system,
 * mcsADP        - information about ADP-compatible stations at large, 
 * sdfADP        - read in and interpret ADP-compatible SDF files, and
 * metabundleADP - read in and work with ADP-compatible MCS metadata.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    