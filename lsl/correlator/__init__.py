# -*- coding: utf-8 -*

"""
lsl.correlator - Modules dealing with correlating LWA data and analyzing
the results.  Modules included are:
 * uvutil   -   get stand positions, baseline uvw coordinates, cable delays, 
                etc.,
 * fx         - calculate spectra from LWA observations and correlate data using 
                and FX-style correlator, and
 * filterbank - calculate spectra using a polyphase filter bank.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    