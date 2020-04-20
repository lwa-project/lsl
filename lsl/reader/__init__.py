"""
lsl.reader - Modular readers for the various LWA data formats:
 * tbw
 * tbn
 * drx
 * drspec
 * vdif
 * tbf
 * cor

A ring buffer for re-ordering TBN data is included in the 'buffer'
module.

Also include as part of this module are the LWA Development
Primities (LDP), a collection of utilites that make reading in
data files fun.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    