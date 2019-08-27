# -*- coding: utf-8 -*-

"""lsl.sim - Simulate various types of LWA data.  The following follow DP
format writers are avaliable:
  * tbn
  * drx

In addition, there are two simulation modules to generate fake data sets::
  * dp  - generate DP-level data sets for basic signals and point source
  * vis - generate visibility data sets for use with the aipy module
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    