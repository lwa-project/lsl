# -*- coding: utf-8 -*-

"""lsl.sim - Simulate various types of LWA data.  The following follow DP
format writers are avaliable:
  * tbw
  * tbn
  * drx

In addition, there are two simulation modules to generate fake data sets::
  * dp  - generate DP-level data sets for basic signals and point source
  * vis - generate visibility data sets for use with the aipy module
"""

import errors

import s60
import tbw
import tbn
import drx

import dp
import vis

import nec_util
