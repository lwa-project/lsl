# -*- coding: utf-8 -*

"""LWA Software Library
Provided packages:
  * lsl.common
  * lsl.reader
  * lsl.writer
  * lsl.correlator
  * lsl.statistics
  * lsl.sim
  * lsl.misc

Provided modules:
  * lsl.astro_array
  * lsl.astro
  * lsl.catalog
  * lsl.skymap
  * lsl.transform

Provided tests:
  * lsl.test()

See the individual package descriptions for more information.
"""

__version__ = '0.3'
__revision__ = '$ Revision: 5 $'
__author__ = "Jayce Dowell"

import common
import reader
import writer
import correlator
import statistics
import sim
import misc

import astro_array
import astro
import catalog
import skymap
import transform

import tests
from tests.test_lsl import main as test
