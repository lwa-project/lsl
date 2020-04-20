"""
LWA Software Library

Provided packages:
  * lsl.common
  * lsl.reader
  * lsl.writer
  * lsl.correlator
  * lsl.statistics
  * lsl.sim
  * lsl.imaging
  * lsl.misc

Provided modules:
  * lsl.astro
  * lsl.catalog
  * lsl.skymap
  * lsl.transform

See the individual package descriptions for more information.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import subprocess

from lsl import version
from lsl.common import paths

__version__ = '0.8'
__author__ = "Jayce Dowell"


def test(): 
    """ 
    Internal test fixture 
    """ 
    
    eggPath = os.path.split(paths.MODULE)[0] 
    testPath = os.path.join(eggPath, 'tests', 'test_lsl.py') 
    subprocess.check_call([sys.executable, testPath, '-v'])

