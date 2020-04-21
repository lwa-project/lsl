"""
Module defining the LSL package test suite.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
__version__   = "0.3"
__author__    = "D. L. Wood"
__maintainer__ = "Jayce Dowell"

# Path adjustment
import os
sys.path.insert(1, os.path.dirname(__file__))

from lsl.misc import telemetry
telemetry.ignore()

from . import test_lsl
