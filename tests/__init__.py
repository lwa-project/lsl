# -*- coding: utf-8 -*-

"""
Module defining the LSL package test suite.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
__revision__  = "$Rev$"
__version__   = "0.3"
__author__    = "D. L. Wood"
__maintainer__ = "Jayce Dowell"

from . import test_lsl
