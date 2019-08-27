# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import absolute_import

"""
Modules defining package tests.
"""

__revision__  = "$Rev$"
__version__   = "0.3"
__author__    = "D. L. Wood"
__maintainer__ = "Jayce Dowell"

from . import test_lsl

from lsl.misc import telemetry
telemetry.track_test_suite()

