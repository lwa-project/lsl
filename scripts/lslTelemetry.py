#!/usr/bin/env python

"""
Utility for updating/changing the LSL telemetry setting.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import argparse

from lsl.misc import telemetry
telemetry.track_script()


if __name__ == "__main__":
    os.system("%s -m lsl.misc.telemetry %s" % (sys.executable, " ".join(sys.argv[1:])))
