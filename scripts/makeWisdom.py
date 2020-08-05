#!/usr/bin/env python

"""
Utility to generate FFTW wisdom for use in LSL.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import sys
import argparse

from lsl.misc import wisdom

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    if args.show:
        wisdom.show()
    else:
        wisdom.make()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='build LSL-specific FFTW wisdom files', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-s', '--show', action='store_true', 
                        help='show information about the avaliable LSL wisdom')
    args = parser.parse_args()
    main(args)
    