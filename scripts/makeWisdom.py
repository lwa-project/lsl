#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from lsl.misc import wisdom


def main(args):
    do_fftw = True if not args.pyfftw_only else False
    do_pyfftw = True if not args.fftw_only else False
    
    if args.show:
        wisdom.show(FFTW=do_fftw, PyFFTW=do_pyfftw)
    else:
        wisdom.make(FFTW=do_fftw, PyFFTW=do_pyfftw)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='build LSL-specific FFTW and PyFFTW wisdom files', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-f', '--fftw-only', action='store_true', 
                        help='build/show only FFTW wisdom only')
    parser.add_argument('-p', '--pyfftw-only', action='store_true', 
                        help='build/show only PyFFTW wisdom only')
    parser.add_argument('-s', '--show', action='store_true', 
                        help='show information about the avaliable LSL wisdom')
    args = parser.parse_args()
    main(args)
    
