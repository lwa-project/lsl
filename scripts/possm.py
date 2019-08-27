#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utility script similar to the AIPS task 'possm' for plotting visibility data
stored in a file.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import sys
import numpy
import argparse

from lsl import astro
from lsl.imaging import utils
from lsl.common.progress import ProgressBar
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


def main(args):
    # Grab the filename
    filename = args.filename
    
    idi = utils.CorrelatedData(filename)
    aa = idi.get_antennaarray()
    lo = idi.get_observer()
    lo.date = idi.date_obs.strftime("%Y/%m/%d %H:%M:%S")
    jd = lo.date + astro.DJD_OFFSET
    lst = str(lo.sidereal_time())   # pylint:disable=no-member

    nStand = len(idi.stands)
    nchan = len(idi.freq)
    freq = idi.freq
    
    print("Raw Stand Count: %i" % nStand)
    print("Final Baseline Count: %i" % (nStand*(nStand-1)/2,))
    print("Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nchan, (freq[1] - freq[0])/1e3))
    print("Polarization Products: %i starting with %i" % (len(idi.pols), idi.pols[0]))
    print("JD: %.3f" % jd)

    print("Reading in FITS IDI data")
    nSets = idi.integration_count
    for set in range(1, nSets+1):
        if args.dataset != -1 and args.dataset != set:
            continue
            
        print("Set #%i of %i" % (set, nSets))
        dataDict = idi.get_data_set(set, min_uv=args.uv_min)
        
        # Prune out what needs to go
        if args.include != 'all' or args.exclude != 'none':
            print("    Processing include/exclude lists")
            dataDict = dataDict.get_antenna_subset(include=args.include, 
                                                   exclude=args.exclude, 
                                                   indicies=False)
            
            ## Report
            for pol in dataDict.pols:
                print("        %s now has %i baselines" % (pol, len(dataDict.baselines)))
                
        # Pull out the right channels
        toWork = numpy.where( (freq >= args.freq_start) & (freq <= args.freq_stop) )[0]
        if len(toWork) == 0:
            raise RuntimeError("Cannot find data between %.2f and %.2f MHz" % (args.freq_start/1e6, args.freq_stop/1e6))
        if args.frequency:
            xValues = freq[toWork]/1e6
            xLabel = 'Frequency [MHz]'
        else:
            xValues = toWork
            xLabel = 'Channel'
            
        # Plot
        print("    Plotting the first 50 baselines")
        pols = dataDict['jd'].keys()
        nBL = len(dataDict['bls'][pols[0]])
        pb = ProgressBar(max=nBL)
        i = 0
        for k in range(2):
            fig = plt.figure()

            for j in range(25):
                try:
                    stnd1, stnd2 = dataDict['bls'][pols[0]][i]
                    stnd1 = idi.stands[stnd1]
                    stnd2 = idi.stands[stnd2]
                    vis = dataDict['vis'][pols[0]][i]
                    i += 1
                except IndexError:
                    plt.draw()
                    break

                amp = numpy.abs(vis)
                if args.log:
                    amp = numpy.log10(amp)
                phs = numpy.angle(vis)*180/numpy.pi

                ax = fig.add_subplot(10, 5, 2*(j//5)*5+j%5+1)
                if ((phs+360)%360).std() < phs.std():
                    ax.plot(xValues, (phs[toWork]+360)%360, linestyle=' ', marker='x')
                    ax.set_ylim([0, 360])
                else:
                    ax.plot(xValues, phs[toWork], linestyle=' ', marker='x')
                    ax.set_ylim([-180, 180])
                ax.set_title('%i - %i' % (stnd1, stnd2))
                if j % 5 == 0:
                    ax.set_ylabel('Phs')
                    
                ax = fig.add_subplot(10, 5, 2*(j//5)*5+j%5+1+5)
                ax.plot(xValues, amp[toWork], linestyle=' ', marker='x', color='green')
                ax.set_title('%i - %i' % (stnd1, stnd2))
                if j % 5 == 0:
                    ax.set_xlabel(xLabel)
                    ax.set_ylabel('%sAmp' % '' 'Log ' if args.log else '')
                    
                pb.inc(amount=1)
                if pb.amount != 0 and pb.amount % 10 == 0:
                    sys.stdout.write(pb.show()+'\r')
                    sys.stdout.flush()
            plt.draw()

        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        plt.show()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='create AIPS POSSM style plots of a FITS IDI file', 
        epilog='NOTE:  If both -i/--include and -e/--exclude are specified the include list has priority.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to plot')
    parser.add_argument('-1', '--freq-start', type=aph.frequency, default='10.0', 
                        help='first frequency to analyze in MHz')
    parser.add_argument('-2', '--freq-stop', type=aph.frequency, default='88.0', 
                        help='last frequency to analyze in MHz')
    parser.add_argument('-s', '--dataset', type=int, default=-1, 
                        help='data set to image')
    parser.add_argument('-m', '--uv-min', type=float, default=0.0, 
                        help='minimun baseline uvw length to include in lambda at the midpoint frequency')
    parser.add_argument('-i', '--include', type=aph.csv_int_list, default='all', 
                        help='comma seperated list of dipoles to include')
    parser.add_argument('-e', '--exclude', type=aph.csv_int_list, default='none', 
                        help='comma seperated list of dipoles to exclude')
    parser.add_argument('-f', '--frequency', action='store_true', 
                        help='label the channels in frequency rather than channel')
    parser.add_argument('-l', '--log', action='store_true', 
                        help='use a log scale for the visbility amplitudes')
    args = parser.parse_args()
    main(args)
