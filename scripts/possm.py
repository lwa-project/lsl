#!/usr/bin/env python3

"""
Utility script similar to the AIPS task 'possm' for plotting visibility data
stored in a file.
"""

import sys
import numpy as np
import argparse

from lsl import astro
from lsl.imaging import utils
from lsl.common.progress import ProgressBar
from lsl.misc import parser as aph

import matplotlib.pyplot as plt

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    # Grab the filename
    filename = args.filename
    
    idi = utils.CorrelatedData(filename)
    lo = idi.get_observer()
    lo.date = idi.date_obs.strftime("%Y/%m/%d %H:%M:%S")
    jd = lo.date + astro.DJD_OFFSET
    
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
        dataDict = idi.get_data_set(set, include_auto=args.include_auto, min_uv=args.uv_min)
        
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
        toWork = np.where( (freq >= args.freq_start) & (freq <= args.freq_stop) )[0]
        if len(toWork) == 0:
            raise RuntimeError("Cannot find data between %.2f and %.2f MHz" % (args.freq_start/1e6, args.freq_stop/1e6))
        if args.frequency:
            xValues = freq[toWork]/1e6
            xLabel = 'Frequency [MHz]'
        else:
            xValues = toWork
            xLabel = 'Channel'
            
        # Plot
        pols = dataDict.pols
        nBL = len(dataDict.baselines)
        args.nbaseline = min([args.nbaseline, nBL])
        nplot = min([25, args.nbaseline])
        nrow = int(round(numpy.sqrt(nplot)))
        while nplot % nrow != 0:
            nrow -= 1
        ncol = nplot // nrow
        print("    Plotting the first %i baselines" % args.nbaseline)
        pb = ProgressBar(max=args.nbaseline)
        i = 0
        for k in range((args.nbaseline-1)//25+1):
            fig = plt.figure()
            gs = fig.add_gridspec(2*nrow, ncol)

            for j in range(min([nplot, 25])):
                try:
                    stnd1, stnd2 = dataDict.baselines[i]
                    stnd1 = idi.stands[stnd1]
                    stnd2 = idi.stands[stnd2]
                    vis = getattr(dataDict, pols[0]).data[i]
                    i += 1
                except IndexError:
                    plt.draw()
                    break

                amp = np.abs(vis)
                if args.log:
                    amp = np.log10(amp)
                phs = np.angle(vis)*180/np.pi

                ax = fig.add_subplot(gs[j//ncol*2+0, j%ncol])
                if ((phs+360)%360).std() < phs.std():
                    ax.plot(xValues, (phs[toWork]+360)%360, linestyle=' ', marker='x')
                    ax.set_ylim([0, 360])
                else:
                    ax.plot(xValues, phs[toWork], linestyle=' ', marker='x')
                    ax.set_ylim([-180, 180])
                ax.set_title('%i - %i' % (stnd1, stnd2))
                if j % nrow == 0:
                    ax.set_ylabel('Phs')
                    
                ax = fig.add_subplot(gs[j//ncol*2+1, j%ncol])
                ax.plot(xValues, amp[toWork], linestyle=' ', marker='x', color='green')
                ax.set_title('%i - %i' % (stnd1, stnd2))
                if j % nrow == 0:
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
    parser.add_argument('-a', '--include-auto', action='store_true',
                        help='also plot auto-correlations')
    parser.add_argument('-i', '--include', type=aph.csv_int_list, default='all', 
                        help='comma seperated list of dipoles to include')
    parser.add_argument('-e', '--exclude', type=aph.csv_int_list, default='none', 
                        help='comma seperated list of dipoles to exclude')
    parser.add_argument('-f', '--frequency', action='store_true', 
                        help='label the channels in frequency rather than channel')
    parser.add_argument('-l', '--log', action='store_true', 
                        help='use a log scale for the visbility amplitudes')
    parser.add_argument('-n', '--nbaseline', type=int, default=50,
                        help='maximum number of baselines to plot')
    args = parser.parse_args()
    main(args)
