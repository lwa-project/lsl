#!/usr/bin/env python

"""
Given a OVRO-LWA triggered voltage buffer dump file, plot the time averaged
spectra for each digitizer input.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import math
import numpy
import argparse

from lsl.reader import ovro, errors
from lsl.common import stations, metabundleADP
from lsl.misc import parser as aph

from matplotlib import pyplot as plt

from lsl.misc import telemetry
telemetry.track_script()


def _best_freq_units(freq):
    """Given a numpy array of frequencies in Hz, return a new array with the
    frequencies in the best units possible (kHz, MHz, etc.)."""
    
    # Figure out how large the data are
    scale = int(math.log10(freq.max()))
    if scale >= 9:
        divis = 1e9
        units = 'GHz'
    elif scale >= 6:
        divis = 1e6
        units = 'MHz'
    elif scale >= 3:
        divis = 1e3
        units = 'kHz'
    else:
        divis = 1
        units = 'Hz'
        
    # Convert the frequency
    newFreq = freq / divis
    
    # Return units and freq
    return (newFreq, units)


def main(args):
    fh = open(args.filename, 'rb')
    
    # Read in the frame and pull out the start time as well as the number of stand/pols.
    frame = ovro.read_frame(fh)
    beginDate = frame.time
    antpols = frame.payload.data.shape[-2]*frame.payload.data.shape[-1]
    
    # Make a dummy list of antennas
    antennas = list(range(antpols))
    
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nFrames = nFramesPerObs = 1
    nchannels = frame.header.nchan
    nSamples = frame.payload.data.shape[0]
    
    # Calculate the frequencies
    freq = frame.channel_freqs
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Frames per Observation: %i" % nFramesPerObs)
    print("Channel Count: %i" % nchannels)
    print("Frames: %i" % nFrames)
    print("===")
    
    # Convert to spectra
    spec = numpy.abs(frame.payload.data)**2
    spec = spec.mean(axis=0)
    
    fh.close()
    
    # Reshape and transpose to get it in to a "normal" order
    spec.shape = (spec.shape[0], spec.shape[1]*spec.shape[2])
    spec = spec.T
    
    # Put the frequencies in the best units possible
    freq, units = _best_freq_units(freq)
    
    # Deal with the `keep` options
    if args.keep == 'all':
        antpolsDisp = int(numpy.ceil(antpols/20))
        js = [i for i in range(antpols)]
    else:
        antpolsDisp = int(numpy.ceil(len(args.keep)*2/20))
        if antpolsDisp < 1:
            antpolsDisp = 1
            
        js = []
        for k in args.keep:
            for i,ant in enumerate(antennas):
                if ant == k:
                    js.append(i)
                    
    nPlot = len(js)
    if nPlot < 16:
        if nPlot % 4 == 0 and nPlot != 4:
            figsY = 4
        else:
            figsY = 2
        figsX = int(numpy.ceil(1.0*nPlot/figsY))
    else:
        figsY = 4
        figsX = 4
    figsN = figsX*figsY
    for i in range(antpolsDisp):
        # Normal plotting
        fig = plt.figure()
        for k in range(i*figsN, i*figsN+figsN):
            try:
                j = js[k]
                currSpectra = numpy.squeeze( numpy.log10(spec[j,:])*10.0 )
            except IndexError:
                break
            ax = fig.add_subplot(figsX, figsY, (k%figsN)+1)
            ax.plot(freq, currSpectra, label='Dig: %i' % (antennas[j]))
            
            ax.set_title('Dig: %i' % (antennas[j],))
            ax.set_xlabel('Frequency [%s]' % units)
            ax.set_ylabel('P.S.D. [dB/RBW]')
            ax.set_ylim([-10, 30])
            
        # Save spectra image if requested
        if args.output is not None:
            base, ext = os.path.splitext(args.output)
            outFigure = "%s-%02i%s" % (base, i+1, ext)
            fig.savefig(outFigure)
            
        plt.draw()
        
    print("RBW: %.4f %s" % ((freq[1]-freq[0]), units))
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='read in a OVRO-LWA triggered voltage buffer dump file and create a collection of time-averaged spectra', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-k', '--keep', type=aph.csv_int_list, default='all', 
                        help='only display the following comma-seperated list of stands')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for spectra image')
    args = parser.parse_args()
    main(args)
    
