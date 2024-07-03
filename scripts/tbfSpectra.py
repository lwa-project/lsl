#!/usr/bin/env python3

"""
Given a TBF file, plot the time averaged spectra for each digitizer input.
"""
    
import os
import sys
import math
import numpy as np
import argparse

from lsl.common import stations, metabundle
from lsl.reader.ldp import LWADataFile, TBFFile
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
    # Setup the LWA station information
    if args.metadata is not None:
        try:
            station = stations.parse_ssmif(args.metadata)
        except ValueError:
            station = metabundle.get_station(args.metadata, apply_sdm=True)
    elif args.lwana:
        station = stations.lwana
    else:
        station = stations.lwasv
    antennas = station.antennas
    
    idf = LWADataFile(args.filename)
    if not isinstance(idf, TBFFile):
        raise RuntimeError("File '%s' does not appear to be a valid TBF file" % os.path.basename(filename))
        
    nFrames = idf.get_info('nframe')
    antpols = len(antennas)
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    beginDate = idf.get_info('start_time').datetime
    
    # Make sure the TBF stand count is consistent with how many antennas we have
    if antpols//2 != junkFrame.nstand:
        raise RuntimeError("Number of stands in the station (%i) does not match what is in the data (%i)" % (antpols//2, junkFrame.nstand))
        
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nchannels = idf.get_info('nchan')
    
    # Figure out how many chunks we need to work with
    nChunks = 1
    
    # Calculate the frequencies
    freq = idf.get_into('freq1')
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Channel Count: %i" % nchannels)
    print("Frames: %i" % nFrames)
    print("===")
    print("Chunks: %i" % nChunks)
    
    for i in range(nChunks):
        print("Working on chunk #%i of %i" % (i+1, nChunks))
        
        try:
            readT, t, data = idf.read()
        except Exception as e:
            print("Error: %s" % str(e))
            continue
            
        # Detect power and integrate
        data = np.abs(data)**2
        spec = data.mean(axis=2)
        
    # Apply the cable loss corrections, if requested
    if args.gain_correct:
        for s in range(spec.shape[0]):
            currGain = antennas[s].cable.gain(freq)
            spec[s,:] /= currGain
            
    # Put the frequencies in the best units possible
    freq, units = _best_freq_units(freq)
    
    # Deal with the `keep` options
    if args.keep == 'all':
        antpolsDisp = int(np.ceil(antpols/20))
        js = [i for i in range(antpols)]
    else:
        antpolsDisp = int(np.ceil(len(args.keep)*2/20))
        if antpolsDisp < 1:
            antpolsDisp = 1
            
        js = []
        for k in args.keep:
            for i,ant in enumerate(antennas):
                if ant.stand.id == k:
                    js.append(i)
                    
    nPlot = len(js)
    if nPlot < 16:
        if nPlot % 4 == 0 and nPlot != 4:
            figsY = 4
        else:
            figsY = 2
        figsX = int(np.ceil(1.0*nPlot/figsY))
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
                currSpectra = np.squeeze( np.log10(spec[j,:])*10.0 )
            except IndexError:
                break
            ax = fig.add_subplot(figsX, figsY, (k%figsN)+1)
            ax.plot(freq, currSpectra, label='Stand: %i, Pol: %i (Dig: %i)' % (antennas[j].stand.id, antennas[j].pol, antennas[j].digitizer))
            
            ax.set_title('Stand: %i (%i); Dig: %i [%i]' % (antennas[j].stand.id, antennas[j].pol, antennas[j].digitizer, antennas[j].combined_status))
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
            description='read in a TBF file and create a collection of time-averaged spectra', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of the SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-n', '--lwana', action='store_true', 
                        help='use LWA-NA instead of LWA-SV')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-g', '--gain-correct', action='store_true',
                        help='correct signals for the cable losses')
    parser.add_argument('-k', '--keep', type=aph.csv_int_list, default='all', 
                        help='only display the following comma-seperated list of stands')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for spectra image')
    args = parser.parse_args()
    main(args)
    
