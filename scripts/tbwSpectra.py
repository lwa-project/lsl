#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBW file, plot the time averaged spectra for each digitizer input.
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import math
import numpy
import ephem
import argparse

from lsl.common import stations, metabundle
from lsl.reader.ldp import LWA1DataFile
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.misc import parser as aph

import matplotlib.pyplot as plt

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    # Setup the LWA station information
    if args.metadata is not None:
        try:
            station = stations.parse_ssmif(args.metadata)
        except ValueError:
            station = metabundle.get_station(args.metadata, apply_sdm=True)
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    # Length of the FFT
    LFFT = args.fft_length
    
    # Make sure that the file chunk size contains is an integer multiple
    # of the FFT length so that no data gets dropped
    maxFrames = int((30000*260)/float(LFFT))*LFFT
    # It seems like that would be a good idea, however...  TBW data comes one
    # capture at a time so doing something like this actually truncates data 
    # from the last set of stands for the first integration.  So, we really 
    # should stick with
    maxFrames = (30000*260)
    
    idf = LWA1DataFile(args.filename)
    
    nFrames = idf.get_info('nframe')
    srate = idf.get_info('sample_rate')
    dataBits = idf.get_info('data_bits')
    # The number of ant/pols in the file is hard coded because I cannot figure out 
    # a way to get this number in a systematic fashion
    antpols = len(antennas)
    nChunks = int(math.ceil(1.0*nFrames/maxFrames))
    if dataBits == 12:
        nSamples = 400
    else:
        nSamples = 1200
        
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    beginDate = ephem.Date(unix_to_utcjd(idf.get_info('start_time')) - DJD_OFFSET)
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Length: %i-bit" % dataBits)
    print("Frames: %i" % nFrames)
    print("Chunks: %i" % nChunks)
    print("===")
    
    # Setup the window function to use
    if args.bartlett:
        window = numpy.bartlett
    elif args.blackman:
        window = numpy.blackman
    elif args.hanning:
        window = numpy.hanning
    else:
        window = fxc.null_window
        
    # Master loop over all of the file chunks
    nChunks = 1
    masterSpectra = numpy.zeros((nChunks, antpols, LFFT))
    masterWeight = numpy.zeros((nChunks, antpols, LFFT))
    
    readT, t, data = idf.read(0.061)
    
    # Calculate the spectra for this block of data and then weight the results by 
    # the total number of frames read.  This is needed to keep the averages correct.
    # NB:  The weighting is the same for the x and y polarizations because of how 
    # the data are packed in TBW
    freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=window, pfb=args.pfb, verbose=args.verbose)
    for stand in xrange(masterSpectra.shape[1]):
        masterSpectra[0,stand,:] = tempSpec[stand,:]
        masterWeight[0,stand,:] = int(readT*srate/LFFT)
        
    # We don't really need the data array anymore, so delete it
    del(data)
    
    # Apply the cable loss corrections, if requested
    if args.gain_correct:
        for s in xrange(masterSpectra.shape[1]):
            currGain = antennas[s].cable.gain(freq)
            for c in xrange(masterSpectra.shape[0]):
                masterSpectra[c,s,:] /= currGain
                
    # Now that we have read through all of the chunks, perform the final averaging by
    # dividing by all of the chunks
    spec = masterSpectra.mean(axis=0)
    
    # The plots:  This is setup for the current configuration of 20 antpols
    if args.gain_correct & args.stack:
        # Stacked spectra - only if cable loss corrections are to be applied
        colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 
            'purple', 'salmon', 'olive', 'maroon', 'saddlebrown', 'yellowgreen', 
            'teal', 'steelblue', 'seagreen', 'slategray', 'mediumorchid', 'lime', 
            'dodgerblue', 'darkorange']
            
        for f in xrange(int(numpy.ceil(antpols/20.))):
            fig = plt.figure()
            ax1 = fig.add_subplot(1, 1, 1)
            for i in xrange(f*20, f*20+20):
                currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
                ax1.plot(freq/1e6, currSpectra, label='%i,%i' % (antennas[i].stand.id, antennas[i].pol), color=colors[i % 20])
                
            ax1.set_xlabel('Frequency [MHz]')
            ax1.set_ylabel('P.S.D. [dB/RBW]')
            ax1.set_xlim([20,88])
            #ax1.set_ylim([10,90])
            leg = ax1.legend(loc=0, ncol=3)
            for l in leg.get_lines():
                l.set_linewidth(1.7)  # the legend line width
    else:
        for f in xrange(int(numpy.ceil(antpols/20))):
            # Normal plotting
            fig = plt.figure()
            figsY = 4
            figsX = 5
            fig.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.20, hspace=0.50)
            for i in xrange(f*20, f*20+20):
                ax = fig.add_subplot(figsX, figsY, (i%20)+1)
                try:
                    currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
                except IndexError:
                    break
                ax.plot(freq/1e6, currSpectra, label='Stand: %i, Pol: %i (Dig: %i)' % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer))
                
                # If there is more than one chunk, plot the difference between the global 
                # average and each chunk
                if nChunks > 1:
                    for j in xrange(nChunks):
                        # Some files are padded by zeros at the end and, thus, carry no 
                        # weight in the average spectra.  Skip over those.
                        if masterWeight[j,i,:].sum() == 0:
                            continue
                            
                        # Calculate the difference between the spectra and plot
                        subspectra = numpy.squeeze( numpy.log10(masterSpectra[j,i,:])*10.0 )
                        diff = subspectra - currSpectra
                        ax.plot(freq/1e6, diff)
                        
                ax.set_title('Stand: %i (%i); Dig: %i [%i]' % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer, antennas[i].combined_status))
                ax.set_xlabel('Frequency [MHz]')
                ax.set_ylabel('P.S.D. [dB/RBW]')
                ax.set_xlim([10,90])
                ax.set_ylim([10,80])
                
            # Save spectra image if requested
            if args.output is not None:
                base, ext = os.path.splitext(args.output)
                outFigure = "%s-%02i%s" % (base, f+1, ext)
                fig.savefig(outFigure)
                
        plt.draw()
        
    print("RBW: %.1f Hz" % (freq[1]-freq[0]))
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='read in a TBW file and create a collection of time-averaged spectra', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of the SSMIF or metadata tarball file to use for mappings')
    wgroup = parser.add_mutually_exclusive_group(required=False)
    wgroup.add_argument('-t', '--bartlett', action='store_true', 
                        help='apply a Bartlett window to the data')
    wgroup.add_argument('-b', '--blackman', action='store_true', 
                        help='apply a Blackman window to the data')
    wgroup.add_argument('-n', '--hanning', action='store_true', 
                        help='apply a Hanning window to the data')
    wgroup.add_argument('-p', '--pfb', action='store_true', 
                        help='enabled the PFB on the F-engine')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=4096, 
                        help='set FFT length')
    parser.add_argument('-g', '--gain-correct', action='store_true', 
                        help='correct signals for the cable losses')
    parser.add_argument('-s', '--stack', action='store_true', 
                        help="stack spectra in groups of 6 (if '-g' is enabled only)")
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for spectra image')
    args = parser.parse_args()
    main(args)

