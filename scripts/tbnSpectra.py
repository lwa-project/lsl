#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBN file, plot the time averaged spectra for each digitizer input.
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

from lsl.common import stations, metabundle, metabundleADP
from lsl.reader.ldp import LWA1DataFile
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.misc import parser as aph

import matplotlib.pyplot as plt

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
            try:
                station = metabundle.get_station(args.metadata, apply_sdm=True)
            except:
                station = metabundleADP.get_station(args.metadata, apply_sdm=True)
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    # Length of the FFT
    LFFT = args.fft_length
    
    idf = LWA1DataFile(args.filename)
    
    nFramesFile = idf.get_info('nframe')
    srate = idf.get_info('sample_rate')
    antpols = len(antennas)
    
    # Offset in frames for beampols beam/tuning/pol. sets
    args.skip = idf.offset(args.skip)
    
    # Make sure that the file chunk size contains is an integer multiple
    # of the FFT length so that no data gets dropped.  This needs to
    # take into account the number of antpols in the data, the FFT length,
    # and the number of samples per frame.
    maxFrames = int((2*260*750)/antpols*512/float(LFFT))*LFFT/512*antpols
    
    # Number of frames to integrate over
    nFrames = int(args.average * srate / 512 * antpols)
    nFrames = int(1.0 * nFrames / antpols*512/float(LFFT))*LFFT/512*antpols
    args.average = 1.0 * nFrames / antpols * 512 / srate
    
    # Number of remaining chunks
    nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    beginDate = ephem.Date(unix_to_utcjd(idf.get_info('start_time')) - DJD_OFFSET)
    central_freq = idf.get_info('freq1')
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz" % central_freq)
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, args.skip*srate*antpols/512))
    print("Integration: %.3f s (%i frames; %i frames per stand/pol)" % (args.average, nFrames, nFrames / antpols))
    print("Chunks: %i" % nChunks)
    
    # Sanity check
    if args.skip*srate*antpols/512 > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - args.skip*srate*antpols/512):
        raise RuntimeError("Requested integration time+offset is greater than file length")
        
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
    masterWeight = numpy.zeros((nChunks, antpols, LFFT))
    masterSpectra = numpy.zeros((nChunks, antpols, LFFT))
    
    for i in xrange(nChunks):
        print("Working on chunk #%i of %i" % (i+1, nChunks))
        
        try:
            readT, t, data = idf.read(args.average/nChunks)
        except Exception, e:
            print("Error: %s" % str(e))
            continue
            
        # Calculate the spectra for this block of data and then weight the results by 
        # the total number of frames read.  This is needed to keep the averages correct.
        
        freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=window, pfb=args.pfb, verbose=args.verbose, sample_rate=srate)
        for stand in xrange(tempSpec.shape[0]):
            masterSpectra[i,stand,:] = tempSpec[stand,:]
            masterWeight[i,stand,:] = int(readT*srate/LFFT)
            
    # Apply the cable loss corrections, if requested
    if False:
        for s in range(masterSpectra.shape[1]):
            currGain = antennas[s].cable.gain(freq)
            for c in range(masterSpectra.shape[0]):
                masterSpectra[c,s,:] /= currGain
                
    # Now that we have read through all of the chunks, perform the final averaging by
    # dividing by all of the chunks
    spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
    
    # Put the frequencies in the best units possible
    freq += central_freq
    freq, units = _best_freq_units(freq)
    
    # Deal with the `keep` options
    if args.keep == 'all':
        antpolsDisp = int(numpy.ceil(antpols/20))
        js = [i for i in xrange(antpols)]
    else:
        antpolsDisp = int(numpy.ceil(len(args.keep)*2/20))
        if antpolsDisp < 1:
            antpolsDisp = 1
            
        js = []
        for k in args.keep:
            for i,ant in enumerate(antennas):
                if ant.stand.id == k:
                    js.append(i)
                    
    nPlot = len(js)
    if nPlot < 20:
        if nPlot % 4 == 0 and nPlot != 4:
            figsY = 4
        else:
            figsY = 2
        figsX = int(numpy.ceil(1.0*nPlot/figsY))
    else:
        figsY = 4
        figsX = 5
    figsN = figsX*figsY
    for i in xrange(antpolsDisp):
        # Normal plotting
        fig = plt.figure()
        for k in xrange(i*figsN, i*figsN+figsN):
            try:
                j = js[k]
                currSpectra = numpy.squeeze( numpy.log10(spec[j,:])*10.0 )
            except IndexError:
                break
            ax = fig.add_subplot(figsX, figsY, (k%figsN)+1)
            ax.plot(freq, currSpectra, label='Stand: %i, Pol: %i (Dig: %i)' % (antennas[j].stand.id, antennas[j].pol, antennas[j].digitizer))
            
            # If there is more than one chunk, plot the difference between the global 
            # average and each chunk
            if nChunks > 1 and not args.disable_chunks:
                for k in xrange(nChunks):
                    # Some files are padded by zeros at the end and, thus, carry no 
                    # weight in the average spectra.  Skip over those.
                    if masterWeight[k,j,:].sum() == 0:
                        continue
                        
                    # Calculate the difference between the spectra and plot
                    subspectra = numpy.squeeze( numpy.log10(masterSpectra[k,j,:])*10.0 )
                    diff = subspectra - currSpectra
                    ax.plot(freq, diff)
                    
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
            description='read in a TBN file and create a collection of time-averaged spectra', 
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
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=10.0, 
                        help='number of seconds of data to average for spectra')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=4096, 
                        help='set FFT length')
    parser.add_argument('-d', '--disable-chunks', action='store_true', 
                        help='disable plotting chunks in addition to the global average')
    parser.add_argument('-k', '--keep', type=aph.csv_int_list, default='all', 
                        help='only display the following comma-seperated list of stands')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for spectra image')
    args = parser.parse_args()
    main(args)
    
