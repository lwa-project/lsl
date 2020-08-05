#!/usr/bin/env python

"""
Given a DR spectrometer file, plot the time averaged spectra for each 
polarization product.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import sys
import math
import numpy
import ephem
import argparse

from lsl.reader.ldp import LWA1DataFile
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
    idf = LWA1DataFile(args.filename)
    
    # Basic file informaiton
    nFramesFile = idf.get_info('nframe')
    srate = idf.get_info('sample_rate')
    beam = idf.get_info('beam')
    beampols = idf.get_info('nbeampol')
    tInt = idf.get_info('tint')
    LFFT = idf.get_info('LFFT')
    products = idf.get_info('data_products')
    
    # Offset in frames for beampols beam/tuning/pol. sets
    args.skip = idf.offset(args.skip)
    
    # Number of frames to integrate over
    maxFrames = 10000
    nFrames = int(args.average / tInt)
    args.average = nFrames * tInt
    
    # Number of remaining chunks
    maxFramesTime = maxFrames*tInt
    nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))
    
    # Date & Central Frequnecy
    beginDate = ephem.Date(unix_to_utcjd(idf.get_info('start_time')) - DJD_OFFSET)
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    freq = numpy.fft.fftfreq(LFFT, d=1.0/srate)
    freq = numpy.fft.fftshift(freq)
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Beam: %i" % beam)
    print("Tune/Pols: %i" % beampols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (central_freq1, central_freq2))
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile*tInt))
    print("---")
    print("Transform Length: %i channels" % LFFT)
    print("Integration Time: %.3f s" % tInt)
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, args.skip*srate*beampols/4096))
    print("Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (args.average, nFrames, nFrames))
    print("Chunks: %i" % nChunks)
    
    # Sanity check
    if args.skip/tInt > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - args.skip/tInt):
        raise RuntimeError("Requested integration time+offset is greater than file length")
        
    # Master loop over all of the file chunks
    masterWeight = numpy.zeros((nChunks, 2*len(products), LFFT))
    masterSpectra = numpy.zeros((nChunks, 2*len(products), LFFT))
    for i in range(nChunks):
        print("Working on chunk #%i of %i" % (i+1, nChunks))
        
        try:
            readT, t, data = idf.read(args.average/nChunks)
        except Exception as e:
            print("Error: %s" % str(e))
            continue
            
        ## Integrate up the chunck
        data = data.mean(axis=1)
        
        ## Save
        for stand in range(data.shape[0]):
            masterSpectra[i,stand,:] = data[stand,:]
            masterWeight[i,stand,:] = int(readT*srate/LFFT)
            
        ## We don't really need the data array anymore, so delete it
        del(data)
        
    # Now that we have read through all of the chunks, perform the final averaging by
    # dividing by all of the chunks
    spec = numpy.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
    
    # Frequencies
    freq1 = freq + central_freq1
    freq2 = freq + central_freq2
    
    # The plots:  This is setup for the current configuration of 20 beampols
    fig = plt.figure()
    figsX = int(round(math.sqrt(2*len(products))))
    figsY = 2*len(products) // figsX
    # Put the frequencies in the best units possible
    freq1, units1 = _best_freq_units(freq1)
    freq2, units2 = _best_freq_units(freq2)
    
    for i in range(masterSpectra.shape[1]):
        if i/len(products) == 0:
            freq = freq1
            units = units1
        else:
            freq = freq2
            units = units2
            
        ax = fig.add_subplot(figsX,figsY,i+1)
        currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
        ax.plot(freq, currSpectra, label='%i (avg)' % (i+1))
        
        # If there is more than one chunk, plot the difference between the global 
        # average and each chunk
        if nChunks > 1 and args.disable_chunks:
            for j in range(nChunks):
                # Some files are padded by zeros at the end and, thus, carry no 
                # weight in the average spectra.  Skip over those.
                if masterWeight[j,i,:].sum() == 0:
                    continue
                    
                # Calculate the difference between the spectra and plot
                subspectra = numpy.squeeze( numpy.log10(masterSpectra[j,i,:])*10.0 )
                diff = subspectra - currSpectra
                ax.plot(freq, diff, label='%i' % j)
                
        ax.set_title('Beam %i, Tune. %i, %s' % (beam, i//len(products), products[i % len(products)]))
        ax.set_xlabel('Frequency [%s]' % units)
        ax.set_ylabel('P.S.D. [dB/RBW]')
        ax.set_xlim([freq.min(), freq.max()])
        ax.legend(loc=0)
        
        print("For beam %i, tune. %i, %s maximum in PSD at %.3f %s" % (beam, i//len(products), products[i % len(products)], freq[numpy.where( spec[i,:] == spec[i,:].max() )][0], units))
        
    print("RBW: %.4f %s" % ((freq[1]-freq[0]), units))
    plt.subplots_adjust(hspace=0.35, wspace=0.30)
    plt.show()
    
    # Save spectra image if requested
    if args.output is not None:
        fig.savefig(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DR spectrometer files and create a collection of time-averaged spectra', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                            help='filename to process')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=10.0, 
                        help='number of seconds of data to average for spectra')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-d', '--disable-chunks', action='store_true', 
                        help='disable plotting chunks in addition to the global average')
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for spectra image')
    args = parser.parse_args()
    main(args)
