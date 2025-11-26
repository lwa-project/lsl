#!/usr/bin/env python3

"""
Given a DR spectrometer file, plot the time averaged spectra for each 
polarization product.
"""

import os
import sys
import math
import numpy as np
import argparse

from lsl.reader.ldp import LWADataFile, DRSpecFile
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
    idf = LWADataFile(args.filename)
    if not isinstance(idf, DRSpecFile):
        raise RuntimeError("File '%s' does not appear to be a valid DR spectrometer file" % os.path.basename(args.filename))
        
    # Basic file information
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
    
    # Date & Central Frequency
    beginDate = idf.get_info('start_time').datetime
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    freq = np.fft.fftfreq(LFFT, d=1.0/srate)
    freq = np.fft.fftshift(freq)
    
    # File summary
    print(f"Filename: {args.filename}")
    print(f"Date of First Frame: {str(beginDate)}")
    print(f"Beam: {beam}")
    print(f"Tune/Pols: {beampols}")
    print(f"Sample Rate: {srate} Hz")
    print(f"Tuning Frequency: {central_freq1:.3f} Hz (1); {central_freq2:.3f} Hz (2)")
    print(f"Frames: {nFramesFile} ({nFramesFile*tInt:.3f} s)")
    print("---")
    print(f"Transform Length: {LFFT} channels")
    print(f"Integration Time: {tInt:.3f} s")
    print("---")
    print(f"Offset: {args.skip:.3f} s ({args.skip*srate*beampols/4096} frames)")
    print(f"Integration: {args.average:.3f} s ({nFrames} frames; {nFrames} frames per beam/tune/pol)")
    print(f"Chunks: {nChunks}")
    
    # Sanity check
    if args.skip/tInt > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - args.skip/tInt):
        raise RuntimeError("Requested integration time+offset is greater than file length")
        
    # Master loop over all of the file chunks
    masterWeight = np.zeros((nChunks, 2*len(products), LFFT))
    masterSpectra = np.zeros((nChunks, 2*len(products), LFFT))
    for i in range(nChunks):
        print(f"Working on chunk #{i+1} of {nChunks}")
        
        try:
            readT, t, data = idf.read(args.average/nChunks)
        except Exception as e:
            print(f"Error: {str(e)}")
            continue
            
        ## Integrate up the chunk
        data = data.mean(axis=1)
        
        ## Save
        for stand in range(data.shape[0]):
            masterSpectra[i,stand,:] = data[stand,:]
            masterWeight[i,stand,:] = int(readT*srate/LFFT)
            
        ## We don't really need the data array anymore, so delete it
        del(data)
        
    # Now that we have read through all of the chunks, perform the final averaging by
    # dividing by all of the chunks
    spec = np.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
    
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
        currSpectra = np.squeeze( np.log10(spec[i,:])*10.0 )
        ax.plot(freq, currSpectra, label=f"{i+1} (avg)")
        
        # If there is more than one chunk, plot the difference between the global 
        # average and each chunk
        if nChunks > 1 and args.disable_chunks:
            for j in range(nChunks):
                # Some files are padded by zeros at the end and, thus, carry no 
                # weight in the average spectra.  Skip over those.
                if masterWeight[j,i,:].sum() == 0:
                    continue
                    
                # Calculate the difference between the spectra and plot
                subspectra = np.squeeze( np.log10(masterSpectra[j,i,:])*10.0 )
                diff = subspectra - currSpectra
                ax.plot(freq, diff, label='%i' % j)
                
        ax.set_title(f"Beam {beam}, Tune. {i//len(products)}, {products[i % len(products)]}")
        ax.set_xlabel(f"Frequency [{units}]")
        ax.set_ylabel('P.S.D. [dB/RBW]')
        ax.set_xlim([freq.min(), freq.max()])
        ax.legend(loc=0)
        
        print(f"For beam {beam}, tune. {i//len(products)}, {products[i % len(products)]} maximum in PSD at {freq[np.where(spec[i,:]==spec[i,:].max())][0]:.3f} {units}")
        
    print(f"RBW: {freq[1]-freq[0]:.4f} {units}")
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
