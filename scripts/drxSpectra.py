#!/usr/bin/env python3

"""
Given a DRX file, plot the time averaged spectra for each beam output.
"""

import os
import sys
import math
import numpy as np
import argparse

import lsl.correlator.fx as fxc
from lsl.reader.ldp import LWADataFile, DRXFile, DRX8File
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
    # Length of the FFT
    LFFT = args.fft_length
    
    idf = LWADataFile(args.filename)
    if not isinstance(idf, (DRXFile, DRX8File)):
        raise RuntimeError(f"File '{os.path.basename(args.filename)}' does not appear to be a valid DRX file")
        
    # Offset in frames for beampols beam/tuning/pol. sets
    args.skip = idf.offset(args.skip)
    
    nFramesFile = idf.get_info('nframe')
    srate = idf.get_info('sample_rate')
    beampols = idf.get_info('nbeampol')
    
    # Make sure that the file chunk size contains is an integer multiple
    # of the FFT length so that no data gets dropped.  This needs to
    # take into account the number of beampols in the data, the FFT length,
    # and the number of samples per frame.
    maxFrames = int(1.0*19144/beampols*4096/float(LFFT))*LFFT/4096*beampols
    
    # Number of frames to integrate over
    nFrames = int(args.average * srate / 4096) * beampols
    nFrames = int(1.0 * (nFrames // beampols)*4096/float(LFFT))*LFFT/4096 * beampols
    args.average = 1.0 * (nFrames // beampols) * 4096 / srate
    
    # Number of remaining chunks
    nChunks = int(math.ceil(1.0*(nFrames)/maxFrames))
    
    # Date & Central Frequency
    beginDate = idf.get_info('start_time').datetime
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    beam = idf.get_info('beam')
    
    # File summary
    print(f"Filename: {args.filename}")
    print(f"Date of First Frame: {str(beginDate)}")
    print(f"Beam: {beam}")
    print(f"Tune/Pols: {beampols}")
    print(f"Sample Rate: {srate} Hz")
    print(f"Bit depth: {'8' if isinstance(idf, DRXFile) else '16'}")
    print(f"Tuning Frequency: {central_freq1:.3f} Hz (1); {central_freq2:.3f} Hz (2)")
    print(f"Frames: {nFramesFile} ({nFramesFile/beampols*4096/srate:.3f} s)")
    print("---")
    print(f"Offset: {args.skip:.3f} s ({args.skip*srate*beampols/4096} frames)")
    print(f"Integration: {args.average:.3f} s ({nFrames} frames; {nFrames//beampols} frames per beam/tune/pol)")
    print(f"Chunks: {nChunks}")
    
    # Sanity check
    if args.skip*srate*beampols/4096 > nFramesFile:
        raise RuntimeError("Requested offset is greater than file length")
    if nFrames > (nFramesFile - args.skip*srate*beampols/4096):
        raise RuntimeError("Requested integration time+offset is greater than file length")
        
    # Setup the window function to use
    if args.bartlett:
        window = np.bartlett
    elif args.blackman:
        window = np.blackman
    elif args.hanning:
        window = np.hanning
    else:
        window = fxc.null_window
        
    # Master loop over all of the file chunks
    standMapper = [4*(beam-1) + i for i in range(4)]
    masterWeight = np.zeros((nChunks, 4, LFFT))
    masterSpectra = np.zeros((nChunks, 4, LFFT))
    for i in range(nChunks):
        print(f"Working on chunk #{i+1} of {nChunks}")
        
        try:
            readT, t, data = idf.read(args.average/nChunks, return_ci8=True)
        except Exception as e:
            print(f"Error: {str(e)}")
            continue
            
        # Calculate the spectra for this block of data and then weight the results by 
        # the total number of frames read.  This is needed to keep the averages correct.
        freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=window, pfb=args.pfb, verbose=args.verbose, sample_rate=srate, clip_level=0)
        for stand in range(tempSpec.shape[0]):
            masterSpectra[i,stand,:] = tempSpec[stand,:]
            masterWeight[i,stand,:] = int(readT*srate/LFFT)
            
        # We don't really need the data array anymore, so delete it
        del(data)
        
    # Now that we have read through all of the chunks, perform the final averaging by
    # dividing by all of the chunks
    spec = np.squeeze( (masterWeight*masterSpectra).sum(axis=0) / masterWeight.sum(axis=0) )
    
    # Frequencies
    freq1 = freq + central_freq1
    freq2 = freq + central_freq2
    
    # The plots:  This is setup for the current configuration of 20 beampols
    fig = plt.figure()
    figsX = int(round(math.sqrt(4)))
    figsY = 4 // figsX
    # Put the frequencies in the best units possible
    freq1, units1 = _best_freq_units(freq1)
    freq2, units2 = _best_freq_units(freq2)
    
    sortedMapper = sorted(standMapper)
    for k, aStand in enumerate(sortedMapper):
        i = standMapper.index(aStand)
        if standMapper[i]%4//2+1 == 1:
            freq = freq1
            units = units1
        else:
            freq = freq2
            units = units2
            
        ax = fig.add_subplot(figsX,figsY,k+1)
        currSpectra = np.squeeze( np.log10(spec[i,:])*10.0 )
        ax.plot(freq, currSpectra, label=f"{i+1} (avg)")
        
        # If there is more than one chunk, plot the difference between the global 
        # average and each chunk
        if nChunks > 1 and not args.disable_chunks:
            for j in range(nChunks):
                # Some files are padded by zeros at the end and, thus, carry no 
                # weight in the average spectra.  Skip over those.
                if masterWeight[j,i,:].sum() == 0:
                    continue
                    
                # Calculate the difference between the spectra and plot
                subspectra = np.squeeze( np.log10(masterSpectra[j,i,:])*10.0 )
                diff = subspectra - currSpectra
                ax.plot(freq, diff, label='%i' % j)
                
        ax.set_title(f"Beam {standMapper[i]//4+1}, Tune. {standMapper[i]%4//2+1}, Pol. {standMapper[i]%2}")
        ax.set_xlabel(f"Frequency [{units}]")
        ax.set_ylabel('P.S.D. [dB/RBW]')
        ax.set_xlim([freq.min(), freq.max()])
        ax.legend(loc=0)
        
        print(f"For beam {standMapper[i]//4+1}, tune. {standMapper[i]%4//2+1}, pol. {standMapper[i]%2} maximum in PSD at {freq[np.where(spec[i,:]==spec[i,:].max())][0]:.3f} {units}")
        
    print(f"RBW: {freq[1]-freq[0]:.4f} {units}")
    plt.subplots_adjust(hspace=0.35, wspace=0.30)
    plt.show()
    
    # Save spectra image if requested
    if args.output is not None:
        fig.savefig(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='read in a DRX file and create a collection of time-averaged spectra', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
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
    parser.add_argument('-o', '--output', type=str, 
                        help='output file name for spectra image')
    args = parser.parse_args()
    main(args)
