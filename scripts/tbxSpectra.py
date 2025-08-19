#!/usr/bin/env python3

"""
Given a TBT/TBS file, plot the time averaged spectra for each digitizer input.
"""
    
import os
import sys
import math
import numpy as np
import argparse

from lsl.common import stations, metabundle
from lsl.reader.ldp import LWADataFile, TBXFile
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
    elif args.lwasv:
        station = stations.lwasv
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    idf = LWADataFile(args.filename)
    if not isinstance(idf, TBXFile):
        raise RuntimeError(f"File '{os.path.basename(args.filename)}' does not appear to be a valid TBT or TBS file")
        
    nFrames = idf.get_info('nframe')
    antpols = len(antennas)
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    beginDate = idf.get_info('start_time').datetime
    
    # Make sure the TBX stand count is consistent with how many antennas we have
    if antpols != idf.get_info('nantenna'):
        raise RuntimeError(f"Number of stands in the station ({antpols//2}) does not match what is in the data ({idf.get_info('nantenna')//2})")
        
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nchannels = idf.get_info('nchan')
    nFpO = nchannels // idf.get_info('frame_channel_count')
    
    # Figure out how many chunks we need to work with
    nChunks = 1
    if args.average > 0.5:
        nChunks = int(np.ceil(args.average/0.5))
        
    # Calculate the frequencies
    freq = idf.get_info('freq1')
    sample_rate = idf.get_info('sample_rate')
    central_freq = freq*1.0
    chan = np.round(central_freq / sample_rate)
    nif = len(np.where(np.diff(chan) > 1)[0]) + 1
    central_freq = central_freq.reshape(nif, -1)
    central_freq = central_freq[:,central_freq.shape[1]//2]
    
    # File summary
    print(f"Filename: {args.filename}")
    print(f"Date of First Frame: {str(beginDate)}")
    print(f"Channel Count: {nchannels}")
    print(f"Sampling rate: {sample_rate} Hz")
    print("Tuning frequency: %s Hz" % (', '.join("%.3f" % v for v in central_freq)))
    print(f"Frames: {nFrames} ({nFrames/nFpO/sample_rate:.3f} s)")
    print("===")
    print(f"Chunks: {nChunks}")
    
    for i in range(nChunks):
        print(f"Working on chunk #{i+1} of {nChunks}")
        
        try:
            readT, t, data = idf.read(args.average/nChunks, return_ci8=True)
        except Exception as e:
            print(f"Error: {str(e)}")
            continue
            
        # Detect power and integrate
        data = data['re']**2 + data['im']**2
        try:
            spec += data.mean(axis=2)
        except NameError:
            spec = data.mean(axis=2)
    spec /= nChunks
    
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
            ax.plot(freq, currSpectra, label=f"Stand: {antennas[j].stand.id}, Pol: {antennas[j].pol} (Dig: {antennas[j].digitizer})")
            
            ax.set_title(f"Stand: {antennas[j].stand.id} ({antennas[j].pol}); Dig: {antennas[j].digitizer} [{antennas[j].combined_status}]")
            ax.set_xlabel(f"Frequency [{units}]")
            ax.set_ylabel('P.S.D. [dB/RBW]')
            ax.set_ylim([-10, 30])
            
        # Save spectra image if requested
        if args.output is not None:
            base, ext = os.path.splitext(args.output)
            outFigure = f"{base}-{i+1:02d}{ext}"
            fig.savefig(outFigure)
            
        plt.draw()
        
    print(f"RBW: {freq[1]-freq[0]:.4f} {units}")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='read in a TBT or TBS file and create a collection of time-averaged spectra', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of the SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-s', '--lwasv', action='store_true', 
                        help='use LWA-SV instead of LWA1')
    parser.add_argument('-n', '--lwana', action='store_true', 
                        help='use LWA-NA instead of LWA1')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=1.0, 
                        help='number of seconds of data to average for spectra')
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
    
