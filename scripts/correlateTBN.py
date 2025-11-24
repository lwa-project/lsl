#!/usr/bin/env python3

"""
Example script that reads in TBN data and runs a cross-correlation on it.  
The results are saved in the FITS IDI format.
"""

import os
import sys
import time
import numpy as np
import argparse

from lsl.reader.ldp import LWADataFile, TBNFile
from lsl.common import stations, metabundle
from lsl.correlator import fx as fxc
from lsl.writer import fitsidi, measurementset
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


def process_chunk(idf, site, good, filename, int_time=5.0, LFFT=64, overlap=1, pfb=False, pols=['xx',], chunk_size=100):
    """
    Given a lsl.reader.ldp.TBNFile instances and various parameters for the 
    cross-correlation, write cross-correlate the data and save it to a file.
    """
    
    # Get antennas
    antennas = site.antennas
    
    # Get the metadata
    sample_rate = idf.get_info('sample_rate')
    central_freq = idf.get_info('freq1')
    
    # Create the list of good digitizers and a digitizer to Antenna instance mapping.  
    # These are:
    #  toKeep  -> mapping of digitizer number to array location
    #  mapper -> mapping of Antenna instance to array location
    toKeep = [antennas[i].digitizer-1 for i in good]
    mapper = [antennas[i] for i in good]
    
    # Create a list of unique stands to know what style of IDI file to create
    stands = set( [antennas[i].stand.id for i in good] )
    
    # Figure out the output mode
    if os.path.splitext(filename)[1].find('.ms_') != -1:
        writer_class = measurementset.Ms
    else:
        if len(stands) > 255:
            writer_class = fitsidi.ExtendedIdi
        else:
            writer_class = fitsidi.Idi
            
    # Main loop over the input file to read in the data and organize it.  Several control 
    # variables are defined for this:
    #  ref_time -> time (in seconds since the UNIX epoch) for the first data set
    #  setTime -> time (in seconds since the UNIX epoch) for the current data set
    ref_time = 0.0
    setTime = 0.0
    wallTime = time.time()
    for s in range(chunk_size):
        try:
            readT, t, data = idf.read(int_time)
        except Exception as e:
            print(f"Error: {str(e)}")
            continue
            
        ## Prune out what we don't want
        data = data[toKeep,:]
        
        setTime = t
        if s == 0:
            ref_time = setTime
            
        # Setup the set time as a python datetime instance so that it can be easily printed
        setDT = setTime.datetime
        print(f"Working on set #{s+1} ({setTime-ref_time:.3f} seconds after set #1 = {setDT.strftime('%Y/%m/%d %H:%M:%S.%f')}")
        
        # Loop over polarization products
        for pol in pols:
            print(f"->  {pol}")
            blList, freq, vis = fxc.FXMaster(data, mapper, LFFT=LFFT, overlap=overlap, pfb=pfb, include_auto=True, verbose=False, sample_rate=sample_rate, central_freq=central_freq, pol=pol, return_baselines=True, gain_correct=True)
            
            # Select the right range of channels to save
            toUse = np.where( (freq>5.0e6) & (freq<93.0e6) )
            toUse = toUse[0]
            
            # If we are in the first polarization product of the first iteration,  setup
            # the FITS IDI file.
            if s  == 0 and pol == pols[0]:
                pol1, pol2 = fxc.pol_to_pols(pol)
                
                fits = writer_class(filename, ref_time=ref_time)
                fits.set_stokes(pols)
                fits.set_frequency(freq[toUse])
                fits.set_geometry(site, [a for a in mapper if a.pol == pol1])
                
            # Convert the setTime to a MJD and save the visibilities to the FITS IDI file
            fits.add_data_set(setTime, readT, blList, vis[:,toUse], pol=pol)
        print(f"->  Cumulative Wall Time: {time.time()-wallTime:.3f} s ({(time.time()-wallTime)/(s+1):.3f} s per integration)")
        
    # Cleanup after everything is done
    fits.write()
    fits.close()
    del(fits)
    del(data)
    del(vis)
    return True


def main(args):
    # Setup the LWA station information
    if args.metadata is not None:
        try:
            station = stations.parse_ssmif(args.metadata)
        except ValueError:
            station = metabundle.get_station(args.metadata, apply_sdm=True)
    elif args.lwasv:
        station = stations.lwasv
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    idf = LWADataFile(args.filename)
    if not isinstance(idf, TBNFile):
        raise RuntimeError(f"File '{os.path.basename(args.filename)}' does not appear to be a valid TBN file")
        
    jd = idf.get_info('start_time').jd
    date = idf.get_info('start_time').datetime
    nFpO = len(antennas)
    sample_rate = idf.get_info('sample_rate')
    nInts = idf.get_info('nframe') // nFpO
    
    # Get valid stands for both polarizations
    goodX = []
    goodY = []
    for i in range(len(antennas)):
        ant = antennas[i]
        if ant.combined_status != 33 and not args.all:
            pass
        else:
            if ant.pol == 0:
                goodX.append(ant)
            else:
                goodY.append(ant)
                
    # Now combine both lists to come up with stands that
    # are in both so we can form the cross-polarization 
    # products if we need to
    good = []
    for antX in goodX:
        for antY in goodY:
            if antX.stand.id == antY.stand.id:
                good.append( antX.digitizer-1 )
                good.append( antY.digitizer-1 )
                
    # Report on the valid stands found.  This is a little verbose,
    # but nice to see.
    print(f"Found {len(good)//2} good stands to use")
    for i in good:
        print(f"{antennas[i].stand.id:3d}, {antennas[i].pol}")
        
    # Number of frames to read in at once and average
    nFrames = int(args.avg_time*sample_rate/512)
    nSets = idf.get_info('nframe') // nFpO // nFrames
    args.offset = idf.offset(args.offset)
    nSets = nSets - int(args.offset*sample_rate/512) // nFrames
    
    central_freq = idf.get_info('freq1')
    
    print(f"Data type: {str(type(idf))}")
    print(f"Samples per observations: {nFpO//2} per pol.")
    print(f"Sampling rate: {sample_rate} Hz")
    print(f"Tuning frequency: {central_freq:.3f} Hz")
    print(f"Captures in file: {nInts} ({nInts*512/sample_rate:.3f} s)")
    print("==")
    print(f"Station: {station.name}")
    print(f"Date observed: {date}")
    print(f"Julian day: {jd:.5f}")
    print(f"Offset: {args.offset:.3f} s ({args.offset*sample_rate/512} frames)")
    print(f"Integration Time: {512*nFrames/sample_rate:.3f} s")
    print(f"Number of integrations in file: {nSets}")
    
    # Make sure we don't try to do too many sets
    if args.samples > nSets:
        args.samples = nSets
        
    # Loop over chunks of 300 integrations to make sure that we don't overflow 
    # the FITS IDI memory buffer
    s = 0
    leftToDo = args.samples
    basename = os.path.split(args.filename)[1]
    basename, ext = os.path.splitext(basename)
    while leftToDo > 0:
        if args.casa:
            fitsFilename = f"{basename}.ms_{s+1}"
        else:
            fitsFilename = f"{basename}.FITS_{s+1}"
            
        if leftToDo > 100:
            chunk = 100
        else:
            chunk = leftToDo
            
        process_chunk(idf, station, good, fitsFilename, int_time=args.avg_time, LFFT=args.fft_length, 
                    overlap=1, pfb=args.pfb, pols=args.products, chunk_size=chunk)
                    
        s += 1
        leftToDo = leftToDo - chunk
        
    idf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='cross-correlate data in a TBN file', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to correlate')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-v', '--lwasv', action='store_true', 
                        help='use LWA-SV instead of LWA1')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=16, 
                        help='set FFT length')
    parser.add_argument('-p', '--pfb', action='store_true', 
                        help='enabled the PFB on the F-engine')
    parser.add_argument('-t', '--avg-time', type=aph.positive_float, default=5.0, 
                        help='time window to average visibilities in seconds')
    parser.add_argument('-s', '--samples', type=aph.positive_int, default=10, 
                        help='number of average visibilities to generate')
    parser.add_argument('-o', '--offset', type=aph.positive_or_zero_float, default=0.0, 
                        help='offset into the file before starting correlation in seconds')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', 
                        help='run %(prog)s in silent mode')
    parser.add_argument('-a', '--all', action='store_true', 
                        help='correlated all dipoles regardless of their status')
    pgroup = parser.add_mutually_exclusive_group(required=True)
    pgroup.add_argument('-x', '--xx', dest='products', action='store_const', const=['xx',], 
                        help='compute only the XX polarization product')
    pgroup.add_argument('-y', '--yy', dest='products', action='store_const', const=['yy',], 
                        help='compute only the YY polarization product')
    pgroup.add_argument('-2', '--two-products', dest='products', action='store_const', const=['xx','yy'], 
                        help='compute only the XX and YY polarization products')
    pgroup.add_argument('-4', '--four-products', dest='products', action='store_const', const=['xx','yy','xy','yx'], 
                        help='compute the XX, XY, YX, and YY polarization products')
    parser.add_argument('--casa', action='store_true',
                        help='write out measurement sets instead of FITS-IDI files')
    args = parser.parse_args()
    main(args)
