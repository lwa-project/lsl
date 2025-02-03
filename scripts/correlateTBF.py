#!/usr/bin/env python3

"""
Example script that reads in TBF data and runs a cross-correlation on it.  
The results are saved in the FITS IDI format.
"""

import os
import sys
import time
import numpy as np
import argparse

from astropy.constants import c as speedOfLight
speedOfLight = speedOfLight.to('m/s').value

from lsl.reader.base import CI8
from lsl.reader.ldp import LWADataFile, TBFFile
from lsl.common import stations, metabundle
from lsl.correlator import uvutils
from lsl.correlator import fx as fxc
from lsl.correlator._core import XEngine2
from lsl.writer import fitsidi, measurementset
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


def process_chunk(idf, site, good, filename, freq_decim=1, int_time=5.0, pols=['xx',], chunk_size=100):
    """
    Given a lsl.reader.ldp.TBFFile instances and various parameters for the
    cross-correlation, write cross-correlate the data and save it to a file.
    """
    
    # Get antennas
    antennas = site.antennas
    
    # Sort out the integration time
    if int_time == 0.0:
        int_time = None
        
    # Get the metadata
    freq = idf.get_info('freq1')
    srate = idf.get_info('sample_rate')
    
    # Break the frequency range into IFs
    chan = np.round(freq / srate)
    nif = len(np.where(np.diff(chan) > 1)[0]) + 1
    freq = freq.reshape(nif, -1)
    
    # Decimate in frequency if requested
    if freq_decim > 1:
        freq = freq.reshape(nif, -1, freq_decim)
        freq = freq.mean(axis=2)
    freq_flat = freq.ravel()
    
    # Create the list of good digitizers and a digitizer to Antenna instance mapping.  
    # These are:
    #  toKeep  -> mapping of digitizer number to array location
    #  mapper -> mapping of Antenna instance to array location
    toKeep = [antennas[i].digitizer-1 for i in good]
    mapper = [antennas[i] for i in good]
    
    # Create a list of unqiue stands to know what style of IDI file to create
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
            readT, t, data = idf.read(int_time, return_ci8=(freq_decim==1))
        except Exception as e:
            print(f"Error: {str(e)}")
            continue
            
        ## Prune out what we don't want
        data = data[toKeep,:,:]
        
        ## Apply frequency decimation
        if freq_decim > 1:
            data = data.reshape(data.shape[0], -1, freq_decim, data.shape[2])
            data = data.mean(axis=2)
            
        ## CI8 data handling
        if data.dtype == CI8:
            data = data.view(np.int8)
            data = data.reshape(data.shape[:-1]+(-1,2))
            
        ## Split the polarizations
        antennasX, antennasY = [a for i,a in enumerate(antennas) if a.pol == 0 and i in toKeep], [a for i,a in enumerate(antennas) if a.pol == 1 and i in toKeep]
        dataX, dataY = data[0::2,...], data[1::2,...]
        validX = np.ones((dataX.shape[0],dataX.shape[2]), dtype=np.uint8)
        validY = np.ones((dataY.shape[0],dataY.shape[2]), dtype=np.uint8)
        
        setTime = t
        if s == 0:
            ref_time = setTime
            
        # Setup the set time as a python datetime instance so that it can be easily printed
        setDT = setTime.datetime
        print(f"Working on set #{s+1} ({setTime-ref_time:.3f} seconds after set #1 = {setDT.strftime('%Y/%m/%d %H:%M:%S.%f')}")
        
        # Loop over polarization products
        for pol in pols:
            print(f"->  {pol}")
            if pol[0] == 'x':
                a1, d1, v1 = antennasX, dataX, validX
            else:
                a1, d1, v1 = antennasY, dataY, validY
            if pol[1] == 'x':
                a2, d2, v2 = antennasX, dataX, validX
            else:
                a2, d2, v2 = antennasY, dataY, validY
                
            ## Get the baselines
            baselines = uvutils.get_baselines(a1, antennas2=a2, include_auto=True)
            
            ## Run the cross multiply and accumulate
            vis = XEngine2(d1, d2, v1, v2)
            
            ## Apply the cable delays as phase rotations
            for k,(ant1,ant2) in enumerate(baselines):
                gain1 = np.sqrt( ant1.cable.gain(freq_flat) )
                phaseRot1 = np.exp(2j*np.pi*freq_flat*(ant1.cable.delay(freq_flat) \
                                                       -ant1.stand.z/speedOfLight))
                gain2 = np.sqrt( ant2.cable.gain(freq_flat) )
                phaseRot2 = np.exp(2j*np.pi*freq_flat*(ant2.cable.delay(freq_flat) \
                                                       -ant2.stand.z/speedOfLight))
                vis[k,:] *= phaseRot2.conj()*phaseRot1 / gain2 / gain1
                
            # If we are in the first polarazation product of the first iteration,  setup
            # the FITS IDI file.
            if s  == 0 and pol == pols[0]:
                pol1, pol2 = fxc.pol_to_pols(pol)
                
                fits = writer_class(filename, ref_time=ref_time)
                fits.set_stokes(pols)
                for f in range(freq.shape[0]):
                    fits.set_frequency(freq[f,:])
                fits.set_geometry(site, [a for a in mapper if a.pol == pol1])
                
            # Convert the setTime to a MJD and save the visibilities to the FITS IDI file
            fits.add_data_set(setTime, readT, baselines, vis, pol=pol)
        print(f"->  Cummulative Wall Time: {time.time()-wallTime:.3f} s ({(time.time()-wallTime)/(s+1):.3f} s per integration)")
        
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
    elif args.lwana:
        station = stations.lwana
    else:
        station = stations.lwasv
    antennas = station.antennas
    
    with LWADataFile(args.filename) as idf:
        if not isinstance(idf, TBFFile):
            raise RuntimeError(f"File '{os.path.basename(args.filename)}' does not appear to be a valid TBF file")
            
        jd = idf.get_info('start_time').jd
        date = idf.get_info('start_time').datetime
        nFpO = idf.get_info('nchan') // 12
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
        if args.avg_time == 0.0:
            args.avg_time = nInts/sample_rate
        nFrames = min([int(args.avg_time*sample_rate), nInts])
        nSets = idf.get_info('nframe') // nFpO // nFrames
        args.offset = idf.offset(args.offset)
        nSets = nSets - int(args.offset*sample_rate) // nFrames
        
        central_freq = idf.get_info('freq1')
        chan = np.round(central_freq / sample_rate)
        nif = len(np.where(np.diff(chan) > 1)[0]) + 1
        central_freq = central_freq.reshape(nif, -1)
        central_freq = central_freq[:,central_freq.shape[1]//2]
        
        print(f"Data type: {str(type(idf))}")
        print(f"Samples per observations: {nFpO}")
        print(f"Sampling rate: {sample_rate} Hz")
        print("Tuning frequency: %s Hz" % (', '.join("%.3f" % v for v in central_freq)))
        print(f"Captures in file: {nInts} ({nInts/sample_rate:.3f} s)")
        print("==")
        print(f"Station: {station.name}")
        print(f"Date observed: {date}")
        print(f"Julian day: {jd:.5f}")
        print(f"Offset: {args.offset:.3f} s ({args.offset*sample_rate} frames)")
        print(f"Integration Time: {nFrames/sample_rate:.3f} s")
        print(f"Number of integrations in file: {nSets}")
        
        # Make sure we don't try to do too many sets
        args.samples = min([args.samples, nSets])
        
        # Loop over junks of 100 integrations to make sure that we don't overflow 
        # the FITS IDI memory buffer
        s = 0
        basename = os.path.split(args.filename)[1]
        basename, ext = os.path.splitext(basename)
        
        leftToDo = args.samples
        while leftToDo > 0:
            if args.casa:
                fitsFilename = f"{basename}.ms_{s+1}"
            else:
                fitsFilename = f"{basename}.FITS_{s+1}"
                
            if leftToDo > 100:
                chunk = 100
            else:
                chunk = leftToDo
            process_chunk(idf, station, good, fitsFilename, int_time=args.avg_time,
                         freq_decim=args.decimate, pols=args.products,
                         chunk_size=chunk)
            
            s += 1
            leftToDo = leftToDo - chunk


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='cross-correlate data in a TBF file', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to correlate')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-n', '--lwana', action='store_true',
                        help='use LWA-NA instead of LWA-SV')
    parser.add_argument('-t', '--avg-time', type=aph.positive_or_zero_float, default=0.0, 
                        help='time window to average visibilities in seconds; 0 = integrate the entire file')
    parser.add_argument('-s', '--samples', type=aph.positive_int, default=1, 
                        help='number of average visibilities to generate')
    parser.add_argument('-o', '--offset', type=aph.positive_or_zero_float, default=0.0, 
                        help='offset into the file before starting correlation in seconds')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', 
                        help='run %(prog)s in silent mode')
    parser.add_argument('-a', '--all', action='store_true', 
                        help='correlated all dipoles regardless of their status')
    parser.add_argument('-d', '--decimate', type=aph.positive_int, default=1,
                        help='frequency decimation factor')
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
    
