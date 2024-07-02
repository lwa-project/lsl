#!/usr/bin/env python3

"""
Example script that reads in TBW data and runs a cross-correlation on it.  
The results are saved in the FITS IDI format.
"""

import os
import sys
import time
import numpy as np
import argparse

from lsl.reader.ldp import LWA1DataFile, TBWFile
from lsl.common import stations, metabundle
from lsl.correlator import fx as fxc
from lsl.writer import fitsidi, measurementset
from lsl.common.progress import ProgressBar
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


def process_chunk(idf, site, good, filename, LFFT=64, overlap=1, pfb=False, pols=['xx','yy']):
    """
    Given an lsl.reader.ldp.TBWFile instances and various parameters for the 
    cross-correlation, write cross-correlate the data and save it to a file.
    """
    
    # Get antennas
    antennas = site.antennas
    
    # Get the metadata
    sample_rate = idf.get_info('sample_rate')
    
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
            
    wallTime = time.time()
    readT, t, data = idf.read()
    setTime = t
    ref_time = t
    
    # Setup the set time as a python datetime instance so that it can be easily printed
    setDT = setTime.datetime
    print("Working on set #1 (%.3f seconds after set #1 = %s)" % ((setTime-ref_time), setDT.strftime("%Y/%m/%d %H:%M:%S.%f")))
    
    # In order for the TBW stuff to actaully run, we need to run in with sub-
    # integrations.  8 sub-integrations (61.2 ms / 8 = 7.7 ms per section) 
    # seems to work ok with a "reasonable" number of channels.
    nSec = 8
    secSize = data.shape[1]//nSec
    
    # Loop over polarizations (there should be only 1)
    for pol in pols:
        print("-> %s" % pol)
        try:
            tempVis *= 0    # pylint:disable=undefined-variable
        except NameError:
            pass
            
        # Set up the progress bar so we can keep up with how the sub-integrations 
        # are progressing
        pb = ProgressBar(max=nSec)
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.flush()
        
        # Loop over sub-integrations (set by nSec)
        for k in range(nSec):
            blList, freq, vis = fxc.FXMaster(data[toKeep,k*secSize:(k+1)*secSize], mapper, LFFT=LFFT, overlap=overlap, pfb=pfb, include_auto=True, verbose=False, sample_rate=sample_rate, central_freq=0.0, pol=pol, return_baselines=True, gain_correct=True)
            
            toUse = np.where( (freq>=5.0e6) & (freq<=93.0e6) )
            toUse = toUse[0]
            
            try:
                tempVis += vis
            except NameError:
                tempVis = vis
                
            pb.inc(amount=1)
            sys.stdout.write(pb.show()+'\r')
            sys.stdout.flush()
            
        # Average the sub-integrations together
        vis = tempVis / float(nSec)
        
        # Set up the FITS IDI file is we need to
        if pol == pols[0]:
            pol1, pol2 = fxc.pol_to_pols(pol)
            
            fits = writer_class(filename, ref_time=ref_time)
            fits.set_stokes(pols)
            fits.set_frequency(freq[toUse])
            fits.set_geometry(site, [a for a in mapper if a.pol == pol1])
            
        # Add the visibilities
        fits.add_data_set(setTime, readT, blList, vis[:,toUse], pol=pol)
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
    print("->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)))
    
    fits.write()
    fits.close()
    del(fits)
    del(data)
    del(vis)
    return True


def main(args):
    # Parse command line options
    filename = args.filename

    # Setup the LWA station information
    if args.metadata is not None:
        try:
            station = stations.parse_ssmif(args.metadata)
        except ValueError:
            station = metabundle.get_station(args.metadata, apply_sdm=True)
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    idf = LWA1DataFile(filename)
    if not isinstance(idf, TBWFile):
        raise RuntimeError("File '%s' does not appear to be a valid TBW file" % os.path.basename(filename))
        
    jd = idf.get_info('start_time').jd
    date = idf.get_info('start_time').datetime
    sample_rate = idf.get_info('sample_rate')
    nInts = idf.get_info('nframe') // (30000 * len(antennas) // 2)
    
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
                
    # Select which polarization to use
    good = []
    for antX in goodX:
        for antY in goodY:
            if antX.stand.id == antY.stand.id:
                good.append(antX.digitizer-1)
                good.append(antY.digitizer-1)
                break
                
    # Report on the valid stands found.  This is a little verbose,
    # but nice to see.
    print("Found %i good stands to use" % (len(good)//2,))
    for i in good:
        print("%3i, %i" % (antennas[i].stand.id, antennas[i].pol))
        
    # Number of frames to read in at once and average
    nFrames = 30000
    nSets = idf.get_info('nframe') // (30000*len(antennas)//2)
    
    print("Data type:  %s" % type(idf))
    print("Captures in file: %i (%.3f s)" % (nInts, nInts*30000*400/sample_rate))
    print("==")
    print("Station: %s" % station.name)
    print("Date observed: %s" % date)
    print("Julian day: %.5f" % jd)
    print("Integration Time: %.3f s" % (400*nFrames/sample_rate))
    print("Number of integrations in file: %i" % nSets)
    print("==")
    
    basename = os.path.split(filename)[1]
    basename, ext = os.path.splitext(basename)
    
    if args.casa:
        fitsFilename = "%s.ms_1" % (basename,)
    else:
        fitsFilename = "%s.FITS_1" % (basename,)
        
    process_chunk(idf, station, good, fitsFilename, LFFT=args.fft_length, overlap=1, pfb=args.pfb, pols=args.products)
    
    idf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='cross-correlate data in a TBW file', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to correlate')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=16, 
                        help='set FFT length')
    parser.add_argument('-p', '--pfb', action='store_true', 
                        help='enabled the PFB on the F-engine')
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
    
