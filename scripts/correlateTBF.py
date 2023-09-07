#!/usr/bin/env python

"""
Example script that reads in TBF data and runs a cross-correlation on it.  
The results are saved in the FITS IDI format.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import time
import numpy
import argparse

from astropy.constants import c as speedOfLight
speedOfLight = speedOfLight.to('m/s').value

from lsl.reader.ldp import LWASVDataFile, LWANADataFile, TBFFile
from lsl.common import stations, metabundleADP, metabundleNDP
from lsl.correlator import uvutils
from lsl.correlator import fx as fxc
from lsl.correlator._core import XEngine2
from lsl.writer import fitsidi, measurementset
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


def process_chunk(idf, site, good, filename, int_time=5.0, pols=['xx',], chunk_size=100):
    """
    Given a lsl.reader.ldp.TBNFile instances and various parameters for the 
    cross-correlation, write cross-correlate the data and save it to a file.
    """
    
    # Get antennas
    antennas = site.antennas
    
    # Get the metadata
    freq = idf.get_info('freq1')
    
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
            readT, t, data = idf.read(int_time)
        except Exception as e:
            print("Error: %s" % str(e))
            continue
            
        ## Prune out what we don't want
        data = data[toKeep,:,:]
        
        ## Split the polarizations
        antennasX, antennasY = [a for i,a in enumerate(antennas) if a.pol == 0 and i in toKeep], [a for i,a in enumerate(antennas) if a.pol == 1 and i in toKeep]
        dataX, dataY = data[0::2,:,:], data[1::2,:,:]
        validX = numpy.ones((dataX.shape[0],dataX.shape[2]), dtype=numpy.uint8)
        validY = numpy.ones((dataY.shape[0],dataY.shape[2]), dtype=numpy.uint8)
        
        ## Apply the cable delays as phase rotations
        for i in range(dataX.shape[0]):
            gain = numpy.sqrt( antennasX[i].cable.gain(freq) )
            phaseRot = numpy.exp(2j*numpy.pi*freq*(antennasX[i].cable.delay(freq) \
                                                   -antennasX[i].stand.z/speedOfLight))
            for j in range(dataX.shape[2]):
                dataX[i,:,j] *= phaseRot / gain
        for i in range(dataY.shape[0]):
            gain = numpy.sqrt( antennasY[i].cable.gain(freq) )
            phaseRot = numpy.exp(2j*numpy.pi*freq*(antennasY[i].cable.delay(freq)\
                                                   -antennasY[i].stand.z/speedOfLight))
            for j in range(dataY.shape[2]):
                dataY[i,:,j] *= phaseRot / gain
                
        setTime = t
        if s == 0:
            ref_time = setTime
            
        # Setup the set time as a python datetime instance so that it can be easily printed
        setDT = setTime.datetime
        print("Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-ref_time), setDT.strftime("%Y/%m/%d %H:%M:%S.%f")))
        
        # Loop over polarization products
        for pol in pols:
            print("->  %s" % pol)
            if pol[0] == 'x':
                a1, d1, v1 = antennasX, dataX, validX
            else:
                a1, d1, v1 = antennasY, dataY, validY
            if pol[1] == 'x':
                a2, d2, v2 = antennasX, dataX, validX
            else:
                a2, d2, v2 = antennasY, dataY, validY
                
            ## Get the baselines
            baselines = uvutils.get_baselines(a1, antennas2=a2, include_auto=True, indicies=True)
            blList = []
            for bl in range(len(baselines)):
                blList.append( (a1[baselines[bl][0]], a2[baselines[bl][1]]) )
                
            ## Run the cross multiply and accumulate
            vis = XEngine2(d1, d2, v1, v2)
            
            # Select the right range of channels to save
            toUse = numpy.where( (freq>5.0e6) & (freq<93.0e6) )
            toUse = toUse[0]
            
            # If we are in the first polarazation product of the first iteration,  setup
            # the FITS IDI file.
            if s  == 0 and pol == pols[0]:
                pol1, pol2 = fxc.pol_to_pols(pol)
                
                fits = writer_class(filename, ref_time=ref_time)
                fits.set_stokes(pols)
                fits.set_frequency(freq[toUse])
                fits.set_geometry(site, [a for a in mapper if a.pol == pol1])
                
            # Convert the setTime to a MJD and save the visibilities to the FITS IDI file
            fits.add_data_set(setTime, readT, blList, vis[:,toUse], pol=pol)
        print("->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)/(s+1)))
        
    # Cleanup after everything is done
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
            try:
                station = metabundleADP.get_station(args.metadata, apply_sdm=True)
            except ValueError:
                station = metabundleNDP.get_station(args.metadata, apply_sdm=True)
    elif args.lwana:
        station = stations.lwana
    else:
        station = stations.lwasv
    antennas = station.antennas
    
    try:
        idf = LWASVDataFile(filename)
    except RuntimeError:
        idf = LWANADataFile(filename)
    if not isinstance(idf, TBFFile):
        raise RuntimeError("File '%s' does not appear to be a valid TBF file" % os.path.basename(filename))
        
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
    print("Found %i good stands to use" % (len(good)//2,))
    for i in good:
        print("%3i, %i" % (antennas[i].stand.id, antennas[i].pol))
        
    # Number of frames to read in at once and average
    nFrames = min([int(args.avg_time*sample_rate), nInts])
    nSets = idf.get_info('nframe') // nFpO // nFrames
    args.offset = idf.offset(args.offset)
    nSets = nSets - int(args.offset*sample_rate) // nFrames
    
    central_freq = idf.get_info('freq1')
    central_freq = central_freq[len(central_freq)//2]
    
    print("Data type:  %s" % type(idf))
    print("Samples per observations: %i" % nFpO)
    print("Sampling rate: %.1f Hz" % sample_rate)
    print("Tuning frequency: %.3f Hz" % central_freq)
    print("Captures in file: %i (%.3f s)" % (nInts, nInts / sample_rate))
    print("==")
    print("Station: %s" % station.name)
    print("Date observed: %s" % date)
    print("Julian day: %.5f" % jd)
    print("Offset: %.3f s (%i frames)" % (args.offset, args.offset*sample_rate))
    print("Integration Time: %.3f s" % (nFrames/sample_rate))
    print("Number of integrations in file: %i" % nSets)
    
    # Make sure we don't try to do too many sets
    if args.samples > nSets:
        args.samples = nSets
        
    # Loop over junks of 100 integrations to make sure that we don't overflow 
    # the FITS IDI memory buffer
    s = 0
    leftToDo = args.samples
    basename = os.path.split(filename)[1]
    basename, ext = os.path.splitext(basename)
    while leftToDo > 0:
        if args.casa:
            fitsFilename = "%s.ms_%i" % (basename, (s+1),)
        else:
            fitsFilename = "%s.FITS_%i" % (basename, (s+1),)
            
        if leftToDo > 100:
            chunk = 100
        else:
            chunk = leftToDo
            
        process_chunk(idf, station, good, fitsFilename, int_time=args.avg_time, 
                     pols=args.products, chunk_size=chunk)
                    
        s += 1
        leftToDo = leftToDo - chunk
        
    idf.close()


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
    parser.add_argument('-t', '--avg-time', type=aph.positive_float, default=1.0, 
                        help='time window to average visibilities in seconds')
    parser.add_argument('-s', '--samples', type=aph.positive_int, default=1, 
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
    
