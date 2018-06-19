#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example script that reads in TBF data and runs a cross-correlation on it.  
The results are saved in the FITS IDI format.
"""

import os
import sys
import time
import ephem
import numpy
import getopt
from datetime import datetime, timedelta, tzinfo

from lsl import astro
from lsl.reader import tbf, errors
from lsl.common import stations, metabundleADP
from lsl.correlator import uvUtils
from lsl.correlator import fx as fxc
from lsl.correlator._core import XEngine2
from lsl.writer import fitsidi


class UTC(tzinfo):
    """tzinfo object for UTC time."""
    
    def utcoffset(self, dt):
        return timedelta(0)
        
    def tzname(self, dt):
        return "UTC"
        
    def dst(self, dt):
        return timedelta(0)


def usage(exitCode=None):
    print """correlateTBF.py - cross-correlate data in a TBF file

Usage: correlateTBF.py [OPTIONS] file

Options:
-h, --help             Display this help information
-m, --metadata         Name of SSMIF or metadata tarball file to use for 
                    mappings
-q, --quiet            Run correlateTBN in silent mode
-a, --all              Correlated all dipoles regardless of their status 
                    (default = no)
-x, --xx               Compute only the XX polarization product (default)
-y, --yy               Compute only the YY polarization product
-2, --two-products     Compute both the XX and YY polarization products
-4, --four-products    Compute all for polariation products:  XX, YY, XY, 
                    and YX.
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseConfig(args):
    config = {}
    # Command line flags - default values
    config['metadata'] = ''
    config['verbose'] = True
    config['all'] = False
    config['products'] = ['xx',]
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "hm:qa24xy", ["help", "metadata=", "quiet", "all", "two-products", "four-products", "xx", "yy"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-m', '--metadata'):
            config['metadata'] = value
        elif opt in ('-q', '--quiet'):
            config['verbose'] = False
        elif opt in ('-a', '--all'):
            config['all'] = True
        elif opt in ('-2', '--two-products'):
            config['products'] = ['xx', 'yy']
        elif opt in ('-4', '--four-products'):
            config['products'] = ['xx', 'yy', 'xy', 'yx']
        elif opt in ('-x', '--xx'):
            config['products'] = ['xx',]
        elif opt in ('-y', '--yy'):
            config['products'] = ['yy',]
        else:
            assert False
            
    # Add in arguments
    config['args'] = arg
    
    # Return configuration
    return config


def processChunk(fh, site, good, filename, intTime=5.0, pols=['xx',], ChunkSize=100):
    """
    Given a lsl.reader.ldp.TBNFile instances and various parameters for the 
    cross-correlation, write cross-correlate the data and save it to a file.
    """
    
    # Get antennas
    antennas = site.getAntennas()
    
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nFrames = os.path.getsize(fh.name) / tbf.FrameSize
    nFramesPerObs = tbf.getFramesPerObs(fh)
    nChannels = tbf.getChannelCount(fh)
    nSamples = 7840
    
    # Figure out how many chunks we need to work with
    nChunks = nFrames / nFramesPerObs
    
    # Pre-load the channel mapper
    cMapper = []
    for i in xrange(2*nFramesPerObs):
        cFrame = tbf.readFrame(fh)
        if cFrame.header.firstChan not in cMapper:
            cMapper.append( cFrame.header.firstChan )
    fh.seek(-2*nFramesPerObs*tbf.FrameSize, 1)
    cMapper.sort()
    
    # Calculate the frequencies
    freq = numpy.zeros(nChannels)
    for i,c in enumerate(cMapper):
        freq[i*12:i*12+12] = c + numpy.arange(12)
    freq *= 25e3
    
    # Get the metadata
    sampleRate = 25e3
    
    # Create the list of good digitizers and a digitizer to Antenna instance mapping.  
    # These are:
    #  toKeep  -> mapping of digitizer number to array location
    #  mapper -> mapping of Antenna instance to array location
    toKeep = [antennas[i].digitizer-1 for i in good]
    mapper = [antennas[i] for i in good]
    
    # Create a list of unqiue stands to know what style of IDI file to create
    stands = set( [antennas[i].stand.id for i in good] )
    
    # Main loop over the input file to read in the data and organize it.  Several control 
    # variables are defined for this:
    #  refTime -> time (in seconds since the UNIX epoch) for the first data set
    #  setTime -> time (in seconds since the UNIX epoch) for the current data set
    refTime = 0.0
    setTime = 0.0
    wallTime = time.time()
    for s in xrange(ChunkSize):
        data = numpy.zeros((256*2,nChannels,nChunks), dtype=numpy.complex64)
        for i in xrange(nChunks):
            # Inner loop that actually reads the frames into the data array
            for j in xrange(nFramesPerObs):
                # Read in the next frame and anticipate any problems that could occur
                try:
                    cFrame = tbf.readFrame(fh)
                except errors.eofError:
                    break
                except errors.syncError:
                    print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbf.FrameSize-1)
                    continue
                if not cFrame.header.isTBF():
                    continue
                    
                firstChan = cFrame.header.firstChan
                
                # Figure out where to map the channel sequence to
                try:
                    aStand = cMapper.index(firstChan)
                except ValueError:
                    cMapper.append(firstChan)
                    aStand = cMapper.index(firstChan)
                
                # Actually load the data.
                if i == 0 and j == 0:
                    t = cFrame.getTime()
                    refCount = cFrame.header.frameCount
                count = cFrame.header.frameCount - refCount
                subData = cFrame.data.fDomain
                subData.shape = (12,512)
                subData = subData.T
                
                data[:,aStand*12:aStand*12+12,count] = subData
        readT = ChunkSize*40e-6
        
        ## Prune out what we don't want
        data = data[toKeep,:,:]
        
        ## Split the polarizations
        antennasX, antennasY = [a for i,a in enumerate(antennas) if a.pol == 0 and i in toKeep], [a for i,a in enumerate(antennas) if a.pol == 1 and i in toKeep]
        dataX, dataY = data[0::2,:,:], data[1::2,:,:]
        validX = numpy.ones((dataX.shape[0],dataX.shape[2]), dtype=numpy.uint8)
        validY = numpy.ones((dataY.shape[0],dataY.shape[2]), dtype=numpy.uint8)
        
        ## Apply the cable delays as phase rotations
        for i in xrange(dataX.shape[0]):
            phaseRot = numpy.exp(2j*numpy.pi*freq*antennasX[i].cable.delay(freq))
            for j in xrange(dataX.shape[2]):
                dataX[i,:,j] *= phaseRot
        for i in xrange(dataY.shape[0]):
            phaseRot = numpy.exp(2j*numpy.pi*freq*antennasY[i].cable.delay(freq))
            for j in xrange(dataY.shape[2]):
                dataY[i,:,j] *= phaseRot
                
        setTime = t
        if s == 0:
            refTime = setTime
            
        # Setup the set time as a python datetime instance so that it can be easily printed
        setDT = datetime.utcfromtimestamp(setTime)
        setDT.replace(tzinfo=UTC())
        print "Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-refTime), setDT.strftime("%Y/%m/%d %H:%M:%S.%f"))
        
        # Loop over polarization products
        for pol in pols:
            print "->  %s" % pol
            if pol[0] == 'x':
                a1, d1, v1 = antennasX, dataX, validX
            else:
                a1, d1, v1 = antennasY, dataY, validY
            if pol[1] == 'x':
                a2, d2, v2 = antennasX, dataX, validX
            else:
                a2, d2, v2 = antennasY, dataY, validY
                
            ## Get the baselines
            baselines = uvUtils.getBaselines(a1, antennas2=a2, IncludeAuto=True, Indicies=True)
            blList = []
            for bl in xrange(len(baselines)):
                blList.append( (a1[baselines[bl][0]], a2[baselines[bl][1]]) )
                
            ## Run the cross multiply and accumulate
            vis = XEngine2(d1, d2, v1, v2)
            
            ## Apply the cable gains
            for bl in xrange(vis.shape[0]):
                cableGain1 = a1[baselines[bl][0]].cable.gain(freq)
                cableGain2 = a2[baselines[bl][1]].cable.gain(freq)
                
                vis[bl,:] /= numpy.sqrt(cableGain1*cableGain2)
            
            # Select the right range of channels to save
            toUse = numpy.where( (freq>5.0e6) & (freq<93.0e6) )
            toUse = toUse[0]
            
            # If we are in the first polarazation product of the first iteration,  setup
            # the FITS IDI file.
            if s  == 0 and pol == pols[0]:
                pol1, pol2 = fxc.pol2pol(pol)
                
                if len(stands) > 255:
                    fits = fitsidi.ExtendedIDI(filename, refTime=refTime)
                else:
                    fits = fitsidi.IDI(filename, refTime=refTime)
                fits.setStokes(pols)
                fits.setFrequency(freq[toUse])
                fits.setGeometry(site, [a for a in mapper if a.pol == pol1])
                
            # Convert the setTime to a MJD and save the visibilities to the FITS IDI file
            obsTime = astro.unix_to_taimjd(setTime)
            fits.addDataSet(obsTime, readT, blList, vis[:,toUse], pol=pol)
        print "->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)/(s+1))
        
    # Cleanup after everything is done
    fits.write()
    fits.close()
    del(fits)
    del(data)
    del(vis)
    return True


def main(args):
    # Parse command line options
    config = parseConfig(args)
    filename = config['args'][0]
    
    # Setup the LWA station information
    if config['metadata'] != '':
        try:
            station = stations.parseSSMIF(config['metadata'])
        except ValueError:
            station = metabundleADP.getStation(config['metadata'], ApplySDM=True)
    else:
        station = stations.lwasv
    antennas = station.getAntennas()
    
    fh = open(config['args'][0], 'rb')
    nFrames = os.path.getsize(config['args'][0]) / tbf.FrameSize
    antpols = len(antennas)
    
    
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbf.readFrame(fh)
    fh.seek(0)
    beginJD = astro.unix_to_utcjd(junkFrame.getTime())
    beginDate = ephem.Date(astro.unix_to_utcjd(junkFrame.getTime()) - astro.DJD_OFFSET)
    
    # Figure out how many frames there are per observation and the number of
    # channels that are in the file
    nFramesPerObs = tbf.getFramesPerObs(fh)
    nChannels = tbf.getChannelCount(fh)
    nSamples = 7840
    
    # Figure out how many chunks we need to work with
    nChunks = nFrames / nFramesPerObs
    
    # Pre-load the channel mapper
    mapper = []
    for i in xrange(2*nFramesPerObs):
        cFrame = tbf.readFrame(fh)
        if cFrame.header.firstChan not in mapper:
            mapper.append( cFrame.header.firstChan )
    fh.seek(-2*nFramesPerObs*tbf.FrameSize, 1)
    mapper.sort()
    
    # Calculate the frequencies
    freq = numpy.zeros(nChannels)
    for i,c in enumerate(mapper):
        freq[i*12:i*12+12] = c + numpy.arange(12)
    freq *= 25e3
    
    nInts = nFrames / nFramesPerObs
    
    # Get valid stands for both polarizations
    goodX = []
    goodY = []
    for i in xrange(len(antennas)):
        ant = antennas[i]
        if ant.getStatus() != 33 and not config['all']:
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
    print "Found %i good stands to use" % (len(good)/2,)
    for i in good:
        print "%3i, %i" % (antennas[i].stand.id, antennas[i].pol)
        
    # Number of frames to read in at once and average
    nFrames = nFrames
    nSets = 1
    
    centralFreq = freq.mean()
    
    print "Samples per observations: %i per pol." % (nFramesPerObs/2)
    print "Tuning frequency: %.3f Hz" % centralFreq
    print "Captures in file: %i (%.1f s)" % (nInts, nInts*40e-6)
    print "=="
    print "Station: %s" % station.name
    print "Date observed: %s" % beginDate
    print "Julian day: %.5f" % beginJD
    print "Integration Time: %.3f s" % (40e-6*nFrames/nFramesPerObs)
    print "Number of integrations in file: %i" % nSets
    
    # Loop over junks of 300 integrations to make sure that we don't overflow 
    # the FITS IDI memory buffer
    s = 0
    leftToDo = 1
    basename = os.path.split(filename)[1]
    basename, ext = os.path.splitext(basename)
    while leftToDo > 0:
        fitsFilename = "%s.FITS_%i" % (basename, (s+1),)
        
        if leftToDo > 100:
            chunk = 100
        else:
            chunk = leftToDo
            
        processChunk(fh, station, good, fitsFilename, pols=config['products'], ChunkSize=chunk)
                    
        s += 1
        leftToDo = leftToDo - chunk
        
    fh.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    