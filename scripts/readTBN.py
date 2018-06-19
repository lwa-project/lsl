#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in TBN data and writing it to a TS-FITS file."""

import sys
import time
import ephem

from lsl.reader.ldp import LWA1DataFile
from lsl.writer import tsfits
from lsl.astro import unix_to_utcjd, DJD_OFFSET


def main(args):
    idf = LWA1DataFile(args[0])
    print idf.getInfo()
    nFramesFile = idf.getInfo('nFrames')
    
    srate = idf.getInfo('sampleRate')
    antpols = idf.getInfo('nAntenna')
    
    # Date
    beginDate = ephem.Date(unix_to_utcjd(idf.getInfo('tStart')) - DJD_OFFSET)
    
    # File summary
    print "Filename: %s" % args[0]
    print "Date of First Frame: %s" % str(beginDate)
    print "Ant/Pols: %i" % antpols
    print "Sample Rate: %i Hz" % srate
    print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
    print "---"
    
    tStart = time.time()
    
    # Create a new FITS file with the name 'tbn-tsfits.fits'
    fitsFile = tsfits.TBN('tbn-tsfits.fits')
    
    nSamples = 2
    
    count = {}
    masterCount = 0
    for i in xrange(nSamples):
        for j in xrange(antpols):
            frame = idf.readFrame()
            stand, pol = frame.parseID()
            try:
                count[stand] += 1
            except KeyError:
                count[stand] = 1
                
            fitsFile.addStandData(frame)
            
        masterCount = masterCount + 1
        
    tEnd = time.time()
    print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (antpols*nSamples, (tEnd-tStart), antpols*nSamples/(tEnd-tStart))
    
    idf.close()
    fitsFile.close()
    fitsFile.info()
    
    # Summary information about the file that was just read in
    print "Summary:"
    for stand in sorted(count.keys()):
        print "Stand: %2i, Frames: %5i" % (stand, count[stand])


if __name__ == "__main__":
    main(sys.argv[1:])
