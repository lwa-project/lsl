#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in DRX data and writing it to a TS-FITS file."""

import sys
import time
import ephem

from lsl.reader.ldp import LWA1DataFile
from lsl.writer import tsfits
from lsl.astro import unix_to_utcjd, DJD_OFFSET


def main(args):
    idf = LWA1DataFile(args[0])
    nFramesFile = idf.get_info('nFrames')
    
    srate = idf.get_info('sample_rate')
    beam = idf.get_info('beam')
    beampols = idf.get_info('beampols')
    
    # Date
    beginDate = ephem.Date(unix_to_utcjd(idf.get_info('tStart')) - DJD_OFFSET)
    
    # File summary
    print "Filename: %s" % args[0]
    print "Date of First Frame: %s" % str(beginDate)
    print "Beam: %i" % beam
    print "Tune/Pols: %ii" % beampols
    print "Sample Rate: %i Hz" % srate
    print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate)
    print "---"

    tStart = time.time()
    
    # Create a new FITS file with the name 'drx-tsfits.fits'
    fitsFile = tsfits.TBN('drx-tsfits.fits')
    
    nSamples = 3400
    
    count = {}
    masterCount = 0
    for i in xrange(nSamples):
        for j in xrange(beampols):
            frame = idf.read_frame()
            beam, pol, tune = frame.parse_id()
            try:
                count[beam] += 1
            except KeyError:
                count[beam] = 1
                
            fitsFile.addStandData(frame)
            
        masterCount = masterCount + 1
        
    tEnd = time.time()
    print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (beampols*nSamples, (tEnd-tStart), beampols*nSamples/(tEnd-tStart))
    
    idf.close()
    fitsFile.close()
    fitsFile.info()
    
    # Summary information about the file that was just read in
    print "Summary:"
    for beam in sorted(count.keys()):
        print "Beam: %2i, Frames: %5i" % (beam, count[beam])


if __name__ == "__main__":
    main(sys.argv[1:])
