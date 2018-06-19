#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for reading in TBW data and writing it to a TS-FITS file."""

import sys
import time
from lsl.reader import tbw
from lsl.reader import errors
from lsl.writer import tsfits

import matplotlib.pyplot as plt


def main(args):
    # Determine the number of samples in the specified file
    nSamples = os.path.getsize(args[0]) / tbw.FrameSize
    print "Samples in file: ", nSamples
    fh = open(args[0], "rb", buffering=tbw.FrameSize)

    # Make sure that the data is TBW and determine the data length
    test = tbw.readFrame(fh)
    print "TBW Data:  %s" % test.header.isTBW()
    if not test.header.isTBW():
        raise errors.notTBWError()
    print "Data Length: %i bits" % test.getDataBits()
    if test.header.getDataBits() == 12:
        nData = 400
    else:
        nData = 1200
    fh.seek(-tbw.FrameSize, 1)
    
    # Due to the size of the FITS files being generated, the number of frames that 
    # can be read in is limited to 300,000, or 30,000 frames for 10 stands.  Getting
    # around this limit will require finding out how to do on-the-fly FITS binary 
    # table resizing.  
    nSamples = 900000
    
    tStart = time.time()
    
    # Create a new FITS file with the name 'tbw.fits'
    fitsFile = tsfits.TBW('tbw-tsfits-test.fits')

    # Read in the data and add it to the FITS file created above
    count = {}
    syncCount = 0
    masterCount = 0
    for i in range(nSamples):
        # Read in the next frame and anticipate any problems that could occur
        try:
            cFrame = tbw.readFrame(fh, Verbose=False)
        except errors.eofError:
            break
        except errors.syncError:
            #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFrameSize-1)
            syncCount = syncCount + 1
            continue
            
        stand = cFrame.header.parseID()
        if cFrame.header.frameCount % 10000 == 0:
            print "%2i  %14i  %6.3f  %5i  %5i" % (stand, cFrame.data.timeTag, cFrame.getTime(), cFrame.header.frameCount, cFrame.header.secondsCount)
        if stand not in count.keys():
            count[stand] = 0
            
        fitsFile.addStandData(cFrame)
        
        count[stand] = count[stand] + 1
        masterCount = masterCount + 1
        
    tEnd = time.time()
    print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (masterCount, (tEnd-tStart), masterCount/(tEnd-tStart))
    
    fh.close()
    fitsFile.close()
    fitsFile.info()
    
    # Summary information about the file that was just read in
    print "Summary:"
    for stand in sorted(count.keys()):
        print "Stand: %2i, Frames: %5i" % (stand, count[stand])
    print "Sync Errors: %5i" % syncCount


if __name__ == "__main__":
    main(sys.argv[1:])
