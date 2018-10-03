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
    nSamples = os.path.getsize(args[0]) / tbw.FRAME_SIZE
    print "Samples in file: ", nSamples
    fh = open(args[0], "rb", buffering=tbw.FRAME_SIZE)

    # Make sure that the data is TBW and determine the data length
    test = tbw.read_frame(fh)
    print "TBW Data:  %s" % test.header.is_tbw()
    if not test.header.is_tbw():
        raise errors.notTBWError()
    print "Data Length: %i bits" % test.get_data_bits()
    if test.header.get_data_bits() == 12:
        nData = 400
    else:
        nData = 1200
    fh.seek(-tbw.FRAME_SIZE, 1)
    
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
            cFrame = tbw.read_frame(fh, Verbose=False)
        except errors.EOFError:
            break
        except errors.SyncError:
            #print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/TBWFRAME_SIZE-1)
            syncCount = syncCount + 1
            continue
            
        stand = cFrame.header.parse_id()
        if cFrame.header.frame_count % 10000 == 0:
            print "%2i  %14i  %6.3f  %5i  %5i" % (stand, cFrame.data.timetag, cFrame.get_time(), cFrame.header.frame_count, cFrame.header.second_count)
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
