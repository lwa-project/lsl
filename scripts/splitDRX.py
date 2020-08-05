#!/usr/bin/env python

"""
Example script for splitting a DRX file into smaller pieces.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import math
import time
import argparse
from datetime import datetime

from lsl.reader import drx, errors
from lsl.common.progress import ProgressBar
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


def split_file(fhIn, fhOut, nCaptures, nBeampols):
    pb = ProgressBar(max=nCaptures)
    
    for c in range(int(nCaptures)):
        for i in range(nBeampols):
            cFrame = fhIn.read(drx.FRAME_SIZE)
            fhOut.write(cFrame)
            
        pb.inc(amount=1)
        if c != 0 and c % 100 == 0:
            sys.stdout.write(pb.show()+'\r')
            sys.stdout.flush()
            
    sys.stdout.write(pb.show()+'\r')
    sys.stdout.write('\n')
    sys.stdout.flush()


def main(args):
    filename = args.filename

    sizeB = os.path.getsize(filename)

    # Open the file and get some basic info about the data contained
    fh = open(filename, 'rb')
    nFramesFile = sizeB / drx.FRAME_SIZE

    while True:
        try:
            junkFrame = drx.read_frame(fh)
            try:
                srate = junkFrame.sample_rate
                t0 = junkFrame.time
                break
            except ZeroDivisionError:
                pass
        except errors.SyncError:
            fh.seek(-drx.FRAME_SIZE+1, 1)
            
    fh.seek(-drx.FRAME_SIZE, 1)
    
    beams = drx.get_beam_count(fh)
    tunepols = drx.get_frames_per_obs(fh)
    tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
    beampols = tunepol

    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(args.offset * srate / 4096 * beampols))
    offset = int(1.0 * offset / beampols) * beampols
    fh.seek(offset*drx.FRAME_SIZE, 1)
    
    # Iterate on the offsets until we reach the right point in the file.  This
    # is needed to deal with files that start with only one tuning and/or a 
    # different sample rate.  
    while True:
        ## Figure out where in the file we are and what the current tuning/sample 
        ## rate is
        junkFrame = drx.read_frame(fh)
        srate = junkFrame.sample_rate
        t1 = junkFrame.time
        tunepols = drx.get_frames_per_obs(fh)
        tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepol
        fh.seek(-drx.FRAME_SIZE, 1)
        
        ## See how far off the current frame is from the target
        tDiff = t1 - (t0 + args.offset)
        
        ## Half that to come up with a new seek parameter
        tCorr = -tDiff / 2.0
        cOffset = int(tCorr * srate / 4096 * beampols)
        cOffset = int(1.0 * cOffset / beampols) * beampols
        offset += cOffset
        
        ## If the offset is zero, we are done.  Otherwise, apply the offset
        ## and check the location in the file again/
        if cOffset is 0:
            break
        fh.seek(cOffset*drx.FRAME_SIZE, 1)
    
    # Update the offset actually used
    args.offset = t1 - t0
    
    nCaptures = nFramesFile/beampols

    print("Filename:     %s" % filename)
    print("Size:         %.1f MB" % (float(sizeB)/1024/1024))
    print("Captures:     %i (%.2f seconds)" % (nCaptures, nCaptures*4096/srate))
    print("Tuning/Pols.: %i " % tunepol)
    print("Sample Rate: %.2f MHz" % (srate/1e6))
    print("===")

    if args.count > 0:
        nCaptures = args.count * srate // 4096
    else:
        nCaptures -= args.offset * srate // 4096
        
        args.count = nCaptures * 4096 // srate

    print("Seconds to Skip:  %.2f (%i captures)" % (args.offset, offset/beampols))
    print("Seconds to Split: %.2f (%i captures)" % (args.count, nCaptures))

    # Make sure that the first frame in the file is the first frame of a capture 
    # (tuning 1, pol 0).  If not, read in as many frames as necessary to get to 
    # the beginning of a complete capture.
    beam, tune, pol = junkFrame.id

    skip = 0
    while (2*(tune-1)+pol) != 0:
        frame = drx.read_frame(fh)
        beam, tune, pol = frame.id
        skip += 1

    nFramesRemaining = (sizeB - fh.tell()) // drx.FRAME_SIZE
    nRecursions = int(nFramesRemaining // (nCaptures*beampols))
    if not args.recursive:
        nRecursions = 1
        
    scale = int(math.log10(nRecursions)) + 1
    ifString = "Working on #%%%ii of %i (%%s)" % (scale, nRecursions)
    
    for r in range(nRecursions):
        if args.date:
            filePos = fh.tell()
            junkFrame = drx.read_frame(fh)
            fh.seek(filePos)

            dt = junkFrame.time.datetime
            captFilename = "%s_%s.dat" % (os.path.splitext(os.path.basename(filename))[0], dt.isoformat())
        else:
            captFilename = "%s_s%04i_p%%0%ii.dat" % (os.path.splitext(os.path.basename(filename))[0], args.count, scale)
            captFilename = captFilename % r
            if not args.recursive:
                captFilename = "%s_s%04i.dat" % (os.path.splitext(os.path.basename(filename))[0], args.count)
            
        print(ifString % (r+1, captFilename))
        
        t0 = time.time()
        fhOut = open(captFilename, 'wb')
        split_file(fh, fhOut, nCaptures, beampols)
        fhOut.close()
        t1 = time.time()
        print("  Copied %i bytes in %.3f s (%.3f MB/s)" % (os.path.getsize(captFilename), t1-t0, os.path.getsize(captFilename)/1024.0**2/(t1-t0)))
        
    fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='split a DRX file into several files', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to split')
    parser.add_argument('-c', '--count', type=aph.positive_float, default=10.0, 
                        help='number of seconds to keep')
    parser.add_argument('-o', '--offset', type=aph.positive_or_zero_float, default=0.0, 
                        help='number of seconds to skip before splitting')
    parser.add_argument('-d', '--date', action='store_true', 
                        help='label the split files with a date rather than a sequence number')
    parser.add_argument('-r', '--recursive', action='store_true', 
                        help='recursively split the file')
    args = parser.parse_args()
    main(args)
    
