#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import time
import getopt
from datetime import datetime

from lsl.reader import drx, errors
from lsl.common.progress import ProgressBar


def usage(exitCode=None):
    print """splitDRX.py - split a DRX file containing multiple seconds into
several files

Usage: splitDRX.py [OPTIONS] file

Options:
-h, --help             Display this help information
-c, --count            Number of seconds to keep
-o, --offset           Number of seconds to skip before splitting
-d, --date             Label the split files with a date rather than a 
                       sequence number
-r, --recurvsive       Recursively split the file
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseConfig(args):
    config = {}
    # Command line flags - default values
    config['offset'] = 0
    config['count'] = 0
    config['date'] = False
    config['recursive'] = False
    
    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "hc:o:dr", ["help", "count=", "offset=", "date", "recursive"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-c', '--count'):
            config['count'] = float(value)
        elif opt in ('-o', '--offset'):
            config['offset'] = float(value)
        elif opt in ('-d', '--date'):
            config['date'] = True
        elif opt in ('-r', '--recursive'):
            config['recursive'] = True
        else:
            assert False
            
    # Add in arguments
    config['args'] = arg
    
    # Return configuration
    return config


def fileSplitFunction(fhIn, fhOut, nCaptures, nBeampols):
    pb = ProgressBar(max=nCaptures)
    
    for c in xrange(int(nCaptures)):
        for i in xrange(nBeampols):
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
    config = parseConfig(args)
    filename = config['args'][0]

    sizeB = os.path.getsize(filename)

    # Open the file and get some basic info about the data contained
    fh = open(filename, 'rb')
    nFramesFile = sizeB / drx.FRAME_SIZE

    while True:
        try:
            junkFrame = drx.read_frame(fh)
            try:
                srate = junkFrame.sample_rate
                ti0, tf0 = junkFrame.time
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
    offset = int(round(config['offset'] * srate / 4096 * beampols))
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
        ti1, tf1 = junkFrame.time
        tunepols = drx.get_frames_per_obs(fh)
        tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepol
        fh.seek(-drx.FRAME_SIZE, 1)
        
        ## See how far off the current frame is from the target
        tDiff = ti1 - (ti0 + config['offset']) + tf1 - tf0
        
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
    config['offset'] = t1 - t0
    
    nCaptures = nFramesFile/beampols

    print "Filename:     %s" % filename
    print "Size:         %.1f MB" % (float(sizeB)/1024/1024)
    print "Captures:     %i (%.2f seconds)" % (nCaptures, nCaptures*4096/srate)
    print "Tuning/Pols.: %i " % tunepol
    print "Sample Rate: %.2f MHz" % (srate/1e6)
    print "==="

    if config['count'] > 0:
        nCaptures = config['count'] * srate / 4096
    else:
        nCaptures -= config['offset'] * srate / 4096
        
        config['count'] = nCaptures * 4096 / srate

    print "Seconds to Skip:  %.2f (%i captures)" % (config['offset'], offset/beampols)
    print "Seconds to Split: %.2f (%i captures)" % (config['count'], nCaptures)

    # Make sure that the first frame in the file is the first frame of a capture 
    # (tuning 1, pol 0).  If not, read in as many frames as necessary to get to 
    # the beginning of a complete capture.
    beam, tune, pol = junkFrame.id

    skip = 0
    while (2*(tune-1)+pol) != 0:
        frame = drx.read_frame(fh)
        beam, tune, pol = frame.id
        skip += 1

    nFramesRemaining = (sizeB - fh.tell()) / drx.FRAME_SIZE
    nRecursions = int(nFramesRemaining / (nCaptures*beampols))
    if not config['recursive']:
        nRecursions = 1
        
    scale = int(math.log10(nRecursions)) + 1
    ifString = "Working on #%%%ii of %i (%%s)" % (scale, nRecursions)
    
    for r in xrange(nRecursions):
        if config['date']:
            filePos = fh.tell()
            junkFrame = drx.read_frame(fh)
            fh.seek(filePos)

            dt = datetime.utcfromtimestamp(sum(junkFrame.time))
            captFilename = "%s_%s.dat" % (os.path.splitext(os.path.basename(filename))[0], dt.isoformat())
        else:
            captFilename = "%s_s%04i_p%%0%ii.dat" % (os.path.splitext(os.path.basename(filename))[0], config['count'], scale)
            captFilename = captFilename % r
            if not config['recursive']:
                captFilename = "%s_s%04i.dat" % (os.path.splitext(os.path.basename(filename))[0], config['count'])
            
        print ifString % (r+1, captFilename)
        
        t0 = time.time()
        fhOut = open(captFilename, 'wb')
        fileSplitFunction(fh, fhOut, nCaptures, beampols)
        fhOut.close()
        t1 = time.time()
        print "  Copied %i bytes in %.3f s (%.3f MB/s)" % (os.path.getsize(captFilename), t1-t0, os.path.getsize(captFilename)/1024.0**2/(t1-t0))
        
    fh.close()
    
    

if __name__ == "__main__":
    main(sys.argv[1:])
