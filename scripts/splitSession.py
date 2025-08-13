#!/usr/bin/env python3

"""
Script for splitting a DRX file based on the information contained in a MCS
metadata tarball.
"""

import os
import sys
import time
import argparse
from datetime import datetime, timedelta

from lsl.common.mcs import mjdmpm_to_datetime
from lsl.common import metabundle
from lsl.reader import tbs, drx, errors

from lsl.misc import telemetry
telemetry.track_script()


def _obs_comp(x, y):
    """
    Function to help sort observations in time.
    """
    
    tX = mjdmpm_to_datetime(x['mjd'], x['mpm'])
    tY = mjdmpm_to_datetime(y['mjd'], y['mpm'])
    if tX < tY:
        return -1
    elif tX > tY:
        return 1
    else:
        return 0


def main(args):
    # Get the file names
    meta = args.metadata
    data = args.filename

    # Get all observations and their start times
    sdf = metabundle.get_sdf(meta)
    ses = metabundle.get_session_spec(meta)
    obs = metabundle.get_observation_spec(meta)
    obs.sort(_obs_comp)
    tStart = []
    oDetails = []
    for i,o in enumerate(obs):
        tStart.append( mjdmpm_to_datetime(o['mjd'], o['mpm']) )
        oDetails.append( {'m': o['mode'], 'd': o['dur'] / 1000.0, 'f': o['bw'], 
                          'p': o['project_id'], 's': o['session_id'], 'o': o['obs_id'], 
                          't': sdf.sessions[0].observations[o['obs_id']-1].target} )

        print(f"Observation #{o['obs_id']}")
        print(f" Start: {o['mjd']}, {o['mpm']} -> {tStart[-1]}")
        print(f" Mode: {o['mode'].name}")
        print(f" BW: {o['bw']}")
        print(f" Target: {sdf.sessions[0].observations[o['obs_id']-1].target}")
    print(" ")

    # Figure out where in the file the various bits are.
    fh = open(data, 'rb')
    lf = drx.read_frame(fh)
    beam, j, k = lf.id
    if beam != obs[0]['drx_beam']:
        print(f"ERROR: Beam mis-match, metadata is for #{obs[0]['drx_beam']}, file is for #{beam}")
        sys.exit()
    firstFrame = lf.time.datetime
    if abs(firstFrame - min(tStart)) > timedelta(seconds=30):
        print(f"ERROR: Time mis-match, metadata is for {min(tStart)}, file is for {firstFrame}")
        sys.exit()
    fh.seek(0)

    for i in range(len(tStart)):
        eof = False

        ## Get observation properties
        oStart = tStart[i]
        oMode = oDetails[i]['m'].name
        oDur  = oDetails[i]['d']
        oBW   = oDetails[i]['f']
        print(f"Seeking {oMode} observation of {oDur:.3f} seconds at {oStart}")

        ## Get the correct reader to use
        if oMode == 'TBS':
            reader = tbs
            reader.FRAME_SIZE = tbs.get_frame_size(fh)
            bwKey = {8: 200e3}
            bwMult = 1.0
            fCount = 8
        else:
            reader = drx
            bwKey = drx.FILTER_CODES
            bwMult = 4.0 / 4096
            fCount = 4096

        ## Jump ahead to where the next frame should be, if needed
        if i != 0:
            pDur  = oDetails[i-1]['d']
            pBW   = oDetails[i-1]['f']

            nFramesSkip = int(pDur*bwKey[pBW]*bwMult)
            fh.seek(nFramesSkip*reader.FRAME_SIZE, 1)
            if fh.tell() >= os.path.getsize(data):
                fh.seek(-10*reader.FRAME_SIZE, 2)
                
        ## Figure out where we are and make sure we line up on a frame
        ## NOTE: This should never be needed
        fail = True
        while fail:
            try:
                frame = reader.read_frame(fh)
                fail = False
            except errors.SyncError:
                fh.seek(1, 1)
            except errors.EOFError:
                break
        fh.seek(-reader.FRAME_SIZE, 1)	

        ## Go in search of the start of the observation
        if frame.time.datetime < oStart:
            ### We aren't at the beginning yet, seek fowards
            print(f"-> At byte {fh.tell()}, time is {frame.time.datetime} < {oStart}")

            while frame.time.datetime < oStart:
                try:
                    frame = reader.read_frame(fh)
                except errors.SyncError:		
                    fh.seek(1, 1)
                except errors.EOFError:
                    break
                #print(frame.time.datetime, oStart)

        elif frame.time.datetime > oStart:
            ### We've gone too far, seek backwards
            print(f"-> At byte {fh.tell()}, time is {frame.time.datetime} > {oStart}")

            while frame.time.datetime > oStart:
                if fh.tell() == 0:
                    break
                fh.seek(-2*reader.FRAME_SIZE, 1)
                try:
                    frame = reader.read_frame(fh)
                except errors.SyncError:		
                    fh.seek(-1, 1)
                except errors.EOFError:
                    break
                #print(frame.time.datetime, oStart)
                
        else:
            ### We're there already
            print(f"-> At byte {fh.tell()}, time is {frame.time.datetime} = {oStart}")
            
        ## Jump back exactly one frame so that the filehandle is in a position 
        ## to read the first frame that is part of the observation
        try:
            frame = reader.read_frame(fh)
            print(f"-> At byte {fh.tell()}, time is {frame.time.datetime} = {oStart}")
            fh.seek(-reader.FRAME_SIZE, 1)
        except errors.EOFError:
            pass
            
        ## Update the bytes ranges
        if fh.tell() < os.path.getsize(data):
            oDetails[i]['b'] = fh.tell()
            oDetails[i]['e'] = -1
        else:
            oDetails[i]['b'] = -1
            oDetails[i]['e'] = -1

        if i != 0:
            oDetails[i-1]['e'] = fh.tell()

        ## Progress report
        if oDetails[i]['b'] >= 0:
            print(f"-> Obs. {oDetails[i]['o']} starts at byte {oDetails[i]['b']}")
        else:
            print(f"-> Obs. {oDetails[i]['o']} starts after the end of the file")
    print(" ")

    # Report
    for i in range(len(tStart)):
        if oDetails[i]['b'] < 0:
            print(f"{oDetails[i]['p']}, Session {oDetails[i]['s']}, Observation {oDetails[i]['o']}: not found")

        else:
            print(f"{oDetails[i]['p']}, Session {oDetails[i]['s']}, Observation {oDetails[i]['o']}: {oDetails[i]['b']} to {oDetails[i]['e']} ({oDetails[i]['e']-oDetails[i]['b']} bytes)")
    print(" ")

    # Split
    if not args.list:
        for i in range(len(tStart)):
            if oDetails[i]['b'] < 0:
                continue
                
            ## Report
            print(f"Working on Observation #{i+1}")
            
            ## Create the output name
            if args.source:
                outname = '%s_%i_%s.dat' % (oDetails[i]['p'], oDetails[i]['s'], oDetails[i]['t'].replace(' ', '').replace('/','').replace('&','and'))
            else:
                outname = f"{oDetails[i]['p']}_{oDetails[i]['s']}_{oDetails[i]['o']}.dat"
                
            oMode = oDetails[i]['m'].name

            ## Get the correct reader to use
            if oMode == 'TBS':
                reader = tbs
            else:
                reader = drx

            ## Get the number of frames
            if oDetails[i]['e'] > 0:
                nFramesRead = (oDetails[i]['e'] - oDetails[i]['b']) // reader.FRAME_SIZE
            else:
                nFramesRead = (os.path.getsize(data) - oDetails[i]['b']) // reader.FRAME_SIZE

            ## Split
            if os.path.exists(outname):
                if not args.force:
                    yn = input(f"WARNING: '{outname}' exists, overwrite? [Y/n] ")
                else:
                    yn = 'y'
                    
                if yn not in ('n', 'N'):
                    os.unlink(outname)
                else:
                    print(f"WARNING: output file '{outname}' already exists, skipping")
                    continue
                    
            fh.seek(oDetails[i]['b'])
            
            t0 = time.time()
            with open(outname, 'wb') as oh:
                for sl in [2**i for i in range(17)[::-1]]:
                    while nFramesRead >= sl:
                        temp = fh.read(sl*reader.FRAME_SIZE)
                        oh.write(temp)
                        nFramesRead -= sl
            t1 = time.time()
            print(f"  Copied {os.path.getsize(outname)} bytes in {t1-t0:.3f} s ({os.path.getsize(outname)/1024.0**2/(t1-t0):.3f} MB/s)")
    print(" ")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='given a MCS metadata tarball and a session data file, split the data file into individual observations', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('metadata', type=str, 
                        help='metadata file for the observation')
    parser.add_argument('filename', type=str, 
                        help='data file for the observation')
    parser.add_argument('-l', '--list', action='store_true', 
                        help='list source names')
    parser.add_argument('-s', '--source', action='store_true', 
                        help='split by source name instead of observation ID')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing split files')
    args = parser.parse_args()

    main(args)
