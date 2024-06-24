#!/usr/bin/env python

"""
Utility to display observation information contained in a MCS metadata
tarball.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import sys
import pytz
import argparse

from datetime import datetime, timedelta

from lsl.common import stations
from lsl.astro import utcjd_to_unix, MJD_OFFSET
from lsl.common import metabundle, metabundleADP
from lsl.common.sdf import get_observation_start_stop

from lsl.misc import telemetry
telemetry.track_script()


__version__ = "0.2"

# Date/time manipulation
_UTC = pytz.utc
_FORMAT_STRING = '%Y/%m/%d %H:%M:%S.%f %Z'


def main(args):
    # Get the site and observer
    site = stations.lwa1
    observer = site.get_observer()
    
    # Filenames in an easier format
    inputTGZ  = args.filename
    
    # Parse the input file and get the dates of the observations.  Be default 
    # this is for LWA1 but we switch over to LWA-SV if an error occurs.
    try:
        # LWA1
        project = metabundle.get_sdf(inputTGZ)
        obsImpl = metabundle.get_observation_spec(inputTGZ)
        fileInfo = metabundle.get_session_metadata(inputTGZ)
        aspConfigB = metabundle.get_asp_configuration_summary(inputTGZ, which='Beginning')
        aspConfigE = metabundle.get_asp_configuration_summary(inputTGZ, which='End')
    except:
        # LWA-SV
        ## Site changes
        site = stations.lwasv
        observer = site.get_observer()
        ## Try again
        project = metabundleADP.get_sdf(inputTGZ)
        obsImpl = metabundleADP.get_observation_spec(inputTGZ)
        fileInfo = metabundleADP.get_session_metadata(inputTGZ)
        aspConfigB = metabundleADP.get_asp_configuration_summary(inputTGZ, which='Beginning')
        aspConfigE = metabundleADP.get_asp_configuration_summary(inputTGZ, which='End')
        
    nObs = len(project.sessions[0].observations)
    tStart = [None,]*nObs
    for i in range(nObs):
        tStart[i]  = utcjd_to_unix(project.sessions[0].observations[i].mjd + MJD_OFFSET)
        tStart[i] += project.sessions[0].observations[i].mpm / 1000.0
        tStart[i]  = datetime.utcfromtimestamp(tStart[i])
        tStart[i]  = _UTC.localize(tStart[i])
    
    # Get the LST at the start
    observer.date = (min(tStart)).strftime('%Y/%m/%d %H:%M:%S')
    lst = observer.sidereal_time()
    
    # Report on the file
    print("Filename: %s" % inputTGZ)
    print(" Project ID: %s" % project.id)
    print(" Session ID: %i" % project.sessions[0].id)
    print(" Observations appear to start at %s" % (min(tStart)).strftime(_FORMAT_STRING))
    print(" -> LST at %s for this date/time is %s" % (site.name, lst))
    
    lastDur = project.sessions[0].observations[nObs-1].dur
    lastDur = timedelta(seconds=int(lastDur/1000), microseconds=(lastDur*1000) % 1000000)
    sessionDur = max(tStart) - min(tStart) + lastDur
    
    print(" ")
    print(" Total Session Duration: %s" % sessionDur)
    print(" -> First observation starts at %s" % min(tStart).strftime(_FORMAT_STRING))
    print(" -> Last observation ends at %s" % (max(tStart) + lastDur).strftime(_FORMAT_STRING))
    if project.sessions[0].observations[0].mode not in ('TBW', 'TBN'):
        drspec = 'No'
        if project.sessions[0].spcSetup[0] != 0 and project.sessions[0].spcSetup[1] != 0:
            drspec = 'Yes'
        drxBeam = project.sessions[0].drx_beam
        if drxBeam < 1:
            drxBeam = "MCS decides"
        else:
            drxBeam = "%i" % drxBeam
        print(" DRX Beam: %s" % drxBeam)
        print(" DR Spectrometer used? %s" % drspec)
        if drspec == 'Yes':
            print(" -> %i channels, %i windows/integration" % tuple(project.sessions[0].spcSetup))
    else:
        tbnCount = 0
        tbwCount = 0
        for obs in project.sessions[0].observations:
            if obs.mode == 'TBW':
                tbwCount += 1
            else:
                tbnCount += 1
        if tbwCount > 0 and tbnCount == 0:
            print(" Transient Buffer Mode: TBW")
        elif tbwCount == 0 and tbnCount > 0:
            print(" Transient Buffer Mode: TBN")
        else:
            print(" Transient Buffer Mode: both TBW and TBN")
    print(" ")
    print("File Information:")
    for obsID in fileInfo.keys():
        print(" Obs. #%i: %s" % (obsID, fileInfo[obsID]['tag']))
    
    print(" ")
    print("ASP Configuration:")
    print('  Beginning')
    for k,v in aspConfigB.items():
        print('    %s: %i' % (k, v))
    print('  End')
    for k,v in aspConfigE.items():
        print('    %s: %i' % (k, v))
        
    print(" ")
    print(" Number of observations: %i" % nObs)
    print(" Observation Detail:")
    for i in range(nObs):
        currDur = project.sessions[0].observations[i].dur
        currDur = timedelta(seconds=int(currDur/1000), microseconds=(currDur*1000) % 1000000)
        
        print("  Observation #%i" % (i+1,))
        currObs = None
        for j in range(len(obsImpl)):
            if obsImpl[j]['obs_id'] == i+1:
                currObs = obsImpl[j]
                break
                
        ## Basic setup
        print("   Target: %s" % project.sessions[0].observations[i].target)
        print("   Mode: %s" % project.sessions[0].observations[i].mode)
        if project.sessions[0].observations[i].mode == 'STEPPED':
            print("    Step Mode: %s" % ('RA/Dec' if project.sessions[0].observations[i].steps[0].is_radec else 'az/alt'))
            print("    Step Count: %i" % len( project.sessions[0].observations[i].steps))
        print("   Start:")
        print("    MJD: %i" % project.sessions[0].observations[i].mjd)
        print("    MPM: %i" % project.sessions[0].observations[i].mpm)
        print("    -> %s" % get_observation_start_stop(project.sessions[0].observations[i])[0].strftime(_FORMAT_STRING))
        print("   Duration: %s" % currDur)
        
        ## DP setup
        if project.sessions[0].observations[i].mode not in ('TBW',):
            print("   Tuning 1: %.3f MHz" % (project.sessions[0].observations[i].frequency1/1e6,))
        if project.sessions[0].observations[i].mode not in ('TBW', 'TBN'):
            print("   Tuning 2: %.3f MHz" % (project.sessions[0].observations[i].frequency2/1e6,))
        if project.sessions[0].observations[i].mode not in ('TBW',):
            print("   Filter code: %i" % project.sessions[0].observations[i].filter)
        if currObs is not None:
            if project.sessions[0].observations[i].mode not in ('TBW',):
                if project.sessions[0].observations[i].mode == 'TBN':
                    print("   Gain setting: %i" % currObs['tbn_gain'])
                else:
                    print("   Gain setting: %i" % currObs['drx_gain'])
        else:
            print("   WARNING: observation specification not found for this observation")
            
        ## Comments/notes
        print("   Observer Comments: %s" % project.sessions[0].observations[i].comments)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='display information about an LWA metadata tarball', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('filename', type=str, 
                        help='metadata file to display')
    args = parser.parse_args()
    main(args)
    
