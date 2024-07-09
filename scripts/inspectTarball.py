#!/usr/bin/env python3

"""
Utility to display observation information contained in a MCS metadata
tarball.
"""

import sys
import pytz
import argparse

from datetime import datetime, timedelta

from lsl.common import stations
from lsl.astro import utcjd_to_unix, MJD_OFFSET
from lsl.common import metabundle
from lsl.common.sdf import get_observation_start_stop

from lsl.misc import telemetry
telemetry.track_script()


__version__ = "0.3"

# Date/time manipulation
_UTC = pytz.utc
_FORMAT_STRING = '%Y/%m/%d %H:%M:%S.%f %Z'


def main(args):
    # Filenames in an easier format
    inputTGZ  = args.filename
    
    # Parse the input file and get the dates of the observations.
    project = metabundle.get_sdf(inputTGZ)
    obsImpl = metabundle.get_observation_spec(inputTGZ)
    fileInfo = metabundle.get_session_metadata(inputTGZ)
    aspConfigB = metabundle.get_asp_configuration_summary(inputTGZ, which='Beginning')
    aspConfigE = metabundle.get_asp_configuration_summary(inputTGZ, which='End')
    mdStyle = metabundle.get_style(inputTGZ)
    site = metabundle.get_station(inputTGZ)
    if site is None:
        if mdStyle.endswith('ADP'):
            # Assume LWA-SV
            site = stations.lwasv
        elif mdStyle.endswith('NDP'):
            # Assume LWA-NA
            site = stations.lwana
        else:
            # Assume LWA1
            site = stations.lwa1
    observer = site.get_observer()
    
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
    print(f"Filename: {inputTGZ}")
    print(f" Project ID: {project.id}")
    print(f" Session ID: {project.sessions[0].id}")
    print(f" Metadata Style: {mdStyle}")
    print(f" Observations appear to start at {min(tStart).strftime(_FORMAT_STRING)}")
    print(f" -> LST at {site.name} for this date/time is {str(lst)}")
    
    lastDur = project.sessions[0].observations[nObs-1].dur
    lastDur = timedelta(seconds=int(lastDur/1000), microseconds=(lastDur*1000) % 1000000)
    sessionDur = max(tStart) - min(tStart) + lastDur
    
    print(" ")
    print(f" Total Session Duration: {sessionDur}")
    print(f" -> First observation starts at {min(tStart).strftime(_FORMAT_STRING)}")
    print(f" -> Last observation ends at {(max(tStart) + lastDur).strftime(_FORMAT_STRING)}")
    if project.sessions[0].observations[0].mode not in ('TBW', 'TBN'):
        drspec = 'No'
        if project.sessions[0].spcSetup[0] != 0 and project.sessions[0].spcSetup[1] != 0:
            drspec = 'Yes'
        drxBeam = project.sessions[0].drx_beam
        if drxBeam < 1:
            drxBeam = "MCS decides"
        else:
            drxBeam = str(drxBeam)
        print(f" DRX Beam: {drxBeam}")
        print(f" DR Spectrometer used? {drspec}")
        if drspec == 'Yes':
            print(f" -> {project.sessions[0].spcSetup[0]} channels, {project.sessions[0].spcSetup[1]} windows/integration")
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
        print(f" Obs. #{obsID}: {fileInfo[obsID]['tag']}")
    
    print(" ")
    print("ASP Configuration:")
    print('  Beginning')
    for k,v in aspConfigB.items():
        print(f"    {k}: {v}")
    print('  End')
    for k,v in aspConfigE.items():
        print(f"    {k}: {v}")
        
    print(" ")
    print(f" Number of observations: {nObs}")
    print(" Observation Detail:")
    for i in range(nObs):
        currDur = project.sessions[0].observations[i].dur
        currDur = timedelta(seconds=int(currDur/1000), microseconds=(currDur*1000) % 1000000)
        
        print(f"  Observation #{i+1}")
        currObs = None
        for j in range(len(obsImpl)):
            if obsImpl[j]['obs_id'] == i+1:
                currObs = obsImpl[j]
                break
                
        ## Basic setup
        print(f"   Target: {project.sessions[0].observations[i].target}")
        print(f"   Mode: {project.sessions[0].observations[i].mode}")
        if project.sessions[0].observations[i].mode == 'STEPPED':
            print("    Step Mode: %s" % ('RA/Dec' if project.sessions[0].observations[i].steps[0].is_radec else 'az/alt'))
            print(f"    Step Count: {len( project.sessions[0].observations[i].steps)}")
        print("   Start:")
        print(f"    MJD: {project.sessions[0].observations[i].mjd}")
        print(f"    MPM: {project.sessions[0].observations[i].mpm}")
        print(f"    -> {get_observation_start_stop(project.sessions[0].observations[i])[0].strftime(_FORMAT_STRING)}")
        print(f"   Duration: {currDur}")
        
        ## DP setup
        if project.sessions[0].observations[i].mode not in ('TBW',):
            print(f"   Tuning 1: {project.sessions[0].observations[i].frequency1/1e6:.3f} MHz")
        if project.sessions[0].observations[i].mode not in ('TBW', 'TBN'):
            print(f"   Tuning 2: {project.sessions[0].observations[i].frequency2/1e6:.3f} MHz")
        if project.sessions[0].observations[i].mode not in ('TBW',):
            print(f"   Filter code: {project.sessions[0].observations[i].filter}")
        if currObs is not None:
            if project.sessions[0].observations[i].mode not in ('TBW',):
                if project.sessions[0].observations[i].mode == 'TBN':
                    print(f"   Gain setting: {currObs['tbn_gain']}")
                else:
                    print(f"   Gain setting: {currObs['drx_gain']}")
        else:
            print("   WARNING: observation specification not found for this observation")
            
        ## Comments/notes
        print(f"   Observer Comments: {project.sessions[0].observations[i].comments}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='display information about an LWA metadata tarball', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('filename', type=str, 
                        help='metadata file to display')
    args = parser.parse_args()
    main(args)
