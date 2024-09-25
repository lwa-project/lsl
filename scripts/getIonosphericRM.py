#!/usr/bin/env python3

"""
Utility for estimating the ionospheric contribution to the DM and RM for
a given point on the sky at a given time.
"""

import sys
import ephem
import numpy as np
import argparse
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.common.mcs import datetime_to_mjdmpm
from lsl.misc import ionosphere
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    # Inputs
    if args.file is not None:
        mjdList = np.loadtxt(args.file)
        mjdList = mjdList.ravel()
        
        mjdList = np.sort(mjdList)
        
    else:
        tStart = f"{args.StartDate} {args.StartTime}"
        tStop = f"{args.StopDate} {args.StopTime}"
        
        # YYYY/MM/DD HH:MM:SS -> datetime instance
        tStart = datetime.strptime(tStart, "%Y/%m/%d %H:%M:%S.%f")
        tStop = datetime.strptime(tStop, "%Y/%m/%d %H:%M:%S.%f")
        
        # datetime instance to MJD
        mjd,mpm = datetime_to_mjdmpm(tStart)
        mjdStart = mjd + mpm/1000.0/86400.0
        mjd,mpm = datetime_to_mjdmpm(tStop)
        mjdStop = mjd + mpm/1000.0/86400.0
        
        mjdList = np.linspace(mjdStart, mjdStop, args.n_samples)
    
    # Setup everything for computing the position of the source
    if args.lwasv:
        site = stations.lwasv
    elif args.lwana:
        site = stations.lwasv
    elif args.ovrolwa:
        site = stations.lwa1
        site.lat, site.lon, site.elev = ('37.23977727', '-118.2816667', 1183.48)
    else:
        site = stations.lwa1
    obs = site.get_observer()
    bdy = ephem.FixedBody()
    bdy._ra = args.RA
    bdy._dec = args.Dec
    
    # Setup the ionospheric model source
    if args.igs:
        mtype = 'IGS'
    elif args.jpl:
        mtype = 'JPL'
    elif args.emr:
        mtype = 'EMR'
    elif args.code:
        mtype = 'CODE'
    elif args.ustec:
        mtype = 'USTEC'
    elif args.uqr:
        mtype = 'UQR'
    elif args.glotec:
        mtype = 'GloTEC'
    else:
        mtype = 'IGS'
        
    # Go!
    print("%-13s  %-6s  %-6s  %-21s  %-15s" % ("MJD", "Az.", "El.", "DM [pc/cm^3]", "RM [1/m^2]"))
    print("-"*(13+2+6+2+6+2+21+2+15))
    for mjd in mjdList:
        # Set the date and compute the location of the target
        obs.date = mjd + astro.MJD_OFFSET - astro.DJD_OFFSET
        bdy.compute(obs)
        az = bdy.az*180/np.pi
        el = bdy.alt*180/np.pi
        
        if el > 0:
            # Get the latitude, longitude, and height of the ionospheric pierce 
            # point in this direction
            ippLat, ippLon, ippElv = ionosphere.get_ionospheric_pierce_point(site, az, el)
            
            # Load in the TEC value and the RMS from the IGS final data product
            tec, rms = ionosphere.get_tec_value(mjd, lat=ippLat, lng=ippLon, include_rms=True, type=mtype)
            tec, rms = tec[0][0], rms[0][0]
            
            # Use the IGRF to compute the ECEF components of the Earth's magnetic field
            Bx, By, Bz = ionosphere.get_magnetic_field(ippLat, ippLon, ippElv, mjd=mjd, ecef=True)
            
            # Rotate the ECEF field into topocentric coordinates so that we can 
            # get the magnetic field along the line of sight
            rot = np.array([[ np.sin(site.lat)*np.cos(site.long), np.sin(site.lat)*np.sin(site.long), -np.cos(site.lat)], 
                            [-np.sin(site.long),                     np.cos(site.long),                      0                  ],
                            [ np.cos(site.lat)*np.cos(site.long), np.cos(site.lat)*np.sin(site.long),  np.sin(site.lat)]])
            ## ECEF -> SEZ
            sez = np.dot(rot, np.array([Bx, By, Bz]))
            ## SEZ -> NEZ
            enz = 1.0*sez[[1,0,2]]
            enz[1] *= -1.0
            
            # Compute the pointing vector for this direction and use that to get
            # B parallel.  Note that we need a negative sign when we dot to get
            # the direction of propagation right.
            pnt = np.array([np.cos(el*np.pi/180)*np.sin(az*np.pi/180),
                            np.cos(el*np.pi/180)*np.cos(az*np.pi/180), 
                            np.sin(el*np.pi/180)])
            Bparallel = -np.dot(pnt, enz)
            
            # Compute the dispersion measure and the RMS
            DM    = 3.24078e-23 * (tec*1e16)
            rmsDM = 3.24078e-23 * (rms*1e16)
            
            # Compute the rotation measure and the RMS
            RM    = 2.62e-13 * (tec*1e16) * (Bparallel*1e-9)
            rmsRM = 2.62e-13 * (rms*1e16) * (Bparallel*1e-9)
            
            # Report
            print("%013.6f  %6.2f  %6.2f  %8.6f +/- %8.6f  %5.3f +/- %5.3f" % (mjd, az, el, DM, rmsDM, RM, rmsRM))
        else:
            # Write out dummy values since the source isn't up
            print("%013.6f  %6.2f  %6.2f  %8s +/- %8s  %5s +/- %5s" % (mjd, az, el, '---', '---', '---', '---'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='estimate the ionospheric contribution to the RM for an observation using the IGS final product and the IGRF', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('RA', type=aph.hours, 
                        help='J2000 right ascension in HH:MM:SS[.SSS]')
    parser.add_argument('Dec', type=aph.degrees, 
                        help='J2000 declination in sDD:MM:SS[.SSS]')
    parser.add_argument('StartDate', type=aph.date, 
                        help='UTC start date in YYYY/MM/DD')
    parser.add_argument('StartTime', type=aph.time, 
                        help='UTC start time in HH:MM:SS')
    parser.add_argument('StopDate', type=aph.date, 
                        help='UTC stop date in YYYY/MM/DD')
    parser.add_argument('StopTime', type=aph.time, 
                        help='UTC stop time in HH:MM:SS')
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-s', '--lwasv', action='store_true', 
                        help='calculate for LWA-SV instead of LWA1')
    sgroup.add_argument('-a', '--lwana', action='store_true', 
                        help='calculate for LWA-NA instead of LWA1')
    sgroup.add_argument('-o', '--ovrolwa', action='store_true',
                        help='calculate for OVRO-LWA instead of LWA1')
    parser.add_argument('-n', '--n-samples', type=aph.positive_int, default=11, 
                        help='number of samples to take between the start and stop times')
    parser.add_argument('-f', '--file', type=str, 
                        help='read MJDs to compute for from a file')
    mgroup = parser.add_mutually_exclusive_group(required=False)
    mgroup.add_argument('-i', '--igs', action='store_true', 
                        help='use the IGS data products')               
    mgroup.add_argument('-j', '--jpl', action='store_true', 
                        help='use the JPL data products')
    mgroup.add_argument('-e', '--emr', action='store_true', 
                        help='use the EMR data products')
    mgroup.add_argument('-c', '--code', action='store_true', 
                        help='use the CODE data products')
    mgroup.add_argument('-u', '--ustec', action='store_true', 
                        help='use the USTEC data products')
    mgroup.add_argument('-q', '--uqr', action='store_true', 
                        help='use the high time resolution UQRG data products')
    mgroup.add_argument('-g', '--glotec', action='store_true', 
                        help='use the experimental GloTEC data products')
    args = parser.parse_args()
    main(args)
