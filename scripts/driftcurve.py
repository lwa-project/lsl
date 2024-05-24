#! /usr/bin/env python

"""
Predict driftcurve for a given site using a given antenna model.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import math
import numpy
import pylab
import argparse

from scipy.interpolate import interp1d

from lsl import skymap, astro
from lsl.common import stations
from lsl.common.data_access import DataAccess
from lsl.misc import parser as aph

from lsl.misc import telemetry
telemetry.track_script()


__version__  = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


def main(args):
    # Validate
    if args.pol not in ('EW', 'NS'):
        raise ValueError("Invalid polarization: %s" % args.pol)
        
    # Get the site information
    if args.lwasv:
        nam = 'lwasv'
        sta = stations.lwasv
    elif args.ovrolwa:
        nam = 'ovro'
        sta = stations.lwa1
        sta.lat, sta.lon, sta.elev = ('37.23977727', '-118.2816667', 1183.48)
    else:
        nam = 'lwa1'
        sta = stations.lwa1
        
    # Read in the skymap (GSM or LF map @ 74 MHz)
    if not args.lfsm:
        smap = skymap.SkyMapGSM(freq_MHz=args.frequency/1e6)
        if args.verbose:
            print("Read in GSM map at %.2f MHz of %s pixels; min=%f, max=%f" % (args.frequency/1e6, len(smap.ra), smap._power.min(), smap._power.max()))
    else:
        smap = skymap.SkyMapLFSM(freq_MHz=args.frequency/1e6)
        if args.verbose:
            print("Read in LFSM map at %.2f MHz of %s pixels; min=%f, max=%f" % (args.frequency/1e6, len(smap.ra), smap._power.min(), smap._power.max()))
    
    # Get the emperical model of the beam and compute it for the correct frequencies
    with DataAccess.open('antenna/lwa1-dipole-emp.npz', 'rb') as fh:
        beamDict = numpy.load(fh)
        if args.pol == 'EW':
            beamCoeff = beamDict['fitX']
        else:
            beamCoeff = beamDict['fitY']
    alphaE = numpy.polyval(beamCoeff[0,0,:], args.frequency)
    betaE =  numpy.polyval(beamCoeff[0,1,:], args.frequency)
    gammaE = numpy.polyval(beamCoeff[0,2,:], args.frequency)
    deltaE = numpy.polyval(beamCoeff[0,3,:], args.frequency)
    alphaH = numpy.polyval(beamCoeff[1,0,:], args.frequency)
    betaH =  numpy.polyval(beamCoeff[1,1,:], args.frequency)
    gammaH = numpy.polyval(beamCoeff[1,2,:], args.frequency)
    deltaH = numpy.polyval(beamCoeff[1,3,:], args.frequency)
    if args.verbose:
        print("Beam Coeffs. X: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaH, betaH, gammaH, deltaH))
        print("Beam Coeffs. Y: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaE, betaE, gammaE, deltaE))
        
    if args.empirical:
        with DataAccess.open('antenna/lwa1-dipole-cor.npz', 'rb') as fh:
            corrDict = numpy.load(fh)
            cFreqs = corrDict['freqs']
            cAlts  = corrDict['alts']
            if corrDict['degrees'].item():
                cAlts *= numpy.pi / 180.0
            cCorrs = corrDict['corrs']
            
        if args.frequency/1e6 < cFreqs.min() or args.frequency/1e6 > cFreqs.max():
            print("WARNING: Input frequency of %.3f MHz is out of range, skipping correction" % (args.frequency/1e6,))
            corrFnc = None
        else:
            fCors = cAlts*0.0
            for i in range(fCors.size):
                ffnc = interp1d(cFreqs, cCorrs[:,i], bounds_error=False)
                fCors[i] = ffnc(args.frequency/1e6)
            corrFnc = interp1d(cAlts, fCors, bounds_error=False)
            
    else:
        corrFnc = None
        
    def compute_beam_pattern(az, alt, corr=corrFnc):
        zaR = numpy.pi/2 - alt*numpy.pi / 180.0 
        azR = az*numpy.pi / 180.0
        
        c = 1.0
        if corrFnc is not None:
            c = corrFnc(alt*numpy.pi / 180.0)
            c = numpy.where(numpy.isfinite(c), c, 1.0)
            
        pE = (1-(2*zaR/numpy.pi)**alphaE)*numpy.cos(zaR)**betaE + gammaE*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaE
        pH = (1-(2*zaR/numpy.pi)**alphaH)*numpy.cos(zaR)**betaH + gammaH*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaH

        return c*numpy.sqrt((pE*numpy.cos(azR))**2 + (pH*numpy.sin(azR))**2)

    if args.do_plot:
        az = numpy.zeros((90,360))
        alt = numpy.zeros((90,360))
        for i in range(360):
            az[:,i] = i
        for i in range(90):
            alt[i,:] = i
        pylab.figure(1)
        pylab.title("Beam Response: %s pol. @ %0.2f MHz" % (args.pol, args.frequency/1e6))
        pylab.imshow(compute_beam_pattern(az, alt), extent=(0,359, 0,89), origin='lower')
        pylab.xlabel("Azimuth [deg]")
        pylab.ylabel("Altitude [deg]")
        pylab.grid(1)
        pylab.draw()
    
    # Calculate times in both site LST and UTC
    t0 = astro.get_julian_from_sys()
    lst = astro.get_local_sidereal_time(sta.long*180.0/math.pi, t0) / 24.0
    t0 -= lst*(23.933/24.0) # Compensate for shorter sidereal days
    times = numpy.arange(0.0, 1.0, args.time_step/1440.0) + t0
    
    lstList = []
    powListAnt = [] 
    
    for t in times:
        # Project skymap to site location and observation time
        pmap = skymap.ProjectedSkyMap(smap, sta.lat*180.0/math.pi, sta.long*180.0/math.pi, t)
        lst = astro.get_local_sidereal_time(sta.long*180.0/math.pi, t)
        lstList.append(lst)
        
        # Convolution of user antenna pattern with visible skymap
        gain = compute_beam_pattern(pmap.visibleAz, pmap.visibleAlt)
        powerAnt = (pmap.visiblePower * gain).sum() / gain.sum()
        powListAnt.append(powerAnt)

        if args.verbose:
            lstH = int(lst)
            lstM = int((lst - lstH)*60.0)
            lstS = ((lst - lstH)*60.0 - lstM)*60.0
            sys.stdout.write("LST: %02i:%02i:%04.1f, Power_ant: %.1f K\r" % (lstH, lstM, lstS, powerAnt))
            sys.stdout.flush()
    sys.stdout.write("\n")
            
    # plot results
    if args.do_plot:
        pylab.figure(2)
        pylab.title("Driftcurve: %s pol. @ %0.2f MHz - %s" % \
            (args.pol, args.frequency/1e6, nam.upper()))
        pylab.plot(lstList, powListAnt, "ro", label="Antenna Pattern")
        pylab.xlabel("LST [hours]")
        pylab.ylabel("Temp. [K]")
        pylab.grid(2)
        pylab.draw()
        pylab.show()
    
    outputFile = "driftcurve_%s_%s_%.2f.txt" % (nam, args.pol, args.frequency/1e6)
    print("Writing driftcurve to file '%s'" % outputFile)
    mf = file(outputFile, "w")
    for lst,pow in zip(lstList, powListAnt):
        mf.write("%f  %f\n" % (lst, pow))
    mf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='simulate a drift curve for a dipole at LWA1 observing at a given frequency in MHz', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-f', '--frequency', type=aph.frequency, default='74.0', 
                    help='frequency of the simulation in MHz')
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-s', '--lwasv', action='store_true', 
                        help='calculate for LWA-SV instead of LWA1')
    sgroup.add_argument('-o', '--ovrolwa', action='store_true',
                        help='calculate for OVRO-LWA instead of LWA1')
    parser.add_argument('-p', '--pol', type=str, default='EW', 
                        help='polarization of the simulations (NS or EW)')
    parser.add_argument('-e', '--empirical', action='store_true', 
                    help='enable empirical corrections to the dipole model (valid from 35 to 80 MHz)')
    parser.add_argument('-l', '--lfsm', action='store_true', 
                        help='use LFSM instead of GSM')
    parser.add_argument('-t', '--time-step', type=float, default=10.0, 
                        help='time step of the simulation in minutes')
    parser.add_argument('-x', '--do-plot', action='store_true', 
                        help='plot the driftcurve data')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='run %(prog)s in verbose mode')
    args = parser.parse_args()
    main(args)
