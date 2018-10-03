#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Predict driftcurve for a given site using a given antenna model."""

import os
import sys
import math
import numpy
import pylab
import getopt

from scipy.interpolate import interp1d

from lsl import skymap, astro
from lsl.common import stations
from lsl.common.paths import data as dataPath

__revision__ = "$Revision: 94 $"
__version__  = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"

def usage(exitCode=None):
    print """driftcurve.py - Generate a drift curve for a dipole at LWA1 
observing at a given frequency in MHz.

Usage: driftcurve.py [OPTIONS]

Options:
-h, --help             Display this help information
-s, --lwasv            Calculate for LWA-SV instead of LWA1
-o, --ovro-lwa         Calculate for OVRO-LWA instead of LWA1
-f, --freq             Frequency of the observations in MHz
                    (default = 74 MHz)
-p, --polarization     Polarization of the observations (NS or EW; 
                    default = EW)
-e, --empirical        Enable empirical corrections to the dipole model
                    (valid from 35 to 80 MHz, default = no)
-l, --lfsm             Use LFSM instead of GSM
-t, --time-step        Time step of simulations in minutes (default = 
                    10)
-x, --do-plot          Plot the driftcurve data
-v, --verbose          Run driftcurve in vebose mode
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['site'] = 'lwa1'
    config['freq'] = 74.0e6
    config['pol'] = 'EW'
    config['corr'] = False
    config['GSM'] = True
    config['tStep'] = 10.0
    config['enableDisplay'] = False
    config['verbose'] = False
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "hvsof:p:elt:x", ["help", "verbose", "lwasv", "ovro-lwa", "freq=", "polarization=", "empirical", "lfsm", "time-step=", "do-plot",])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-v', '--verbose'):
            config['verbose'] = True
        elif opt in ('-s', '--lwasv'):
            config['site'] = 'lwasv'
        elif opt in ('-o', '--ovro-lwa'):
            config['site'] = 'ovro'
        elif opt in ('-f', '--freq'):
            config['freq'] = float(value)*1e6
        elif opt in ('-p', '--polarization'):
            config['pol'] = value.upper()
        elif opt in ('-e', '--empirical'):
            config['corr'] = True
        elif opt in ('-l', '--lfsm'):
            config['GSM'] = False
        elif opt in ('-t', '--time-step'):
            config['tStep'] = float(value)
        elif opt in ('-x', '--do-plot'):
            config['enableDisplay'] = True
        else:
            assert False
    
    # Add in arguments
    config['args'] = arg

    # Check the validity of arguments
    if config['pol'] not in ('NS', 'EW'):
        print "Invalid polarization: '%s'" % config['pol']
        usage(exitCode=2)

    # Return configuration
    return config


def main(args):
    # Parse command line
    config = parseOptions(args)
    
    # Get the site information
    if config['site'] == 'lwa1':
        sta = stations.lwa1
    elif config['site'] == 'lwasv':
        sta = stations.lwasv
    elif config['site'] == 'ovro':
        sta = stations.lwa1
        sta.lat, sta.lon, sta.elev = ('37.2397808', '-118.2816819', 1183.4839)
    else:
        raise RuntimeError("Unknown site: %s" % config['site'])
        
    # Read in the skymap (GSM or LF map @ 74 MHz)
    if config['GSM']:
        smap = skymap.SkyMapGSM(freq_MHz=config['freq']/1e6)
        if config['verbose']:
            print "Read in GSM map at %.2f MHz of %s pixels; min=%f, max=%f" % (config['freq']/1e6, len(smap.ra), smap._power.min(), smap._power.max())
    else:
        smap = skymap.SkyMapLFSM(freq_MHz=config['freq']/1e6)
        if config['verbose']:
            print "Read in LFSM map at %.2f MHz of %d x %d pixels; min=%f, max=%f" % (config['freq']/1e6, smap.numPixelsX, smap.numPixelsY, smap._power.min(), smap._power.max())
    
    # Get the emperical model of the beam and compute it for the correct frequencies
    beamDict = numpy.load(os.path.join(dataPath, 'lwa1-dipole-emp.npz'))
    if config['pol'] == 'EW':
        beamCoeff = beamDict['fitX']
    else:
        beamCoeff = beamDict['fitY']
    try:
        beamDict.close()
    except AttributeError:
        pass
    alphaE = numpy.polyval(beamCoeff[0,0,:], config['freq'])
    betaE =  numpy.polyval(beamCoeff[0,1,:], config['freq'])
    gammaE = numpy.polyval(beamCoeff[0,2,:], config['freq'])
    deltaE = numpy.polyval(beamCoeff[0,3,:], config['freq'])
    alphaH = numpy.polyval(beamCoeff[1,0,:], config['freq'])
    betaH =  numpy.polyval(beamCoeff[1,1,:], config['freq'])
    gammaH = numpy.polyval(beamCoeff[1,2,:], config['freq'])
    deltaH = numpy.polyval(beamCoeff[1,3,:], config['freq'])
    if config['verbose']:
        print "Beam Coeffs. X: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaH, betaH, gammaH, deltaH)
        print "Beam Coeffs. Y: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaE, betaE, gammaE, deltaE)
        
    if config['corr']:
        corrDict = numpy.load(os.path.join(dataPath, 'lwa1-dipole-cor.npz'))
        cFreqs = corrDict['freqs']
        cAlts  = corrDict['alts']
        if corrDict['degrees'].item():
            cAlts *= numpy.pi / 180.0
        cCorrs = corrDict['corrs']
        corrDict.close()
        
        if config['freq']/1e6 < cFreqs.min() or config['freq']/1e6 > cFreqs.max():
            print "WARNING: Input frequency of %.3f MHz is out of range, skipping correction"
            corrFnc = None
        else:
            fCors = cAlts*0.0
            for i in xrange(fCors.size):
                ffnc = interp1d(cFreqs, cCorrs[:,i], bounds_error=False)
                fCors[i] = ffnc(config['freq']/1e6)
            corrFnc = interp1d(cAlts, fCors, bounds_error=False)
            
    else:
        corrFnc = None
        
    def BeamPattern(az, alt, corr=corrFnc):
        zaR = numpy.pi/2 - alt*numpy.pi / 180.0 
        azR = az*numpy.pi / 180.0
        
        c = 1.0
        if corrFnc is not None:
            c = corrFnc(alt*numpy.pi / 180.0)
            c = numpy.where(numpy.isfinite(c), c, 1.0)
            
        pE = (1-(2*zaR/numpy.pi)**alphaE)*numpy.cos(zaR)**betaE + gammaE*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaE
        pH = (1-(2*zaR/numpy.pi)**alphaH)*numpy.cos(zaR)**betaH + gammaH*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaH

        return c*numpy.sqrt((pE*numpy.cos(azR))**2 + (pH*numpy.sin(azR))**2)

    if config['enableDisplay']:
        az = numpy.zeros((90,360))
        alt = numpy.zeros((90,360))
        for i in range(360):
            az[:,i] = i
        for i in range(90):
            alt[i,:] = i
        pylab.figure(1)
        pylab.title("Beam Response: %s pol. @ %0.2f MHz" % (config['pol'], config['freq']/1e6))
        pylab.imshow(BeamPattern(az, alt), extent=(0,359, 0,89), origin='lower')
        pylab.xlabel("Azimuth [deg]")
        pylab.ylabel("Altitude [deg]")
        pylab.grid(1)
        pylab.draw()
    
    # Calculate times in both site LST and UTC
    t0 = astro.get_julian_from_sys()
    lst = astro.get_local_sidereal_time(sta.long*180.0/math.pi, t0) / 24.0
    t0 -= lst*(23.933/24.0) # Compensate for shorter sidereal days
    times = numpy.arange(0.0, 1.0, config['tStep']/1440.0) + t0
    
    lstList = []
    powListAnt = [] 
    
    for t in times:
        # Project skymap to site location and observation time
        pmap = skymap.ProjectedSkyMap(smap, sta.lat*180.0/math.pi, sta.long*180.0/math.pi, t)
        lst = astro.get_local_sidereal_time(sta.long*180.0/math.pi, t)
        lstList.append(lst)
        
        if config['GSM']:
            cdec = numpy.ones_like(pmap.visibleDec)
        else:
            cdec = numpy.cos(pmap.visibleDec * smap.degToRad)
                
        # Convolution of user antenna pattern with visible skymap
        gain = BeamPattern(pmap.visibleAz, pmap.visibleAlt)
        powerAnt = (pmap.visiblePower * gain * cdec).sum() / (gain * cdec).sum()
        powListAnt.append(powerAnt)

        if config['verbose']:
            lstH = int(lst)
            lstM = int((lst - lstH)*60.0)
            lstS = ((lst - lstH)*60.0 - lstM)*60.0
            sys.stdout.write("LST: %02i:%02i:%04.1f, Power_ant: %.1f K\r" % (lstH, lstM, lstS, powerAnt))
            sys.stdout.flush()
    sys.stdout.write("\n")
            
    # plot results
    if config['enableDisplay']:
        pylab.figure(2)
        pylab.title("Driftcurve: %s pol. @ %0.2f MHz - %s" % \
            (config['pol'], config['freq']/1e6, config['site'].upper()))
        pylab.plot(lstList, powListAnt, "ro", label="Antenna Pattern")
        pylab.xlabel("LST [hours]")
        pylab.ylabel("Temp. [K]")
        pylab.grid(2)
        pylab.draw()
        pylab.show()
    
    outputFile = "driftcurve_%s_%s_%.2f.txt" % (config['site'], config['pol'], config['freq']/1e6)
    print "Writing driftcurve to file '%s'" % outputFile
    mf = file(outputFile, "w")
    for lst,pow in zip(lstList, powListAnt):
        mf.write("%f  %f\n" % (lst, pow))
    mf.close()


if __name__ == '__main__':
    main(sys.argv[1:])
