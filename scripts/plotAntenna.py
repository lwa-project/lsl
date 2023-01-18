#!/usr/bin/env python

"""
Example script to plot the relative response of an isolated LWA antenna 
as a function of azimuth and elevation.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import numpy
import argparse

from scipy.interpolate import interp1d

from lsl.common.paths import DATA as DATA_PATH
from lsl.misc import parser as aph

import matplotlib.pyplot as plt

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    # Build the grid of azimuth and elevations to plot
    az = numpy.zeros((90,360))
    alt = numpy.zeros((90,360))
    for i in range(360):
        az[:,i] = i
    for i in range(90):
        alt[i,:] = i
        
    # Build the figure and the axes
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    
    # Get the emperical model of the beam and compute it for the correct frequencies
    i = 0
    beamDict = numpy.load(os.path.join(DATA_PATH, 'lwa1-dipole-emp.npz'))
    for beamCoeff in (beamDict['fitX'], beamDict['fitY']):
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
            corrDict = numpy.load(os.path.join(DATA_PATH, 'lwa1-dipole-cor.npz'))
            cFreqs = corrDict['freqs']
            cAlts  = corrDict['alts']
            if corrDict['degrees'].item():
                cAlts *= numpy.pi / 180.0
            cCorrs = corrDict['corrs']
            corrDict.close()
            
            if args.frequency/1e6 < cFreqs.min() or args.frequency/1e6 > cFreqs.max():
                print("WARNING: Input frequency of %.3f MHz is out of range, skipping correction" % (args.frequency/1e6,))
                corrFnc = None
            else:
                fCors = cAlts*0.0
                for j in range(fCors.size):
                    ffnc = interp1d(cFreqs, cCorrs[:,j], bounds_error=False)
                    fCors[j] = ffnc(args.frequency/1e6)
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
    
        # Calculate the beam
        pattern = compute_beam_pattern(az, alt)

        if i == 0:
            p = ax1.imshow(pattern, origin='lower', vmin=0, vmax=1)
            ax1.set_title('X pol. @ %.2f MHz' % (args.frequency/1e6))
            ax1.set_xlabel('Azimuth [deg.]')
            ax1.set_ylabel('Elevation [deg.]')
            
            cb = fig.colorbar(p, ax=ax1)
            cb.ax.set_ylabel('Relative Response')
            
        else:
            p = ax2.imshow(pattern, origin='lower', vmin=0, vmax=1)
            ax2.set_title('Y pol. @ %.2f MHz' % (args.frequency/1e6))
            ax2.set_xlabel('Azimuth [deg.]')
            ax2.set_ylabel('Elevation [deg.]')
    
            cb = fig.colorbar(p, ax=ax2)
            cb.ax.set_ylabel('Relative Response')
            
        i += 1
    beamDict.close()
    
    # Display
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='plot the relative dipole response for both polarizations of an isolated LWA antenna at a particular frequency', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-f', '--frequency', type=aph.frequency, default='50.0', 
                        help='frequency in MHz to generate the dipole model for')
    parser.add_argument('-e', '--empirical', action='store_true', 
                        help='enable empirical corrections to the dipole model (valid from 35 to 80 MHz)')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='run %(prog)s in verbose mode')
    args = parser.parse_args()
    main(args)
    
