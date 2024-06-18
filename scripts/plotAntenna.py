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

from lsl.sim.beam import get_avaliable_models, beam_response
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
    
    # Get the beam response
    model = 'empirical'
    if args.empirical:
        model = 'llfss'
    for i,pol in enumerate(('XX', 'YY')):
        pattern = beam_response(model, pol, az, alt, frequency=args.frequency)
        
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
    
