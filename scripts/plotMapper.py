#!/usr/bin/env python3

"""
Read and plot the NOSTA_MAPPER table in a FITS IDI file writen by 
lsl.writer.fitsidi if it exists.
"""

import sys
import numpy as np
from astropy.io import fits as astrofits

from lsl.common import stations

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    # Set the station
    station = stations.lwa1
    antennas = station.antennas[0::2]
    stands = []
    for ant in antennas:
        stands.append(ant.stand.id)

    # Open the FITS file for reading
    hdu = astrofits.open(args[0])

    # This try...except block catches files that don't use the stand mapper 
    # feature in lsl.writer.fitsidi.  If a NOSTA_MAPPER table doesn't exist,
    # the script fall back on the ARRAY_GEOMETRY table to gather the stand 
    # numbers.
    try:
        mapper = hdu['NOSTA_MAPPER']
        
        nosta = mapper.data.field('NOSTA')
        noact = mapper.data.field('NOACT')
        anname = mapper.data.field('ANNAME')
    except KeyError:
        ag = hdu['ARRAY_GEOMETRY']

        nosta = ag.data.field('NOSTA')
        noact = ag.data.field('NOSTA')
        anname = ag.data.field('ANNAME')

    # Print the stand mapping out
    print("Stand Mapping:")
    for sta,act,name in zip(nosta, noact, anname):
        print("  %3i -> %s (stand %3i)" % (sta, name, act))
    
    # Load in the positions of all of the stands
    xyz = np.zeros((len(antennas), 3))
    i = 0
    for i,ant in enumerate(antennas):
        xyz[i,:] = ant.stand.xyz

    # Begin the plot
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = plt.axes([0.75, 0.75, 0.1, 0.1])
    
    ## Part 1 - Station between -80 and 80 (the section inside the fence)
    ax1.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], s=40.0, alpha=0.50)
    for mSta,sAct in zip(nosta, noact):
        i = stands.index(sAct)
        ax1.plot(xyz[i,0], xyz[i,1], marker='x', linestyle=' ', alpha=1.0, color='k')
        ax1.text(xyz[i,0], xyz[i,1], ' %i' % mSta)
    ax1.set_xlabel('$\Delta$X [E-W; m]')
    ax1.set_xlim([-80, 80])
    ax1.set_ylabel('$\Delta$Y [N-S; m]')
    ax1.set_ylim([-80, 80])
    
    ## Part 2 - The outlier as inset axes
    ax2.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], s=40.0, alpha=0.50)
    for mSta,sAct in zip(nosta, noact):
        i = stands.index(sAct)
        ax2.plot(xyz[i,0], xyz[i,1], marker='x', linestyle=' ', alpha=1.0, color='k')
        ax2.text(xyz[i,0], xyz[i,1], ' %i' % mSta)
    ax2.set_title('RTA (Outlier)')
    ax2.set_xlim([335, 345])
    ax2.set_ylim([10, 20])
    ax2.xaxis.set_major_formatter( NullFormatter() )
    ax2.yaxis.set_major_formatter( NullFormatter() )
    
    ### Part 3 - Outlier direction indicator
    #x0 = 62.0
    #y0 = 62.0
    #i = stands.index(258)
    #dx = xyz[i,0] - x0
    #dy = xyz[i,1] - y0
    #ax1.arrow(x0, y0, dx*8/dx, dy*8/dx, linewidth=2.0)
    #ax1.text(x0+5, y0, '(%i m)' % np.sqrt(dx**2 + dy**2))

    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
