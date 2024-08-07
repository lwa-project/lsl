#!/usr/bin/env python3

"""
Plot the (u,v)-plane coverage of LWA1 for a zenith snapshot and the expected 
beam.
"""

import sys
import math
import numpy as np
import argparse

from lsl.common import stations, metabundle
from lsl.correlator import uvutils
from lsl.misc import parser as aph

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,  MaxNLocator

from lsl.misc import telemetry
telemetry.track_script()


def main(args):
    # Setup the LWA station information
    if args.metadata is not None:
        try:
            station = stations.parse_ssmif(args.metadata)
        except ValueError:
            station = metabundle.get_station(args.metadata, apply_sdm=True)
    elif args.lwasv:
        station = stations.lwasv
    elif args.lwana:
        station = stations.lwana
    else:
        station = stations.lwa1
        
    antennas = []
    for ant in station.antennas[0::2]:
        if ant.combined_status == 33:
            antennas.append(ant)
    print(f"Displaying (u,v) coverage for {len(antennas)} good stands")
    
    HA = 0.0
    dec = station.lat*180.0/math.pi
    
    uvw = uvutils.compute_uvw(antennas, HA=HA, dec=dec, freq=args.frequency)
    uvw = np.squeeze(uvw[:,:,0])
    
    # Coursely grid the uv data to come up with a rough beam
    grid = np.zeros((1*240,1*240))
    for i in range(uvw.shape[0]):
        u = round((uvw[i,0]+120)*1)
        v = round((uvw[i,1]+120)*1)
        try:
            grid[u,v] += 1
        except IndexError:
            pass
            
    # Plot
    # Part 1 - Setup
    fig = plt.figure(figsize=(8,8))
    ax1 = plt.axes([0.30, 0.30, 0.60, 0.60])
    ax2 = plt.axes([0.30, 0.05, 0.60, 0.15])
    ax3 = plt.axes([0.05, 0.30, 0.15, 0.60])
    ax4 = plt.axes([0.08, 0.08, 0.15, 0.15])
    ax5 = plt.axes([0.32, 0.32, 0.15, 0.15])
    
    # Part 2 - Beam response (in dB)
    beam = np.fft.fft2(grid)
    beam = np.fft.fftshift(beam)
    beam = np.abs(beam)**2
    beam = np.log10(beam)*10.0
    ax5.imshow(beam[40:200,40:200], interpolation="nearest", vmin=np.median(beam), vmax=beam.max())
    ax5.xaxis.set_major_formatter( NullFormatter() )
    ax5.yaxis.set_major_formatter( NullFormatter() )
    
    # Part 3 - uv plane plot
    ax1.scatter(uvw[:,0], uvw[:,1], c=uvw[:,2], s=10.0, alpha=0.75)
    ax1.scatter(-uvw[:,0], -uvw[:,1], c=-uvw[:,2], s=10.0, alpha=0.75)
    ax1.set_xlabel('u [$\\lambda$]')
    ax1.set_ylabel('v [$\\lambda$]')
    ax1.set_title(f"(u,v) Coverage for HA={HA:+.3f}$^h$, $\delta$={dec:+.3f}$^\circ$ at {station.name}")
    
    # Part 4 - uw plane plot
    ax2.scatter(uvw[:,0], uvw[:,2], c=uvw[:,2], s=10.0)
    ax2.scatter(-uvw[:,0], -uvw[:,2], c=-uvw[:,2], s=10.0)
    ax2.xaxis.set_major_formatter( NullFormatter() )
    ax2.set_ylabel('w [$\\lambda$]')
    
    # Part 5 - wv plane plot
    ax3.scatter(uvw[:,2], uvw[:,1], c=uvw[:,2], s=10.0)
    ax3.scatter(-uvw[:,2], -uvw[:,1], c=-uvw[:,2], s=10.0)
    ax3.yaxis.set_major_formatter( NullFormatter() )
    ax3.set_xlabel('w [$\\lambda$]')
    
    # Part 6 - Histogram of uvw distances in lambda
    rad = np.zeros(uvw.shape[0])
    for i in range(rad.shape[0]):
        rad[i] = math.sqrt( uvw[i,0]**2.0 + uvw[i,1]**2.0 + uvw[i,2]**2.0 )
    try:
        ax4.hist(rad, 20)
    except TypeError:
        ## I don't know why this happens
        pass
    ax4.set_xlabel('uvw Radius [$\lambda$]')
    ax4.set_ylabel('Baselines')
    
    # Plot adjustment
    xlim = ax1.get_xlim()
    ylim = ax1.get_ylim()
    ax1.set_xlim([np.floor(xlim[0]/25.0)*25.0, np.ceil(xlim[1]/25.0)*25.0])
    ax1.set_ylim([np.floor(ylim[0]/25.0)*25.0, np.ceil(ylim[1]/25.0)*25.0])
    
    ax2.set_xlim( ax1.get_xlim() )
    ax2.yaxis.set_major_locator( MaxNLocator(nbins=4) )
    
    ax3.set_ylim( ax1.get_ylim() )
    ax3.xaxis.set_major_locator( MaxNLocator(nbins=4) )
    
    xlim = ax4.get_xlim()
    ylim = ax4.get_ylim()
    ax4.set_xlim([np.floor(xlim[0]/25.0)*25.0, np.ceil(xlim[1]/25.0)*25.0])
    ax4.set_ylim([np.floor(ylim[0]/5.e3)*5.e3, np.ceil(ylim[1]/5.e3)*5.e3])
    ax4.xaxis.set_major_locator( MaxNLocator(nbins=4) )
    ax4.yaxis.set_major_locator( MaxNLocator(nbins=4) )
    
    # Show n' save
    plt.show()
    if args.output is not None:
        fig.savefig(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='plot the UV-plane converage of an LWA station', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-s', '--lwasv', action='store_true', 
                        help='use LWA-SV instead of LWA1')
    sgroup.add_argument('-n', '--lwana', action='store_true', 
                        help='use LWA-NA instead of LWA1')
    parser.add_argument('-f', '--frequency', type=aph.frequency, default='50.0', 
                        help='frequency in MHz to compute the (u,v) coverage')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of the SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-o', '--output', type=str, 
                        help='filename to save the plot to')
    args = parser.parse_args()
    main(args)
