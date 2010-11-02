#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import numpy

import lsl.correlator.uvUtils as uvUtils
import lsl.common.stations as lwa_common

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

def main(args):
	# Set the LWA Station
	station = lwa_common.lwa1()

	# Load in the stand position data
	data = uvUtils.getXYZ(numpy.arange(1,257))

	# Color-code the stands by their elevation
	color = data[:,2]

	# Plot the stands as colored circles
	fig = plt.figure(figsize=(8,8))

	ax1 = plt.axes([0.30, 0.30, 0.60, 0.60])
	ax2 = plt.axes([0.30, 0.05, 0.60, 0.15])
	ax3 = plt.axes([0.05, 0.30, 0.15, 0.60])
	ax4 = plt.axes([0.05, 0.05, 0.15, 0.15])
	c = ax1.scatter(data[:,0], data[:,1], c=color, s=40.0, alpha=0.50)
	ax1.set_xlabel('$\Delta$X [E-W; m]')
	ax1.set_xlim([-80, 80])
	ax1.set_ylabel('$\Delta$Y [N-S; m]')
	ax1.set_ylim([-80, 80])
	ax1.set_title('%s Site:  %.3f$^\circ$N, %.3f$^\circ$W' % (station.name, station.lat*180.0/math.pi, -station.long*180.0/math.pi))

	ax2.scatter(data[:,0], data[:,2], c=color, s=40.0)
	ax2.xaxis.set_major_formatter( NullFormatter() )
	ax2.set_ylabel('$\Delta$Z [m]')
	ax3.scatter(data[:,2], data[:,1], c=color, s=40.0)
	ax3.yaxis.set_major_formatter( NullFormatter() )
	ax3.set_xlabel('$\Delta$Z [m]')

	# Add and elevation colorbar to the right-hand side of the figure
	cb = plt.colorbar(c, cax=ax4, orientation='vertical', ticks=[-2, -1, 0, 1, 2])

	# Set the axis limits
	ax1.set_xlim([-60, 60])
	ax1.set_ylim([-60, 60])
	ax2.set_xlim( ax1.get_xlim() )
	ax3.set_ylim( ax1.get_ylim() )

	plt.show()
	fig.savefig('stands.png')


if __name__ == "__main__":
	main(sys.argv[1:])