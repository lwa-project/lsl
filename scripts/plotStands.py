#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script to read in the positions of stands at LWA-1 and make a plot
of the site."""

import os
import sys
import numpy
import getopt

import lsl.correlator.uvUtils as uvUtils
import lsl.common.stations as lwa_common

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def usage(exitCode=None):
	print """plotStands.py - Plot the x, y, and z locations of stands at 
LWA-1.  Also, mark and label particular stands, if requested.

Usage: plotStands.py [OPTIONS] [stand1 [stand2 [...]]]

Options:
-h, --help             Display this help information
-l, --label            Label the stands with their ID numbers
                       (default = No)
-v, --verbose          Run plotStands in vebose mode
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	
# Command line flags - default values
	config['label'] = False
	config['verbose'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hvl", ["help", "verbose", "label"])
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
		elif opt in ('-l', '--label'):
			config['label'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = [int(i) for i in arg]

	# Return configuration
	return config


def main(args):
	# Parse command line
	config = parseOptions(args)
	toMark = numpy.array(config['args'])-1
	
	# Set the LWA Station
	station = lwa_common.lwa1()
	stands = station.getStands()

	# Load in the stand position data
	data = numpy.zeros((len(stands),3))
	
	i = 0
	for stand in stands:
		data[i,0] = stand.x
		data[i,1] = stand.y
		data[i,2] = stand.z
		i += 1

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
	ax1.set_title('%s Site:  %.3f$^\circ$N, %.3f$^\circ$W' % (station.name, station.lat*180.0/numpy.pi, -station.long*180.0/numpy.pi))

	ax2.scatter(data[:,0], data[:,2], c=color, s=40.0)
	ax2.xaxis.set_major_formatter( NullFormatter() )
	ax2.set_ylabel('$\Delta$Z [m]')
	ax3.scatter(data[:,2], data[:,1], c=color, s=40.0)
	ax3.yaxis.set_major_formatter( NullFormatter() )
	ax3.set_xlabel('$\Delta$Z [m]')
	
	# Explicitly mark those that need to be marked
	if toMark.size != 0:
		for i in xrange(toMark.size):
			ax1.plot(data[toMark[i],0], data[toMark[i],1], marker='x', linestyle='x', color='black')
			ax2.plot(data[toMark[i],0], data[toMark[i],2], marker='x', linestyle='x', color='black')
			ax3.plot(data[toMark[i],2], data[toMark[i],1], marker='x', linestyle='x', color='black')
			
			if config['label']:
				ax1.annotate('%i' % (toMark[i]+1), xy=(data[toMark[i],0], data[toMark[i],1]), xytext=(data[toMark[i],0]+1, data[toMark[i],1]+1))

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
