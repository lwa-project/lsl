#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Plot the uv-plane coverage of LWA1 for a zenith snapshot and the expected 
beam.
"""

import os
import sys
import math
import numpy
import getopt

from lsl.common import stations, metabundle, metabundleADP
from lsl.correlator import uvUtils

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,  MaxNLocator


def usage(exitCode=None):
	print """plotUVCoverage.py - Plot the UV-plane converage of LWA1.

Usage: plotUVCoverage.py [OPTIONS]

Options:
-h, --help             Display this help information
-s, --lwasv            Use LWA-SV instead of LWA1
-f, --frequency        Frequency in MHz to compute the uv coverage (default 
                       50 MHz)
-m, --metadata         Name of SSMIF or metadata tarball file to use for 
                       mappings
-o, --output           Filename to save the plot to (default = do not save)
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['site'] = 'lwa1'
	config['metadata'] = ''
	config['freq'] = 50e6
	config['output'] = None
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hsf:m:o:", ["help", "lwasv", "frequency=", "metadata=", "output="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-s', '--lwasv'):
			config['site'] = 'lwasv'
		elif opt in ('-f', '--frequency'):
			config['freq'] = float(value)*1e6
		elif opt in ('-m', '--metadata'):
			config['metadata'] = value
		elif opt in ('-o', '--output'):
			config['output'] = value
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	
	# Setup the LWA station information
	if config['metadata'] != '':
		try:
			station = stations.parseSSMIF(config['metadata'])
		except ValueError:
			try:
				station = metabundle.getStation(config['metadata'], ApplySDM=True)
			except:
				station = metabundleADP.getStation(config['metadata'], ApplySDM=True)
	elif config['site'] == 'lwa1':
		station = stations.lwa1
	elif config['site'] == 'lwasv':
		station = stations.lwasv
	else:
		raise RuntimeError("Unknown site name: %s" % config['site'])
		
	antennas = []
	for ant in station.getAntennas()[0::2]:
		if ant.getStatus() == 33:
			antennas.append(ant)
	print "Displaying uv coverage for %i good stands" % len(antennas)
	
	HA = 0.0
	dec = station.lat*180.0/math.pi
	
	uvw = uvUtils.computeUVW(antennas, HA=HA, dec=dec, freq=config['freq'])
	uvw = numpy.squeeze(uvw[:,:,0])
	
	# Coursely grid the uv data to come up with a rough beam
	grid = numpy.zeros((1*240,1*240))
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
	beam = numpy.fft.fft2(grid)
	beam = numpy.fft.fftshift(beam)
	beam = numpy.abs(beam)**2
	beam = numpy.log10(beam)*10.0
	ax5.imshow(beam[40:200,40:200], interpolation="nearest", vmin=numpy.median(beam), vmax=beam.max())
	ax5.xaxis.set_major_formatter( NullFormatter() )
	ax5.yaxis.set_major_formatter( NullFormatter() )
	
	# Part 3 - uv plane plot
	c = ax1.scatter(uvw[:,0], uvw[:,1], c=uvw[:,2], s=10.0, alpha=0.75)
	d = ax1.scatter(-uvw[:,0], -uvw[:,1], c=-uvw[:,2], s=10.0, alpha=0.75)
	ax1.set_xlabel('u [$\\lambda$]')
	ax1.set_ylabel('v [$\\lambda$]')
	ax1.set_title('UV Coverage for HA=%+.3f$^h$, $\delta$=%+.3f$^\circ$ at %s' % (HA, dec, station.name))
	
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
	rad = numpy.zeros(uvw.shape[0])
	for i in range(rad.shape[0]):
		rad[i] = math.sqrt( uvw[i,0]**2.0 + uvw[i,1]**2.0 + uvw[i,2]**2.0 )
	ax4.hist(rad, 20)
	ax4.set_xlabel('uvw Radius [$\lambda$]')
	ax4.set_ylabel('Baselines')
	
	# Plot adjustment
	xlim = ax1.get_xlim()
	ylim = ax1.get_ylim()
	ax1.set_xlim([numpy.floor(xlim[0]/25.0)*25.0, numpy.ceil(xlim[1]/25.0)*25.0])
	ax1.set_ylim([numpy.floor(ylim[0]/25.0)*25.0, numpy.ceil(ylim[1]/25.0)*25.0])
	
	ax2.set_xlim( ax1.get_xlim() )
	ax2.yaxis.set_major_locator( MaxNLocator(nbins=4) )
	
	ax3.set_ylim( ax1.get_ylim() )
	ax3.xaxis.set_major_locator( MaxNLocator(nbins=4) )
	
	xlim = ax4.get_xlim()
	ylim = ax4.get_ylim()
	ax4.set_xlim([numpy.floor(xlim[0]/25.0)*25.0, numpy.ceil(xlim[1]/25.0)*25.0])
	ax4.set_ylim([numpy.floor(ylim[0]/5.e3)*5.e3, numpy.ceil(ylim[1]/5.e3)*5.e3])
	ax4.xaxis.set_major_locator( MaxNLocator(nbins=4) )
	ax4.yaxis.set_major_locator( MaxNLocator(nbins=4) )
	
	# Show n' save
	plt.show()
	if config['output'] is not None:
		fig.savefig(config['output'])


if __name__ == "__main__":
	main(sys.argv[1:])
