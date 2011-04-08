#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Randomly select 20 antennae from LWA-1 and plot the uv-plane coverage for 
a zenith snapshot and the expected beam."""

import os
import sys
import math
import numpy

import lsl.correlator.uvUtils as uvUtils
import lsl.common.stations as lwa_common

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

allStands = numpy.arange(1,257)


def shift(data, xy=[0,0]):
	data1 = numpy.concatenate([data[xy[0]:], data[:xy[0]]], axis=0)
	data2 = numpy.concatenate([data1[:,xy[1]:], data1[:,:xy[1]]], axis=1)

	return data2


def randomSelection(N, readyStands, allStands, Fraction=0.75, IncludeOutlier=False):
	Nready = int(Fraction*N)
	Nall = N - Nready
	
	standList = []
	for i in range(Nready):
		randInt = numpy.random.randint(0, readyStands.shape[0])
		while readyStands[randInt] in standList or readyStands[randInt] in excludeStands:
			randInt = numpy.random.randint(0, readyStands.shape[0])
		standList.append( readyStands[randInt] )
	
	for i in range(Nall):
		randInt = numpy.random.randint(0, allStands.shape[0])
		while allStands[randInt] in standList or allStands[randInt] in excludeStands or allStands[randInt] in readyStands:
			randInt = numpy.random.randint(0, allStands.shape[0])
		standList.append( allStands[randInt] )

	if IncludeOutlier:
		standList[0] = -1

	return numpy.array(standList)


def main(args):
	# Set the LWA Station
	station = lwa_common.lwa1()

	stands = allStands

	HA = 0.0
	dec = station.lat*180.0/math.pi
	freq = 50.0e6

	uvw = uvUtils.computeUVW(stands, HA=HA, dec=dec, freq=freq)
	uvw = numpy.squeeze(uvw[:,:,0])
	uvw = uvw * 2.9979245800e8 / freq

	# Coursely grid the uv data to come up with a rough beam
	grid = numpy.zeros((1*240,1*240))
	for i in range(uvw.shape[0]):
		u = round((uvw[i,0]+120)*1)
		v = round((uvw[i,1]+120)*1)
		grid[u,v] = grid[u,v]+1

	# Plot
	fig = plt.figure(figsize=(8,8))
	ax1 = plt.axes([0.30, 0.30, 0.60, 0.60])
	ax2 = plt.axes([0.30, 0.05, 0.60, 0.15])
	ax3 = plt.axes([0.05, 0.30, 0.15, 0.60])
	#ax4 = plt.axes([0.05, 0.05, 0.10, 0.15])
	ax4 = plt.axes([0.08, 0.08, 0.15, 0.15])
	ax5 = plt.axes([0.32, 0.32, 0.15, 0.15])

	beam = numpy.fft.fft2(grid)
	beam = shift(beam, xy=(grid.shape[0]/2,grid.shape[1]/2))
	beam = beam.real - beam.real.min()
	beam = numpy.where( beam > 0, beam, 1e-6 )
	beam = numpy.log10(beam.real)*10.0
	ax5.imshow(beam[40:200,40:200], interpolation="nearest", vmin=numpy.median(beam), vmax=beam.max())
	ax5.xaxis.set_major_formatter( NullFormatter() )
	ax5.yaxis.set_major_formatter( NullFormatter() )

	c = ax1.scatter(uvw[:,0], uvw[:,1], c=uvw[:,2], s=10.0, alpha=0.75)
	d = ax1.scatter(-uvw[:,0], -uvw[:,1], c=-uvw[:,2], s=10.0, alpha=0.75)
	ax1.set_xlabel('u [m]')
	ax1.set_ylabel('v [m]')
	ax1.set_title('UV Coverage for HA=%+.3f$^h$, $\delta$=%+.3f$^\circ$ at %s' % (HA, dec, station.name))

	ax2.scatter(uvw[:,0], uvw[:,2], c=uvw[:,2], s=10.0)
	ax2.scatter(-uvw[:,0], -uvw[:,2], c=-uvw[:,2], s=10.0)
	ax2.xaxis.set_major_formatter( NullFormatter() )
	ax2.set_ylabel('w')

	ax3.scatter(uvw[:,2], uvw[:,1], c=uvw[:,2], s=10.0)
	ax3.scatter(-uvw[:,2], -uvw[:,1], c=-uvw[:,2], s=10.0)
	ax3.yaxis.set_major_formatter( NullFormatter() )
	ax3.set_xlabel('w')

	#cb = plt.colorbar(c, cax=ax4, orientation='vertical')

	rad = numpy.zeros(uvw.shape[0])
	for i in range(rad.shape[0]):
		rad[i] = math.sqrt( uvw[i,0]**2.0 + uvw[i,1]**2.0 + uvw[i,2]**2.0 ) * freq / 2.9979245800e8
	ax4.hist(rad, 20)
	ax4.set_xlabel('uv Radius [$\lambda$]')
	ax4.set_ylabel('Baselines')
	ax4.set_xticks([0, 4, 8, 12, 16, 20])
	ax4.set_yticks([0, 500, 1000, 1500, 2000, 2500, 3000])

	ax1.set_xlim([-120, 120])
	ax1.set_ylim([-120, 120])
	ax2.set_xlim( ax1.get_xlim() )
	ax3.set_ylim( ax1.get_ylim() )

	plt.show()
	fig.savefig('uv-plane.png')


if __name__ == "__main__":
	main(sys.argv[1:])
