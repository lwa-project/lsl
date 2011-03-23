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
# List of installed stands from Joe
installedStands = numpy.array([13, 78, 9, 38, 34, 70, 106, 10, 170, 142, 32, 8, 71, 33, 69, 168, 108, 
						140, 110, 76, 161, 101, 99, 31, 67, 159, 185, 97, 127, 157, 125, 208, 
						123, 181, 206, 183, 153, 174, 225, 176, 172, 204, 178, 200, 202, 228, 
						210, 212, 214, 187, 15, 80, 14, 148, 112, 46, 118])
# List of stands we *do not* want to use.  204, 202, 178, 148, 176, 225, and 
# 200 are close to the shelter and may have some self-RFI problems.
excludeStands = numpy.array([204, 202, 178, 148, 176, 225, 200])

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

	if len(args) == 0:
		fraction = 0.75
		stands = randomSelection(20, installedStands, allStands, Fraction=fraction, IncludeOutlier=True)
	elif len(args) > 1:
		stands = numpy.array([int(i) for i in args])
	else:
		fraction = float(args[0])
		stands = randomSelection(20, installedStands, allStands, Fraction=fraction, IncludeOutlier=True)

	print stands.shape[0], stands
	print "New Stands:"
	for stand in stands:
		if stand not in installedStands:
			print " %3i" % stand

	HA = 0.0
	dec = station.lat*180.0/math.pi
	freq = 50.0e6

	uvw = uvUtils.computeUVW(stands, HA=HA, dec=dec, freq=freq)
	uvw = numpy.squeeze(uvw[:,:,0])
	uvw = uvw * 2.9979245800e8 / freq

	# Coursely grid the uv data to come up with a rough beam
	if min(stands) < 0:
		grid = numpy.zeros((900,900))
		for i in range(uvw.shape[0]):
			u = round(uvw[i,0])+450
			v = round(uvw[i,1])+450
			grid[u,v] = grid[u,v]+1
	else:
		grid = numpy.zeros((240,240))
		for i in range(uvw.shape[0]):
			u = round(uvw[i,0])+120
			v = round(uvw[i,1])+120
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
	if min(stands) < 0:
		beam = shift(beam, xy=(450,450))
		ax5.imshow((numpy.log10(beam[375:525, 375:525])*10.0).real, interpolation="nearest")
		ax5.xaxis.set_major_formatter( NullFormatter() )
		ax5.yaxis.set_major_formatter( NullFormatter() )
	else:
		beam = shift(beam, xy=(120,120))
		ax5.imshow((numpy.log10(beam.real[100:140, 100:140])*10.0).real, interpolation="nearest")
		ax5.xaxis.set_major_formatter( NullFormatter() )
		ax5.yaxis.set_major_formatter( NullFormatter() )

	c = ax1.scatter(uvw[:,0], uvw[:,1], c=uvw[:,2], s=40.0, alpha=0.75)
	d = ax1.scatter(-uvw[:,0], -uvw[:,1], c=-uvw[:,2], s=40.0, alpha=0.75)
	ax1.set_xlabel('u [m]')
	ax1.set_ylabel('v [m]')
	ax1.set_title('UV Coverage for HA=%+.3f$^h$, $\delta$=%+.3f$^\circ$ at %s' % (HA, dec, station.name))

	ax2.scatter(uvw[:,0], uvw[:,2], c=uvw[:,2], s=40.0)
	ax2.scatter(-uvw[:,0], -uvw[:,2], c=-uvw[:,2], s=40.0)
	ax2.xaxis.set_major_formatter( NullFormatter() )
	ax2.set_ylabel('w')

	ax3.scatter(uvw[:,2], uvw[:,1], c=uvw[:,2], s=40.0)
	ax3.scatter(-uvw[:,2], -uvw[:,1], c=-uvw[:,2], s=40.0)
	ax3.yaxis.set_major_formatter( NullFormatter() )
	ax3.set_xlabel('w')

	#cb = plt.colorbar(c, cax=ax4, orientation='vertical')

	rad = numpy.zeros(uvw.shape[0])
	for i in range(rad.shape[0]):
		rad[i] = math.sqrt( uvw[i,0]**2.0 + uvw[i,1]**2.0 + uvw[i,2]**2.0 ) * freq / 2.9979245800e8
	ax4.hist(rad, 20)
	ax4.set_xlabel('uv Radius [$\lambda$]')
	ax4.set_ylabel('Baselines')
	if min(stands) < 0:
		ax4.set_xticks([0, 15, 30, 45, 60, 75])
		ax4.set_yticks([0, 10, 20, 30, 40, 50])
	else:
		ax4.set_xticks([0, 4, 8, 12, 16, 20])
		ax4.set_yticks([0, 4, 8, 12, 16, 20])

	if min(stands) < 0:
		ax1.set_xlim([-450, 450])
		ax1.set_ylim([-450, 450])
	else:
		ax1.set_xlim([-120, 120])
		ax1.set_ylim([-120, 120])
	ax2.set_xlim( ax1.get_xlim() )
	ax3.set_ylim( ax1.get_ylim() )

	if min(stands) < 0:
		ax1.text(465, 375, 'Stands:')
		for i in range(stands.shape[0]):
			if stands[i] > 0:
				if stands[i] in installedStands:
					ax1.text(476.25, 337.5-i*37.5, "%i" % stands[i])
				else:
					ax1.text(476.25, 337.5-i*37.5, "%i (*)" % stands[i])
			else:
				ax1.text(476.25, 337.5-i*37.5, 'outlier')
	else:
		ax1.text(124, 100, 'Stands:')
		for i in range(stands.shape[0]):
			ax1.text(127, 90-i*10, str(stands[i]))

	plt.show()
	fig.savefig('uv-plane.png')


if __name__ == "__main__":
	main(sys.argv[1:])
