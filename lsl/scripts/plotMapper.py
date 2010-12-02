#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read and plot the NOSTA_MAPPER table in a FITS IDI file writen by 
lsl.writer.fitsidi if it exists."""

import sys
import numpy
import pyfits

from lsl.correlator import uvUtils

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main(args):
	# Open the FITS file for reading
	hdu = pyfits.open(args[0])

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
	print "Stand Mapping:"
	for sta,act,name in zip(nosta, noact, anname):
		print "  %3i -> %s (stand %3i)" % (sta, name, act)
	
	# Load in the positions of all of the stands
	xyz = uvUtils.getXYZ(numpy.arange(1,259))

	# Begin the plot
	fig = plt.figure()
	ax1 = fig.add_subplot(1, 1, 1)
	ax2 = plt.axes([0.75, 0.75, 0.1, 0.1])
	## Part 1 - Station between -80 and 80 (the section inside the fence)
	c = ax1.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], s=40.0, alpha=0.50)
	ax1.plot(xyz[noact-1,0], xyz[noact-1,1], marker='x', linestyle=' ', alpha=1.0, color='k')
	for m,a in zip(nosta, noact):
		ax1.text(xyz[a-1,0], xyz[a-1,1], '%i' % m)
	ax1.set_xlabel('$\Delta$X [E-W; m]')
	ax1.set_xlim([-80, 80])
	ax1.set_ylabel('$\Delta$Y [N-S; m]')
	ax1.set_ylim([-80, 80])
	## Part 2 - The outlier as inset axes
	ax2.scatter(xyz[:,0], xyz[:,1], c=xyz[:,2], s=40.0, alpha=0.50)
	ax2.plot(xyz[noact-1,0], xyz[noact-1,1], marker='x', linestyle=' ', alpha=1.0, color='k')
	for m,a in zip(nosta, noact):
		ax2.text(xyz[a-1,0], xyz[a-1,1], '%i' % m)
	ax2.set_title('Outlier')
	ax2.set_xlim([335, 345])
	ax2.set_ylim([10, 20])
	ax2.xaxis.set_major_formatter( NullFormatter() )
	ax2.yaxis.set_major_formatter( NullFormatter() )

	plt.show()

if __name__ == "__main__":
	main(sys.argv[1:])
