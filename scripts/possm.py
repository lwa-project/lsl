#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy
import pyfits
from datetime import datetime

from lsl import astro
from lsl.common.progress import ProgressBar

import matplotlib.pyplot as plt


def main(args):
	# Grab the filename and open the FITS file using PyFits
	filename = args[0]
	hdulist = pyfits.open(filename)
	uvData = hdulist['UV_DATA']

	# Extract observation date/time information
	refDate = datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S")
	refJD = astro.unix_to_utcjd(long(refDate.strftime("%s")))
	print "Date Observered:", hdulist[0].header['DATE-OBS']

	# Read in the stand mapping table if it exists
	try:
		mapper = hdu['NOSTA_MAPPER']
		
		nosta = mapper.data.field('NOSTA')
		noact = mapper.data.field('NOACT')
	except KeyError:
		ag = hdu['ARRAY_GEOMETRY']

		nosta = ag.data.field('NOSTA')
		noact = ag.data.field('NOSTA')
	standMap = {}
	stands = []
	for s, a in zip(nosta, noact):
		standMap[s] = a
		stands.append(a)

	# Create output array
	nStand = len(stands)
	nBL = nStand*(nStand-1)/2
	nChan = uvData.header['NO_CHAN']
	baselines = []
	visibilities = numpy.zeros((nBL, nChan), dtype=numpy.complex64)
	visCount = numpy.zeros(nBL)

	print "Reading in FITS IDI data"
	pb = ProgressBar(max=len(uvData.data))
	for row in uvData.data:
		bl = row['BASELINE']
		i = ((bl >> 8) & 255)
		j = (bl & 255)
		if i == j:
			pb.inc(amount=1)
			continue

		if bl not in baselines:
			baselines.append(bl)
		blIndex = baselines.index(bl)
		
		flux = row['FLUX']
		vis = numpy.zeros(len(flux)/2, dtype=numpy.complex64)
		vis.real = flux[0::2]
		vis.imag = flux[1::2]

		#visMean = robustMean(numpy.abs(vis[200:300]))
		#visStd = robustSigma(numpy.abs(vis[200:300]))
		#wgt = numpy.where( (numpy.abs(vis)-visMean) > 2.0*visStd, 0, 1)

		visibilities[blIndex,:] += vis
		visCount[blIndex] += 1

		pb.inc(amount=1)
		if pb.amount != 0 and pb.amount % 100 == 0:
			sys.stdout.write(pb.show()+'\r')
			sys.stdout.flush()
	sys.stdout.write(pb.show()+'\r')
	sys.stdout.write('\n')
	sys.stdout.flush()

	baselines = numpy.array(baselines)
	order = baselines.argsort()
	for i in range(visCount.size):
		visibilities[i,:] /= visCount[i]
	
	print "Plotting"
	pb = ProgressBar(max=nBL)
	for k in range(int(numpy.ceil(nBL/25.0))):
		fig = plt.figure()

		for j in range(25):
			try:
				i = order[j+k*25]
			except IndexError:
				plt.draw()
				break

			amp = numpy.abs(visibilities[i,:])
			phs = numpy.angle(visibilities[i,:])*180/numpy.pi

			ax = fig.add_subplot(10, 5, 2*(j/5)*5+j%5+1)
			if ((phs+360)%360).std() < phs.std():
				ax.plot((phs+360)%360, linestyle=' ', marker='x')
				ax.set_ylim([0, 360])
			else:
				ax.plot(phs, linestyle=' ', marker='x')
				ax.set_ylim([-180, 180])
			stnd1 =  standMap[(baselines[i] >> 8) & 255]
			stnd2 = standMap[baselines[i] & 255]

			ax.set_title('%i - %i' % (stnd1, stnd2))
			ax.set_ylabel('Phs')

			ax = fig.add_subplot(10, 5, 2*(j/5)*5+j%5+1+5)
			ax.plot(amp, linestyle=' ', marker='x', color='green')
			ax.set_title('%i - %i' % (stnd1, stnd2))
			ax.set_ylabel('Amp')

			pb.inc(amount=1)
			if pb.amount != 0 and pb.amount % 10 == 0:
				sys.stdout.write(pb.show()+'\r')
				sys.stdout.flush()
		plt.draw()

	sys.stdout.write(pb.show()+'\r')
	sys.stdout.write('\n')
	sys.stdout.flush()
	plt.show()
	

if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])