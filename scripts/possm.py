#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy
import pyfits
from datetime import datetime

from lsl import astro
from lsl.imaging import utils
from lsl.common.progress import ProgressBar

import matplotlib.pyplot as plt


def main(args):
	# Grab the filename and open the FITS file using PyFits
	filename = args[0]
	
	idi = utils.CorrelatedData(filename)
	aa = idi.getAntennaArray()
	lo = idi.getObserver()
	lo.date = idi.dateObs.strftime("%Y/%m/%d %H:%M:%S")
	jd = lo.date + astro.DJD_OFFSET
	lst = str(lo.sidereal_time())

	nStand = len(idi.stands)
	nChan = len(idi.freq)
	freq = idi.freq
	
	print "Raw Stand Count: %i" % nStand
	print "Final Baseline Count: %i" % (nStand*(nStand-1)/2,)
	print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nChan, (freq[-1] - freq[0])/1e3/nChan)
	print "Polarization Products: %i starting with %i" % (len(idi.pols), idi.pols[0])
	print "JD: %.3f" % jd

	print "Reading in FITS IDI data"
	nSets = idi.integrationCount
	for set in range(1, nSets+1):
		print "Set #%i of %i" % (set, nSets)
		dataDict = idi.getDataSet(set)
	
		print "Plotting the first 50 baselines"
		pols = dataDict['jd'].keys()
		nBL = len(dataDict['bls'][pols[0]])
		pb = ProgressBar(max=nBL)
		i = 0
		for k in range(2):
			fig = plt.figure()

			for j in range(25):
				try:
					stnd1, stnd2 = dataDict['bls'][pols[0]][i]
					vis = dataDict['vis'][pols[0]][i]
					i += 1
				except IndexError:
					plt.draw()
					break

				amp = numpy.abs(vis)
				phs = numpy.angle(vis)*180/numpy.pi

				ax = fig.add_subplot(10, 5, 2*(j/5)*5+j%5+1)
				if ((phs+360)%360).std() < phs.std():
					ax.plot((phs+360)%360, linestyle=' ', marker='x')
					ax.set_ylim([0, 360])
				else:
					ax.plot(phs, linestyle=' ', marker='x')
					ax.set_ylim([-180, 180])
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
