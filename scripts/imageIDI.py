#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import aipy
import pytz
import numpy
import pyfits
from calendar import timegm
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.statistics.robust import *
from lsl.correlator import uvUtils
from lsl.writer.fitsidi import NumericStokes

from lsl.sim import vis as simVis

from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


MST = pytz.timezone('US/Mountain')
UTC = pytz.UTC

def baselineOrder(bls):
	"""
	Like numpy.argsort(), but for a list of two-element tuples of baseline 
	pairs.  The resulting lists can then be used to sort a data dictionary
	a la sortDataDict().
	"""
	
	def __cmpBaseline(bl):
		return 1024*bl[0] + bl[1]
	
	return [i for (v, i) in sorted((v, i) for (i, v) in enumerate([__cmpBaseline(bl) for bl in bls]))]


def sortDataDict(dataDict, order=None):
	"""
	Sort a data dictionary by the specified order.  If no order is supplied, 
	the data dictionary is sorted by baseline using baselineOrder().
	"""
	
	if order is None:
		for pol in ['xx', 'yy']:
			try:
				if len(dataDict['bls'][pol]) == 0:
					continue
				order = baselineOrder(dataDict['bls'][pol])
				break
			except KeyError:
				pass
	
	for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
		for pol in ['xx', 'yy', 'xy', 'yx']:
			try:
				newList = [dataDict[key][pol][i] for i in order]
				dataDict[key][pol] = newList
			except (KeyError, IndexError):
				pass
			
	return dataDict

def graticle(ax, lst, lat, label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot lines of
	constant declinate and RA.  Declinations are spaced at 20 degree intervals
	and RAs are spaced at 2 hour intervals.

	.. note::
		LST and latitude values should be passed as radians.  This is the default
		for lwa1.getObserver.sidereal_time() and lwa1.getObserver().lat.
	"""

	# Lines of constant declination first
	decs = range(-80, 90, 20)
	ras = numpy.linspace(0, 360, 800)

	x = numpy.zeros(ras.size)
	x = numpy.ma.array(x, mask=numpy.zeros(ras.size))
	y = numpy.zeros(ras.size)
	y = numpy.ma.array(y, mask=numpy.zeros(ras.size))
	
	for dec in decs:
		x *= 0
		y *= 0

		# Loop over RA to compute the topocentric coordinates (used by the image) for
		# the lines.  Also, figure out the elevation for each point on the line so
		# we can mask those below the horizon
		for i,ra in enumerate(ras):
			eq = aipy.coord.radec2eq((-lst + ra*numpy.pi/180,dec*numpy.pi/180))
			xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
			az,alt = aipy.coord.top2azalt(xyz)
					
			x[i] = xyz[0]
			y[i] = xyz[1]
			if alt <= 0:
				x.mask[i] = 1
				y.mask[i] = 1
			else:
				x.mask[i] = 0
				y.mask[i] = 0
	
		ax.plot(x, y, color='white', alpha=0.75)
			
		eq = aipy.coord.radec2eq((-lst + lst,(dec+5)*numpy.pi/180))
		xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
		az,alt = aipy.coord.top2azalt(xyz)
			
		if alt > 15*numpy.pi/180 and label:
			ax.text(xyz[0], xyz[1], '%+i$^\circ$' % dec, color='white')

	# Lines of constant RA			
	decs = numpy.linspace(-80, 80, 400)
	ras = range(0,360,30)

	x = numpy.zeros(decs.size)
	x = numpy.ma.array(x, mask=numpy.zeros(decs.size))
	y = numpy.zeros(decs.size)
	y = numpy.ma.array(y, mask=numpy.zeros(decs.size))

	for ra in ras:
		x *= 0
		y *= 0
		
		# Loop over dec to compute the topocentric coordinates (used by the image) for
		# the lines.  Also, figure out the elevation for each point on the line so
		# we can mask those below the horizon
		for i,dec in enumerate(decs):
			eq = aipy.coord.radec2eq((-lst + ra*numpy.pi/180,dec*numpy.pi/180))
			xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
			az,alt = aipy.coord.top2azalt(xyz)
			
			x[i] = xyz[0]
			y[i] = xyz[1]
			if alt <= 0:
				x.mask[i] = 1
				y.mask[i] = 1
			else:
				x.mask[i] = 0
				y.mask[i] = 0
		
		ax.plot(x, y, color='white', alpha=0.75)

		eq = aipy.coord.radec2eq((-lst + ra*numpy.pi/180,0))
		xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
		az,alt = aipy.coord.top2azalt(xyz)

		if alt > 20*numpy.pi/180 and label:
			ax.text(xyz[0], xyz[1], '%i$^h$' % (ra/15,), color='white')



def main(args):
	filename = args[0]
	seqNum = 0
	hdulist = pyfits.open(filename)
	uvData = hdulist['UV_DATA']

	refDate = UTC.localize(datetime.strptime(hdulist[0].header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S"))
	refJD = astro.unix_to_utcjd(timegm(refDate.timetuple()))
	print "Date Observered:", hdulist[0].header['DATE-OBS']

	station = stations.lwa1

	# Figure out how the stands are mapped in the FITS IDI file
	try:
		mapper = hdulist['NOSTA_MAPPER']
		
		nosta = mapper.data.field('NOSTA')
		noact = mapper.data.field('NOACT')
		anname = mapper.data.field('ANNAME')
	except KeyError:
		ag = hdulist['ARRAY_GEOMETRY']

		nosta = ag.data.field('NOSTA')
		noact = ag.data.field('NOSTA')
		anname = ag.data.field('ANNAME')
	
	standMap = {}
	stands = []
	for nosta, noact in zip(nosta, noact):
		standMap[nosta] = noact
		stands.append(noact)
	
	antennaMap = {}
	antennas = []
	for ant in station.getAntennas():
		if ant.stand.id in stands and ant.pol == 0:
			antennas.append(ant)
			antennaMap[ant.stand.id] = ant
	
	# Determine the contents (number of polarization, frequencies, etc.) in the FITS IDI file
	nPol = uvData.header['MAXIS2']
	firstStokes = uvData.header['CRVAL2']
	freq = numpy.arange(0, uvData.header['NO_CHAN'], 1)*uvData.header['CHAN_BW'] + uvData.header['REF_FREQ']
	aa = simVis.buildSimArray(stations.lwa1, antennas, freq/1e9, jd=refJD)
	lo = stations.lwa1.getObserver()
	lo.date = refDate.strftime("%Y/%m/%d %H:%M:%S")
	lst = str(lo.sidereal_time())

	# Report
	nStand = len(stands)
	nBaseline = nStand*(nStand-1)/2
	print "Raw Stand Count: %i" % len(stands)
	print "Final Baseline Count: %i" % nBaseline
	print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, len(freq), (freq[-1]/1e3 - freq[0]/1e3)/len(freq))
	print "Polarization Products: %i starting with %s" % (nPol, NumericStokes[firstStokes])
	print " "
	
	print "Reading in FITS IDI data..."
	nSets = len(uvData.data['BASELINE']) / (nStand*(nStand+1)/2)
	for set in range(1, nSets+1):
		print "  Set #%i of %i" % (set, nSets)
		dataDict = {'freq': freq, 
					'isMasked': False, 
					'uvw': {'xx':[], 'yy': [], 'xy':[], 'yx':[]},
					'vis': {'xx':[], 'yy': [], 'xy':[], 'yx':[]},
					'wgt': {'xx':[], 'yy': [], 'xy':[], 'yx':[]},
					'msk': {'xx':[], 'yy': [], 'xy':[], 'yx':[]},
					'bls': {'xx':[], 'yy': [], 'xy':[], 'yx':[]}, 
					'jd':  {'xx':[], 'yy': [], 'xy':[], 'yx':[]}
				}

		sourceID = range(set,set+1)
		print "    Loading data"
		for row in uvData.data:
			if row['SOURCE'] in sourceID:
				bl = row['BASELINE']
				i = standMap[(bl >> 8) & 255]
				j = standMap[bl & 255]
				if i == j:
					## Skip auto-correlations
					continue
				ri = numpy.where(stands == i)[0][0]
				rj = numpy.where(stands == j)[0][0]

				uvw = numpy.array([row['UU'], row['VV'], row['WW']])
				
				jd = row['DATE'] + row['TIME']
				uvw = numpy.array([numpy.dot(uvw[0], freq), numpy.dot(uvw[1], freq), numpy.dot(uvw[2], freq)])
				flux = row['FLUX']
				
				## Unravel the data
				if nPol == 1:
					if firstStokes == -5:
						visXX = numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
						visXX.real = flux[0::(2*nPol)]
						visXX.imag = flux[1::(2*nPol)]
						wgtXX = numpy.ones(visXX.size)
						
						dataDict['uvw']['xx'].append( uvw ) 
						dataDict['vis']['xx'].append( visXX )
						dataDict['wgt']['xx'].append( wgtXX )
						dataDict['msk']['xx'].append( numpy.zeros(len(visXX), dtype=numpy.int16) )
						dataDict['bls']['xx'].append( (ri,rj) )
						dataDict['jd' ]['xx'].append( jd )
					else:
						visYY = numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
						visYY.real = flux[0::(2*nPol)]
						visYY.imag = flux[1::(2*nPol)]
						wgtYY = numpy.ones(visYY.size)
						
						dataDict['uvw']['yy'].append( uvw ) 
						dataDict['vis']['yy'].append( visYY )
						dataDict['wgt']['yy'].append( wgtYY )
						dataDict['msk']['yy'].append( numpy.zeros(len(visYY), dtype=numpy.int16) )
						dataDict['bls']['yy'].append( (ri,rj) )
						dataDict['jd' ]['yy'].append( jd )
						
				else:
					visXX = numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
					visXX.real = flux[0::(2*nPol)]
					visXX.imag = flux[1::(2*nPol)]
					wgtXX = numpy.ones(visXX.size)
					
					dataDict['uvw']['xx'].append( uvw ) 
					dataDict['vis']['xx'].append( visXX )
					dataDict['wgt']['xx'].append( wgtXX )
					dataDict['msk']['xx'].append( numpy.zeros(len(visXX), dtype=numpy.int16) )
					dataDict['bls']['xx'].append( (ri,rj) )
					dataDict['jd' ]['xx'].append( jd )
				
				if nPol > 1:
					visYY = numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
					visYY.real = flux[2::(2*nPol)]
					visYY.imag = flux[3::(2*nPol)]
					wgtYY = numpy.ones(visYY.size)
					
					dataDict['uvw']['yy'].append( uvw ) 
					dataDict['vis']['yy'].append( visYY )
					dataDict['wgt']['yy'].append( wgtYY )
					dataDict['msk']['yy'].append( numpy.zeros(len(visYY), dtype=numpy.int16) )
					dataDict['bls']['yy'].append( (ri,rj) )
					dataDict['jd' ]['yy'].append( jd )
				
				if nPol > 2:
					visXY =  numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
					visXY.real = flux[4::(2*nPol)]
					visXY.imag = flux[5::(2*nPol)]
					wgtXY = numpy.ones(visXY.size)
					
					dataDict['uvw']['xy'].append( uvw ) 
					dataDict['vis']['xy'].append( visXY )
					dataDict['wgt']['xy'].append( wgtXY )
					dataDict['msk']['xy'].append( numpy.zeros(len(visXY), dtype=numpy.int16) )
					dataDict['bls']['xy'].append( (ri,rj) )
					dataDict['jd' ]['xy'].append( jd )
				
				if nPol > 3:
					visYX =  numpy.zeros(len(flux)/2/nPol, dtype=numpy.complex64)
					visYX.real = flux[6::(2*nPol)]
					visYX.imag = flux[7::(2*nPol)]
					wgtYX = numpy.ones(visYX.size)
					
					dataDict['uvw']['yx'].append( uvw ) 
					dataDict['vis']['yx'].append( visYX )
					dataDict['wgt']['yx'].append( wgtYX )
					dataDict['msk']['yx'].append( numpy.zeros(len(visYX), dtype=numpy.int16) )
					dataDict['bls']['yx'].append( (ri,rj) )
					dataDict['jd' ]['yx'].append( jd )
				
		# Build a list of unique JDs for the data
		jdList = []
		for jd in dataDict['jd']['xx']:
			if jd not in jdList:
				jdList.append(jd)
		
		# Sort and pull out the middle channels (inner 2/3 of the band)
		dataDict = sortDataDict(dataDict)
		toWork = range(freq.size/6, 5*freq.size/6)

		# Build up the images for each polarization
		print "    Gridding"
		try:
			imgXX = simVis.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='xx', chan=toWork)
		except:
			imgXX = None
			
		try:
			imgYY = simVis.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='yy', chan=toWork)
		except:
			imgYY = None
			
		try:
			imgXY = simVis.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='xy', chan=toWork)
		except:
			imgXY = None
			
		try:
			imgYX = simVis.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='yx', chan=toWork)
		except:
			imgYX = None
		
		# Plots
		print "    Plotting"
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		ax2 = fig.add_subplot(2, 2, 2)
		ax3 = fig.add_subplot(2, 2, 3)
		ax4 = fig.add_subplot(2, 2, 4)
		for ax, img, pol in zip([ax1, ax2, ax3, ax4], [imgXX, imgYY, imgXY, imgYX], ['XX', 'YY', 'XY', 'YX']):
			# Skip missing images
			if img is None:
				ax.text(0.5, 0.5, 'Not found in file', color='black', size=12, horizontalalignment='center')

				ax.xaxis.set_major_formatter( NullFormatter() )
				ax.yaxis.set_major_formatter( NullFormatter() )

				ax.set_title("%s @ %s LST" % (pol, lst))
				continue
			
			# Display the image and label with the polarization/LST
			cb = ax.imshow(img.image(center=(80,80)), extent=(1,-1,-1,1), origin='lower', 
					vmin=img.image().min(), vmax=img.image().max())
			fig.colorbar(cb, ax=ax)
			ax.set_title("%s @ %s LST" % (pol, lst))

			# Turn off tick marks
			ax.xaxis.set_major_formatter( NullFormatter() )
			ax.yaxis.set_major_formatter( NullFormatter() )

			# Compute the positions of major sources and label the images
			compSrc = {}
			for name,src in simVis.srcs.iteritems():
				src.compute(aa)
				top = src.get_crds(crdsys='top', ncrd=3)
				az, alt = aipy.coord.top2azalt(top)
				compSrc[name] = [az, alt]
				if alt <= 0:
					continue
				ax.plot(top[0], top[1], marker='x', markerfacecolor='None', markeredgecolor='w', 
						linewidth=10.0, markersize=10)
				ax.text(top[0], top[1], name, color='white', size=12)
				
			# Add lines of constant RA and dec.
			graticle(ax, lo.sidereal_time(), lo.lat)

		plt.show()

	print "...Done"
	hdulist.close()


if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])
