#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import aipy
import pytz
import numpy
import getopt
import pyfits
from calendar import timegm
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.statistics.robust import *
from lsl.correlator import uvUtils
from lsl.writer.fitsidi import NumericStokes

from lsl.imaging import utils
from lsl.sim import vis as simVis

from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


MST = pytz.timezone('US/Mountain')
UTC = pytz.UTC


def usage(exitCode=None):
	print """imageIDI.py - Create images from a FITS-IDI file

Usage: imageIDI.py [OPTIONS] file

Options:
-h, --help             Display this help information
-1, --freq-start       First frequency to image in MHz (Default = 10 MHz)
-2, --freq-stop        Last frequency to image in MHz (Default = 88 MHz)
-s, --dataset          Data set to image (Default = All)
-m, --uv-min           Minimun baseline uvw length to include 
                       (Default = 0 lambda at midpoint frequency)
-n, --no-labels        Disable source and grid labels
-g, --no-grid          Disable the RA/Dec grid
-f, --fits             Save the images to the specified FITS image file
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['freq1'] = 10e6
	config['freq2'] = 88e6
	config['dataset'] = 0
	config['uvMin'] = 0.0
	config['label'] = True
	config['grid'] = True
	config['fits'] = None
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "h1:2:s:m:ngf:", ["help", "freq-start=", "freq-stop=", "dataset=", "uv-min=", "no-labels", "no-grid", "fits="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-1', '--freq-start'):
			config['freq1'] = float(value)*1e6
		elif opt in ('-2', '--freq-stop'):
			config['freq2'] = float(value)*1e6
		elif opt in ('-s', '--dataset'):
			config['dataset'] = int(value)
		elif opt in ('-m', '--uv-min'):
			config['uvMin'] = float(value)
		elif opt in ('-n', '--no-labels'):
			config['label'] = False
		elif opt in ('-g', '--no-grid'):
			config['grid'] = False
		elif opt in ('-f', '--fits'):
			config['fits'] = str(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg
	
	# Return configuration
	return config


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
	config = parseConfig(args)
	filename = config['args'][0]
	
	idi = utils.CorrelatedData(filename)
	aa = idi.getAntennaArray()
	lo = idi.getObserver()
	
	nStand = len(idi.stands)
	nChan = len(idi.freq)
	freq = idi.freq
	
	print "Raw Stand Count: %i" % nStand
	print "Final Baseline Count: %i" % (nStand*(nStand-1)/2,)
	print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nChan, (freq[-1] - freq[0])/1e3/nChan)
	print "Polarization Products: %i starting with %i" % (len(idi.pols), idi.pols[0])
	
	print "Reading in FITS IDI data"
	nSets = idi.integrationCount
	for set in range(1, nSets+1):
		if config['dataset'] != 0 and config['dataset'] != set:
			continue
			
		print "Set #%i of %i" % (set, nSets)
		dataDict = idi.getDataSet(set, uvMin=config['uvMin'])
		
		# Build a list of unique JDs for the data
		jdList = []
		for jd in dataDict['jd']['xx']:
			if jd not in jdList:
				jdList.append(jd)
				
		# Find the LST
		lo.date = jdList[0] - astro.DJD_OFFSET
		lst = str(lo.sidereal_time())
		
		# Pull out the right channels
		toWork = numpy.where( (freq >= config['freq1']) & (freq <= config['freq2']) )[0]
		if len(toWork) == 0:
			raise RuntimeError("Cannot find data between %.2f and %.2f MHz" % (config['freq1']/1e6, config['freq2']/1e6))

		# Build up the images for each polarization
		print "    Gridding"
		try:
			imgXX = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='xx', chan=toWork)
		except:
			imgXX = None
			
		try:
			imgYY = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='yy', chan=toWork)
		except:
			imgYY = None
			
		try:
			imgXY = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='xy', chan=toWork)
		except:
			imgXY = None
			
		try:
			imgYX = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol='yx', chan=toWork)
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
			
			junk = img.image(center=(80,80))
			print "%s: image is %.4f to %.4f with mean %.4f" % (pol, junk.min(), junk.max(), junk.mean())
			
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
				if config['label']:
					ax.text(top[0], top[1], name, color='white', size=12)
					
			# Add in the horizon
			x = numpy.zeros(361)
			y = numpy.zeros(361)
			for i in xrange(361):
				xyz = aipy.coord.azalt2top([i*numpy.pi/180.0, 0])
				x[i] = xyz[0]
				y[i] = xyz[1]
			ax.plot(x, y, color='white')
			
			# Add lines of constant RA and dec.
			if config['grid']:
				graticle(ax, lo.sidereal_time(), lo.lat, label=config['label'])
				
		plt.show()
		
		if config['fits'] is not None:
			## Loop over the images to build up the FITS file
			hdulist = [pyfits.PrimaryHDU(),]
			for img,pol in zip((imgXX,imgYY,imgXY,imgXY), ('XX','YY','XY','YX')):
				if img is None:
					continue
					
				### Create the HDU
				try:
					hdu = pyfits.ImageHDU(data=img.image(center=(80,80)), name=pol)
				except AttributeError:
					hdu = pyfits.ImageHDU(data=img, name=pol)
					
				### Add in the coordinate information
				hdu.header['CTYPE1'] = 'RA---SIN'
				hdu.header['CRPIX1'] = 80
				hdu.header['CDELT1'] = -360.0/160.0/numpy.pi
				hdu.header['CRVAL1'] = lo.sidereal_time()*180/numpy.pi
				hdu.header['CTYPE2'] = 'DEC--SIN'
				hdu.header['CRPIX2'] = 80
				hdu.header['CDELT2'] = 360.0/160.0/numpy.pi
				hdu.header['CRVAL2'] = lo.lat*180/numpy.pi
				hdu.header['LONPOLE'] = 180.0
				hdu.header['LATPOLE'] = 90.0
				
				### Add the HDU to the list
				hdulist.append(hdu)
				
		## Save the FITS file to disk
		hdulist = pyfits.HDUList(hdulist)
		clobber = False
		if os.path.exists(config['fits']):
			yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n]" % config['fits'])
			if yn not in ('n', 'N'):
				clobber = True
		try:
			hdulist.writeto(config['fits'], clobber=clobber)
		except IOError, e:
			print "WARNING: FITS image file not saved"
			
	print "...Done"


if __name__ == "__main__":
	numpy.seterr(all='ignore')
	main(sys.argv[1:])
