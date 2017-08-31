#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import aipy
import pytz
import numpy
import getopt
import pyfits
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.statistics.robust import *
from lsl.correlator import uvUtils
from lsl.writer.fitsidi import NumericStokes

from lsl.imaging import utils, overlay
from lsl.sim import vis as simVis
from lsl.writer.fitsidi import NumericStokes

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
-i, --include          Comma seperated list of dipoles to include 
                       (Default = All)
-e, --exclude          Comma seperated list of dipoles to exclude
                       (Default = None)
-t, --topo             Display an az/el grid instead of a RA/Dec grid
-u, --utc              Display the time as UTC, instead of LST
-n, --no-labels        Disable source and grid labels
-g, --no-grid          Disable the grid
-f, --fits             Save the images to the specified FITS image file

NOTE:  If both -i/--include and -e/--exclude are specified the include list
       has priority.
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
	config['include'] = None
	config['exclude'] = None
	config['label'] = True
	config['grid'] = True
	config['coord'] = 'RADec'
	config['time'] = 'LST'
	config['fits'] = None
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "h1:2:s:m:i:e:tungf:", ["help", "freq-start=", "freq-stop=", "dataset=", "uv-min=", "include=", "exclude=", "topo", "utc", "no-labels", "no-grid", "fits="])
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
		elif opt in ('-i', '--include'):
			config['include'] = [int(v) for v in value.split(',')]
		elif opt in ('-e', '--exclude'):
			config['exclude'] = [int(v) for v in value.split(',')]
		elif opt in ('-t', '--topo'):
			config['coord'] = 'AzEl'
		elif opt in ('-u', '--utc'):
			config['time'] = 'UTC'
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
	print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nChan, (freq[1] - freq[0])/1e3)
	try:
		print "Polarization Products: %s" % ' '.join([NumericStokes[p] for p in idi.pols])
	except KeyError:
		# Catch for CASA MS that use a different numbering scheme
		NumericStokesMS = {1:'I', 2:'Q', 3:'U', 4:'V', 
					    9:'XX', 10:'XY', 11:'YX', 12:'YY'}
		print "Polarization Products: %s" % ' '.join([NumericStokesMS[p] for p in idi.pols])
		
	print "Reading in FITS IDI data"
	nSets = idi.integrationCount
	for set in range(1, nSets+1):
		if config['dataset'] != 0 and config['dataset'] != set:
			continue
			
		print "Set #%i of %i" % (set, nSets)
		dataDict = idi.getDataSet(set, uvMin=config['uvMin'])
		
		# Build a list of unique JDs for the data
		pols = dataDict['bls'].keys()
		jdList = []
		for jd in dataDict['jd'][pols[0]]:
			if jd not in jdList:
				jdList.append(jd)
				
		# Find the LST
		lo.date = jdList[0] - astro.DJD_OFFSET
		utc = str(lo.date)
		lst = str(lo.sidereal_time())
		
		# Pull out the right channels
		toWork = numpy.where( (freq >= config['freq1']) & (freq <= config['freq2']) )[0]
		if len(toWork) == 0:
			raise RuntimeError("Cannot find data between %.2f and %.2f MHz" % (config['freq1']/1e6, config['freq2']/1e6))
			
		# Integration report
		print "    Date Observed: %s" % utc
		print "    LST: %s" % lst
		print "    Selected Frequencies: %.3f to %.3f MHz" % (freq[toWork[0]]/1e6, freq[toWork[-1]]/1e6)
		
		# Prune out what needs to go
		if config['include'] is not None or config['exclude'] is not None:
			print "    Processing include/exclude lists"
			
			## Create an empty output data dictionary
			newDict = {}
			for key in dataDict.keys():
				if key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
					newDict[key] = {}
					for pol in dataDict[key].keys():
						newDict[key][pol] = []
				else:
					newDict[key] = dataDict[key]
					
			## Fill it
			for pol in dataDict['bls'].keys():
				for b in xrange(len(dataDict['bls'][pol])):
					a0,a1 = dataDict['bls'][pol][b]
					
					if config['include'] is not None:
						if idi.stands[a0] in config['include'] and idi.stands[a1] in config['include']:
							for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
								newDict[key][pol].append( dataDict[key][pol][b] )
						else:
							continue
							
					if config['exclude'] is not None:
						if idi.stands[a0] in config['exclude'] or idi.stands[a1] in config['exclude']:
							continue
						else:
							for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
								newDict[key][pol].append( dataDict[key][pol][b] )
								
			## Make the substitution so that we are working with the right data now
			dataDict = newDict
			
			## Report
			for pol in dataDict['bls'].keys():
				print "        %s now has %i baselines" % (pol, len(dataDict['bls'][pol]))
				
		# Build up the images for each polarization
		print "    Gridding"
		img1 = None
		lbl1 = 'XX'
		for p in ('xx', 'rr', 'I'):
			try:
				img1 = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol=p, chan=toWork)
				lbl1 = p.upper()
			except:
				pass
				
		img2 = None
		lbl2 = 'YY'
		for p in ('yy', 'll', 'Q'):
			try:
				img2 = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol=p, chan=toWork)
				lbl2 = p.upper()
			except:
				pass
				
		img3 = None
		lbl3 = 'XY'
		for p in ('xy', 'rl', 'U'):
			try:
				img3 = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol=p, chan=toWork)
				lbl3 = p.upper()
			except:
				pass
				
		img4 = None
		lbl4 = 'YY'
		for p in ('yx', 'lr', 'V'):
			try:
				img4 = utils.buildGriddedImage(dataDict, MapSize=80, MapRes=0.5, pol=p, chan=toWork)
				lbl4 = p.upper()
			except:
				pass
				
		# Plots
		print "    Plotting"
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		ax2 = fig.add_subplot(2, 2, 2)
		ax3 = fig.add_subplot(2, 2, 3)
		ax4 = fig.add_subplot(2, 2, 4)
		for ax, img, pol in zip([ax1, ax2, ax3, ax4], [img1, img2, img3, img4], [lbl1, lbl2, lbl3, lbl4]):
			# Skip missing images
			if img is None:
				ax.text(0.5, 0.5, 'Not found in file', color='black', size=12, horizontalalignment='center')
				
				ax.xaxis.set_major_formatter( NullFormatter() )
				ax.yaxis.set_major_formatter( NullFormatter() )
				
				if config['time'] == 'LST':
					ax.set_title("%s @ %s LST" % (pol, lst))
				else:
					ax.set_title("%s @ %s UTC" % (pol, utc))
				continue
				
			# Display the image and label with the polarization/LST
			cb = ax.imshow(img.image(center=(80,80)), extent=(1,-1,-1,1), origin='lower', 
					vmin=img.image().min(), vmax=img.image().max())
			fig.colorbar(cb, ax=ax)
			if config['time'] == 'LST':
				ax.set_title("%s @ %s LST" % (pol, lst))
			else:
				ax.set_title("%s @ %s UTC" % (pol, utc))
				
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
				if config['coord'] == 'RADec':
					overlay.graticleRADec(ax, aa)
				else:
					overlay.graticleAzEl(ax, aa)
					
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
