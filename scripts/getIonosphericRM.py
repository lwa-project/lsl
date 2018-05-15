#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import ephem
import numpy
import getopt
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.common.mcs import datetime2mjdmpm
from lsl.misc import ionosphere


def usage(exitCode=None):
	print """getIonosphericRM - Estimate the ionospheric contribution to the RM
for an observation using the IGS final product and the IGRF.

Usage: getIonosphericRM.py [OPTIONS] RA Dec Start Stop

RA:     J2000 right ascension in HH:MM:SS[.SSS]
Dec:    J2000 declination in sDD:MM:SS[.SSS]
Start:  YYYY/MM/DD HH:MM:SS start time in UTC
Stop:   YYYY/MM/DD HH:MM:SS stop time in UTC

Options:
-h, --help             Display this help information
-s, --lwasv            Calculate for LWA-SV instead of LWA1
-o, --ovro-lwa         Calculate for OVRO-LWA instead of LWA1
-n, --n-samples        Number of samples to take between the start and stop
                       times (default = 11)
-f, --file             Read MJDs to compute for from a file
-i, --igs              Use the IGS data products (default)
-j, --jpl              Use the JPL data prodcuts
-c, --code             Use the CODE data products
-u, --ustec            Use the USTEC data products
-q, --uqr              Use the high time resolution UQRG data products
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['site'] = 'lwa1'
	config['nSamples'] = 11
	config['mjdFile'] = None
	config['type'] = 'IGS'
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hson:f:ijcuq", ["help", "lwasv", "ovro-lwa", "n-samples", "file=", "igs", "jpl", "code", "ustec", "uqr"])
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
		elif opt in ('-o', '--ovro-lwa'):
			config['site'] = 'ovro'
		elif opt in ('-n', '--n-samples'):
			config['nSamples'] = int(value)
		elif opt in ('-f', '--file'):
			config['mjdFile'] = value
		elif opt in ('-i', '--igs'):
			config['type'] = 'IGS'
		elif opt in ('-j', '--jpl'):
			config['type'] = 'JPL'
		elif opt in ('-c', '--code'):
			config['type'] = 'CODE'
		elif opt in ('-u', '--ustec'):
			config['type'] = 'USTEC'
		elif opt in ('-q', '--uqr'):
			config['type'] = 'UQR'
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg
	
	# Validate
	if config['mjdFile'] is None:
		if len(config['args']) != 6:
			raise RuntimeError("Invalid number of arguments")
	else:
		if len(config['args']) != 2:
			raise RuntimeError("Invalid number of arguments")
	if config['nSamples'] < 1:
		raise RuntimeError("Invalid number of samples to generate")
		
	# Return configuration
	return config


def main(args):
	# Parse the command line
	config = parseConfig(args)
	
	# Inputs
	RA = ephem.hours( config['args'][0] )
	dec = ephem.degrees( config['args'][1] )
	if config['mjdFile'] is not None:
		mjdList = numpy.loadtxt(config['mjdFile'])
		mjdList = mjdList.ravel()
		
		mjdList = numpy.sort(mjdList)
		
	else:
		tStart = "%s %s" % (config['args'][2], config['args'][3])
		tStop = "%s %s" % (config['args'][4], config['args'][5])
		
		# YYYY/MM/DD HH:MM:SS -> datetime instance
		tStart = datetime.strptime(tStart, "%Y/%m/%d %H:%M:%S")
		tStop = datetime.strptime(tStop, "%Y/%m/%d %H:%M:%S")
		
		# datetime instance to MJD
		mjd,mpm = datetime2mjdmpm(tStart)
		mjdStart = mjd + mpm/1000.0/86400.0
		mjd,mpm = datetime2mjdmpm(tStop)
		mjdStop = mjd + mpm/1000.0/86400.0
		
		mjdList = numpy.linspace(mjdStart, mjdStop, config['nSamples'])
	
	# Setup everthing for computing the position of the source
	if config['site'] == 'lwa1':
		site = stations.lwa1
	elif config['site'] == 'lwasv':
		site = stations.lwasv
	elif config['site'] == 'ovro':
		site = stations.lwa1
		site.lat, site.lon, site.elev = ('37.2397808', '-118.2816819', 1183.4839)
	else:
		raise RuntimeError("Unknown site: %s" % config['site'])
	obs = site.getObserver()
	bdy = ephem.FixedBody()
	bdy._ra = RA
	bdy._dec = dec
	
	# Go!
	print "%-13s  %-6s  %-6s  %-21s  %-15s" % ("MJD", "Az.", "El.", "DM [pc/cm^3]", "RM [1/m^2]")
	print "-"*(13+2+6+2+6+2+21+2+15)
	for mjd in mjdList:
		# Set the date and compute the location of the target
		obs.date = mjd + astro.MJD_OFFSET - astro.DJD_OFFSET
		bdy.compute(obs)
		az = bdy.az*180/numpy.pi
		el = bdy.alt*180/numpy.pi
		
		if el > 0:
			# Get the latitude, longitude, and height of the ionospheric pierce 
			# point in this direction
			ippLat, ippLon, ippElv = ionosphere.getIonosphericPiercePoint(site, az, el)
			
			# Load in the TEC value and the RMS from the IGS final data product
			tec, rms = ionosphere.getTECValue(mjd, lat=ippLat, lng=ippLon, includeRMS=True, type=config['type'])
			tec, rms = tec[0][0], rms[0][0]
			
			# Use the IGRF to compute the ECEF components of the Earth's magnetic field
			Bx, By, Bz = ionosphere.getMagneticField(ippLat, ippLon, ippElv, mjd=mjd, ecef=True)
			
			# Rotate the ECEF field into topocentric coordinates so that we can 
			# get the magnetic field along the line of sight
			rot = numpy.array([[ numpy.sin(site.lat)*numpy.cos(site.long), numpy.sin(site.lat)*numpy.sin(site.long), -numpy.cos(site.lat)], 
						[-numpy.sin(site.long),                     numpy.cos(site.long),                      0                  ],
						[ numpy.cos(site.lat)*numpy.cos(site.long), numpy.cos(site.lat)*numpy.sin(site.long),  numpy.sin(site.lat)]])
			## ECEF -> SEZ
			sez = numpy.dot(rot, numpy.array([Bx, By, Bz]))
			## SEZ -> NEZ
			enz = 1.0*sez[[1,0,2]]
			enz[1] *= -1.0
			
			# Compute the pointing vector for this direction and use that to get
			# B parallel.  Note that we need a negative sign when we dot to get
			# the direction of propagation right.
			pnt = numpy.array([numpy.cos(el*numpy.pi/180)*numpy.sin(az*numpy.pi/180),
						numpy.cos(el*numpy.pi/180)*numpy.cos(az*numpy.pi/180), 
						numpy.sin(el*numpy.pi/180)])
			Bparallel = -numpy.dot(pnt, enz)
			
			# Compute the dispersion measure and the RMS
			DM    = 3.24078e-23 * (tec*1e16)
			rmsDM = 3.24078e-23 * (rms*1e16)
			
			# Compute the rotation measure and the RMS
			RM    = 2.62e-13 * (tec*1e16) * (Bparallel*1e-9)
			rmsRM = 2.62e-13 * (rms*1e16) * (Bparallel*1e-9)
			
			# Report
			print "%013.6f  %6.2f  %6.2f  %8.6f +/- %8.6f  %5.3f +/- %5.3f" % (mjd, az, el, DM, rmsDM, RM, rmsRM)
		else:
			# Write out dummy values since the source isn't up
			print "%013.6f  %6.2f  %6.2f  %8s +/- %8s  %5s +/- %5s" % (mjd, az, el, '---', '---', '---', '---')


if __name__ == "__main__":
	main(sys.argv[1:])
	
