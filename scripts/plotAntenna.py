#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script to plot the relative response of an isolated LWA antenna 
as a function of azimuth and elevation."""

import os
import sys
import numpy
import getopt

from scipy.interpolate import interp1d

from lsl.common.paths import data as dataPath

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """plotAntenna.py - Plot the relative dipole response for both
polarizations of an isolated LWA-1 antenna at a particular frequency.

Usage: plotAntenna.py [OPTIONS]

Options:
-h, --help             Display this help information
-f, --freq             Frequency of the observations in MHz
                       (default = 74 MHz)
-e, --empirical        Enable empirical corrections to the dipole model
                       (valid from 35 to 80 MHz, default = no)
-v, --verbose          Run plotAntenna in vebose mode
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['freq'] = 74.0e6
	config['corr'] = False
	config['verbose'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "hvf:e", ["help", "verbose", "freq=", "empirical"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-v', '--verbose'):
			config['verbose'] = True
		elif opt in ('-f', '--freq'):
			config['freq'] = float(value)*1e6
		elif opt in ('-e', '--empirical'):
			config['corr'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = arg

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	
	# Build the grid of azimuth and elevations to plot
	az = numpy.zeros((90,360))
	alt = numpy.zeros((90,360))
	for i in range(360):
		az[:,i] = i
	for i in range(90):
		alt[i,:] = i
		
	# Build the figure and the axes
	fig = plt.figure()
	ax1 = fig.add_subplot(2, 1, 1)
	ax2 = fig.add_subplot(2, 1, 2)
	
	# Get the emperical model of the beam and compute it for the correct frequencies
	i = 0
	beamDict = numpy.load(os.path.join(dataPath, 'lwa1-dipole-emp.npz'))
	for beamCoeff in (beamDict['fitX'], beamDict['fitY']):
		alphaE = numpy.polyval(beamCoeff[0,0,:], config['freq'])
		betaE =  numpy.polyval(beamCoeff[0,1,:], config['freq'])
		gammaE = numpy.polyval(beamCoeff[0,2,:], config['freq'])
		deltaE = numpy.polyval(beamCoeff[0,3,:], config['freq'])
		alphaH = numpy.polyval(beamCoeff[1,0,:], config['freq'])
		betaH =  numpy.polyval(beamCoeff[1,1,:], config['freq'])
		gammaH = numpy.polyval(beamCoeff[1,2,:], config['freq'])
		deltaH = numpy.polyval(beamCoeff[1,3,:], config['freq'])
		if config['verbose']:
			print "Beam Coeffs. X: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaH, betaH, gammaH, deltaH)
			print "Beam Coeffs. Y: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaE, betaE, gammaE, deltaE)
			
		if config['corr']:
			corrDict = numpy.load(os.path.join(dataPath, 'lwa1-dipole-cor.npz'))
			cFreqs = corrDict['freqs']
			cAlts  = corrDict['alts']
			if corrDict['degrees'].item():
				cAlts *= numpy.pi / 180.0
			cCorrs = corrDict['corrs']
			corrDict.close()
			
			if config['freq']/1e6 < cFreqs.min() or config['freq']/1e6 > cFreqs.max():
				print "WARNING: Input frequency of %.3f MHz is out of range, skipping correction"
				corrFnc = None
			else:
				fCors = cAlts*0.0
				for j in xrange(fCors.size):
					ffnc = interp1d(cFreqs, cCorrs[:,j], bounds_error=False)
					fCors[j] = ffnc(config['freq']/1e6)
				corrFnc = interp1d(cAlts, fCors, bounds_error=False)
				
		else:
			corrFnc = None
			
		def BeamPattern(az, alt, corr=corrFnc):
			zaR = numpy.pi/2 - alt*numpy.pi / 180.0 
			azR = az*numpy.pi / 180.0
			
			c = 1.0
			if corrFnc is not None:
				c = corrFnc(alt*numpy.pi / 180.0)
				c = numpy.where(numpy.isfinite(c), c, 1.0)
				
			pE = (1-(2*zaR/numpy.pi)**alphaE)*numpy.cos(zaR)**betaE + gammaE*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaE
			pH = (1-(2*zaR/numpy.pi)**alphaH)*numpy.cos(zaR)**betaH + gammaH*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaH

			return c*numpy.sqrt((pE*numpy.cos(azR))**2 + (pH*numpy.sin(azR))**2)
	
		# Calculate the beam
		pattern = BeamPattern(az, alt)

		if i == 0:
			p = ax1.imshow(pattern, origin='lower', vmin=0, vmax=1)
			ax1.set_title('X pol. @ %.2f MHz' % (config['freq']/1e6))
			ax1.set_xlabel('Azimuth [deg.]')
			ax1.set_ylabel('Elevation [deg.]')
			
			cb = fig.colorbar(p, ax=ax1)
			cb.ax.set_ylabel('Relative Response')
			
		else:
			p = ax2.imshow(pattern, origin='lower', vmin=0, vmax=1)
			ax2.set_title('Y pol. @ %.2f MHz' % (config['freq']/1e6))
			ax2.set_xlabel('Azimuth [deg.]')
			ax2.set_ylabel('Elevation [deg.]')
	
			cb = fig.colorbar(p, ax=ax2)
			cb.ax.set_ylabel('Relative Response')
			
		i += 1
	beamDict.close()
	
	# Display
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
