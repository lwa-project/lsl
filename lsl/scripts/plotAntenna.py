#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script to plot the relative response of an isolated LWA antenna 
as a function of azimuth and elevation."""

import os
import sys
import shutil
import tempfile
import numpy as n
from scipy.special import sph_harm

from lsl.common.paths import data as dataPath
from lsl.sim import nec_util

import matplotlib.pyplot as plt

def main(args):
	if len(args) > 0:
		freq = float(args[0])
	else:
		freq = 50.0
	
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
		alphaE = numpy.polyval(beamCoeff[0,0,:], freq)
		betaE =  numpy.polyval(beamCoeff[0,1,:], freq)
		gammaE = numpy.polyval(beamCoeff[0,2,:], freq)
		deltaE = numpy.polyval(beamCoeff[0,3,:], freq)
		alphaH = numpy.polyval(beamCoeff[1,0,:], freq)
		betaH =  numpy.polyval(beamCoeff[1,1,:], freq)
		gammaH = numpy.polyval(beamCoeff[1,2,:], freq)
		deltaH = numpy.polyval(beamCoeff[1,3,:], freq)
		print "Beam Coeffs. X: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaH, betaH, gammaH, deltaH)
		print "Beam Coeffs. Y: a=%.2f, b=%.2f, g=%.2f, d=%.2f" % (alphaE, betaE, gammaE, deltaE)
		
		def BeamPattern(az, alt):
			zaR = numpy.pi/2 - alt*numpy.pi / 180.0 
			azR = az*numpy.pi / 180.0

			pE = (1-(2*zaR/numpy.pi)**alphaE)*numpy.cos(zaR)**betaE + gammaE*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaE
			pH = (1-(2*zaR/numpy.pi)**alphaH)*numpy.cos(zaR)**betaH + gammaH*(2*zaR/numpy.pi)*numpy.cos(zaR)**deltaH

			return numpy.sqrt((pE*numpy.cos(azR))**2 + (pH*numpy.sin(azR))**2)
	
		# Calculate the beam
		pattern = BeamPattern(az, alt)

		if i == 0:
			p = ax1.imshow(pattern.T, origin='lower')
			ax1.set_title('X pol. @ %.2f MHz' % freq)
			ax1.set_xlabel('Azimuth [deg.]')
			ax1.set_ylabel('Elevation [deg.]')
		else:
			p = ax2.imshow(pattern.T, origin='lower')
			ax2.set_title('Y pol. @ %.2f MHz' % freq)
			ax2.set_xlabel('Azimuth [deg.]')
			ax2.set_ylabel('Elevation [deg.]')
	
			cb = fig.colorbar(p, ax=ax)
			cb.ax.set_ylabel('Relative Response')

	# Display
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
