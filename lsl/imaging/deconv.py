# -*- coding: utf-8 -*-

"""
Deconvolution support for images made with :func:`lsl.imaging.utils.buildGriddedImage`.
"""

import numpy
from aipy.img import ImgW
from aipy.coord import eq2radec, top2azalt
from aipy.fit import RadioFixedBody
from scipy.signal import fftconvolve as convolve
from scipy.signal import convolve2d

from lsl.sim import vis as simVis
from lsl.imaging import utils
from lsl.common import stations
from lsl.correlator import uvUtils
from lsl.astro import deg_to_dms, deg_to_hms
from lsl.misc.mathutil import gaussparams, gaussian2d

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['estimateBeam', 'deconvolve', '__version__', '__revision__', '__all__']


def estimateBeam(aa, HA, dec, MapSize=80, MapRes=0.50, MapWRes=0.10, freq=49e6):
	"""
	Compute the beam shape for a specified pointing and array configuration.  To 
	get the scale of the beam and the gridding correct, the MapSize, MapRes, and 
	MapWRes for the image all need to be specified.
	
	Returns a numpy array of the beam response.
	"""
	
	# Build the point source
	src = {'pnt': RadioFixedBody(HA*15/180.0*numpy.pi, dec/180.0*numpy.pi, jys=1e3, mfreq=freq/1e9, index=0)}
	
	# Simulate the source - the JD value shouldn't matter
	simDict = simVis.buildSimData(aa, src, jd=2455659.25544, pols=['xx',])
	
	# Break out what is needed
	beamImg = utils.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, pol='xx', chan=[0,])
	
	return beamImg.image(center=(MapSize,MapSize))


def deconvolve(aa, aipyImg, MapSize=80, MapRes=0.50, MapWRes=0.10, lat=34.070, freq=49e6, gain=0.3, maxIter=150, verbose=True):
	"""
	Given a AIPY antenna array instance and an AIPY ImgW instance filled 
	with data, return a deconvolved image.  This function uses a CLEAN-like
	method that computes the array beam for each peak in the flux.  Thus the
	CLEAN loop becomes:
	  1.  Find the peak flux in the residual image
	  2.  Compute the systems response to a point source at that location
	  3.  Remove the scaled porition of this beam from the residuals
	  4.  Go to 1.
	
	CLEAN tuning parameters:
	  * gain - CLEAN loop gain (default 0.1)
	  * maxIter - Maximum number of iteration (default 150)
	"""
	
	# Get a grid of hour angle and dec values for the image we are working with
	xyz = aipyImg.get_eq(0.0, lat*numpy.pi/180.0, center=(MapSize,MapSize))
	HA, dec = eq2radec(xyz)
	top = aipyImg.get_top(center=(MapSize,MapSize))
	az,alt = top2azalt(top)
	
	# Get the actual image out of the ImgW instance
	img = aipyImg.image(center=(MapSize,MapSize))
	
	# Setup the arrays to hold the point sources and the residual.
	cleaned = numpy.zeros_like(img)
	working = numpy.zeros_like(img)
	working += img
	
	# Setup the dictionary that will hold the beams as they are computed
	prevBeam = {}
	
	# Go!
	for i in xrange(maxIter):
		# Find the location of the peak in the flux density
		peak = numpy.where( working == working.max() )
		peakX = peak[0][0]
		peakY = peak[1][0]
		peakV = working[peakX,peakY]
		
		# Pixel coordinates to hour angle, dec.
		##peakHA = (HA[peakX-2:peakX+2,peakY-2:peakY+2]*working[peakX-2:peakX+2,peakY-2:peakY+2]).sum()
		##peakHA /= working[peakX-2:peakX+2,peakY-2:peakY+2].sum()
		##peakHA *= 180/numpy.pi / 15.0
		peakHA = HA[peakX, peakY] * 180/numpy.pi / 15.0
		
		##peakDec = (dec[peakX-2:peakX+2,peakY-2:peakY+2]*working[peakX-2:peakX+2,peakY-2:peakY+2]).sum()
		##peakDec /= working[peakX-2:peakX+2,peakY-2:peakY+2].sum()
		##peakDec *= 180/numpy.pi
		peakDec = dec[peakX,peakY] * 180/numpy.pi
		
		if verbose:
			currHA  = deg_to_hms(peakHA*15.0)
			currDec = deg_to_dms(peakDec)
			
			print "Iteration %i:  Log peak of %.2f at row: %i, column: %i" % (i+1, numpy.log10(peakV), peakX, peakY)
			print "               -> HA: %s, Dec: %s" % (currHA, currDec)
			print "               -> %.1f, %.1f" % (az[peakX, peakY]*180/numpy.pi, alt[peakX, peakY]*180/numpy.pi)
		
		# Check for the exit criteria
		if peakV < 0:
			break
		
		# Find the beam index and see if we need to compute the beam or not
		beamIndex = (peakX,peakY)
		try:
			beam = prevBeam[beamIndex]
		except KeyError:
			if verbose:
				print "               -> Computing beam"
				
			beam = estimateBeam(aa, peakHA, peakDec, 
						MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, freq=freq)
			beam /= beam.max()
			print beam.mean(), beam.min(), beam.max(), beam.sum()
			prevBeam[beamIndex] = beam
		
		# Calculate how much signal needs to be removed...
		toRemove = gain*peakV*beam
		
		if i == 0:
			globalMax = working.max()
			calib = working.max()/beam.max()
		
		# And then remove it and add it into list of CLEAN components
		working -= toRemove
		cleaned[peakX,peakY] += gain*peakV*calib/globalMax
	
	# Calculate what the restore beam should look like
	beam = estimateBeam(aa, 0.0, lat, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, freq=freq)
	beamClean = numpy.zeros_like(beam)
	beamClean[40:120,40:120] += beam[40:120,40:120]
	beamClean /= beamClean.max()
	
	# Fit a Guassian to the zenith beam response and use that for the restore beam
	params = gaussparams(beamClean)
	gauGen = gaussian2d(1.0, params[1], params[2], 0.8, 1.5)
	beamClean *= 0
	for i in xrange(beamClean.shape[0]):
		for j in xrange(beamClean.shape[1]):
			beamClean[i,j] = gauGen(i,j)
	beamClean /= beamClean.max()
	
	# Restore
	conv = convolve(cleaned, beamClean, mode='same')
	
	if verbose:
		# Make an image for comparison purposes if we are verbose
		from matplotlib import pyplot as plt
		
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		ax2 = fig.add_subplot(2, 2, 2)
		ax3 = fig.add_subplot(2, 2, 3)
		ax4 = fig.add_subplot(2, 2, 4)
		
		c = ax1.imshow(img, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(c, ax=ax1)
		ax1.set_title('Input')
		
		d = ax2.imshow(cleaned, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(d, ax=ax2)
		ax2.set_title('CLEAN Comps.')
		
		e = ax3.imshow(working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest', vmin=-working.max(), vmax=working.max())
		fig.colorbar(e, ax=ax3)
		ax3.set_title('Residuals')
		
		f = ax4.imshow(conv + working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(f, ax=ax4)
		ax4.set_title('Final')
		
		plt.show()
	
	# Return
	return conv + working