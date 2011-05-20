# -*- coding: utf-8 -*-

"""
Deconvolution support for images made with lsl.sim.vis.
"""

import numpy
from aipy.img import ImgW
from aipy.coord import eq2radec
from aipy.fit import RadioFixedBody
from scipy.signal import fftconvolve as convolve

from lsl.sim import vis
from lsl.common import stations
from lsl.correlator import uvUtils
from lsl.misc.mathutil import gaussparams, gaussian2d

def estimate_beam(antennas, HA, dec, MapSize=30, MapRes=0.50, MapWRes=0.10, freq=49e6):
	"""
	Compute the beam shape for a specified pointing and array configuration.  To 
	get the scale of the beam and the gridding correct, the MapSize, MapRes, and 
	MapWRes for the image all need to be specified.
	
	Returns a numpy array of the beam response.
	"""
	
	# Build the simulated array and the point source
	aa = vis.buildSimArray(stations.lwa1, antennas, numpy.array([freq,freq,])/1e9, jd=2455659.25544)
	src = {'pnt': RadioFixedBody(HA*15/180.0*numpy.pi, dec/180.0*numpy.pi, jys=5e4, mfreq=freq/1e9, index=0)}
	
	# Simulate the source
	simDict = vis.buildSimData(aa, src, jd=2455659.25544, pols=['xx',])
	
	# Break out what is needed
	beamImg = vis.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, pol='xx', chan=[0,])
	
	return beamImg.image(center=(MapSize,MapSize))


def deconvolve(antennas, aipyImg, MapSize=30, MapRes=0.50, MapWRes=0.10, lat=34.070, freq=49e6):
	"""
	Given an AIPY ImgW instance filled with data, return a deconvolved image.
	"""
	
	import matplotlib.pyplot as plt
	
	gain = 0.1
	
	xyz = aipyImg.get_eq(0.0, lat*numpy.pi/180.0, center=(MapSize,MapSize))
	HA, dec = eq2radec(xyz)
	img = aipyImg.image(center=(MapSize,MapSize))
	
	cleaned = numpy.zeros_like(img)
	working = numpy.zeros_like(img)
	working += img
	
	prevBeam = {}
	
	for i in xrange(250):
		# Find the location of the peak in the flux density
		peak = numpy.where( working == working.max() )
		peakX = peak[0][0]
		peakY = peak[1][0]
		peakV = working[peakX,peakY]
		
		peakHA = HA[peakX, peakY]*180/numpy.pi / 15.0
		peakDec = dec[peakX, peakY]*180/numpy.pi
		
		print i ,peakX, peakY, peakV, working.max(), peakHA, peakDec
		if peakV < 0:
			break
			
		beamIndex = peakX * 10000 + peakY
		try:
			beam = prevBeam[beamIndex]
		except KeyError:
			beam = estimate_beam(antennas, peakHA, peakDec, 
						MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, freq=freq)
			beam /= beam.max()
			prevBeam[beamIndex] = beam
			junk = numpy.where( beam == beam.max() )
			print 'beam ->', junk[0][0], junk[1][0], beam.min(), beam.max(), beam.sum()
		
		toRemove = gain*peakV*beam
		
		working -= toRemove
		cleaned[peakX,peakY] += gain*peakV
		#print 'sum ->', toRemove.sum(), toRemove.max(), toRemove.min()
	
	beam = estimate_beam(antennas, 0.0, lat, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, freq=freq)
	beamClean = numpy.zeros_like(beam)
	beamClean[40:120,40:120] += beam[40:120,40:120]
	beamClean /= beamClean.max()
	
	params = gaussparams(beamClean)
	gauGen = gaussian2d(1.0, params[1], params[2], 0.8, 1.5)
	beamClean *= 0
	for i in xrange(beamClean.shape[0]):
		for j in xrange(beamClean.shape[1]):
			beamClean[i,j] = gauGen(i,j)
	beamClean /= beamClean.max()
	
	conv = convolve(cleaned, beamClean, mode='same')
	
	fig = plt.figure()
	ax1 = fig.add_subplot(2, 2, 1)
	ax2 = fig.add_subplot(2, 2, 2)
	ax3 = fig.add_subplot(2, 2, 3)
	ax4 = fig.add_subplot(2, 2, 4)
	c = ax1.imshow(img)
	fig.colorbar(c, ax=ax1)
	d = ax2.imshow(cleaned)
	fig.colorbar(d, ax=ax2)
	e = ax3.imshow(working)
	fig.colorbar(e, ax=ax3)
	f = ax4.imshow(conv + working)
	fig.colorbar(f, ax=ax4)
	plt.show()
	
	
	