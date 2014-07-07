# -*- coding: utf-8 -*-

"""
Deconvolution support for images made with :func:`lsl.imaging.utils.buildGriddedImage`.
"""

import numpy
from aipy.img import ImgW
from aipy.coord import eq2radec, top2azalt
from aipy.fit import RadioFixedBody
from scipy.signal import fftconvolve as convolve, convolve2d

from lsl.sim import vis as simVis
from lsl.imaging import utils
from lsl.common import stations
from lsl.correlator import uvUtils
from lsl.astro import deg_to_dms, deg_to_hms
from lsl.statistics.robust import std as rStd
from lsl.misc.mathutil import gaussparams, gaussian2d

__version__ = '0.4'
__revision__ = '$Rev$'
__all__ = ['clean', 'cleanSources', 'lsq', '__version__', '__revision__', '__all__']


def _interpolateValues(data, peakX, peakY):
	x1 = int(peakX)
	x2 = x1 + 1
	y1 = int(peakY)
	y2 = y1 + 1
	
	x = peakX
	if x < 0:
		raise IndexError("Invalid x")
	y = peakY
	if y < 0:
		raise IndexError("Invalid y")
		
	dataPrime  = data[x1,y1]*(x2-x)*(y2-y)
	dataPrime += data[x2,y1]*(x-x1)*(y2-y)
	dataPrime += data[x1,y2]*(x2-x)*(y-y1)
	dataPrime += data[x2,y2]*(x-x1)*(y-y1)
	dataPrime /= (x2-x1)*(y2-y1)
	return dataPrime


def _fit2DGaussian(data):
	"""
	Fit a 2D Gaussian to the provided data.  This function returns a
	five-element tuple of:
	  * height
	  * center - X and Y
	  * sigma - X and Y
	"""
	
	from scipy.optimize import leastsq
	
	def gaussian(height, center_x, center_y, width_x, width_y):
		"""
		Returns a gaussian function with the given parameters
		
		From:  http://wiki.scipy.org/Cookbook/FittingData
		"""
		
		width_x = float(width_x)
		width_y = float(width_y)
		return lambda x,y: height*numpy.exp(
					-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

	def moments(data):
		"""
		Returns (height, x, y, width_x, width_y) the gaussian parameters 
		of a 2D distribution by calculating its moments
		
		From:  http://wiki.scipy.org/Cookbook/FittingData
		"""
		
		total = data.sum()
		X, Y = numpy.indices(data.shape)
		x = (X*data).sum()/total
		y = (Y*data).sum()/total
		col = data[:, int(y)]
		width_x = numpy.sqrt(abs((numpy.arange(col.size)-y)**2*col).sum()/col.sum())
		row = data[int(x), :]
		width_y = numpy.sqrt(abs((numpy.arange(row.size)-x)**2*row).sum()/row.sum())
		height = data.max()
		return height, x, y, width_x, width_y

	def fitgaussian(data):
		"""
		Returns (height, x, y, width_x, width_y) the gaussian parameters 
		of a 2D distribution found by a fit
		
		From:  http://wiki.scipy.org/Cookbook/FittingData
		"""
		
		params = moments(data)
		errorfunction = lambda p: numpy.ravel(gaussian(*p)(*numpy.indices(data.shape)) -
								data)
		p, success = leastsq(errorfunction, params)
		return p
		
	params = fitgaussian(data)
	fit = gaussian(*params)
	
	#pylab.matshow(data, cmap=pylab.cm.gist_earth_r)
	#pylab.contour(fit(*numpy.indices(data.shape)), cmap=pylab.cm.copper)
	#pylab.show()
	
	return params


def clean(aa, dataDict, aipyImg, imageInput=None, MapSize=80, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None, gain=0.2, maxIter=150, sigma=3.0, verbose=True, plot=False):
	"""
	Given a AIPY antenna array instance, a data dictionary, and an AIPY ImgW 
	instance filled with data, return a deconvolved image.  This function 
	uses a CLEAN-like method that computes the array beam for each peak in 
	the flux.  Thus the CLEAN loop becomes:
	  1.  Find the peak flux in the residual image
	  2.  Compute the systems response to a point source at that location
	  3.  Remove the scaled porition of this beam from the residuals
	  4.  Go to 1.
	
	CLEAN tuning parameters:
	  * gain - CLEAN loop gain (default 0.2)
	  * maxIter - Maximum number of iterations (default 150)
	  * sigma - Threshold in sigma to stop cleaning (default 3.0)
	"""
	
	# Sort out the channels to work on
	if chan is None:
		chan = range(dataDict['freq'].size)
		
	# Get a grid of right ascensions and dec values for the image we are working with
	xyz = aipyImg.get_eq(0.0, aa.lat, center=(MapSize,MapSize))
	RA, dec = eq2radec(xyz)
	RA += aa.sidereal_time()
	RA %= (2*numpy.pi)
	top = aipyImg.get_top(center=(MapSize,MapSize))
	az,alt = top2azalt(top)
	
	# Get the list of baselines to generate visibilites for
	baselines = dataDict['bls'][pol]
	
	# Get the actual image out of the ImgW instance
	if imageInput is None:
		img = aipyImg.image(center=(MapSize,MapSize))
	else:
		img = imageInput*1.0
		
	# Setup the arrays to hold the point sources and the residual.
	cleaned = numpy.zeros_like(img)
	working = numpy.zeros_like(img)
	working += img
	
	# Setup the dictionary that will hold the beams as they are computed
	prevBeam = {}
	
	# Estimate the zenith beam response
	psfSrc = {'z': RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0)}
	psfDict = simVis.buildSimData(aa, psfSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flatResponse=True)
	psf = utils.buildGriddedImage(psfDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, pol=pol, verbose=verbose)
	psf = psf.image(center=(MapSize,MapSize))
	psf /= psf.max()
	
	# Fit a Guassian to the zenith beam response and use that for the restore beam
	beamCutout = psf[MapSize/2:3*MapSize/2, MapSize/2:3*MapSize/2]
	beamCutout = numpy.where( beamCutout > 0.0, beamCutout, 0.0 )
	h, cx, cy, sx, sy = _fit2DGaussian( beamCutout )
	gauGen = gaussian2d(1.0, MapSize/2+cx, MapSize/2+cy, sx, sy)
	FWHM = int( round( (sx+sy)/2.0 * 2.0*numpy.sqrt(2.0*numpy.log(2.0)) ) )
	beamClean = psf * 0.0
	for i in xrange(beamClean.shape[0]):
		for j in xrange(beamClean.shape[1]):
			beamClean[i,j] = gauGen(i,j)
	beamClean /= beamClean.sum()
	convMask = xyz.mask[0,:,:]
	
	# Go!
	if plot:
		import pylab
		from matplotlib import pyplot as plt
		
		pylab.ion()
		
	exitStatus = 'iteration limit'
	for i in xrange(maxIter):
		# Find the location of the peak in the flux density
		peak = numpy.where( working == working.max() )
		peakX = peak[0][0]
		peakY = peak[1][0]
		peakV = working[peakX,peakY]
		
		# Optimize the location
		peakParams = _fit2DGaussian(working[peakX-FWHM/2:peakX+FWHM/2+1, peakY-FWHM/2:peakY+FWHM/2+1])
		peakVO = peakParams[0]
		peakXO = peakX - FWHM/2 + peakParams[1]
		peakYO = peakY - FWHM/2 + peakParams[2]
		
		# Quantize to try and keep the computation down without over-simplifiying things
		subpixelationLevel = 5
		peakXO = round(peakXO*subpixelationLevel)/float(subpixelationLevel)
		peakYO = round(peakYO*subpixelationLevel)/float(subpixelationLevel)
		#print 'X', peakX, peakXO
		#print 'Y', peakY, peakYO
		
		# Pixel coordinates to right ascension, dec.
		try:
			peakRA = _interpolateValues(RA, peakXO, peakYO)
		except IndexError:
			peakRA = RA[peakX, peakY]
		try:
			peakDec = _interpolateValues(dec, peakXO, peakYO)
		except IndexError:
			peakDec = dec[peakX, peakY]
		#print 'R', peakRA, RA[peakX, peakY], peakRA - RA[peakX, peakY]
		#print 'D', peakDec, dec[peakX, peakY], peakDec - dec[peakX, peakY]
		
		# Pixel coordinates to az, el
		try:
			peakAz = _interpolateValues(az, peakXO, peakYO)
		except IndexError:
			peakAz = az[peakX, peakY]
		try:
			peakEl = _interpolateValues(alt, peakX, peakY)
		except IndexError:
			peakEl = alt[peakX, peakY]
			
		if verbose:
			currRA  = deg_to_hms(peakRA * 180/numpy.pi)
			currDec = deg_to_dms(peakDec * 180/numpy.pi)
			currAz  = deg_to_dms(peakAz * 180/numpy.pi)
			currEl  = deg_to_dms(peakEl * 180/numpy.pi)
			
			print "Iteration %i:  Log peak of %.2f at row: %i, column: %i" % (i+1, numpy.log10(peakV), peakX, peakY)
			print "               -> RA: %s, Dec: %s" % (currRA, currDec)
			print "               -> az: %s, el: %s" % (currAz, currEl)
			
		# Check for the exit criteria
		if peakV < 0:
			exitStatus = 'peak value is negative'
			
			break
		
		# Find the beam index and see if we need to compute the beam or not
		beamIndex = (int(peakXO*subpixelationLevel), int(peakYO*subpixelationLevel))
		try:
			beam = prevBeam[beamIndex]
			
		except KeyError:
			if verbose:
				print "               -> Computing beam(s)"
				
			beamSrc = {'Beam': RadioFixedBody(peakRA, peakDec, jys=1.0, index=0)}
			beamDict = simVis.buildSimData(aa, beamSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flatResponse=True)
			beam = utils.buildGriddedImage(beamDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, pol=pol, verbose=verbose)
			beam = beam.image(center=(MapSize,MapSize))
			beam /= beam.max()
			print "                  ", beam.mean(), beam.min(), beam.max(), beam.sum()
			
			prevBeam[beamIndex] = beam
			if verbose:
				print "               -> Beam cache contains %i entries" % len(prevBeam.keys())
				
		# Calculate how much signal needs to be removed...
		toRemove = gain*peakV*beam
		working -= toRemove
		asum = 0.0
		for l in xrange(int(peakXO), int(peakXO)+2):
			if l > peakXO:
				side1 = (peakXO+0.5) - (l-0.5)
			else:
				side1 = (l+0.5) - (peakXO-0.5)
				
			for m in xrange(int(peakYO), int(peakYO)+2):
				if m > peakYO:
					side2 = (peakYO+0.5) - (m-0.5)
				else:
					side2 = (m+0.5) - (peakYO-0.5)
					
				area = side1*side2
				asum += area
				#print 'II', l, m, area, asum
				cleaned[l,m] += gain*area*peakV
				
		#print working.min(), working.max(), working.std(), working.max()/working.std(), toRemove.sum()
		
		if plot:
			pylab.subplot(2, 2, 1)
			pylab.imshow(working+toRemove, origin='lower')
			pylab.title('Before')
			
			pylab.subplot(2, 2, 2)
			pylab.imshow(working, origin='lower')
			pylab.title('After')
			
			pylab.subplot(2, 2, 3)
			pylab.imshow(toRemove, origin='lower')
			pylab.title('CLEAN Comps.')
			
			pylab.subplot(2, 2, 4)
			pylab.imshow(convolve(cleaned, beamClean, mode='same'), origin='lower')
			
			pylab.draw()
			
		if working.max()/working.std() < sigma:
			exitStatus = 'peak is less than %.3f-sigma' % sigma
			
			break
			
	# Summary
	print "Exited after %i iterations with status '%s'" % (i+1, exitStatus)
	
	# Restore
	conv = convolve(cleaned, beamClean, mode='same')
	conv = numpy.ma.array(conv, mask=convMask)
	conv *= ((img-working).max() / conv.max())
	
	if plot:
		# Make an image for comparison purposes if we are verbose
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		ax2 = fig.add_subplot(2, 2, 2)
		ax3 = fig.add_subplot(2, 2, 3)
		ax4 = fig.add_subplot(2, 2, 4)
		
		c = ax1.imshow(img, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(c, ax=ax1)
		ax1.set_title('Input')
		
		d = ax2.imshow(conv, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(d, ax=ax2)
		ax2.set_title('CLEAN Comps.')
		
		e = ax3.imshow(working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(e, ax=ax3)
		ax3.set_title('Residuals')
		
		f = ax4.imshow(conv + working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(f, ax=ax4)
		ax4.set_title('Final')
		
		plt.show()
		
	if plot:
		pylab.ioff()
		
	# Return
	conv = conv + working
	return conv


def cleanSources(aa, dataDict, aipyImg, srcs, imageInput=None, MapSize=80, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None, gain=0.1, maxIter=150, sigma=2.0, verbose=True, plot=False):
	"""
	Given a AIPY antenna array instance, a data dictionary, an AIPY ImgW 
	instance filled with data, and a dictionary of sources, return the CLEAN
	components and the residuals map.  This function uses a CLEAN-like method
	that computes the array beam for each peak in the flux.  Thus the CLEAN 
	loop becomes: 
	  1.  Find the peak flux in the residual image
	  2.  Compute the systems response to a point source at that location
	  3.  Remove the scaled porition of this beam from the residuals
	  4.  Go to 1.
	  
	This function differs from clean() in that it only cleans localized 
	regions around each source rather than the whole image.  This is
	intended to help the mem() function along.
	
	CLEAN tuning parameters:
	  * gain - CLEAN loop gain (default 0.1)
	  * maxIter - Maximum number of iterations (default 150)
	  * sigma - Threshold in sigma to stop cleaning (default 2.0)
	"""
	
	# Sort out the channels to work on
	if chan is None:
		chan = range(dataDict['freq'].size)
		
	# Get a grid of right ascensions and dec values for the image we are working with
	xyz = aipyImg.get_eq(0.0, aa.lat, center=(MapSize,MapSize))
	RA, dec = eq2radec(xyz)
	RA += aa.sidereal_time()
	RA %= (2*numpy.pi)
	top = aipyImg.get_top(center=(MapSize,MapSize))
	az,alt = top2azalt(top)
	
	# Get the list of baselines to generate visibilites for
	baselines = dataDict['bls'][pol]
	
	# Get the actual image out of the ImgW instance
	if imageInput is None:
		img = aipyImg.image(center=(MapSize,MapSize))
	else:
		img = imageInput*1.0
		
	# Setup the arrays to hold the point sources and the residual.
	cleaned = numpy.zeros_like(img)
	working = numpy.zeros_like(img)
	working += img
	
	# Setup the dictionary that will hold the beams as they are computed
	prevBeam = {}
	
	# Estimate the zenith beam response
	psfSrc = {'z': RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0)}
	psfDict = simVis.buildSimData(aa, psfSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flatResponse=True)
	psf = utils.buildGriddedImage(psfDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, pol=pol, verbose=verbose)
	psf = psf.image(center=(MapSize,MapSize))
	psf /= psf.max()
	
	# Fit a Guassian to the zenith beam response and use that for the restore beam
	beamCutout = psf[MapSize/2:3*MapSize/2, MapSize/2:3*MapSize/2]
	beamCutout = numpy.where( beamCutout > 0.0, beamCutout, 0.0 )
	h, cx, cy, sx, sy = _fit2DGaussian( beamCutout )
	gauGen = gaussian2d(1.0, MapSize/2+cx, MapSize/2+cy, sx, sy)
	FWHM = int( round( (sx+sy)/2.0 * 2.0*numpy.sqrt(2.0*numpy.log(2.0)) ) )
	beamClean = psf * 0.0
	for i in xrange(beamClean.shape[0]):
		for j in xrange(beamClean.shape[1]):
			beamClean[i,j] = gauGen(i,j)
	beamClean /= beamClean.sum()
	convMask = xyz.mask[0,:,:]
	
	# Go!
	if plot:
		import pylab
		from matplotlib import pyplot as plt
		
		pylab.ion()
		
	for name in srcs.keys():
		src = srcs[name]
		
		# Make sure the source is up
		src.compute(aa)
		print 'Source: %s @ %s degrees elevation' % (name, src.alt)
		if src.alt <= 10*numpy.pi/180.0:
			continue
			
		# Locate the approximate position of the source
		srcDist = (src.ra-RA)**2 + (src.dec-dec)**2
		srcPeak = numpy.where( srcDist == srcDist.min() )
		
		# Define the clean box - this is fixed at 2*FWHM in width on each side
		rx0 = srcPeak[0][0] - FWHM/2
		rx1 = rx0 + FWHM + 1
		ry0 = srcPeak[1][0] - FWHM/2
		ry1 = ry0 + FWHM + 1
		#print rx0, rx1, ry0, ry1, '@', FWHM
		
		# Define the background box - this lies outside the clean box and serves
		# as a reference for the background
		X, Y = numpy.indices(working.shape)
		R = numpy.sqrt( (X-srcPeak[0][0])**2 + (Y-srcPeak[1][0])**2 )
		background = numpy.where( (R <= 2*FWHM+3) & (R > 2*FWHM) )
		
		px0 = min(background[0])-1
		px1 = max(background[0])+2
		py0 = min(background[1])-1
		py1 = max(background[1])+2
		
		exitStatus = 'iteration'
		for i in xrange(maxIter):
			# Find the location of the peak in the flux density
			peak = numpy.where( working[rx0:rx1,ry0:ry1] == working[rx0:rx1,ry0:ry1].max() )
			peakX = peak[0][0] + rx0
			peakY = peak[1][0] + ry0
			peakV = working[peakX,peakY]
			
			# Optimize the location
			try:
				peakParams = _fit2DGaussian(working[peakX-FWHM/2:peakX+FWHM/2+1, peakY-FWHM/2:peakY+FWHM/2+1])
			except IndexError:
				peakParams = [peakV, peakX, peakY]
			peakVO = peakParams[0]
			peakXO = peakX - FWHM/2 + peakParams[1]
			peakYO = peakY - FWHM/2 + peakParams[2]
			
			# Quantize to try and keep the computation down without over-simplifiying things
			subpixelationLevel = 5
			peakXO = round(peakXO*subpixelationLevel)/float(subpixelationLevel)
			peakYO = round(peakYO*subpixelationLevel)/float(subpixelationLevel)
			#print 'X', peakX, peakXO, srcPeak[0][0]
			#print 'Y', peakY, peakYO, srcPeak[1][0]
			
			# Pixel coordinates to right ascension, dec.
			try:
				peakRA = _interpolateValues(RA, peakXO, peakYO)
			except IndexError:
				peakXO, peakY0 = peakX, peakY
				peakRA = RA[peakX, peakY]
			try:
				peakDec = _interpolateValues(dec, peakXO, peakYO)
			except IndexError:
				peakDec = dec[peakX, peakY]
			#print 'R', peakRA, RA[peakX, peakY], peakRA - RA[peakX, peakY]
			#print 'D', peakDec, dec[peakX, peakY], peakDec - dec[peakX, peakY]
			
			# Pixel coordinates to az, el
			try:
				peakAz = _interpolateValues(az, peakXO, peakYO)
			except IndexError:
				peakXO, peakYO = peakX, peakY
				peakAz = az[peakX, peakY]
			try:
				peakEl = _interpolateValues(alt, peakX, peakY)
			except IndexError:
				peakEl = alt[peakX, peakY]
				
			if verbose:
				currRA  = deg_to_hms(peakRA * 180/numpy.pi)
				currDec = deg_to_dms(peakDec * 180/numpy.pi)
				currAz  = deg_to_dms(peakAz * 180/numpy.pi)
				currEl  = deg_to_dms(peakEl * 180/numpy.pi)
				
				print "%s - Iteration %i:  Log peak of %.2f at row: %i, column: %i" % (name, i+1, numpy.log10(peakV), peakX, peakY)
				print "               -> RA: %s, Dec: %s" % (currRA, currDec)
				print "               -> az: %s, el: %s" % (currAz, currEl)
				
			# Check for the exit criteria
			if peakV < 0:
				exitStatus = 'peak value is negative'
				
				break
			
			# Find the beam index and see if we need to compute the beam or not
			beamIndex = (int(peakXO*subpixelationLevel), int(peakYO*subpixelationLevel))
			try:
				beam = prevBeam[beamIndex]
				
			except KeyError:
				if verbose:
					print "               -> Computing beam(s)"
					
				beamSrc = {'Beam': RadioFixedBody(peakRA, peakDec, jys=1.0, index=0)}
				beamDict = simVis.buildSimData(aa, beamSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flatResponse=True)
				beam = utils.buildGriddedImage(beamDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, pol=pol, verbose=verbose)
				beam = beam.image(center=(MapSize,MapSize))
				beam /= beam.max()
				print "                  ", beam.mean(), beam.min(), beam.max(), beam.sum()
				
				prevBeam[beamIndex] = beam
				if verbose:
					print "               -> Beam cache contains %i entries" % len(prevBeam.keys())
					
			# Calculate how much signal needs to be removed...
			toRemove = gain*peakV*beam
			working -= toRemove
			asum = 0.0
			for l in xrange(int(peakXO), int(peakXO)+2):
				if l > peakXO:
					side1 = (peakXO+0.5) - (l-0.5)
				else:
					side1 = (l+0.5) - (peakXO-0.5)
					
				for m in xrange(int(peakYO), int(peakYO)+2):
					if m > peakYO:
						side2 = (peakYO+0.5) - (m-0.5)
					else:
						side2 = (m+0.5) - (peakYO-0.5)
						
					area = side1*side2
					asum += area
					#print 'II', l, m, area, asum
					cleaned[l,m] += gain*area*peakV
					
			#print 'S1', numpy.max(working[rx0:rx1,ry0:ry1]), numpy.median(working[background]), numpy.std(working[background])
			#print 'S2', numpy.abs(numpy.max(working[rx0:rx1,ry0:ry1])-numpy.median(working[background]))/rStd(working[background])
			
			if plot:
				try:
					pylab.subplot(2, 2, 1)
					pylab.imshow((working+toRemove)[px0:px1,py0:py1], origin='lower', interpolation='nearest')
					pylab.title('Before')
					
					pylab.subplot(2, 2, 2)
					pylab.imshow(working[px0:px1,py0:py1], origin='lower', interpolation='nearest')
					pylab.title('After')
					
					pylab.subplot(2, 2, 3)
					pylab.imshow(toRemove[px0:px1,py0:py1], origin='lower', interpolation='nearest')
					pylab.title('Removed')
					
					pylab.subplot(2, 2, 4)
					pylab.imshow(convolve(cleaned, beamClean, mode='same')[px0:px1,py0:py1], origin='lower', interpolation='nearest')
					pylab.title('CLEAN Comps.')
				except:
					pass
					
				try:
					st.set_text('%s @ %i' % (name, i+1))
				except NameError:
					st = pylab.suptitle('%s @ %i' % (name, i+1))
				pylab.draw()
				
			if numpy.abs(numpy.max(working[rx0:rx1,ry0:ry1])-numpy.median(working[background]))/rStd(working[background]) <= sigma:
				exitStatus = 'peak is less than %.3f-sigma' % sigma
				
				break
				
		# Summary
		print "Exited after %i iterations with status '%s'" % (i+1, exitStatus)
		
	# Restore
	conv = convolve(cleaned, beamClean, mode='same')
	conv = numpy.ma.array(conv, mask=convMask)
	conv *= ((img-working).max() / conv.max())
	
	if plot:
		# Make an image for comparison purposes if we are verbose
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		ax2 = fig.add_subplot(2, 2, 2)
		ax3 = fig.add_subplot(2, 2, 3)
		ax4 = fig.add_subplot(2, 2, 4)
		
		c = ax1.imshow(img, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(c, ax=ax1)
		ax1.set_title('Input')
		
		d = ax2.imshow(conv, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(d, ax=ax2)
		ax2.set_title('CLEAN Comps.')
		
		e = ax3.imshow(working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(e, ax=ax3)
		ax3.set_title('Residuals')
		
		f = ax4.imshow(conv + working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(f, ax=ax4)
		ax4.set_title('Final')
		
		plt.show()
		
	if plot:
		pylab.ioff()
		
	# Return
	return conv, working


def lsq(aa, dataDict, aipyImg, imageInput=None, MapSize=80, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None, gain=0.2, maxIter=150, verbose=True, plot=False):
	"""
	Given a AIPY antenna array instance, a data dictionary, and an AIPY ImgW 
	instance filled with data, return a deconvolved image.  This function 
	implements a least squares deconvolution.
	
	Least squares tuning parameters:
	  * gain - least squares loop gain (default 0.2)
	  * maxIter - Maximum number of iteration (default 150)
	"""
	
	# Sort out the channels to work on
	if chan is None:
		chan = range(dataDict['freq'].size)
		
	# Get a grid of right ascensions and dec values for the image we are working with
	xyz = aipyImg.get_eq(aa.sidereal_time(), aa.lat, center=(MapSize,MapSize))
	top = aipyImg.get_top(center=(MapSize,MapSize))
	ra, dec = eq2radec(xyz)
	
	# Get the list of baselines to generate visibilites for
	baselines = dataDict['bls'][pol]
	
	# Estimate the zenith beam response
	psfSrc = {'z': RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0)}
	psfDict = simVis.buildSimData(aa, psfSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flatResponse=True)
	psf = utils.buildGriddedImage(psfDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, pol=pol, verbose=verbose)
	psf = psf.image(center=(MapSize,MapSize))
	psf /= psf.max()
	
	# Fit a Guassian to the zenith beam response and use that for the restore beam
	beamCutout = psf[MapSize/2:3*MapSize/2, MapSize/2:3*MapSize/2]
	beamCutout = numpy.where( beamCutout > 0.0, beamCutout, 0.0 )
	h, cx, cy, sx, sy = _fit2DGaussian( beamCutout )
	gauGen = gaussian2d(1.0, MapSize/2+cx, MapSize/2+cy, sx, sy)
	FWHM = int( round( (sx+sy)/2.0 * 2.0*numpy.sqrt(2.0*numpy.log(2.0)) ) )
	beamClean = psf * 0.0
	for i in xrange(beamClean.shape[0]):
		for j in xrange(beamClean.shape[1]):
			beamClean[i,j] = gauGen(i,j)
	beamClean /= beamClean.sum()
	convMask = xyz.mask[0,:,:]
	
	# Get the actual image out of the ImgW instance
	if imageInput is None:
		img = aipyImg.image(center=(MapSize,MapSize))
	else:
		img = imageInput*1.0
		
	# Build the initial model
	mdl = img*40
	mdl[numpy.where(mdl < 0)] = 0
	mdl[numpy.where(ra.mask == 1)] = 0
	
	# Go!
	if plot:
		import pylab
		from matplotlib import pyplot as plt
		
		pylab.ion()
		
	rChan = [chan[0], chan[-1]]
	diff = img - mdl
	diffScaled = diff/gain
	oldRMS = diff.std()
	rHist = []
	exitStatus = 'iteration'
	for k in xrange(maxIter):
		## Update the model image but don't allow negative flux
		mdl += diffScaled * gain
		mdl[numpy.where( mdl <= 0 )] = 0.0
		
		## Convert the model image to an ensemble of point sources for forward 
		## modeling
		bSrcs = {}
		for i in xrange(mdl.shape[0]):
			for j in xrange(mdl.shape[1]):
				if dec.mask[i,j]:
					continue
				if mdl[i,j] <= 0:
					continue
					
				nm = '%i-%i' % (i,j)
				bSrcs[nm] = RadioFixedBody(ra[i,j], dec[i,j], name=nm, jys=mdl[i,j], index=0)
				
		## Model the visibilities
		simDict = simVis.buildSimData(aa, bSrcs, jd=aa.get_jultime(), pols=[pol,], chan=rChan, baselines=baselines, flatResponse=True)
		
		## Form the simulated image
		simImg = utils.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=rChan, pol=pol, verbose=verbose)
		simImg = simImg.image(center=(MapSize,MapSize))
		
		## Difference the image and the simulated image and scale it to the 
		## model's peak flux
		diff = img - simImg
		diffScaled = diff * (mdl.max() / img.max())
		
		## Compute the RMS and create an appropriately scaled version of the model
		RMS = diff.std()
		mdl2 = mdl*(img-diff).max()/mdl.max()
		
		## Status report
		if verbose:
			print "Iteration %i:  %i sources used, RMS is %.4e" % (k+1, len(bSrcs.keys()), RMS)
			print "               -> maximum residual: %.4e (%.2f%% of peak)" % (max, 100.0*max/img.max())
			print "               -> delta RMS: %.4e" % (RMS-oldRMS,)
			print "               -> delta max residual: %.4e" % (max-oldMax)
			
		## Has the RMS gone up?  If so, it is time to exit.  But first, restore 
		## the previous iteration
		if RMS-oldRMS > 0:
			mdl = oldModel
			diff = oldDiff
			
			exitStatus = 'residuals'
			
			break
			
		## Save the current iteration as the previous state
		rHist.append(RMS)
		oldRMS = RMS
		oldModel = mdl
		oldDiff = diff
		oldMax = max
		
		if plot:
			pylab.subplot(3, 2, 1)
			pylab.imshow(img, origin='lower', interpolation='nearest', vmin=img.min(), vmax=img.max())
			pylab.subplot(3, 2, 2)
			pylab.imshow(simImg, origin='lower', interpolation='nearest', vmin=img.min(), vmax=img.max())
			pylab.subplot(3, 2, 3)
			pylab.imshow(diff, origin='lower', interpolation='nearest')
			pylab.subplot(3, 2, 4)
			pylab.imshow(mdl, origin='lower', interpolation='nearest')
			pylab.subplot(3, 1, 3)
			pylab.cla()
			pylab.semilogy(rHist)
			pylab.draw()
			
		## Jump start the next iteration if this is our first pass through the data
		if k == 0:
			scaleRatio = 0.75 * img.max() / simImg.max()
			diffScaled *= scaleRatio / gain
			
	# Summary
	print "Exited after %i iterations with status '%s'" % (k+1, exitStatus)
	
	# Restore
	conv = convolve(mdl2, beamClean, mode='same')
	conv = numpy.ma.array(conv, mask=convMask)
	conv *= ((img-diff).max() / conv.max())
	
	if plot:
		# Make an image for comparison purposes if we are verbose
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		ax2 = fig.add_subplot(2, 2, 2)
		ax3 = fig.add_subplot(2, 2, 3)
		ax4 = fig.add_subplot(2, 2, 4)
		
		c = ax1.imshow(img, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(c, ax=ax1)
		ax1.set_title('Input')
		
		d = ax2.imshow(simImg, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(d, ax=ax2)
		ax2.set_title('Realized Model')
		
		e = ax3.imshow(diff, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(e, ax=ax3)
		ax3.set_title('Residuals')
		
		f = ax4.imshow(conv + diff, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(f, ax=ax4)
		ax4.set_title('Final')
		
		plt.show()
		
	if plot:
		pylab.ioff()
		
	return conv + diff