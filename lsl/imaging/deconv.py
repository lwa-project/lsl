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


def estimateBeam(aa, HA, dec, MapSize=80, MapRes=0.50, MapWRes=0.10, chan=None, baselines=None, scale=0.0):
	"""
	Compute the beam shape for a specified pointing and array configuration.  To 
	get the scale of the beam and the gridding correct, the MapSize, MapRes, and 
	MapWRes for the image all need to be specified.
	
	Returns a numpy array of the beam response.
	"""
	
	# Build the point source
	src = {'pnt': RadioFixedBody(HA*15/180.0*numpy.pi, dec/180.0*numpy.pi, jys=1e3, mfreq=0.050, index=0)}
	
	# Simulate the source - the JD value shouldn't matter
	simDict = simVis.buildSimData(aa, src, jd=2455659.25544, pols=['xx',], baselines=baselines)
	
	# Break out what is needed
	beamImg = utils.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, pol='xx', chan=chan)
	
	return beamImg.image(center=(MapSize,MapSize))


def estimateBeam2(aa, HA, dec, MapSize=80, MapRes=0.50, MapWRes=0.10, chan=None, baselines=None, scale=0.0):
	"""
	Compute the beam shape for a specified pointing and array configuration.  To 
	get the scale of the beam and the gridding correct, the MapSize, MapRes, and 
	MapWRes for the image all need to be specified.
	
	Returns a numpy array of the beam response.
	"""
	
	if chan is None:
		chan = range( len(aa.get_afreqs()) )
	freq = aa.get_afreqs().ravel()
	print freq.shape
	freq = freq[chan]
		
	# Build the point source
	srcshape = (scale*numpy.pi/180.0, scale*numpy.pi/180.0, 0.0)
	src = {'pnt': RadioFixedBody(HA*15/180.0*numpy.pi, dec/180.0*numpy.pi, jys=1e3, mfreq=50e6/1e9, index=0, srcshape=srcshape)}
	
	# Simulate the source - the JD value shouldn't matter
	#simDict = simVis.buildSimData(aa, src, jd=2455659.25544, pols=['xx',], baselines=baselines)
	antennas = []
	for a in stations.lwa1.getAntennas()[::2]:
		if a.stand.id in aa.get_stands():
			antennas.append(a)
	uvws  = uvUtils.computeUVW(antennas, -HA, dec, freq=freq.mean()*1e9)
	uvwsZ = uvUtils.computeUVW(antennas, freq=freq.mean()*1e9)
	simDict = {}
	simDict['bls'] = {'xx':[]}
	simDict['uvw'] = {'xx':[]}
	simDict['vis'] = {'xx':[]}
	simDict['wgt'] = {'xx':[]}
	simDict['msk'] = {'xx':[]}
	print uvws.shape
	for i in xrange(uvws.shape[0]):
		simDict['bls']['xx'].append((i,i))
		simDict['uvw']['xx'].append(uvwsZ[i,:,:])
		simDict['vis']['xx'].append(numpy.exp(2j*numpy.pi*uvws[i,2,:])*numpy.exp(-2j*numpy.pi*uvwsZ[i,2,:]))
		simDict['wgt']['xx'].append(numpy.array([1,]))
		simDict['msk']['xx'].append(numpy.array([0,]))
	
	# Break out what is needed
	beamImg = utils.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, pol='xx')
	
	im = beamImg.image(center=(MapSize,MapSize))
	bm = beamImg.bm_image(center=(MapSize,MapSize), term=0)
	print im.min(), im.mean(), im.max(), im.sum()
	print bm.min(), bm.mean(), bm.max(), bm.sum()
	
	#from matplotlib import pyplot as plt
	#fig = plt.figure()
	#ax1 = fig.add_subplot(1, 2, 1)
	#ax2 = fig.add_subplot(1, 2, 2)
	#ax1.imshow(im)
	#ax2.imshow(bm)
	#plt.show()
	
	return beamImg.image(center=(MapSize,MapSize))


def deconvolve(aa, dataDict, aipyImg, MapSize=80, MapRes=0.50, MapWRes=0.10, chan=None, gain=0.3, maxIter=150, scales=[0,], verbose=True):
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
	xyz = aipyImg.get_eq(0.0, aa.lat*numpy.pi/180.0, center=(MapSize,MapSize))
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
	
	beamClean = estimateBeam(aa, 0.0, lat, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, baselines=dataDict['bls']['xx'])
	beamClean /= beamClean.max()
	
	#import pylab
	#pylab.imshow(beamClean)
	#pylab.show()
	
	# Fit a Guassian to the zenith beam response and use that for the restore beam
	profileX = beamClean[MapSize/2:3*MapSize/2,MapSize]
	profileY = beamClean[MapSize,MapSize/2:3*MapSize/2]
	goodX = numpy.where( profileX >= 0.5 )[0]
	goodY = numpy.where( profileY >= 0.5 )[0]
	print goodX.max()-goodX.min(), goodY.max()-goodY.min()
	width = ((goodX.max() - goodX.min()) + (goodY.max() - goodY.min())) / 2.0
	print 'width', width
	
	# Go!
	import pylab
	pylab.ion()
	for i in xrange(maxIter):
		# Find the location of the peak in the flux density
		peak = numpy.where( working == working.max() )
		peakX = peak[0][0]
		peakY = peak[1][0]
		peakV = working[peakX,peakY]
		print i, working.max(), peakV
		
		# Pixel coordinates to hour angle, dec.
		##peakHA = (HA[peakX-2:peakX+2,peakY-2:peakY+2]*working[peakX-2:peakX+2,peakY-2:peakY+2]).sum()
		##peakHA /= working[peakX-2:peakX+2,peakY-2:peakY+2].sum()
		##peakHA *= 180/numpy.pi / 15.0
		xStart = peakX - width/2
		xStop  = peakX + width/2
		yStart = peakY - width/2
		yStop  = peakY + width/2
		peakHA = (HA[xStart:xStop,yStart:yStop]*working[xStart:xStop,yStart:yStop]).sum() / working[xStart:xStop,yStart:yStop].sum() * 180/numpy.pi / 15.0
		print peakHA
		
		peakHA = HA[peakX, peakY] * 180/numpy.pi / 15.0
		print peakHA
		
		#peakDec = (dec[xStart:xStop,yStart:yStop]*working[xStart:xStop,yStart:yStop]).sum() / working[xStart:xStop,yStart:yStop].sum() * 180/numpy.pi
		
		##peakDec *= 180/numpy.pi
		peakDec = dec[peakX,peakY] * 180/numpy.pi
		#print peakDec
		#peakDec = (dec[xStart:xStop,yStart:yStop]*working[xStart:xStop,yStart:yStop]).sum() / working[xStart:xStop,yStart:yStop].sum() * 180/numpy.pi
		#print peakDec
		
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
		#template = beamClean*1.0
		#template = numpy.roll(template, -MapSize+peakX, axis=0)
		#template = numpy.roll(template, -MapSize+peakY, axis=1)
		#beams = [template,]
		try:
			beams = prevBeam[beamIndex]
		except KeyError:
			if verbose:
				print "               -> Computing beam(s)"
				
			beams = []
			for scale in scales:
				beam = estimateBeam(aa, peakHA, peakDec, 
								MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=chan, baselines=dataDict['bls']['xx'], scale=scale)
				beam /= beam.max()
				print "                  ", scale, beam.mean(), beam.min(), beam.max(), beam.sum()
				beams.append(beam)
			prevBeam[beamIndex] = beams
		
		# Calculate how much signal needs to be removed...
		bestScale = 0
		bestValue = 1e9
		for scale in xrange(len(scales)):
			temp = 1.0*working
			toRemove = gain*peakV*beams[scale]
			
			# And then remove it and add it into list of CLEAN components
			temp -= toRemove
			value = temp.std()
			print "-> ", scale, value
			if value < bestValue:
				bestScale = scale
				bestValue = value
				
		print " @ ", bestScale, scales[bestScale], bestValue
		toRemove = gain*peakV*beams[bestScale]
		working -= toRemove
		cleaned[peakX,peakY] += toRemove.sum()
		
		print working.min(), working.max(), working.std(), working.max()/working.std(), toRemove.sum()
		pylab.imshow(working, origin='lower')
		pylab.draw()
		
		if working.max()/working.std() < 3:
			break
			
	# Calculate what the restore beam should look like
	gauGen = gaussian2d(1.0, MapSize, MapSize, width/numpy.sqrt(8*numpy.log(2)), width/numpy.sqrt(8*numpy.log(2)))
	beamClean *= 0
	for i in xrange(beamClean.shape[0]):
		for j in xrange(beamClean.shape[1]):
			beamClean[i,j] = gauGen(i,j)
	beamClean /= beamClean.sum()
	
	# Restore
	conv = convolve(cleaned, beamClean, mode='same')
	conv = numpy.ma.array(conv, mask=mask)
	conv /= conv.max()
	conv *= (img.sum() - working.sum())
	
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
		#ax1.plot(profileX/profileX.max())
		#ax1.plot(profileY/profileY.max())
		#ax1.plot(beamClean[MapSize,MapSize/2:3*MapSize/2]/beamClean.max())
		print img.sum()
		
		d = ax2.imshow(conv, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(d, ax=ax2)
		ax2.set_title('CLEAN Comps.')
		print cleaned.sum(), conv.sum()
		
		e = ax3.imshow(working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(e, ax=ax3)
		ax3.set_title('Residuals')
		print working.sum()
		
		f = ax4.imshow(conv + working, extent=(1,-1,-1,1), origin='lower', interpolation='nearest')
		fig.colorbar(f, ax=ax4)
		ax4.set_title('Final')
		print (conv + working).sum()
		
		plt.show()
	
	# Return
	pylab.ioff()
	return conv + working


def buildSky(aa, jd, MapSize=80, MapRes=0.50, MapWRes=0.10):
	import aipy
	import math
	from lsl import skymap
	
	f = aa.get_afreqs()*1e9/1e6
	smap = skymap.SkyMapGSM(freqMHz=f.mean())
	pmap = skymap.ProjectedSkyMap(smap, aa.lat*180.0/math.pi, aa.long*180.0/math.pi, aa.get_jultime())
	
	az  = pmap.visibleAz*math.pi/180.0
	alt = pmap.visibleAlt*math.pi/180.0
	pwr = pmap.visiblePower
	
	img = aipy.img.ImgW(MapSize, MapRes, MapWRes)
	iTop = img.get_top(center=(MapSize,MapSize))
	iAz,iAlt = aipy.coord.top2azalt(iTop)
	
	output = numpy.zeros((MapSize*2, MapSize*2))
	for i in xrange(iAz.shape[0]):
		for j in xrange(iAlt.shape[1]):
			if iAlt.mask[i,j]:
				continue
				
			diff = numpy.sqrt( (iAz[i,j]-az)**2 + (iAlt[i,j]-alt)**2 )
			output[i,j] += pwr[diff.argmin()] * math.sin(iAlt[i,j])**1.5
			
	return output


def mem(aa, dataDict, aipyImg, MapSize=80, MapRes=0.50, MapWRes=0.10, chan=None, maxIter=150, verbose=True):
	import aipy
	
	image = aipyImg.image(center=(MapSize,MapSize))
	xyz = aipyImg.get_eq(aa.sidereal_time(), aa.lat, center=(MapSize,MapSize))
	top = aipyImg.get_top(center=(MapSize,MapSize))
	ra, dec = eq2radec(xyz)
	az,alt = aipy.coord.top2azalt(top)
	
	baselines = dataDict['bls']['xx']
	psfSrcs = {'z': aipy.amp.RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0)}
	simDict = simVis.buildSimData(aa, psfSrcs, jd=aa.get_jultime(), pols=['xx',], baselines=baselines)
	psf = utils.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, chan=[0,], pol='xx')
	psf = psf.image(center=(MapSize,MapSize))
	psf /= psf.max()
	
	# Fit a Guassian to the zenith beam response and use that for the restore beam
	profileX = psf[MapSize/2:3*MapSize/2,MapSize]
	profileY = psf[MapSize,MapSize/2:3*MapSize/2]
	goodX = numpy.where( profileX >= 0.5 )[0]
	goodY = numpy.where( profileY >= 0.5 )[0]
	print goodX.max()-goodX.min(), goodY.max()-goodY.min()
	width = ((goodX.max() - goodX.min()) + (goodY.max() - goodY.min())) / 2.0
	print 'width', width
	
	# Calculate what the restore beam should look like
	gauGen = gaussian2d(1.0, MapSize, MapSize, width/numpy.sqrt(8*numpy.log(2)), width/numpy.sqrt(8*numpy.log(2)))
	psf *= 0
	for i in xrange(psf.shape[0]):
		for j in xrange(psf.shape[1]):
			psf[i,j] = gauGen(i,j)
	psf /= psf.sum()
	
	d_i = image.flatten()
	q = numpy.sqrt( (psf**2).sum() )
	minus_two_q = -2*q
	two_q_sq = 2*q**2
	
	#model = buildSky(aa, aa.get_jultime(), MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes)
	#model *= image.sum() / model.sum()
	model = numpy.ones(image.shape, dtype=image.dtype)
	for n,src in simVis.srcs.iteritems():
		src.compute(aa)
		if src.alt > 0:
			print n, src.az, src.alt, src.jys.mean()
			d = numpy.sqrt( (az-src.az)**2 + (alt-src.alt)**2 )
			best = numpy.where( d == d.min() )
			newPSF = numpy.roll(psf, best[0][0]-MapSize, axis=0)
			newPSF = numpy.roll(newPSF, best[1][0]-MapSize, axis=1)
			model += src.jys.mean()/1e2 * numpy.sin(src.alt) * newPSF
	model[numpy.where(ra.mask == 1)] = 0
	model *= image.mean() / psf.sum()
	
	#import pylab
	#pylab.imshow(model, origin='lower')
	#pylab.show()
	
	valid = numpy.where( ra.mask.flatten() == 0 )[0]
	valid2 = numpy.where( ra.mask == 0 )
	print 'V', valid.shape
	
	var0 = 0
	Nvar0 = d_i.size * var0
	m_i = model.flatten()
	print 'I', model.shape, m_i.shape, image.shape
	def next_step(b_i, alpha, verbose=True):
		b_i.shape = image.shape
		bSrcs = {}
		for i in xrange(b_i.shape[0]):
			for j in xrange(b_i.shape[1]):
				if dec.mask[i,j]:
					continue
				nm = '%i-%i' % (i,j)
				bSrcs[nm] = aipy.amp.RadioFixedBody(ra[i,j], dec[i,j], name=nm, jys=b_i[i,j], index=0)
				
		import time
		t0 = time.time()
		simDict = simVis.buildSimData(aa, bSrcs, jd=aa.get_jultime(), pols=['xx',], baselines=baselines)
		b_i_conv_ker = utils.buildGriddedImage(simDict, MapSize=MapSize, MapRes=MapRes, MapWRes=MapWRes, pol='xx')
		b_i_conv_ker = b_i_conv_ker.image(center=(MapSize,MapSize))
		b_i_conv_ker *= image.sum()/b_i_conv_ker.sum()
		print 'Sim. Time:  %.1f s' % (time.time() - t0)
		
		b_i.shape = m_i.shape
		diff = (image - b_i_conv_ker).flatten()
		chi2 = numpy.dot(diff[valid],diff[valid]) - Nvar0
		g_chi2 = minus_two_q*(diff[valid])
		g_J = (-numpy.log(b_i[valid]/m_i[valid]) - 1) - alpha * g_chi2
		g_J = numpy.where( numpy.isfinite(g_J), g_J, 100 )
		gg_J = (-1/b_i[valid]) - alpha * two_q_sq
		# Define dot product using the metric of -gg_J^-1
		def dot(x, y):
			return (x*y/-gg_J).sum()
		score = dot(g_J, g_J) / dot(1,1)
		d_alpha = (chi2 + dot(g_chi2,g_J)) / dot(g_chi2,g_chi2)
		print chi2, dot(g_chi2,g_J), dot(g_chi2,g_chi2)
		
		# For especially clean images, gg_J in some components can go to 0.
		# This check makes lsq a little slower for most images, though...
		d_b_i = numpy.where(abs(gg_J) > 0, -1/gg_J * (g_J - d_alpha*g_chi2), 0)
		if verbose:
			print '    score', score, 'fit', numpy.dot(diff,diff)
			print '    alpha', alpha, 'd_alpha', d_alpha
			
		return d_b_i, d_alpha, score, b_i_conv_ker
		
	alpha = 0.
	b_i = m_i.copy()
	info = {'success':True, 'term':'maxiter', 'var0':var0, 'tol':1e-3}
	for i in range(maxIter):
		if verbose:
			print 'Step %d:' % i
		d_b_i, d_alpha, score, b_i_conv_ker = next_step(b_i, alpha, verbose=verbose)
		
		import pylab
		pylab.ion()
		pylab.subplot(2, 2, 1)
		pylab.imshow(image, origin='lower')
		pylab.subplot(2, 2, 2)
		pylab.imshow(b_i_conv_ker, origin='lower')
		pylab.subplot(2, 2, 3)
		pylab.imshow(image-b_i_conv_ker, origin='lower')
		pylab.subplot(2, 2, 4)
		temp = b_i.copy()
		temp.shape = image.shape
		pylab.imshow(temp, origin='lower')
		pylab.draw()
		
		if score < 1e-3 and score > 0:
			info['term'] = 'tol'
			break
		elif score > 1e10 or numpy.isnan(score) or score <= 0:
			info.update({'term':'divergence', 'success':False})
			break
		print b_i.size, d_b_i.size
		b_i[valid] = numpy.clip(b_i[valid] + 0.1 * d_b_i, numpy.finfo(numpy.float).tiny , numpy.Inf)
		alpha += 0.1 * d_alpha
		
	b_i.shape = image.shape
	info.update({'res':b_i_conv_ker, 'score': score, 'alpha': alpha, 'iter':i+1})
	
	import pylab
	pylab.ioff()
	pylab.subplot(1, 3, 1)
	pylab.imshow(image, origin='lower')
	pylab.subplot(1, 3, 2)
	pylab.imshow(b_i_conv_ker, origin='lower')
	pylab.subplot(1, 3, 3)
	pylab.imshow(b_i, origin='lower')
	pylab.draw()
	
	return b_i, info