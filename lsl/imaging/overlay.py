# -*- coding: utf-8 -*-

"""
Module that provides a variety of overlays for all-sky images.  These overlays include:
 * the locations and names of sources, 
 * the horizon, 
 * a graticle showing lines of constant RA and dec., and
 * a graticle showing lines of constant azimuth and elevation.

All of the functions in this module accept a matplotlib axes instances that is used for
plotting.

.. versionadded:: 1.0.1
"""

import aipy
import numpy


__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ["sources", "horizon", "graticleRADec", "graticleAzEl", 
		 "__version__", "__revision__", "__all__"]


def sources(ax, aa, srcs, label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot the
	locations of the srcs given in the 'srcs' dictionary.
	"""
	
	# Compute the positions of major sources and label the images
	compSrc = {}
	for name,src in srcs.iteritems():
		src.compute(aa)
		top = src.get_crds(crdsys='top', ncrd=3)
		az, alt = aipy.coord.top2azalt(top)
		compSrc[name] = [az, alt]
		if alt <= 0:
			continue
		ax.plot(top[0], top[1], marker='x', markerfacecolor='None', markeredgecolor='w', 
				linewidth=10.0, markersize=10)
		if label:
			ax.text(top[0], top[1], name, color='white', size=12)


def horizon(ax):
	"""
	For a matplotlib axis instance showing an image of the sky, plot the horizon.
	"""
	
	# Add in the horizon
	x = numpy.zeros(361)
	y = numpy.zeros(361)
	for i in xrange(361):
		xyz = aipy.coord.azalt2top([i*numpy.pi/180.0, 0])
		x[i] = xyz[0]
		y[i] = xyz[1]
	ax.plot(x, y, color='white')


def graticleRADec(ax, aa, label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot lines of
	constant declinate and RA.  Declinations are spaced at 20 degree intervals
	and RAs are spaced at 2 hour intervals.
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
			eq = aipy.coord.radec2eq((-aa.sidereal_time() + ra*numpy.pi/180,dec*numpy.pi/180))
			xyz = numpy.dot(aipy.coord.eq2top_m(0, aa.lat), eq)
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
		
		eq = aipy.coord.radec2eq((-aa.sidereal_time() + aa.sidereal_time(),(dec+5)*numpy.pi/180))
		xyz = numpy.dot(aipy.coord.eq2top_m(0, aa.lat), eq)
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
			eq = aipy.coord.radec2eq((-aa.sidereal_time() + ra*numpy.pi/180,dec*numpy.pi/180))
			xyz = numpy.dot(aipy.coord.eq2top_m(0, aa.lat), eq)
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
		
		eq = aipy.coord.radec2eq((-aa.sidereal_time() + ra*numpy.pi/180,0))
		xyz = numpy.dot(aipy.coord.eq2top_m(0, aa.lat), eq)
		az,alt = aipy.coord.top2azalt(xyz)
		
		if alt > 20*numpy.pi/180 and label:
			ax.text(xyz[0], xyz[1], '%i$^h$' % (ra/15,), color='white')


def graticleAzEl(ax, aa, label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot lines of
	constant azimuth and elevation.  Elevations are spaced at 20 degree intervals
	and azimuths are spaced at 45 degree intervals
	"""
	
	# Lines of constant elevation
	els = range(0, 90, 20)
	
	for el in els:
		x = numpy.zeros(361)
		y = numpy.zeros(361)
		for i in xrange(361):
			xyz = aipy.coord.azalt2top([i*numpy.pi/180.0, el*numpy.pi/180.0])
			x[i] = xyz[0]
			y[i] = xyz[1]
			
		ax.plot(x, y, color='white')
		if el > 0:
			ax.text(x[160], y[160], '%i$^\circ$' % el, color='white')
			
	# Lines of constant azimuth
	azs = range(0, 360, 45)
	
	for az in azs:
		x = numpy.zeros(81)
		y = numpy.zeros(81)
		
		for i in xrange(81):
			xyz = aipy.coord.azalt2top([az*numpy.pi/180.0, i*numpy.pi/180.0])
			x[i] = xyz[0]
			y[i] = xyz[1]
			
		ax.plot(x, y, color='white')
		ax.text(x[45], y[45], '%i$^\circ$' % az, color='white')