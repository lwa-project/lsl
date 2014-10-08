# -*- coding: utf-8 -*-

"""
Module that provides a variety of overlays for all-sky images.  These overlays
include:
 * the locations and names of sources, 
 * the horizon, 
 * a graticle showing lines of constant RA and dec., and
 * a graticle showing lines of constant azimuth and elevation.

All of the functions in this module accept a matplotlib axes instances that 
is used for plotting.

.. versionadded:: 1.0.1

.. versionchanged:: 1.1.0
	Added support for overlaying on images with non-zenith phase centers
"""

import aipy
import ephem
import numpy


__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ["sources", "horizon", "graticleRADec", "graticleAzEl", 
		 "__version__", "__revision__", "__all__"]


def _radec_of(aa, az, alt):
	sinDec = numpy.sin(alt)*numpy.sin(aa.lat) + numpy.cos(alt)*numpy.cos(aa.lat)*numpy.cos(az)
	dec = numpy.arcsin( numpy.clip(sinDec, -1, 1) )
	cosHA = (numpy.sin(alt) - numpy.sin(dec)*numpy.sin(aa.lat)) / (numpy.cos(dec)*numpy.cos(aa.lat))
	HA = numpy.arccos( numpy.clip(cosHA, -1, 1) )
	try:
		if numpy.sin(az) > 0:
			HA = 2*numpy.pi - HA
	except ValueError:
		toShift = numpy.where( numpy.sin(az) > 0 )
		HA[toShift] = 2*numpy.pi - HA[toShift]	
	RA = aa.sidereal_time() - HA
	return RA, dec 


def sources(ax, aa, srcs, phaseCenter='z', label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot the
	locations of the srcs given in the 'srcs' dictionary.
	"""
	
	# Get the phase center
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)
		pcRA, pcDec = phaseCenter.ra, phaseCenter.dec
	else:
		pcRA, pcDec = aa.sidereal_time(), aa.lat
	rot = aipy.coord.eq2top_m(0, pcDec)
		
	# Compute the positions of major sources and label the images
	for name,src in srcs.iteritems():
		src.compute(aa)
		eq = aipy.coord.radec2eq((src.ra-pcRA, src.dec))
		top = numpy.dot(rot, eq)
		junk,alt = aipy.coord.top2azalt(top)
		if alt >= 0:
			ax.plot(top[0], top[1], marker='x', markerfacecolor='None', markeredgecolor='w', 
				linewidth=10.0, markersize=10)
			if label:
				ax.text(top[0], top[1], name, color='white', size=12)


def horizon(ax, aa, phaseCenter='z'):
	"""
	For a matplotlib axis instance showing an image of the sky, plot the horizon.
	
	.. versionchanged:: 1.1.0
		Added a new argument for the AntennaArray instance to provide a 
		uniform default call for all functions.
	"""
		
	# Get the phase center
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)
		pcRA, pcDec = phaseCenter.ra, phaseCenter.dec
	else:
		pcRA, pcDec = aa.sidereal_time(), aa.lat
	rot = aipy.coord.eq2top_m(0, pcDec)
	
	# Add in the horizon
	x = numpy.zeros(361) + numpy.nan
	y = numpy.zeros(361) + numpy.nan
	for i in xrange(361):
		ra, dec = _radec_of(aa, i*numpy.pi/180.0, 0.0)
		eq = aipy.coord.radec2eq((ra-pcRA,dec))
		top = numpy.dot(rot, eq)
		junk,alt = aipy.coord.top2azalt(top)
		if alt >= -1e-5:
			x[i] = top[0]
			y[i] = top[1]
	ax.plot(x, y, color='white')


def graticleRADec(ax, aa, phaseCenter='z', label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot lines of
	constant declinate and RA.  Declinations are spaced at 20 degree intervals
	and RAs are spaced at 2 hour intervals.
	"""
	
	# Get the phase center
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)
		pcRA, pcDec = phaseCenter.ra, phaseCenter.dec
	else:
		pcRA, pcDec = aa.sidereal_time(), aa.lat
	rot = aipy.coord.eq2top_m(0, pcDec)
	
	# Lines of constant declination first
	decs = range(-80, 90, 20)
	ras = numpy.linspace(0, 360, 800)
	
	x = numpy.zeros(ras.size) + numpy.nan
	y = numpy.zeros(ras.size) + numpy.nan
	for dec in decs:
		x *= numpy.nan
		y *= numpy.nan
		
		# Loop over RA to compute the topocentric coordinates (used by the image) for
		# the lines.  Also, figure out the elevation for each point on the line so
		# we can mask those below the horizon
		for i,ra in enumerate(ras):
			eq = aipy.coord.radec2eq((ra*numpy.pi/180-pcRA, dec*numpy.pi/180))
			top = numpy.dot(rot, eq)
			junk,alt = aipy.coord.top2azalt(top)
			if alt >= -1e-5:
				x[i] = top[0]
				y[i] = top[1]
				
		ax.plot(x, y, color='white', alpha=0.75)
		
		eq = aipy.coord.radec2eq((pcRA-pcRA, (dec+5)*numpy.pi/180))
		top = numpy.dot(rot, eq)
		az,alt = aipy.coord.top2azalt(top)
		if alt > 15*numpy.pi/180 and label:
			ax.text(top[0], top[1], '%+i$^\circ$' % dec, color='white')
			
	# Lines of constant RA			
	decs = numpy.linspace(-80, 80, 400)
	ras = range(0,360,30)
	
	x = numpy.zeros(decs.size) + numpy.nan
	y = numpy.zeros(decs.size) + numpy.nan
	for ra in ras:
		x *= numpy.nan
		y *= numpy.nan
		
		# Loop over dec to compute the topocentric coordinates (used by the image) for
		# the lines.  Also, figure out the elevation for each point on the line so
		# we can mask those below the horizon
		for i,dec in enumerate(decs):
			eq = aipy.coord.radec2eq((ra*numpy.pi/180-pcRA, dec*numpy.pi/180))
			top = numpy.dot(rot, eq)
			junk,alt = aipy.coord.top2azalt(top)
			if alt >= -1e-5:
				x[i] = top[0]
				y[i] = top[1]
				
		ax.plot(x, y, color='white', alpha=0.75)
		
		eq = aipy.coord.radec2eq((ra*numpy.pi/180-pcRA, 0))
		top = numpy.dot(rot, eq)
		az,alt = aipy.coord.top2azalt(top)
		if alt > 20*numpy.pi/180 and label:
			ax.text(top[0], top[1], '%i$^h$' % (ra/15,), color='white')


def graticleAzEl(ax, aa, phaseCenter='z', label=True):
	"""
	For a matplotlib axis instance showing an image of the sky, plot lines of
	constant azimuth and elevation.  Elevations are spaced at 20 degree intervals
	and azimuths are spaced at 45 degree intervals
	"""
	
	# Get the phase center
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)
		pcRA, pcDec = phaseCenter.ra, phaseCenter.dec
	else:
		pcRA, pcDec = aa.sidereal_time(), aa.lat
	rot = aipy.coord.eq2top_m(0, pcDec)
	
	# Lines of constant elevation
	els = range(0, 90, 20)
	
	x = numpy.zeros(361) + numpy.nan
	y = numpy.zeros(361) + numpy.nan
	for el in els:
		x *= numpy.nan
		y *= numpy.nan
		
		for i in xrange(361):
			ra, dec = _radec_of(aa, i*numpy.pi/180.0, el*numpy.pi/180.0)
			eq = aipy.coord.radec2eq((ra-pcRA,dec))
			top = numpy.dot(rot, eq)
			junk,alt = aipy.coord.top2azalt(top)
			if alt >= -1e-5:
				x[i] = top[0]
				y[i] = top[1]
				
		ax.plot(x, y, color='white')
		
		if el > 0 or phaseCenter is not 'z':
			valid = numpy.where( numpy.isfinite(x) & numpy.isfinite(y) )[0]
			pos = valid.size / 2 - valid.size / 5
			if valid.size > 10:
				ax.text(x[valid[pos]], y[valid[pos]], '%i$^\circ$' % el, color='white')
			
	# Lines of constant azimuth
	azs = range(0, 360, 45)
	
	x = numpy.zeros(81) + numpy.nan
	y = numpy.zeros(81) + numpy.nan
	for az in azs:
		x *= numpy.nan
		y *= numpy.nan
		
		for i in xrange(81):
			ra, dec = _radec_of(aa, az*numpy.pi/180.0, i*numpy.pi/180.0)
			eq = aipy.coord.radec2eq((ra-pcRA,dec))
			top = numpy.dot(rot, eq)
			junk,alt = aipy.coord.top2azalt(top)
			if alt >= -1e-5:
				x[i] = top[0]
				y[i] = top[1]
				
		ax.plot(x, y, color='white')
		
		valid = numpy.where( numpy.isfinite(x) & numpy.isfinite(y) )[0]
		pos = valid.size / 2 - valid.size / 5
		if valid.size > 10:
			ax.text(x[valid[pos]], y[valid[pos]], '%i$^\circ$' % az, color='white')
