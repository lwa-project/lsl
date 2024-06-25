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
import numpy as np

from lsl import astro

from lsl.misc import telemetry
telemetry.track_module()


__version__ = "0.3"
__all__ = ["sources", "horizon", "graticule_radec", "graticule_azalt"]


def _radec_of(antennaarray, az, alt, degrees=False):
    """
    Given a lsl.sim.vis.AntennaArray instance and an azimuth/elevation pair,
    find the RA/dec of that point in the FK5 frame (equinox=J2000).
    
    .. note:: If 'degrees' is True then the input azimuth/elevation is taken to
              be in degrees and the returned RA/dec is also in degrees.
    """
    
    if not degrees:
        # radians -> degrees
        az = az * 180/np.pi
        alt = alt * 180/np.pi
        
    hrz = astro.hrz_posn()
    hrz.alt = alt
    hrz.az = az % 360
    
    geo = astro.geo_posn()
    # radians -> degrees
    geo.lat = antennaarray.lat*180/np.pi
    geo.lng = antennaarray.lon*180/np.pi
    geo.elv = antennaarray.elev
    
    # az/el -> RA/dec
    equ = astro.get_equ_from_hrz(hrz, geo, antennaarray.date+astro.DJD_OFFSET)
    
    RA = equ.ra
    dec = equ.dec
    if not degrees:
        # degrees -> radians
        RA = RA * np.pi/180.0
        dec = dec * np.pi/180.0
        
    return RA, dec 


def sources(ax, antennaarray, srcs, phase_center='z', label=True, marker='x', color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot the
    locations of the srcs given in the 'srcs' dictionary.
    """
    
    # Get the phase center
    if phase_center != 'z':
        phase_center.compute(antennaarray)
        pcRA, pcDec = phase_center.ra, phase_center.dec
    else:
        pcRA, pcDec = antennaarray.sidereal_time(), antennaarray.lat
    rot = aipy.coord.eq2top_m(0, pcDec)
        
    # Compute the positions of major sources and label the images
    for name,src in srcs.items():
        src.compute(antennaarray)
        eq = aipy.coord.radec2eq((src.ra-pcRA, src.dec))
        top = np.dot(rot, eq)
        junk,alt = aipy.coord.top2azalt(top)
        if alt >= 0:
            ax.plot(top[0], top[1], marker=marker, markerfacecolor='None', markeredgecolor=color, 
                linewidth=10.0, markersize=10)
            if label:
                ax.text(top[0], top[1], name, color=color, size=12)


def horizon(ax, antennaarray, phase_center='z', color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot the horizon.
    
    .. versionchanged:: 1.1.0
        Added a new argument for the AntennaArray instance to provide a 
        uniform default call for all functions.
    """
        
    # Get the phase center
    if phase_center != 'z':
        phase_center.compute(antennaarray)
        pcRA, pcDec = phase_center.ra, phase_center.dec
    else:
        pcRA, pcDec = antennaarray.sidereal_time(), antennaarray.lat
    rot = aipy.coord.eq2top_m(0, pcDec)
    
    # Add in the horizon
    x = np.zeros(361) + np.nan
    y = np.zeros(361) + np.nan
    for i in range(361):
        ra, dec = _radec_of(antennaarray, i*np.pi/180.0, 0.0, degrees=False)
        eq = aipy.coord.radec2eq((ra-pcRA,dec))
        top = np.dot(rot, eq)
        junk,alt = aipy.coord.top2azalt(top)
        if alt >= -0.01:
            x[i] = top[0]
            y[i] = top[1]
    ax.plot(x, y, color=color)


def graticule_radec(ax, antennaarray, phase_center='z', label=True, color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot lines of
    constant declinate and RA.  Declinations are spaced at 20 degree intervals
    and RAs are spaced at 2 hour intervals.
    """
    
    # Get the phase center
    if phase_center != 'z':
        phase_center.compute(antennaarray)
        pcRA, pcDec = phase_center.ra, phase_center.dec
    else:
        pcRA, pcDec = antennaarray.sidereal_time(), antennaarray.lat
    rot = aipy.coord.eq2top_m(0, pcDec)
    
    # Lines of constant declination first
    decs = range(-80, 90, 20)
    ras = np.linspace(0, 360, 800)
    
    x = np.zeros(ras.size) + np.nan
    y = np.zeros(ras.size) + np.nan
    for dec in decs:
        x *= np.nan
        y *= np.nan
        
        # Loop over RA to compute the topocentric coordinates (used by the image) for
        # the lines.  Also, figure out the elevation for each point on the line so
        # we can mask those below the horizon
        for i,ra in enumerate(ras):
            eq = aipy.coord.radec2eq((ra*np.pi/180-pcRA, dec*np.pi/180))
            top = np.dot(rot, eq)
            junk,alt = aipy.coord.top2azalt(top)
            if alt >= -1e-5:
                x[i] = top[0]
                y[i] = top[1]
                
        ax.plot(x, y, color=color, alpha=0.75)
        
        eq = aipy.coord.radec2eq((pcRA-pcRA, (dec+5)*np.pi/180))
        top = np.dot(rot, eq)
        az,alt = aipy.coord.top2azalt(top)
        if alt > 15*np.pi/180 and label:
            ax.text(top[0], top[1], r'%+i$^\circ$' % dec, color=color)
            
    # Lines of constant RA			
    decs = np.linspace(-80, 80, 400)
    ras = range(0, 360, 30)
    
    x = np.zeros(decs.size) + np.nan
    y = np.zeros(decs.size) + np.nan
    for ra in ras:
        x *= np.nan
        y *= np.nan
        
        # Loop over dec to compute the topocentric coordinates (used by the image) for
        # the lines.  Also, figure out the elevation for each point on the line so
        # we can mask those below the horizon
        for i,dec in enumerate(decs):
            eq = aipy.coord.radec2eq((ra*np.pi/180-pcRA, dec*np.pi/180))
            top = np.dot(rot, eq)
            junk,alt = aipy.coord.top2azalt(top)
            if alt >= -1e-5:
                x[i] = top[0]
                y[i] = top[1]
                
        ax.plot(x, y, color=color, alpha=0.75)
        
        eq = aipy.coord.radec2eq((ra*np.pi/180-pcRA, 0))
        top = np.dot(rot, eq)
        az,alt = aipy.coord.top2azalt(top)
        if alt > 20*np.pi/180 and label:
            ax.text(top[0], top[1], '%i$^h$' % (ra/15,), color=color)


def graticule_azalt(ax, antennaarray, phase_center='z', label=True, color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot lines of
    constant azimuth and elevation.  Elevations are spaced at 20 degree intervals
    and azimuths are spaced at 45 degree intervals
    """
    
    # Get the phase center
    if phase_center != 'z':
        phase_center.compute(antennaarray)
        pcRA, pcDec = phase_center.ra, phase_center.dec
    else:
        pcRA, pcDec = antennaarray.sidereal_time(), antennaarray.lat
    rot = aipy.coord.eq2top_m(0, pcDec)
    
    # Lines of constant elevation
    els = range(0, 90, 20)
    
    x = np.zeros(361) + np.nan
    y = np.zeros(361) + np.nan
    for el in els:
        x *= np.nan
        y *= np.nan
        
        for i in range(361):
            ra, dec = _radec_of(antennaarray, i*np.pi/180.0, el*np.pi/180.0, degrees=False)
            eq = aipy.coord.radec2eq((ra-pcRA,dec))
            top = np.dot(rot, eq)
            junk,alt = aipy.coord.top2azalt(top)
            if alt >= -1e-5:
                x[i] = top[0]
                y[i] = top[1]
                
        ax.plot(x, y, color=color)
        
        if el > 0 or phase_center != 'z':
            valid = np.where( np.isfinite(x) & np.isfinite(y) )[0]
            pos = valid.size // 2 - valid.size // 5
            if valid.size > 10:
                ax.text(x[valid[pos]], y[valid[pos]], r'%i$^\circ$' % el, color=color)
            
    # Lines of constant azimuth
    azs = range(0, 360, 45)
    
    x = np.zeros(81) + np.nan
    y = np.zeros(81) + np.nan
    for az in azs:
        x *= np.nan
        y *= np.nan
        
        for i in range(81):
            ra, dec = _radec_of(antennaarray, az*np.pi/180.0, i*np.pi/180.0, degrees=False)
            eq = aipy.coord.radec2eq((ra-pcRA,dec))
            top = np.dot(rot, eq)
            junk,alt = aipy.coord.top2azalt(top)
            if alt >= -1e-5:
                x[i] = top[0]
                y[i] = top[1]
                
        ax.plot(x, y, color=color)
        
        valid = np.where( np.isfinite(x) & np.isfinite(y) )[0]
        pos = valid.size // 2 - valid.size // 5
        if valid.size > 10:
            ax.text(x[valid[pos]], y[valid[pos]], r'%i$^\circ$' % az, color=color)
