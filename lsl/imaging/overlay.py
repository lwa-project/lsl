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

.. versionchanged:: 3.0.0
    Switch to using imaging.utils.ImgWPlus for all image coordinate info.
"""

import aipy
import ephem
import numpy as np

from astropy import units as astrounits
from astropy.time import Time as AstroTime
from astropy.coordinates import FK5, AltAz

from lsl import astro
from lsl.imaging.utils import ImgWPlus

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


def sources(ax, gimg, srcs, phase_center='z', label=True, marker='x', color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot the
    locations of the srcs given in the 'srcs' dictionary.
    
    .. versionchanged:: 3.0.0
        Switch to using imaging.utils.ImgWPlus for all image coordinate info.
    """
    
    # Setup
    mjd = gimg.mjd
    wcs = gimg.wcs
    antennaarray = gimg.antennaarray
    
    # Compute the positions of major sources and label the images
    old_jultime = antennaarray.get_jultime()*1.0
    antennaarray.set_jultime(mjd + astro.MJD_OFFSET)
    for name,src in srcs.items():
        src.compute(antennaarray)
        x, y = wcs.all_world2pix(src.ra*180/np.pi, src.dec*180/np.pi, 0)
        
        if src.alt >= 0:
            ax.plot(x, y, marker=marker, markerfacecolor='None', markeredgecolor=color, 
                linewidth=10.0, markersize=10)
            if label:
                ax.text(x, y, name, color=color, size=12)
    antennaarray.set_jultime(old_jultime)


def horizon(ax, gimg, elevation_cut=1e-3, color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot the horizon.
    
    .. versionchanged:: 3.0.0
        Switch to using imaging.utils.ImgWPlus for all image coordinate info.
    
    .. versionchanged:: 1.1.0
        Added a new argument for the AntennaArray instance to provide a 
        uniform default call for all functions.
    """
    
    # Setup
    mjd = gimg.mjd
    wcs = gimg.wcs
    el = gimg.antennaarray.earth_location
    
    # Find the horizon (well elevation of 0.001 deg)
    ot = AstroTime(mjd, format='mjd', scale='utc')
    tc = AltAz(np.arange(361)*astrounits.deg, np.ones(361)*elevation_cut*astrounits.deg,
               location=el, obstime=ot)
    eq = tc.transform_to(FK5(equinox=ot))
    
    x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
    ax.plot(x, y, color=color)


def graticule_radec(ax, gimg, label=True, color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot lines of
    constant declinate and RA.  Declinations are spaced at 20 degree intervals
    and RAs are spaced at 2 hour intervals.
    
    .. versionchanged:: 3.0.0
        Switch to using imaging.utils.ImgWPlus for all image coordinate info.
    """
    
    # Setup
    mjd = gimg.mjd
    wcs = gimg.wcs
    el = gimg.antennaarray.earth_location
    
    # Lines of constant dec.
    ot = AstroTime(mjd, format='mjd', scale='utc')
    for dec in range(-80, 90, 20):
        eq = FK5(np.arange(361)*astrounits.deg, np.ones(361)*dec*astrounits.deg,
                 equinox=ot)
        
        x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
        ax.plot(x, y, color=color, alpha=0.75)
        
        eq = FK5(gimg.wcs.wcs.crval[0]*astrounits.deg, (dec+5)*astrounits.deg,
                 equinox=ot)
        tc = eq.transform_to(AltAz(location=el, obstime=ot))
        
        if tc.alt > 15*astrounits.deg and label:
            x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
            ax.text(x, y, r'%+i$^\circ$' % dec, color=color)
            
    # Lines of constant RA
    for ra in range(0, 360, 30):
        eq = FK5(np.ones(161)*ra*astrounits.deg, (np.arange(161)-80)*astrounits.deg,
                 equinox=ot)
        
        x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
        ax.plot(x, y, color=color, alpha=0.75)
        
        eq = FK5((ra-5)*astrounits.deg, '0deg',
                 equinox=ot)
        tc = eq.transform_to(AltAz(location=el, obstime=ot))
        
        if tc.alt > 20*astrounits.deg and label:
            x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
            ax.text(x, y, '%i$^h$' % (ra/15,), color=color)


def graticule_azalt(ax, gimg, label=True, color='white'):
    """
    For a matplotlib axis instance showing an image of the sky, plot lines of
    constant azimuth and elevation.  Elevations are spaced at 20 degree intervals
    and azimuths are spaced at 45 degree intervals.

    .. versionchanged:: 3.0.0
        Switch to using imaging.utils.ImgWPlus for all image coordinate info.
    """
    
    # Setup
    mjd = gimg.mjd
    wcs = gimg.wcs
    el = gimg.antennaarray.earth_location
    
    # Lines of constant elevation
    ot = AstroTime(mjd, format='mjd', scale='utc')
    for alt in range(0, 90, 20):
        tc = AltAz(np.arange(361)*astrounits.deg, np.ones(361)*alt*astrounits.deg,
                   location=el, obstime=ot)
        eq = tc.transform_to(FK5(equinox=ot))
        
        x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
        ax.plot(x, y, color=color, alpha=0.75)
        
        if label:
            valid = np.where( np.isfinite(x) & np.isfinite(y) )[0]
            if valid.size > 10:
                pos = valid.size // 2 - valid.size // 5
                ax.text(x[valid[pos]], y[valid[pos]], r'%i$^\circ$' % alt, color=color)
        
    # Lines of constant azimuth
    for az in range(0, 360, 45):
        tc = AltAz(np.ones(161)*az*astrounits.deg, (np.arange(161)-80)*astrounits.deg,
                   location=el, obstime=ot)
        eq = tc.transform_to(FK5(equinox=ot))
        
        x, y = wcs.all_world2pix(eq.ra.deg, eq.dec.deg, 0)
        ax.plot(x, y, color=color, alpha=0.75)
        
        if label:
            valid = np.where( np.isfinite(x) & np.isfinite(y) )[0]
            if valid.size > 10:
                pos = valid.size // 2 - valid.size // 5
                ax.text(x[valid[pos]], y[valid[pos]], r'%i$^\circ$' % az, color=color)
