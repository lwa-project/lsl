"""
Deconvolution support for images made with :func:`lsl.imaging.utils.build_gridded_image`.

.. versionchanged:: 3.0.0
    Switched from AntennaArray to imaging.utils.ImgWPlus for all image coordinate info
"""

import numpy as np
import logging
from aipy.fit import RadioFixedBody
from scipy.signal import fftconvolve as convolve

from lsl.sim.vis import build_sim_data
from lsl.imaging import utils
from lsl.astro import MJD_OFFSET, deg_to_dms, deg_to_hms
from lsl.statistics.robust import std as rStd
from lsl.misc.mathutils import gaussian2d

from lsl.logger import LSL_LOGGER

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.6'
__all__ = ['clean', 'clean_sources', 'lsq']


def _interpolate(data, peak_x, peak_y):
    x1 = int(peak_x)
    x2 = x1 + 1
    y1 = int(peak_y)
    y2 = y1 + 1
    
    x = peak_x
    if x < 0:
        raise IndexError("Invalid x")
    y = peak_y
    if y < 0:
        raise IndexError("Invalid y")
        
    dataPrime  = data[x1,y1]*(x2-x)*(y2-y)
    dataPrime += data[x2,y1]*(x-x1)*(y2-y)
    dataPrime += data[x1,y2]*(x2-x)*(y-y1)
    dataPrime += data[x2,y2]*(x-x1)*(y-y1)
    dataPrime /= (x2-x1)*(y2-y1)
    return dataPrime


def _fit_gaussian(data):
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
        return lambda x,y: height*np.exp(
                    -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

    def moments(data):
        """
        Returns (height, x, y, width_x, width_y) the gaussian parameters 
        of a 2D distribution by calculating its moments
        
        From:  http://wiki.scipy.org/Cookbook/FittingData
        """
        
        total = data.sum()
        X, Y = np.indices(data.shape)
        x = (X*data).sum()/total
        y = (Y*data).sum()/total
        col = data[:, int(y)]
        width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
        row = data[int(x), :]
        width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
        height = data.max()
        return height, x, y, width_x, width_y

    def fitgaussian(data, params=None):
        """
        Returns (height, x, y, width_x, width_y) the gaussian parameters 
        of a 2D distribution found by a fit
        
        From:  http://wiki.scipy.org/Cookbook/FittingData
        """
        
        if params is None:
            params = moments(data)
        errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                data)
        p, success = leastsq(errorfunction, params)
        return p
        
    level = 0.5
    params = None
    while level > 0.1:
        data_clipped = np.where(data/data.max() > level, data, 0)
        params = fitgaussian(data_clipped, params=params)
        level /= 2.0
        
    #fit = gaussian(*params)
    #import pylab
    #pylab.matshow(data, cmap=pylab.cm.gist_earth_r)
    #pylab.contour(fit(*np.indices(data.shape)), cmap=pylab.cm.copper)
    #pylab.show()
    
    return params


def clean(dataDict, gimg, input_image=None, size=80, res=0.50, wres=0.10, pol='XX', chan=None, gain=0.2, max_iter=150, sigma=3.0, verbose=True, plot=False):
    """
    Given a AIPY antenna array instance, a data dictionary, and an AIPY ImgW 
    instance filled with data, return a deconvolved image.  This function 
    uses a CLEAN-like method that computes the array beam for each peak in 
    the flux.  Thus the CLEAN loop becomes:
      1.  Find the peak flux in the residual image
      2.  Compute the systems response to a point source at that location
      3.  Remove the scaled portion of this beam from the residuals
      4.  Go to 1.
    
    CLEAN tuning parameters:
      * gain - CLEAN loop gain (default 0.2)
      * max_iter - Maximum number of iterations (default 150)
      * sigma - Threshold in sigma to stop cleaning (default 3.0)
    
    .. versionchanged:: 3.0.0
        Switched from AntennaArray to imaging.utils.ImgWPlus for all image coordinate info
    """
    
    # Setup
    mjd = gimg.mjd
    aa = gimg.antennaarray
    aa.set_jultime(gimg.mjd + MJD_OFFSET)
    
    # Sort out the channels to work on
    if chan is None:
        chan = range(dataDict.freq.size)
        
    # Get a grid of right ascensions and dec values for the image we are working with
    ra, dec = utils.get_image_radec(gimg)
    az, alt = utils.get_image_azalt(gimg)
    
    # Get the list of baselines to generate visibilites for
    baselines = dataDict.baselines
    
    # Get the actual image out of the ImgW instance
    if input_image is None:
        img = gimg.image()
        imgSize = img.shape[0]	# should be square
        
        img = np.roll(img, imgSize//2, axis=0)
        img = np.roll(img, imgSize//2, axis=1)
    else:
        img = input_image*1.0
        imgSize = img.shape[0]	# should be square
        
    # Setup the arrays to hold the point sources and the residual.
    cleaned = np.zeros_like(img)
    working = np.zeros_like(img)
    working += img
    
    # Setup the dictionary that will hold the beams as they are computed
    prevBeam = {}
    
    # Estimate the zenith beam response
    psfSrc = {'z': RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0, epoch=aa.date)}
    psfDict = build_sim_data(aa, psfSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flat_response=True)
    psf = utils.build_gridded_image(psfDict, size=size, res=res, wres=wres, chan=chan, pol=pol, verbose=verbose)
    psf = psf.image(center=(imgSize//2,imgSize//2))
    psf /= psf.max()
    
    # Fit a Guassian to the zenith beam response and use that for the restore beam
    beamCutout = psf[size//2:3*size//2, size//2:3*size//2]
    beamCutout = np.where( beamCutout > 0.0, beamCutout, 0.0 )
    h, cx, cy, sx, sy = _fit_gaussian( beamCutout )
    gauGen = gaussian2d(1.0, size/2+cx, size/2+cy, sx, sy)
    FWHM = int( round( (sx+sy)/2.0 * 2.0*np.sqrt(2.0*np.log(2.0)) ) )
    beamClean = psf * 0.0
    for i in range(beamClean.shape[0]):
        for j in range(beamClean.shape[1]):
            beamClean[i,j] = gauGen(i,j)
    beamClean /= beamClean.sum()
    convMask = np.where(np.isfinite(ra), False, True)
    
    # Go!
    if plot:
        import pylab
        from matplotlib import pyplot as plt
        
        pylab.ion()
        
    exitStatus = 'iteration limit'
    for i in range(max_iter):
        # Find the location of the peak in the flux density
        peak = np.where( working == working.max() )
        peak_x = peak[0][0]
        peak_y = peak[1][0]
        peakV = working[peak_x,peak_y]
        
        # Optimize the location
        peakParams = _fit_gaussian(working[peak_x-FWHM//2:peak_x+FWHM//2+1, peak_y-FWHM//2:peak_y+FWHM//2+1])
        peakVO = peakParams[0]
        peak_xO = peak_x - FWHM//2 + peakParams[1]
        peak_yO = peak_y - FWHM//2 + peakParams[2]
        
        # Quantize to try and keep the computation down without over-simplifiying things
        subpixelationLevel = 5
        peak_xO = round(peak_xO*subpixelationLevel)/float(subpixelationLevel)
        peak_yO = round(peak_yO*subpixelationLevel)/float(subpixelationLevel)
        
        # Pixel coordinates to right ascension, dec.
        try:
            peakRA = _interpolate(ra, peak_xO, peak_yO)
        except IndexError:
            peakRA = ra[peak_x, peak_y]
        try:
            peakDec = _interpolate(dec, peak_xO, peak_yO)
        except IndexError:
            peakDec = dec[peak_x, peak_y]
            
        # Pixel coordinates to az, el
        try:
            peakAz = _interpolate(az, peak_xO, peak_yO)
        except IndexError:
            peakAz = az[peak_x, peak_y]
        try:
            peakAlt = _interpolate(alt, peak_x, peak_y)
        except IndexError:
            peakAlt = alt[peak_x, peak_y]
            
        if verbose:
            currRA  = deg_to_hms(peakRA * 180/np.pi)
            currDec = deg_to_dms(peakDec * 180/np.pi)
            currAz  = deg_to_dms(peakAz * 180/np.pi)
            currAlt = deg_to_dms(peakAlt * 180/np.pi)
            
            print("Iteration %i:  Log peak of %.3f at row: %i, column: %i" % (i+1, np.log10(peakV), peak_x, peak_y))
            print("               -> RA: %s, Dec: %s" % (currRA, currDec))
            print("               -> az: %s, el: %s" % (currAz, currAlt))
        LSL_LOGGER.info("Iteration %i:  Log peak of %.3f at row: %i, column: %i" % (i+1, np.log10(peakV), peak_x, peak_y))
        LSL_LOGGER.info("               -> RA: %s, Dec: %s" % (currRA, currDec))
        LSL_LOGGER.info("               -> az: %s, el: %s" % (currAz, currAlt))
        
        # Check for the exit criteria
        if peakV < 0:
            exitStatus = 'peak value is negative'
            
            break
        
        # Find the beam index and see if we need to compute the beam or not
        beamIndex = (int(peak_xO*subpixelationLevel), int(peak_yO*subpixelationLevel))
        try:
            beam = prevBeam[beamIndex]
            
        except KeyError:
            if verbose:
                print("               -> Computing beam(s)")
            LSL_LOGGER.debug("               -> Computing beam(s)")
            
            beamSrc = {'Beam': RadioFixedBody(peakRA, peakDec, jys=1.0, index=0, epoch=aa.date)}
            beamDict = build_sim_data(aa, beamSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flat_response=True)
            beam = utils.build_gridded_image(beamDict, size=size, res=res, wres=wres, chan=chan, pol=pol, verbose=verbose)
            beam = beam.image(center=(imgSize//2,imgSize//2))
            beam /= beam.max()
            if verbose:
                print("                  ", beam.mean(), beam.min(), beam.max(), beam.sum())
            LSL_LOGGER.debug("                  ", beam.mean(), beam.min(), beam.max(), beam.sum())
            
            prevBeam[beamIndex] = beam
            if verbose:
                print("               -> Beam cache contains %i entries" % len(prevBeam.keys()))
            LSL_LOGGER.debug("               -> Beam cache contains %i entries" % len(prevBeam.keys()))
            
        # Calculate how much signal needs to be removed...
        toRemove = gain*peakV*beam
        working -= toRemove
        asum = 0.0
        for l in range(int(peak_xO), int(peak_xO)+2):
            if l > peak_xO:
                side1 = (peak_xO+0.5) - (l-0.5)
            else:
                side1 = (l+0.5) - (peak_xO-0.5)
                
            for m in range(int(peak_yO), int(peak_yO)+2):
                if m > peak_yO:
                    side2 = (peak_yO+0.5) - (m-0.5)
                else:
                    side2 = (m+0.5) - (peak_yO-0.5)
                    
                area = side1*side2
                asum += area
                #print('II', l, m, area, asum)
                cleaned[l,m] += gain*area*peakV
                
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
    print("Exited after %i iterations with status '%s'" % (i+1, exitStatus))
    LSL_LOGGER.info("Exited after %i iterations with status '%s'" % (i+1, exitStatus))
    
    # Restore
    conv = convolve(cleaned, beamClean, mode='same')
    conv = np.ma.masked_array(conv, mask=convMask, dtype=conv.dtype)
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


def clean_sources(dataDict, gimg, srcs, input_image=None, size=80, res=0.50, wres=0.10, pol='XX', chan=None, gain=0.1, max_iter=150, sigma=2.0, verbose=True, plot=False):
    """
    Given a AIPY antenna array instance, a data dictionary, an AIPY ImgW 
    instance filled with data, and a dictionary of sources, return the CLEAN
    components and the residuals map.  This function uses a CLEAN-like method
    that computes the array beam for each peak in the flux.  Thus the CLEAN 
    loop becomes: 
      1.  Find the peak flux in the residual image
      2.  Compute the systems response to a point source at that location
      3.  Remove the scaled portion of this beam from the residuals
      4.  Go to 1.
    
    This function differs from clean() in that it only cleans localized 
    regions around each source rather than the whole image.  This is
    intended to help the mem() function along.
    
    CLEAN tuning parameters:
      * gain - CLEAN loop gain (default 0.1)
      * max_iter - Maximum number of iterations (default 150)
      * sigma - Threshold in sigma to stop cleaning (default 2.0)
    
    .. versionchanged:: 3.0.0
        Switched from AntennaArray to imaging.utils.ImgWPlus for all image coordinate info
    """
    
    
    # Setup
    mjd = gimg.mjd
    aa = gimg.antennaarray
    aa.set_jultime(gimg.mjd + MJD_OFFSET)
    
    # Sort out the channels to work on
    if chan is None:
        chan = range(dataDict.freq.size)
        
    # Get a grid of right ascensions and dec values for the image we are working with
    ra, dec = utils.get_image_radec(gimg)
    az, alt = utils.get_image_azalt(gimg)
    
    # Get the list of baselines to generate visibilites for
    baselines = dataDict.baselines
    
    # Get the actual image out of the ImgW instance
    if input_image is None:
        img = gimg.image()
        imgSize = img.shape[0]	# should be square
        
        img = np.roll(img, imgSize//2, axis=0)
        img = np.roll(img, imgSize//2, axis=1)
    else:
        img = input_image*1.0
        imgSize = img.shape[0]	# should be square
        
    # Setup the arrays to hold the point sources and the residual.
    cleaned = np.zeros_like(img)
    working = np.zeros_like(img)
    working += img
    
    # Setup the dictionary that will hold the beams as they are computed
    prevBeam = {}
    
    # Estimate the zenith beam response
    psfSrc = {'z': RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0, epoch=aa.date)}
    psfDict = build_sim_data(aa, psfSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flat_response=True)
    psf = utils.build_gridded_image(psfDict, size=size, res=res, wres=wres, chan=chan, pol=pol, verbose=verbose)
    psf = psf.image(center=(imgSize//2,imgSize//2))
    psf /= psf.max()
    
    # Fit a Guassian to the zenith beam response and use that for the restore beam
    beamCutout = psf[size//2:3*size//2, size//2:3*size//2]
    beamCutout = np.where( beamCutout > 0.0, beamCutout, 0.0 )
    h, cx, cy, sx, sy = _fit_gaussian( beamCutout )
    gauGen = gaussian2d(1.0, size/2+cx, size/2+cy, sx, sy)
    FWHM = int( round( (sx+sy)/2.0 * 2.0*np.sqrt(2.0*np.log(2.0)) ) )
    beamClean = psf * 0.0
    for i in range(beamClean.shape[0]):
        for j in range(beamClean.shape[1]):
            beamClean[i,j] = gauGen(i,j)
    beamClean /= beamClean.sum()
    convMask = np.where(np.isfinite(ra), False, True)
    
    # Go!
    if plot:
        import pylab
        from matplotlib import pyplot as plt
        
        pylab.ion()
        
    for name,src in srcs.items():
        # Make sure the source is up
        src.compute(aa)
        if verbose:
            print('Source: %s @ %s degrees altitude' % (name, src.alt))
        LSL_LOGGER.info('Source: %s @ %s degrees altitude' % (name, src.alt))
        if src.alt <= 10*np.pi/180.0:
            continue
            
        # Locate the approximate position of the source
        srcDist = (src.ra-ra)**2 + (src.dec-dec)**2
        srcPeak = np.where( srcDist == np.nanmin(srcDist) )
        
        # Define the clean box - this is fixed at 2*FWHM in width on each side
        rx0 = max([0, srcPeak[0][0] - FWHM//2])
        rx1 = min([rx0 + FWHM + 1, img.shape[0]])
        ry0 = max([0, srcPeak[1][0] - FWHM//2])
        ry1 = min([ry0 + FWHM + 1, img.shape[1]])
        
        # Define the background box - this lies outside the clean box and serves
        # as a reference for the background
        X, Y = np.indices(working.shape)
        R = np.sqrt( (X-srcPeak[0][0])**2 + (Y-srcPeak[1][0])**2 )
        bpad = 3
        background = np.where( (R <= FWHM+bpad) & (R > FWHM) )
        while len(background[0]) == 0 and bpad < img.shape[0]:
            bpad += 1
            background = np.where( (R <= FWHM+bpad) & (R > FWHM) )
            
        px0 = min(background[0])-1
        px1 = max(background[0])+2
        py0 = min(background[1])-1
        py1 = max(background[1])+2
        
        exitStatus = 'iteration'
        for i in range(max_iter):
            # Find the location of the peak in the flux density
            peak = np.where( working[rx0:rx1,ry0:ry1] == working[rx0:rx1,ry0:ry1].max() )
            peak_x = peak[0][0] + rx0
            peak_y = peak[1][0] + ry0
            peakV = working[peak_x,peak_y]
            
            # Optimize the location
            try:
                peakParams = _fit_gaussian(working[peak_x-FWHM//2:peak_x+FWHM//2+1, peak_y-FWHM//2:peak_y+FWHM//2+1])
            except IndexError:
                peakParams = [peakV, FWHM//2, FWHM//2]
            peakVO = peakParams[0]
            peak_xO = peak_x - FWHM//2 + peakParams[1]
            peak_yO = peak_y - FWHM//2 + peakParams[2]
            
            # Quantize to try and keep the computation down without over-simplifiying things
            subpixelationLevel = 5
            peak_xO = round(peak_xO*subpixelationLevel)/float(subpixelationLevel)
            peak_yO = round(peak_yO*subpixelationLevel)/float(subpixelationLevel)
            
            # Pixel coordinates to right ascension, dec.
            try:
                peakRA = _interpolate(ra, peak_xO, peak_yO)
            except IndexError:
                peak_xO, peak_yO = peak_x, peak_y
                peakRA = ra[peak_x, peak_y]
            try:
                peakDec = _interpolate(dec, peak_xO, peak_yO)
            except IndexError:
                peakDec = dec[peak_x, peak_y]
                
            # Pixel coordinates to az, el
            try:
                peakAz = _interpolate(az, peak_xO, peak_yO)
            except IndexError:
                peak_xO, peak_yO = peak_x, peak_y
                peakAz = az[peak_x, peak_y]
            try:
                peakAlt = _interpolate(alt, peak_x, peak_y)
            except IndexError:
                peakAlt = alt[peak_x, peak_y]
                
            if verbose:
                currRA  = deg_to_hms(peakRA * 180/np.pi)
                currDec = deg_to_dms(peakDec * 180/np.pi)
                currAz  = deg_to_dms(peakAz * 180/np.pi)
                currAlt = deg_to_dms(peakAlt * 180/np.pi)
                
                print("%s - Iteration %i:  Log peak of %.3f at row: %i, column: %i" % (name, i+1, np.log10(peakV), peak_x, peak_y))
                print("               -> RA: %s, Dec: %s" % (currRA, currDec))
                print("               -> az: %s, el: %s" % (currAz, currAlt))
            LSL_LOGGER.info("%s - Iteration %i:  Log peak of %.3f at row: %i, column: %i" % (name, i+1, np.log10(peakV), peak_x, peak_y))
            LSL_LOGGER.info("               -> RA: %s, Dec: %s" % (currRA, currDec))
            LSL_LOGGER.info("               -> az: %s, el: %s" % (currAz, currAlt))
            
            # Check for the exit criteria
            if peakV < 0:
                exitStatus = 'peak value is negative'
                
                break
            
            # Find the beam index and see if we need to compute the beam or not
            beamIndex = (int(peak_xO*subpixelationLevel), int(peak_yO*subpixelationLevel))
            try:
                beam = prevBeam[beamIndex]
                
            except KeyError:
                if verbose:
                    print("               -> Computing beam(s)")
                LSL_LOGGER.debug("               -> Computing beam(s)")
                
                beamSrc = {'Beam': RadioFixedBody(peakRA, peakDec, jys=1.0, index=0, epoch=aa.date)}
                beamDict = build_sim_data(aa, beamSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flat_response=True)
                beam = utils.build_gridded_image(beamDict, size=size, res=res, wres=wres, chan=chan, pol=pol, verbose=verbose)
                beam = beam.image(center=(imgSize//2,imgSize//2))
                beam /= beam.max()
                if verbose:
                    print("                  ", beam.mean(), beam.min(), beam.max(), beam.sum())
                LSL_LOGGER.debug("                  ", beam.mean(), beam.min(), beam.max(), beam.sum())
                
                prevBeam[beamIndex] = beam
                if verbose:
                    print("               -> Beam cache contains %i entries" % len(prevBeam.keys()))
                LSL_LOGGER.debug("               -> Beam cache contains %i entries" % len(prevBeam.keys()))
                
            # Calculate how much signal needs to be removed...
            toRemove = gain*peakV*beam
            working -= toRemove
            asum = 0.0
            for l in range(int(peak_xO), int(peak_xO)+2):
                if l > peak_xO:
                    side1 = (peak_xO+0.5) - (l-0.5)
                else:
                    side1 = (l+0.5) - (peak_xO-0.5)
                    
                for m in range(int(peak_yO), int(peak_yO)+2):
                    if m > peak_yO:
                        side2 = (peak_yO+0.5) - (m-0.5)
                    else:
                        side2 = (m+0.5) - (peak_yO-0.5)
                        
                    area = side1*side2
                    asum += area
                    #print('II', l, m, area, asum)
                    cleaned[l,m] += gain*area*peakV
                    
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
                
            if np.abs(np.max(working[rx0:rx1,ry0:ry1])-np.median(working[background]))/rStd(working[background]) <= sigma:
                exitStatus = 'peak is less than %.3f-sigma' % sigma
                
                break
                
        # Summary
        print("Exited after %i iterations with status '%s'" % (i+1, exitStatus))
        LSL_LOGGER.info("Exited after %i iterations with status '%s'" % (i+1, exitStatus))
        
    # Restore
    conv = convolve(cleaned, beamClean, mode='same')
    conv = np.ma.masked_array(conv, mask=convMask, dtype=conv.dtype)
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


def _minor_cycle(img, beam, gain=0.2, max_iter=150):
    """
    Function for performing the minor cycle part of lsq().
    """
    
    cleaned = img*0.0
    working = img*1.0
    
    for i in range(max_iter):
        # Find the location of the peak in the flux density
        aw = np.abs( working )
        peak = np.where( aw == aw.max() )
        peak_x = peak[0][0]
        peak_y = peak[1][0]
        peakV = working[peak_x,peak_y]
        
        if np.abs(peakV - working.mean()) < 2*working.std():
            break
            
        # Build the beam
        beam2 = np.roll(beam,  peak_x-beam.shape[0]//2, axis=0)
        beam2 = np.roll(beam2, peak_y-beam.shape[1]//2, axis=1) 
        
        # Calculate how much signal needs to be removed...
        toRemove = gain*peakV*beam2
        working -= toRemove
        cleaned[peak_x,peak_y] += toRemove.sum()
        
    return  cleaned + working


def lsq(dataDict, gimg, input_image=None, size=80, res=0.50, wres=0.10, pol='XX', chan=None, gain=0.05, max_iter=150, rtol=1e-9, verbose=True, plot=False):
    """
    Given a AIPY antenna array instance, a data dictionary, and an AIPY ImgW 
    instance filled with data, return a deconvolved image.  This function 
    implements a least squares deconvolution.
    
    Least squares tuning parameters:
      * gain - least squares loop gain (default 0.05)
      * max_iter - Maximum number of iteration (default 150)
      * rtol - Minimum change in the residual RMS between iterations
               (default 1e-9)
    
    .. versionchanged:: 3.0.0
        Switched from AntennaArray to imaging.utils.ImgWPlus for all image coordinate info
    """
    
    # Setup
    mjd = gimg.mjd
    aa = gimg.antennaarray
    aa.set_jultime(gimg.mjd + MJD_OFFSET)
    
    # Sort out the channels to work on
    if chan is None:
        chan = range(dataDict.freq.size)
        
    # Get a grid of right ascensions and dec values for the image we are working with
    ra, dec = utils.get_image_radec(gimg)
    
    # Get the list of baselines to generate visibilites for
    baselines = dataDict.baselines
    
    # Get the actual image out of the ImgW instance
    if input_image is None:
        img = gimg.image()
        imgSize = img.shape[0]	# should be square
        
        img = np.roll(img, imgSize//2, axis=0)
        img = np.roll(img, imgSize//2, axis=1)
    else:
        img = input_image*1.0
        imgSize = img.shape[0]	# should be square
        
    # Estimate the zenith beam response
    psfSrc = {'z': RadioFixedBody(aa.sidereal_time(), aa.lat, jys=1.0, index=0, epoch=aa.date)}
    psfDict = build_sim_data(aa, psfSrc, jd=aa.get_jultime(), pols=[pol,], chan=chan, baselines=baselines, flat_response=True)
    psf = utils.build_gridded_image(psfDict, size=size, res=res, wres=wres, chan=chan, pol=pol, verbose=verbose)
    psf = psf.image(center=(imgSize//2,imgSize//2))
    psf /= psf.max()
    
    # Fit a Guassian to the zenith beam response and use that for the restore beam
    beamCutout = psf[size//2:3*size//2, size//2:3*size//2]
    beamCutout = np.where( beamCutout > 0.0, beamCutout, 0.0 )
    h, cx, cy, sx, sy = _fit_gaussian( beamCutout )
    gauGen = gaussian2d(1.0, size/2+cx, size/2+cy, sx, sy)
    FWHM = int( round( (sx+sy)/2.0 * 2.0*np.sqrt(2.0*np.log(2.0)) ) )
    beamClean = psf * 0.0
    for i in range(beamClean.shape[0]):
        for j in range(beamClean.shape[1]):
            beamClean[i,j] = gauGen(i,j)
    beamClean /= beamClean.sum()
    convMask = np.where(np.isfinite(ra), False, True)
    
    # Build the initial model
    mdl = img*0 + img.max()
    mdl[np.where(mdl < 0)] = 0
    mdl[np.where(convMask)] = 0
    
    # Determine the overall image->model scale factor
    bSrcs = {}
    rChan = [chan[0], chan[-1]]
    bSrcs['zenith'] = RadioFixedBody(aa.sidereal_time(), aa.lat, name='zenith', jys=1, index=0)
    simDict = build_sim_data(aa, bSrcs, jd=aa.get_jultime(), pols=[pol,], chan=rChan, baselines=baselines, flat_response=True)
    simImg = utils.build_gridded_image(simDict, size=size, res=res, wres=wres, chan=rChan, pol=pol, verbose=verbose)
    simImg = simImg.image(center=(imgSize//2,imgSize//2))
    
    simToModel = 1.0 / simImg.max()
    modelToSim = simImg.max() / 1.0
    
    # Go!
    if plot:
        import pylab
        from matplotlib import pyplot as plt
        
        pylab.ion()
        
    rChan = [chan[0], chan[-1]]
    mdl *= modelToSim
    diff = img - mdl
    diffScaled = 0.0*diff/gain
    oldModel = mdl
    oldRMS = diff.std()*1e6
    oldDiff = diff*0.0
    rHist = []
    exitStatus = 'iteration'
    for k in range(max_iter):
        ## Update the model image but don't allow negative flux
        mdl += diffScaled * gain
        mdl[np.where( mdl <= 0 )] = 0.0
        
        ## Convert the model image to an ensemble of point sources for forward 
        ## modeling
        bSrcs = {}
        for i in range(mdl.shape[0]):
            for j in range(mdl.shape[1]):
                if convMask[i,j]:
                    continue
                if mdl[i,j] <= 0:
                    continue
                    
                nm = '%i-%i' % (i,j)
                bSrcs[nm] = RadioFixedBody(ra[i,j], dec[i,j], name=nm, jys=mdl[i,j], index=0, epoch=aa.date)
                
        ## Model the visibilities
        simDict = build_sim_data(aa, bSrcs, jd=aa.get_jultime(), pols=[pol,], chan=rChan, baselines=baselines, flat_response=True)
        
        ## Form the simulated image
        simImg = utils.build_gridded_image(simDict, size=size, res=res, wres=wres, chan=rChan, pol=pol, verbose=verbose)
        simImg = simImg.image(center=(imgSize//2,imgSize//2))
        
        ## Difference the image and the simulated image and scale it to the 
        ## model's peak flux
        diff = img - simImg
        diff2 = _minor_cycle(diff, beamClean, gain=0.1, max_iter=2000)
        
        ## Compute the RMS and create an appropriately scaled version of the model
        RMS = diff.std()
        mdl2 = mdl*modelToSim
        
        ## Status report
        if verbose:
            print("Iteration %i:  %i sources used, RMS is %.4e" % (k+1, len(bSrcs.keys()), RMS))
            print("               -> maximum residual: %.4e (%.3f%% of peak)" % (diff.max(), 100.0*diff.max()/img.max()))
            print("               -> minimum residual: %.4e (%.3f%% of peak)" % (diff.min(), 100.0*diff.min()/img.max()))
            print("               -> delta RMS: %.4e (%.3f%%)" % (RMS-oldRMS, 100.0*(RMS-oldRMS)/RMS))
        LSL_LOGGER.info("Iteration %i:  %i sources used, RMS is %.4e" % (k+1, len(bSrcs.keys()), RMS))
        LSL_LOGGER.info("               -> maximum residual: %.4e (%.3f%% of peak)" % (diff.max(), 100.0*diff.max()/img.max()))
        LSL_LOGGER.info("               -> minimum residual: %.4e (%.3f%% of peak)" % (diff.min(), 100.0*diff.min()/img.max()))
        LSL_LOGGER.info("               -> delta RMS: %.4e (%.3f%%)" % (RMS-oldRMS, 100.0*(RMS-oldRMS)/RMS))
        
        ## Make the cleaned residuals map ready for updating the model
        diff = diff2
        diffScaled = diff * simToModel
        
        ## Has the RMS gone up?  If so, it is time to exit.  But first, restore 
        ## the previous iteration
        if RMS-oldRMS > 0:
            mdl = oldModel
            diff = oldDiff
            
            exitStatus = 'residuals'
            
            break
            
        ## Is the RMS still changing in an acceptable manner?
        if abs(RMS-oldRMS) < rtol:
            # No need to go back a step			
            #mdl = oldModel
            #diff = oldDiff
            
            exitStatus = 'tolerance'
            
            break
            
        ## Save the current iteration as the previous state
        rHist.append(RMS)
        oldRMS = RMS
        oldModel = mdl
        oldDiff = diff
        
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
            pylab.plot(rHist)
            pylab.draw()
            
    # Summary
    print("Exited after %i iterations with status '%s'" % (k+1, exitStatus))
    LSL_LOGGER.info("Exited after %i iterations with status '%s'" % (k+1, exitStatus))
    
    # Restore
    conv = convolve(mdl2, beamClean, mode='same')
    conv = np.ma.masked_array(conv, mask=convMask, dtype=conv.dtype)
    
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
