# -*- coding: utf-8 -*-

"""
Module for analyzing images.  Currently, this module supports:
  * estimating the position-dependent background level in the image
  * finding point sources

.. versionadded:: 1.1.0
"""

# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info > (3,):
    xrange = range
    
import math
import numpy
from scipy.signal import convolve, medfilt
from scipy.interpolate import bisplrep, bisplev

from lsl.misc import telemetry
telemetry.track_module()


__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ['estimate_background', 'find_point_sources']


def estimate_background(image, window=32):
    """
    Given a 2-D image, estimate and return the background a la SExtractor.
    
    This works by:
      1.  Dividing the image into a number of half-overlapped tiles that are
          'window' by 'window' pixels in size.
      2.  Computing the mean and standard deviation within the tile both
          with and without 3*sigma clipping applied.
      3.  Using the tile statistics, determine if the tile is empty, i.e.,
          the standard deviation has changed by less than 20%, or full.  
            a. If the tile is empty, use the clipped mean for the background.
            b. Otherwise, use 2.5*median - 1.5*mean.
      4.  Once all of the tiles have been processes, median filter them to 
          remove those that are anomalously high.
      5.  Build a bicubic spline to interpolate the background back to the
          original image size.
      6.  Evaluate the spline and return.
    
    For more information on the SExtractor method, see:
      * Bertin & Arnouts (1996, AApS, 317, 393)
      * https://www.astromatic.net/pubsvn/software/sextractor/trunk/doc/sextractor.pdf
      * http://www.astr.tohoku.ac.jp/~akhlaghi/sextractor_notes.html
    """
    
    # Figure out how large the background grid should be assuming a half-
    # tile overlap
    nX = 2 * image.shape[0] / window
    nY = 2 * image.shape[1] / window
    
    # Process the background grid
    backgroundX = numpy.zeros((nX,nY), dtype=numpy.float64)
    backgroundY = numpy.zeros((nX,nY), dtype=numpy.float64)
    backgroundBasis = numpy.zeros((nX,nY), dtype=image.dtype)
    for i in xrange(nX):
        for j in xrange(nY):
            ## Extract the tile area
            tile = image[i*window/2:(i+1)*window/2,j*window/2:(j+1)*window/2]
            
            ## Compute the mean and standard deviation - both with and
            ## without progressive clipping
            m0, s0 = numpy.mean(tile), numpy.std(tile)
            valid = numpy.where( numpy.abs(tile-m0) < 3*s0 )
            m1, s1 = numpy.mean(tile[valid]), numpy.std(tile[valid])
            
            ## If the standard deviation hasn't changed by more than 20% use
            ## it's an uncrowded field and we can use the clipped mean.  Other-
            ## wise estimate the mode.
            if abs(s1-s0)/s0 <= 0.2:
                backgroundBasis[i,j] = m1
            else:
                ### This is SExtractor specific and is a correction for the 
                ### clipped distributions that are being used in this step.
                backgroundBasis[i,j] = 2.5*numpy.median(tile[valid]) - 1.5*m1
                
            ## Save the center of this tile for interpolation
            backgroundX[i,j] = i*window/2 + window/4.0
            backgroundY[i,j] = j*window/2 + window/4.0
            
    # Deal with NaNs in the image
    good = numpy.where( numpy.isfinite(backgroundBasis) )
    bad = numpy.where( ~numpy.isfinite(backgroundBasis) )
    backgroundBasis[bad] = numpy.median(backgroundBasis[good])
    
    # Use a median filter to get rid of unusually high regions
    backgroundBasis = medfilt(backgroundBasis, 3)
    
    # Build up the interpolation function
    backgroundX, backgroundY = backgroundX.ravel(), backgroundY.ravel()
    backgroundBasis = backgroundBasis.ravel()
    backgroundInterp = bisplrep(backgroundX, backgroundY, backgroundBasis, 
                            xb=0, xe=image.shape[0], yb=0, ye=image.shape[1], kx=3, ky=3)
        
    # Evaluate
    x, y = numpy.arange(image.shape[0]), numpy.arange(image.shape[1])
    background = bisplev(x, y, backgroundInterp)
    
    # Done
    return background


def find_point_sources(image, threshold=4.0, fwhm=1.0, sharp=[0.2,1.0], round=[-1.0,1.0], background_size=16, verbose=True):
    """
    Given a 2-D image, find all of the point sources in it that meet the
    selection criteria provided via the keywords.  These are:
      * threshold: detection threshold in counts about the background
      * fwhm: source full width half max. in pixel
      * sharp: two-element array that defines the lower and upper bounds
               on the source sharpness
      * round: two-element array that defines the lower and upper bounds
               on the source roundness
        
    Background estimation and removal is handled by the estimate_background()
    function in this module that implements a SExtractor-like method.  This
    can be disabled by setting the 'background_size' keyword to 0.  For 
    details see :func:`lsl.imaging.analysis.estimate_background`.
        
    The output of this function is a five-element tuple of 1-D NumPy arrays
    that store information for each source.  The elements are:
      * x: intensity-weighted center - x coordinate
      * y: intensity-weighted center - y coordinate
      * flux: peak count value - input image units
      * sharpness: sharpness statistic
      * roundness: roundness statistic
    
    This function is based on the FIND procedure from the AstroIDL User's
    Library that was written by W. Landsman (Hughes STX).
    
    For additional information about the original IDL routines, see:
    http://idlastro.gsfc.nasa.gov/contents.html#C2
    """
    
    # Maximum size of convolution box in pixels 
    maxbox = 13
    
    # Make sure that the FWHM is not too small and compute the corresponding 
    # value for sigma
    if fwhm < 1.0:
        fwhm = 1.0
    sigma = fwhm / (2*math.sqrt(2*math.log(2)))
    
    # Setup the source detection boxes
    ## The detection radius is the larger of 1.5 sigma or 2 pixels
    radius = max([1.5*sigma, 2.0])			# Radius is the
    ## The detection box size is the smaller of the radius or half the
    ## maximum size of the convolution box
    nHalf = min([int(radius),  (maxbox-1)/2])   	
    ## The size of the convolution box
    nBox = 2*nHalf + 1
    
    # Remove the diffuse background.
    if background_size > 0:
        background = estimate_background(image, background_size)
        cleanImage = image - background
    else:
        cleanImage = image
        
    # Create the raw convolution kernel and a mask that excludes pixels beyond
    # the detection radius
    ## Coordinate setup
    row = numpy.arange(nBox, dtype=numpy.float64) - nHalf
    rx, ry = numpy.meshgrid(row, row)
    ## Build the mask
    mask = (rx**2 + ry**2) <= radius**2
    valid = numpy.where( mask )
    pixels = int(mask.sum())
    ## Build the kernel
    kernel = numpy.exp(-(rx**2 + ry**2)/2.0/sigma**2)
    ## Trim the kernel to take into account the mask and normalize it
    kernel *= mask  
    kernel /= kernel[valid].sum()
    
    # Build a cross-section of the raw kernel for use in computing the 
    # roundness
    c1 = numpy.exp(-row**2/2.0/sigma**2)
    c1 /= c1.sum()
    
    # Convolve the image with the trimmed and normalized kernel
    convImage = convolve(cleanImage, kernel, mode='same')
    
    # Clean out the image edges where the convolution is lacking
    minConvValue = convImage[nHalf:-nHalf,nHalf:-nHalf].min()
    convImage[:nHalf,:] = minConvValue
    convImage[-nHalf:,:] = minConvValue
    convImage[:,:nHalf] = minConvValue
    convImage[:,-nHalf:] = minConvValue
    
    # Exclude the central pixel from the mask and create the arrays that store
    # the relative pixel coordinates for the detection box
    mask[nHalf,nHalf] = 0
    pixels -= 1
    xx, yy = numpy.where(mask)
    xx, yy = xx - nHalf, yy - nHalf
    
    # Find all of the pixels above the threshold and run a simple "de-duplication"
    # to get the maximum inside the detection box
    ix, iy = numpy.where( convImage >= threshold )
    for i in xrange(pixels):
        ox = ix + xx[i]
        oy = iy + yy[i]
        stars, = numpy.where( convImage[ix,iy] >= convImage[ox,oy] )
        ix, iy = ix[stars], iy[stars]
    nGood = len(ix)       
    
    # Setup the output arrays
    ## Position
    x = numpy.zeros(nGood)
    y = numpy.zeros(nGood)
    ## Peak flux
    flux = numpy.zeros(nGood)
    ## Shape statistics
    sharpness = numpy.zeros(nGood)
    roundness = numpy.zeros(nGood)
    
    # Loop over the source positions and compute the various statistics
    nStar = 0
    bad = {'round': 0, 'sharp': 0}
    for i in xrange(nGood): 
        ## Sources are valid until they are not
        validSharpness = True
        validRoundness = True
    
        ## Extract a "postage" stamp of the source and get it's peak
        temp = image[ix[i]-nHalf:ix[i]+nHalf+1,iy[i]-nHalf:iy[i]+nHalf+1]
        cPeak = image[ix[i],iy[i]]			# "d" is actual pixel intensity        
        
        ## Compute the sharpness statistic
        cSharpness = (temp[nHalf,nHalf] - ((mask*temp).sum())/pixels)/cPeak
        if cSharpness < sharp[0] or cSharpness > sharp[1]:
            validSharpness = False
            bad['sharp'] += 1
            
        ## Compute the roundness statistic
        dx = (temp.sum(axis=1)*c1).sum()   
        dy = (temp.sum(axis=0)*c1).sum()
        cRoundness = (dx-dy) / (dx+dy)		# Roundness statistic
        if (dx <= 0 or dy <= 0) or (cRoundness < round[0] or cRoundness > round[1]):
            validRoundness = False
            bad['round'] += 1
            
        ## Flux-weighted center
        xcen = ix[i] + (temp*rx).sum()/temp.sum()
        ycen = iy[i] + (temp*ry).sum()/temp.sum()
        
        ## Save the results if they are valid
        if validSharpness and validRoundness:
            x[nStar] = xcen
            y[nStar] = ycen
            flux[nStar] = cPeak
            sharpness[nStar] = cSharpness
            roundness[nStar] = cRoundness
            nStar += 1
            
        ## Print out information, if requested
        if verbose:
            print("Source #%i" % (i+1,))
            print("  Center:  %.3f, %.3f" % (xcen, ycen))
            print("  Peak: %.3f" % cPeak)
            print("  Sharpness: %.3f%s" % (cSharpness, " (rejected)" if not validSharpness else ""))
            print("  Roundness: %.3f%s" % (cRoundness, " (rejected)" if not validRoundness else ""))
            
    # Print out a summary, if requested
    if verbose:
        print(" ")
        print("Summary")
        print("  Detections: %i" % nGood)
        print("  Valid Detections: %i" % nStar)
        print("  Number Rejected:")
        print("    Sharpness: %i" % bad['sharp'])
        print("    Roundness: %i" % bad['round'])
        
    # Trim the output arrays for the actual number of stars found
    x = x[:nStar]
    y = y[:nStar]
    flux = flux[:nStar]
    sharpness = sharpness[:nStar]
    roundness = roundness[:nStar]
    
    # Done
    return x, y, flux, sharpness, roundness
