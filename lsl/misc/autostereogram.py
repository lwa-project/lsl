# -*- coding: utf-8 -*-

# Python3 compatiability
import sys
if sys.version_info > (3,):
    xrange = range
    
"""
Module for taking a 2-D numpy array and turning it into an autostereogram
(http://en.wikipedia.org/wiki/Autostereogram).  This module also takes auto-
stereograms and returns the depth maps.

.. versionadded:: 0.5.3
"""

import numpy

__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ['getASG', 'getDepthMap', '__version__', '__revision__', '__all__']


def getASG(data, patternWidth=100, patternPad=1, maxDepth=30):
    """
    Given a 2-D numpy array of data (image, elevation, etc.), convert the 
    data into a random dot autostereogram.  The options control the size of
    the pattern in the final image, the amount of padding around the image, 
    and how the input data is converted to a depth map.  The output is a 
    three-element tuple of the depth map, input random dot image, and auto-
    stereogram.  
    
    .. note::
        The output autostereogram is a 3-D numpy array with width, depth, and
        color channel.
    """
    
    # Autostereogram final size
    x = data.shape[0] + patternWidth*(2*patternPad)
    y = data.shape[1] + patternWidth*(2*patternPad)
    
    # Pattern to repeat (random colored dots)
    pattern = 255*numpy.random.rand(x, patternWidth, 3)
    pattern = pattern.astype(numpy.int16)
    
    # Number of tiles of the pattern needed to cover the autostereogram
    nTiles = 1.0*y/patternWidth
    
    # Base image and fill
    orig = numpy.zeros((x,y,3), dtype=numpy.int16)
    for i in xrange(int(nTiles)):
        start = i*pattern.shape[1]
        stop = start + pattern.shape[1]
        orig[:,start:stop,:] = pattern
    if nTiles > int(nTiles):
        start = (i+1)*pattern.shape[1]
        stop = orig.shape[1]
        orig[:,start:stop,:] = pattern[:,:(stop-start),:]
    
    # Compute the depth map and convert it to an integer
    dmap = data - data.min()
    dmap = numpy.round( float(maxDepth) * dmap/dmap.max() )
    dmap = dmap.astype(numpy.int16)

    # Shift it base image to make the autostereogram
    final = 1*orig
    for i in xrange(dmap.shape[0]):
        l = i + (patternWidth*patternPad)
        
        for j in xrange(dmap.shape[1]):
            m = j + (patternWidth*patternPad)
            d = patternWidth + dmap[i,j]
            
            final[l,m,:] = 1*final[l,m-d,:]
            
    # Return a three-element tuple of the map, base image, and autostereogram
    return dmap, orig, final


def getDepthMap(final):
    """
    Given an autostereogram image (width x heigth x color channel), determine the 
    depth map and return it.  This function uses the padding area around the input 
    image to determine the pattern width and then uses that width to look for 
    relative shifts.
    """
    
    # Find the pattern width by looking at the correlation of different 
    # columns of the first image channel
    l = numpy.arange(9, final.shape[1]/2)
    p = []
    for r in l:
        part1 = final[:,0:r,0].ravel().astype(numpy.float32)
        part2 = final[:,r:(2*r),0].ravel().astype(numpy.float32)
        p.append( numpy.correlate(part1, part2)[0] / r )
    p = numpy.array(p)
    p = numpy.abs(p)**2

    ## The pattern width occurs at the maximum
    patternWidth = l[numpy.where( p == p.max() )][0]
    
    # Find the depth map by 
    shift = numpy.zeros((final.shape[0]-2*patternWidth, final.shape[1]-2*patternWidth))
    for i in xrange(patternWidth, final.shape[0]-patternWidth):
        for j in xrange(patternWidth, final.shape[1]-patternWidth):
            d = numpy.zeros(patternWidth)
            for k in xrange(patternWidth):
                d[k] = ((final[i,j,:] - final[i,j-(k+patternWidth),:])**2).sum()
            d = numpy.sqrt(d)
            
            try:
                s = numpy.where( d == 0 )[0][0]
            except IndexError:
                s = 0
                
            shift[i-patternWidth,j-patternWidth] = s
            
    # Return the depth map
    return shift
