"""
Useful math functions for LWA work.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import math
import warnings

import numpy
from scipy.special import sph_harm

from lsl.misc import telemetry
telemetry.track_module()


__version__   = '0.6'
__all__ = ['regrid', 'downsample', 'smooth', 'cmagnitude', 'cphase', 'cpolar', 'crect', 
           'creal', 'cimag', 'to_dB', 'from_dB', 'ndft',  'gaussian1d', 'gaussian2d', 
           'gaussparams', 'sphfit', 'sphval']
__author__    = 'P. S. Ray'
__maintainer__ = 'Jayce Dowell'


def regrid(x, y, newx, allow_extrapolation = False, method = 'spline'):
    """
    Regrid data from x,y onto newx. If allow_extrapolation is True,
    extrapolation is attempted if the method supports it.  Supported
    methods are:
      * linear
      * spline
    Use of this function may require the scipy extension package.
    """
    
    if method == 'linear':
        return _regrid_linear(x, y, newx, allow_extrapolation)
    elif method == 'spline':
        return _regrid_spline(x, y, newx, allow_extrapolation)
    else:
        raise ValueError("interp method '%s' not recognized" % method)


def _regrid_linear(x, y, newx, allow_extrapolation=False):
    """
    Implement regrid() function using linear interpolation.
    """
    
    if allow_extrapolation:
        warnings.warn("allow_extrapolation=True not honored for regrid_linear", RuntimeWarning)
        
    if newx.min() < x.min():
        raise ValueError('x.min(%f) must be smaller than newx.min(%f)' % (x.min(), newx.min()))
    if newx.max() > x.max():
        raise ValueError('x.max(%f) must be larger than newx.max(%f)' % (x.max(), newx.max()))
        
    return numpy.interp(newx, x, y)


def _regrid_spline(x, y, newx, allow_extrapolation=False):
    """
    Implement regrid() function using spline fit and interpolation.
    """
    
    from scipy import interpolate

    if not allow_extrapolation:
        if newx.min() < x.min():
            raise ValueError('x.min(%f) must be smaller than newx.min(%f)' % (x.min(), newx.min()))
        if newx.max() > x.max():
            raise ValueError('x.max(%f) must be larger than newx.max(%f)' % (x.max(), newx.max()))
    spl = interpolate.splrep(x, y)
    newy = interpolate.splev(newx, spl)
    
    return newy


def downsample(vector, factor, rescale=True):
    """
    Downsample (i.e. co-add consecutive numbers) a vector by an integer 
    factor.  Trims the input timeseries to be a multiple of the downsample 
    factor, if needed.  If rescale == True, then divides each sum by factor 
    to produce a mean value, otherwise just adds the values in the vector.
    """

    if (len(vector) % factor):
        warnings.warn("Length of 'vector' is not divisible by 'factor'=%d, clipping!" % factor, RuntimeWarning)
        newlen = (len(vector)//factor)*factor
        warnings.warn("Oldlen %d, newlen %d" % (len(vector), newlen), RuntimeWarning)
        vector = vector[:newlen]
    if rescale:
        newvector = numpy.reshape(vector, (len(vector)//factor, factor))/float(factor)
    else:
        newvector = numpy.reshape(vector, (len(vector)//factor, factor))
        
    return numpy.add.reduce(newvector, 1)


def smooth(x, window_len=10, window='hanning'):
    """
    Smooth the data using a window with requested size.  Stolen from SciPy 
    Cookbook at http://www.scipy.org/Cookbook/SignalSmooth
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    Input:
      * x: the input signal 
      * window_len: the dimension of the smoothing window
      * window: the type of window from 'flat', 'hanning', 'hamming', 
                'bartlett', 'blackman' flat window will produce a moving 
                average smoothing.

    Output:
      * the smoothed signal
      
    Example:
        >>> from numpy import *
        >>> t=linspace(-2,2,0.1)
        >>> x=sin(t)+randn(len(t))*0.1
        >>> y=smooth(x)
    
    .. seealso:: 
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
    
    TODO: the window parameter could be the window itself if an array instead of a string
    """
    
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
        
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
        
    if window_len < 3:
        return x
        
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
        
    s = numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w = numpy.ones(window_len, dtype=s.dtype)
    else:
        w = eval('numpy.'+window+'(window_len)')
        
    y = numpy.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]


def cmagnitude(cmplx):
    """
    Return the polar magnitudes of complex values.
    """
    
    return abs(cmplx)


def cphase(cmplx):
    """
    Return the polar phases of complex values as radians.
    """
    
    return numpy.angle(cmplx)


def cpolar(cmplx):
    """
    Return the polar (magnitude, phase) representation of complex
    values (real, imaginary).  The return value is an array of shape (N,2),
    where N is the length of the cmplx input array.
    """
    
    if isinstance(cmplx, (numpy.ndarray, list, tuple)):
        return numpy.array(list(zip(cmagnitude(cmplx), cphase(cmplx))))
    else:
        return (cmagnitude(cmplx), cphase(cmplx))


def creal(cmplx):
    """
    Return the real rectilinear component from complex values
    expressed in polar form (magnitude, phase).
    """
    
    if isinstance(cmplx, numpy.ndarray):
        return (cmplx[...,0] * numpy.cos(cmplx[...,1]))
    else:
        return (cmplx[0] * math.cos(cmplx[1]))


def cimag(cmplx):
    """
    Return the imaginary rectilinear component from complex values
    expressed in polar form (magnitude, phase).
    """
    
    if isinstance(cmplx, numpy.ndarray):
        return (cmplx[...,0] * numpy.sin(cmplx[...,1]))
    else:
        return (cmplx[0] * math.sin(cmplx[1]))


def crect(cmplx):
    """
    Return the rectilinear (real, imaginary) representation of complex
    values (magnitude, phase).
    """
    
    if isinstance(cmplx, numpy.ndarray):
        ret = numpy.empty((len(cmplx),), numpy.complex64 if cmplx.dtype == numpy.float32 else numpy.complex128)
        ret.real = creal(cmplx)
        ret.imag = cimag(cmplx)
        return ret         
    else:
        return complex(creal(cmplx), cimag(cmplx))


def to_dB(factor):
    """
    Convert from linear units to decibels.
    """
    
    return 10.0 * numpy.log10(factor)


def from_dB(dB):
    """
    Convert from decibels to linear units.
    """
    
    return numpy.power(10.0, (dB/10.0))


def ndft(t, x):
    """
    Given a list of times and a list of data values, compute a non-uniform 
    discrete Fourier transform (NDFT) of the data.  Returns a two element 
    tuple of frequency and the complex NDFT result.
    """
    
    # Setup the output dtype to make sure that we have enough precision
    if x.dtype in (numpy.complex128, numpy.float64):
        dtype = numpy.complex128
    else:
        dtype = numpy.complex64
        
    # Create the output NDFT array and fill it
    out = numpy.zeros(x.shape, dtype=dtype)
    for m in range(out.size):
        mPrime = out.size//2 - m
        s = 0.0j
        for n in range(out.size):
            s += x[n]*numpy.exp(-2j*numpy.pi*t[n]*mPrime / (t.max() - t.min()))
        out[m] = s
        
    # Create the output frequency array and fill it
    f = numpy.zeros_like(t)
    for m in range(out.size):
        mPrime = out.size//2 - m
        f[m] = mPrime / (t.max() - t.min())
        
    # Done
    return f, out


def gaussian1d(height, center, width):
    """
    Return a function that generates a 1-D gaussian with the specified
    height, mean, and standard deviation.  
    
    Example:
    >>> height = 1
    >>> center = 5.0
    >>> width = 2.1
    >>> gauFnc = guassian1d(height, center, width)
    >>> value = gauFnc(numpy.arange(0, 100))
    
    Based on: http://code.google.com/p/agpy/source/browse/trunk/agpy/gaussfitter.py
    """
    
    width = float(width)
    return lambda x: height*numpy.exp(-(center-x)**2/2.0/width**2)


def gaussian2d(height, centerX, centerY, widthMaj, widthMin, angle=0.0):
    """
    Return a function that generates a 2-D gaussian with the specified 
    height, mean (for both X and Y), standard deviation (for both major and 
    minor axes), and rotation angle from the X axis in degrees.

    Based on: http://code.google.com/p/agpy/source/browse/trunk/agpy/gaussfitter.py
    """

    widthMaj = float(widthMaj)
    widthMin = float(widthMin)
    pa = angle*numpy.pi/180.0
    return lambda x,y: height*numpy.exp(-((((x-centerX)*numpy.cos(pa)+(y-centerY)*numpy.sin(pa))/widthMaj)**2 + ((-(x-centerX)*numpy.sin(pa)+(y-centerY)*numpy.cos(pa))/widthMin)**2)/2.0)


def gaussparams(data, x=None, y=None):
    """
    Estimate the parameters (height, center, width) for a gaussian.  The 
    return order is:
      1-D: height, center, width
      2-D: height, center x, center y, width x, width y, position angle
    
    .. note:: 
        The 2-D fits always return a position angle of zero since the
        routine decomposes the process into two 1-D fits.
    """
    
    total = data.sum()
    height = data.max()
    
    if len(data.shape) == 1:
        # 1-D Data
        if x is None:
            x = numpy.arange(data.size)
        center = (x*data).sum() / total
        width = numpy.sqrt(abs(numpy.sum((x-center)**2*data)/total))
        return height, center, width
        
    elif len(data.shape) == 2:
        # 2-D Data
        if x is None or y is None:
            x, y = numpy.indices(data.shape)
            
        # Break the problem down into two 1-D filts
        profileX = data.sum(axis=1)
        profileX *= height/profileX.max()
        profileY = data.sum(axis=0)
        profileY *= height/profileY.max()
        
        hx, centerX, widthX = gaussparams(profileX, x[:,0])
        hy, centerY, widthY = gaussparams(profileY, y[0,:])
        return height, centerX, centerY, widthX, widthY, 0.0
        
    else:
        # N-D data
        raise ValueError("Cannot estimate parameters for %i-D" % (len(data.shape),))        


def sphfit(az, alt, data, lmax=5, degrees=False, real_only=False):
    """
    Decompose a spherical or semi-spherical data set into spherical harmonics.  

    Inputs:
      * az: 2-D numpy array of azimuth coordinates in radians or degrees if the 
            `degrees` keyword is set
      * alt: 2-D numpy array of altitude coordinates in radian or degrees if the 
             `degrees` keyword is set
      * data: 2-D numpy array of the data to be fit.  If the data array is purely
             real, then the `real_only` keyword can be set which speeds up the 
             decomposition
      * lmax: integer setting the maximum order harmonic to fit

    Keywords:
      * degrees: boolean of whether or not the input azimuth and altitude coordinates
         are in degrees or not
      * real_only: boolean of whether or not the input data is purely real or not.  If
        the data are real, only coefficients for modes >=0 are computed.

    Returned is a 1-D complex numpy array with the spherical harmonic coefficients 
    packed packed in order of increasing harmonic order and increasing mode, i.e.,
    (0,0), (1,-1), (1,0), (1,1), (2,-2), etc.  If the `real_only` keyword has been 
    set, the negative coefficients for the negative modes are excluded from the 
    output array.
    
    .. note::
        sphfit was designed to fit the LWA dipole response pattern as a function of
        azimuth and elevation.  Elevation angles are mapped to theta angles by adding
        pi/2 so that an elevation of 90 degrees corresponds to a theta of 180 degrees.
        To fit in terms of spherical coordianates, subtract pi/2 from the theta values
        before running.
    """
    
    if degrees:
        rAz = az*numpy.pi/180.0
        rAlt = alt*numpy.pi/180.0
    else:
        rAz = 1.0*az
        rAlt = 1.0*alt
    rAlt += numpy.pi/2
    sinAlt = numpy.sin(rAlt)
    
    if real_only:
        nTerms = (lmax*(lmax+3)+2)/2
        terms = numpy.zeros(nTerms, dtype=numpy.complex64)
        
        t = 0
        for l in range(lmax+1):
            for m in range(0, l+1):
                Ylm = sph_harm(m, l, rAz, rAlt)
                terms[t] = (data*sinAlt*Ylm.conj()).sum() * (rAz[1,0]-rAz[0,0])*(rAlt[0,1]-rAlt[0,0])
                t += 1
                
    else:
        nTerms = (lmax+1)**2
        terms = numpy.zeros(nTerms, dtype=numpy.complex64)
        
        t = 0
        for l in range(lmax+1):
            for m in range(-l, l+1):
                Ylm = sph_harm(m, l, rAz, rAlt)
                terms[t] = (data*sinAlt*Ylm.conj()).sum()
                t += 1
                
    return terms


def sphval(terms, az, alt, degrees=False, real_only=False):
    """
    Evaluate a set of spherical harmonic coefficents at a specified set of
    azimuth and altitude coordinates.

    Inputs:
      * terms: 1-D complex numpy array, typically from sphfit
      * az: 2-D numpy array of azimuth coordinates in radians or degrees if the 
            `degrees` keyword is set
      * alt: 2-D numpy array of altitude coordinates in radian or degrees if the 
             `degrees` keyword is set

    Keywords:
      * degrees: boolean of whether or not the input azimuth and altitude coordinates
                 are in degrees or not
      * real_only: boolean of whether or not the input data is purely real or not.  If
                  the data are real, only coefficients for modes >=0 are computed.

    Returns a 2-D numpy array of the harmoics evalated and summed at the given 
    coordinates.
    
    .. note::
        sphfit was designed to fit the LWA dipole response pattern as a function of
        azimuth and elevation.  Elevation angles are mapped to theta angles by adding
        pi/2 so that an elevation of 90 degrees corresponds to a theta of 180 degrees.
        To spherical harmonics in terms of spherical coordianates, subtract pi/2 from 
        the theta values before running.
    """
    
    if degrees:
        rAz = az*numpy.pi/180.0
        rAlt = alt*numpy.pi/180.0
    else:
        rAz = 1.0*az
        rAlt = 1.0*alt
    rAlt += numpy.pi/2
    
    nTerms = terms.size
    if real_only:
        lmax = int((numpy.sqrt(1+8*nTerms)-3)/2)

        t = 0
        out = numpy.zeros(az.shape, dtype=numpy.float32)
        for l in range(lmax+1):
            Ylm = sph_harm(0, l, rAz, rAlt)
            out += numpy.real(terms[t]*Ylm)
            t += 1
            for m in range(1, l+1):
                Ylm = sph_harm(m, l, rAz, rAlt)
                out += numpy.real(terms[t]*Ylm)
                out += numpy.real(terms[t]*Ylm.conj()/(-1)**m)
                t += 1
                
    else:
        lmax = int(numpy.sqrt(nTerms)-1)
        
        t = 0
        out = numpy.zeros(az.shape, dtype=numpy.complex64)
        for l in range(lmax+1):
            for m in range(-l, l+1):
                Ylm = sph_harm(m, l, rAz, rAlt)
                out += terms[t]*Ylm
                t += 1
                
    return out
