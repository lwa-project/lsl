# -*- coding: utf-8 -*-

"""Useful math functions for LWA work"""

import logging
import math

import numpy

__version__   = "0.1"
__revision__  = "$Revision: 92 $"
__all__ = ['regrid', 'downsample', 'smooth', 'cmagnitude', 'cphase', 'cpolar', 'crect', 'creal', 'cimag', 'to_dB', 'from_dB', '__version__', '__revision__', '__all__']
__author__    = "P.S.Ray"
__maintainer__ = "Jayce Dowell"

_MATHUTIL_LOG = logging.getLogger('mathutil')


def regrid(x, y, newx, allow_extrapolation = False, method = 'spline'):
	"""Regrid data from x,y onto newx. If allow_extrapolation is True,
	extrapolation is attempted if the method supports it.  Supported
	methods are:
		linear
		spline
	Use of this function may require the scipy extension package."""
	
	if method == 'linear':
		return _regrid_linear(x, y, newx, allow_extrapolation)
	elif method == 'spline':
		return _regrid_spline(x, y, newx, allow_extrapolation)
	else:
		raise ValueError("interp method '%s' not recognized" % method)


def _regrid_linear(x, y, newx, allow_extrapolation=False):
	"""Implement regrid() function using linear interpolation."""
	
	if allow_extrapolation:
		_MATHUTIL_LOG.warning("allow_extrapolation=True not honored for regrid_linear")
	
	if newx.min() < x.min():
		raise ValueError('x.min(%f) must be smaller than newx.min(%f)' % (x.min(), newx.min()))
	if newx.max() > x.max():
		raise ValueError('x.max(%f) must be larger than newx.max(%f)' % (x.max(), newx.max()))

	return numpy.interp(newx, x, y)


def _regrid_spline(x, y, newx, allow_extrapolation=False):
	"""Implement regrid() function using spline fit and interpolation."""
	
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
	"""downsample(vector, factor):
	Downsample (i.e. co-add consecutive numbers) a vector
	by an integer factor.  Trims the input timeseries to be 
	a multiple of the downsample factor, if needed.
	If rescale == True, then divides each sum by factor to produce a mean value,
	otherwise just adds the values in the vector."""

	if (len(vector) % factor):
		_MATHUTIL_LOG.warning("Length of 'vector' is not divisible by 'factor'=%d, clipping!", factor)
		newlen = (len(vector)/factor)*factor
		_MATHUTIL_LOG.debug("Oldlen %d, newlen %d", len(vector), newlen)
		vector = vector[:newlen]
	if rescale:
		newvector = numpy.reshape(vector, (len(vector)/factor, factor))/float(factor)
	else:
		newvector = numpy.reshape(vector, (len(vector)/factor, factor))

	return numpy.add.reduce(newvector, 1)
    
    
def smooth(x,window_len=10,window='hanning'):
	"""Smooth the data using a window with requested size.
	
	Stolen from SciPy Cookbook
	http://www.scipy.org/Cookbook/SignalSmooth
	
	This method is based on the convolution of a scaled window with the signal.
	The signal is prepared by introducing reflected copies of the signal 
	(with the window size) in both ends so that transient parts are minimized
	in the begining and end part of the output signal.
	
	input:
		x: the input signal 
		window_len: the dimension of the smoothing window
		window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
			flat window will produce a moving average smoothing.

	output:
		the smoothed signal
		
	example:

	t=linspace(-2,2,0.1)
	x=sin(t)+randn(len(t))*0.1
	y=smooth(x)
	
	see also: 
	
	numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
	scipy.signal.lfilter
	
	TODO: the window parameter could be the window itself if an array instead of a string"""

	if x.ndim != 1:
		raise ValueError("smooth only accepts 1 dimension arrays.")

	if x.size < window_len:
		raise ValueError("Input vector needs to be bigger than window size.")


	if window_len<3:
		return x


	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


	s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
	#print(len(s))
	if window == 'flat': #moving average
		w=numpy.ones(window_len,'d')
	else:
		w=eval('numpy.'+window+'(window_len)')

	y=numpy.convolve(w/w.sum(),s,mode='same')
	return y[window_len-1:-window_len+1]


def cmagnitude(cmplx):
	"""Return the polar magnitudes of complex values."""
	
	return abs(cmplx)
    

def cphase(cmplx):
	"""Return the polar phases of complex values as radians."""
	
	return numpy.angle(cmplx)
    
    
def cpolar(cmplx):
	"""Return the polar (magnitude, phase) representation of complex
	values (real, imaginary).  The return value is an array of shape (N,2),
	where N is the length of the cmplx input array."""
	
	if isinstance(cmplx, (numpy.ndarray, list, tuple)):
		return numpy.array(list(zip(cmagnitude(cmplx), cphase(cmplx))))
	else:
		return (cmagnitude(cmplx), cphase(cmplx))
    

def creal(cmplx):
	"""Return the real rectilinear component from complex values
	expressed in polar form (magnitude, phase)."""
	
	if isinstance(cmplx, numpy.ndarray):
		return (cmplx[...,0] * numpy.cos(cmplx[...,1]))
	else:
		return (cmplx[0] * math.cos(cmplx[1]))   

  
def cimag(cmplx):
	"""Return the imaginary rectilinear component from complex values
	expressed in polar form (magnitude, phase)."""
	
	if isinstance(cmplx, numpy.ndarray):
		return (cmplx[...,0] * numpy.sin(cmplx[...,1]))
	else:
		return (cmplx[0] * math.sin(cmplx[1]))    
    

def crect(cmplx):
	"""Return the rectilinear (real, imaginary) representation of complex
	values (magnitude, phase)."""
	
	if isinstance(cmplx, numpy.ndarray):
		ret = numpy.empty((len(cmplx),), numpy.complex_)
		ret.real = creal(cmplx)
		ret.imag = cimag(cmplx)
		return ret         
	else:
		return complex(creal(cmplx), cimag(cmplx))
        
    
def to_dB(factor):
	"""Convert from linear units to decibels."""
	
	return 10.0 * numpy.log10(factor)
    
    
def from_dB(dB):
	"""Convert from decibels to linear units."""
	
	return numpy.power(10.0, (dB/10.0))
    
    
def robustmean(arr):
	"""Take the robust mean of an array, normally a small section of a 
	spectrum, over which the mean can be assumed to be constant.  Makes two 
	passes discarding outliers >3 sigma ABOVE (not below) the mean."""
	
	# First pass discarding points >3 sigma above mean
	mean = arr.mean()
	sd = arr.std()
	thresh = mean + 3.0*sd
	idx = numpy.where(arr < thresh)
	newarr = arr[idx]
	
	if len(newarr) == 0:
		# Warning, all points discarded.  Just return array mean
		_MATHUTIL_LOG.warning("All points discarded!, %f, %f", mean, sd)
		finalmean = mean
	else:
		# Second pass discarding points >3 sigma above mean
		newmean = newarr.mean()
		newsd = newarr.std()
		newthresh = newmean+3.0*newsd
		newidx = numpy.where(newarr < newthresh)
		finalarr = newarr[newidx]
		if len(finalarr) == 0:
			finalmean = newmean
		else:
			# Final mean of good points
			finalmean = finalarr.mean()

	return finalmean
