# -*- coding: utf-8 -*-

"""Small collection of robust statistical estimators based on functions from
the AstroIDL User's Library.  Function included are:
  * robustMean - robust estimator of the mean of a data set and
  * robustSigma - robust estimator of the standard deviation of a data set.
"""

import math
import numpy

__version__ = '0.1'
__revision__ = '$ Revision: 2 $'
__all__ = ['robustMean', 'robustSigma', '__version__', '__revision__', '__all__']

__epsilon = 1.0e-20

def robustMean(inputData, Cut=3.0):
	"""Robust estimator of the mean of a data set.  Based on the 
	resistant_mean function from the AstroIDL User's Library.

	.. seealso::
		:func:`lsl.misc.mathutil.robustmean`
	"""

	data = inputData.ravel()

	data0 = numpy.median(data)
	maxAbsDev = numpy.median(numpy.abs(data-data0)) / 0.6745
	if maxAbsDev < __epsilon:
		maxAbsDev = (numpy.abs(data-data0)).mean() / 0.8000

	cutOff = Cut*maxAbsDev
	good = numpy.where( numpy.abs(data-data0) <= cutOff )
	good = good[0]
	dataMean = data[good].mean()
	dataSigma = math.sqrt( ((data[good]-dataMean)**2.0).sum() / len(good) )

	if Cut > 1.0:
		sigmaCut = Cut
	else:
		sigmaCut = 1.0
	if sigmaCut <= 4.5:
		dataSigma = dataSigma / (-0.15405 + 0.90723*sigmaCut - 0.23584*sigmaCut**2.0 + 0.020142*sigmaCut**3.0)

	cutOff = Cut*dataSigma
	good = numpy.where(  numpy.abs(data-data0) <= cutOff )
	good = good[0]
	dataMean = data[good].mean()
	if len(good) > 3:
		dataSigma = math.sqrt( ((data[good]-dataMean)**2.0).sum() / len(good) )

	if Cut > 1.0:
		sigmaCut = Cut
	else:
		sigmaCut = 1.0
	if sigmaCut <= 4.5:
		dataSigma = dataSigma / (-0.15405 + 0.90723*sigmaCut - 0.23584*sigmaCut**2.0 + 0.020142*sigmaCut**3.0)

	dataSigma = dataSigma / math.sqrt(len(good)-1)

	return dataMean
	

def robustSigma(inputData):
	"""Robust estimator of the standard deviation of a data set.  Based on 
	robust_sigma function from the AstroIDL User's Library."""

	data = inputData.ravel()

	data0 = numpy.median(data)
	maxAbsDev = numpy.median(numpy.abs(data-data0)) / 0.6745
	if maxAbsDev < __epsilon:
		maxAbsDev = (numpy.abs(data-data0)).mean() / 0.8000
	if maxAbsDev < __epsilon:
		sigma = 0.0
		return sigma

	u = (data-data0) / 6.0 / maxAbsDev
	u2 = u**2.0
	good = numpy.where( u2 <= 1.0 )
	good = good[0]
	if len(good) < 3:
		print "WARNING:  Distribution is too strange to compute standard deviation"
		sigma = -1.0
		return sigma

	numerator = ((data[good]-data0)**2.0 * (1.0-u2[good])**2.0).sum()
	nElements = (data.ravel()).shape[0]
	denominator = ((1.0-u2[good])*(1.0-5.0*u2[good])).sum()
	sigma = nElements*numerator / (denominator*(denominator-1.0))
	if sigma > 0:
		sigma = math.sqrt(sigma)
	else:
		sigma = 0.0

	return sigma
