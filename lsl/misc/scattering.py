# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import print_function

"""
Module for removing multi-path scattering effects in pulsar profiles.  This 
is based on the CLEAN-like deconvolution method presented in Bhat, N., 
Cordes, J., & Chatterjee, S.  2003, ApJ, 584, 782.

http://iopscience.iop.org/0004-637X/584/2/782/fulltext/56392.text.html

.. note::
    All functions assume that the time and scattering time values are in
    seconds.
"""

import numpy
from lsl.statistics import robust


__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ['thin', 'thick', 'uniform', 'unscatter', '__version__', '__revision__', '__all__']


def thin(t, tauScatter):
    """
    Pulsar broadening function for multi-path scattering through a
    thin screen.
    """

    g  = 1.0/tauScatter * numpy.exp(-t/tauScatter)
    g  = numpy.where(t >= 0, g, 0)
    g /= g.sum()

    return g


def thick(t, tauScatter):
    """
    Pulse broadening function for multi-path scattering through a
    thick screen.

    From:
    Williamson, I. P. 1972, MNRAS, 157, 55; first equation on page 62 
    """

    tPrime = t + 1e-15
    
    g  = numpy.sqrt(numpy.pi*tauScatter/4/tPrime**3)
    g *= numpy.exp(-numpy.pi**2*tauScatter/16/tPrime)
    g  = numpy.where(t > 0, g, 0)
    g /= g.sum()

    return g


def uniform(t, tauScatter):
    """
    Pulsr broadening function for multi-path scattering through a 
    uniform screen.

    From:
    Williamson, I. P. 1972, MNRAS, 157, 55; last equation on page 65
    """

    tPrime = t + 1e-15
    
    g  = numpy.sqrt(numpy.pi**5*tauScatter**3/8/tPrime**5)
    g *= numpy.exp(-numpy.pi**2*tauScatter/4/tPrime)
    g  = numpy.where(t > 0, g, 0)
    g /= g.sum()

    return g


def _positivity(t, raw, resids, cc):
    """
    Compute the positivity of the residual profile following Equation 10 of
    Bhat, N., Cordes, J., & Chatterjee, S.  2003, ApJ, 
    584, 782.
    """
    
    # Weight factor
    m = 1.0

    # Tuning factor to find when the CLEANed profile goes negative
    x = 3./2.

    # Get the off-pulsar standard deviation
    sigma = robust.std(raw)
    
    temp = -resids - x*sigma
    f  = (resids**2 * numpy.where(temp>=0, 1, 0)).mean()
    f *= m/sigma**2

    return f

def _skewness(t, raw, resids, cc):
    """
    Compute the skewness for a collection of clean components following 
    Equation 12 of Bhat, N., Cordes, J., & Chatterjee, S.  2003, ApJ, 
    584, 782.
    """

    tM = (t*cc).sum() / cc.sum()
    t2 = ((t - tM)**2 * cc).sum()
    t3 = ((t - tM)**3 * cc).sum()
    
    return t3 / t2**(3./2.)


def _figureOfMerit(t, raw, resids, cc):
    """
    Figure of merit for deconvolution that combines the positivity of the 
    residuals, the skewness of the clean components, the RMS of the 
    residuals, and the number of noise-like points
    
    From:
    Bhat, N., Cordes, J., & Chatterjee, S.  2003, ApJ, 584, 782.
    """

    # Use robust methods to estimate the off-pulse mean and standard deviation
    m = robust.mean(raw)
    s = robust.std(raw)
    
    # Find the fraction of noise-like points in the residuals
    n = len( numpy.where( numpy.abs(resids - m) <= 3*s )[0] )
    n = float(n) / raw.size

    # Compute the positivity and skewness
    f = _positivity(t, raw, resids, cc)
    g = _skewness(t, raw, resids, cc)
    
    # Get the standard deviation of all of the residuals relative to the estimated
    # off-pulse value
    r  = resids.std()
    r /= s
    
    # The figure-of-metric
    return (f + g)/2.0 + r - n


def unscatter(t, raw, tScatMin, tScatMax, tScatStep, gain=0.05, iterMax=10000, broadeningFunction=thin, verbose=True):
    """
    Multi-path scattering deconvolution method based on the method 
    presented in Bhat, N., Cordes, J., & Chatterjee, S.  2003, ApJ, 
    584, 782.
    
    Inputs:
    1) t: List of times in seconds the pulse profile corresponds to
    2) raw: the raw (scattered) pulse profile over time
    3) tScatMin: minimum scattering time to search
    4) tScatMax: maximum scattering time to search
    5) tScatStep: time step for tScat search
        
    Options:
    * gain: CLEAN loop gain (default is 0.05)
    * iterMax: maximum number of iterations to use (default is 10000)
    * broadeningFunction: pulse broadening function (default is thin)
        
    Outputs (as a tuple):
    1) tScat: best-fit scattering time in seconds
    2) merit: figure of merit for the deconvolution/scattering time fit
    3) cleand: CLEANed pulsar profile as a function of time
    """
    
    iList = {}
    meritList = {}
    ccList = {}
    residList = {}

    # Off-pulsar standard deviation
    sigma = robust.std(raw)

    # Loop over tScat values
    best = 1e9
    bestTau = None
    for tScat in numpy.arange(tScatMin, tScatMax, tScatStep):
        ## Setup the temporary variables
        working = raw*1.0
        cc = raw*0.0

        ## Iterate...
        i = 0
        while i < iterMax:
            ### Find the peak and make sure it is really a peak
            peak = numpy.where( working == working.max() )[0][0]
            if working.max() < 3./2.*sigma:
                break

            ### Generate the clean component
            tPeak = t[peak]
            tRel = t - tPeak
            toRemove  = broadeningFunction(tRel, tScat)
            toRemove /= toRemove.sum()
            toRemove *= gain*working.max()

            ### Remove and continue
            cc[peak] += toRemove.sum()
            working -= toRemove
            i += 1
            
        ## Evaluate and save
        iList[tScat] = i
        meritList[tScat] = _figureOfMerit(t, raw, working, cc)
        ccList[tScat] = cc
        residList[tScat] = working
        
        ## Compare
        if meritList[tScat] < best:
            best = meritList[tScat]
            bestTau = tScat
            
    # Make sure we have something to report
    if tScat is None:
        raise RuntimeError("No good solution found, consider changing the search range.")
        
    # Select out what we need, i.e., the best
    tScat = bestTau
    i = iList[bestTau]
    merit = meritList[bestTau]
    cc = ccList[bestTau]
    resids = residList[bestTau]
    
    # Report on the findings
    if verbose:
        print("Multi-path Scattering Results:")
        print("  Iterations Used: %i of %i" % (i, iterMax))
        print("  Best-fit Scattering time: %.3f ms" % (tScat*1000.0,))
        print("  Figure-of-merit:  %.5f" % merit)
        
    # Restore the profile using a Gaussian with a sigma value of 5 time steps
    sigmaRestore = 5*(t[1]-t[0])
    restoreFunction = numpy.exp(-t**2 / (2*sigmaRestore**2))
    restoreFunction /= restoreFunction.sum()
    out = numpy.convolve(cc, restoreFunction, 'same')
    
    # Add back in the residuals
    out += resids
    
    # Done!
    return tScat, merit, out
    