"""
Module for removing multi-path scattering effects in pulsar profiles.  This 
is based on the CLEAN-like deconvolution method presented in Bhat, N., 
Cordes, J., & Chatterjee, S.  2003, ApJ, 584, 782.

http://iopscience.iop.org/0004-637X/584/2/782/fulltext/56392.text.html

.. note::
    All functions assume that the time and scattering time values are in
    seconds.
"""

import numpy as np

from lsl.statistics import robust

from lsl.misc import telemetry
telemetry.track_module()

from typing import Callable, Tuple


__version__ = "0.1"
__all__ = ['thin', 'thick', 'uniform', 'unscatter']


def thin(t: np.ndarray, tau: float) -> np.ndarray:
    """
    Pulsar broadening function for multi-path scattering through a
    thin screen.
    """

    g  = 1.0/tau * np.exp(-t/tau)
    g  = np.where(t >= 0, g, 0)
    g /= g.sum()

    return g


def thick(t: np.ndarray, tau: float) -> np.ndarray:
    """
    Pulse broadening function for multi-path scattering through a
    thick screen.

    From:
    Williamson, I. P. 1972, MNRAS, 157, 55; first equation on page 62 
    """

    tPrime = t + 1e-15
    
    g  = np.sqrt(np.pi*tau/4/tPrime**3)
    g *= np.exp(-np.pi**2*tau/16/tPrime)
    g  = np.where(t > 0, g, 0)
    g /= g.sum()
    
    return g


def uniform(t: np.ndarray, tau: float) -> np.ndarray:
    """
    Pulsr broadening function for multi-path scattering through a 
    uniform screen.

    From:
    Williamson, I. P. 1972, MNRAS, 157, 55; last equation on page 65
    """

    tPrime = t + 1e-15
    
    g  = np.sqrt(np.pi**5*tau**3/8/tPrime**5)
    g *= np.exp(-np.pi**2*tau/4/tPrime)
    g  = np.where(t > 0, g, 0)
    g /= g.sum()

    return g


def _positivity(t: np.ndarray, raw: np.ndarray, resids: np.ndarray, cc: np.ndarray) -> float:
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
    f  = (resids**2 * np.where(temp>=0, 1, 0)).mean()
    f *= m/sigma**2

    return f

def _skewness(t: np.ndarray, raw: np.ndarray, resids: np.ndarray, cc: np.ndarray) -> float:
    """
    Compute the skewness for a collection of clean components following 
    Equation 12 of Bhat, N., Cordes, J., & Chatterjee, S.  2003, ApJ, 
    584, 782.
    """

    tM = (t*cc).sum() / cc.sum()
    t2 = ((t - tM)**2 * cc).sum()
    t3 = ((t - tM)**3 * cc).sum()
    
    return t3 / (t2 + 1e-15)**(3./2.)


def _figure_of_merit(t: np.ndarray, raw: np.ndarray, resids: np.ndarray, cc: np.ndarray) -> float:
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
    n = len( np.where( np.abs(resids - m) <= 3*s )[0] ) / raw.size

    # Compute the positivity and skewness
    f = _positivity(t, raw, resids, cc)
    g = _skewness(t, raw, resids, cc)
    
    # Get the standard deviation of all of the residuals relative to the estimated
    # off-pulse value
    r  = resids.std()
    r /= s
    
    # The figure-of-metric
    return (f + g)/2.0 + r - n


def unscatter(t: np.ndarray, raw: np.ndarray, tScatMin: float, tScatMax: float, tScatStep: float, gain: float=0.05, max_iter: int=10000, screen: Callable[[np.ndarray,float],np.ndarray]=thin, verbose: bool=True) -> Tuple[float,float,np.ndarray]:
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
     * max_iter: maximum number of iterations to use (default is 10000)
     * screen: pulse broadening function (default is thin)
        
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
    best = np.inf
    bestTau = None
    for tScat in np.arange(tScatMin, tScatMax, tScatStep):
        ## Setup the temporary variables
        working = raw*1.0
        cc = raw*0.0

        ## Iterate...
        i = 0
        while i < max_iter:
            ### Find the peak and make sure it is really a peak
            peak = np.where( working == working.max() )[0][0]
            if working.max() < 3./2.*sigma:
                break

            ### Generate the clean component
            tPeak = t[peak]
            tRel = t - tPeak
            toRemove  = screen(tRel, tScat)
            toRemove /= toRemove.sum()
            toRemove *= gain*working.max()
            
            ### Remove and continue
            cc[peak] += toRemove.sum()
            working -= toRemove
            i += 1
            
        ## Evaluate and save
        iList[tScat] = i
        meritList[tScat] = _figure_of_merit(t, raw, working, cc)
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
        print("  Iterations Used: %i of %i" % (i, max_iter))
        print("  Best-fit Scattering time: %.3f ms" % (tScat*1000.0,))
        print("  Figure-of-merit:  %.5f" % merit)
        
    # Restore the profile using a Gaussian with a sigma value of 5 time steps
    sigmaRestore = 5*(t[1]-t[0])
    restoreFunction = np.exp(-t**2 / (2*sigmaRestore**2))
    restoreFunction /= restoreFunction.sum()
    out = np.convolve(cc, restoreFunction, 'same')
    
    # Add back in the residuals
    out += resids
    
    # Done!
    return tScat, merit, out
    
