# -*- coding: utf-8 -*-

"""
Simple self-calibration module for correlated TBW and TBN data.  The 
supported self-calibration methods are:
* phase-only
* amplitude and phase
* delay-only
* amplitude and delay
* delay/phase offset
* amplitude and delay/phase offset

..versionchanged:: 0.6.3
    Reworked the module to make it more versatile
    
..versionadded:: 0.5.5
"""

import numpy

from lsl.statistics import robust

__version__ = '0.2'
__revision__ = '$Rev$'
__all__ = ['phaseOnly', 'delayOnly', 'delayAndPhase', '__version__', '__revision__', '__all__']


def _scaleData(dataDict, amps, delays, phaseOffsets):
    """
    Apply a set of antenna-based real gain value, phase delays in ns, phase 
    offsets in radians to a data dictionary.  Returned the new scaled and 
    delayed dictionary.
    """

    import copy

    # Build the data dictionary to hold the scaled and delayed data
    sclUVData = {'freq': (dataDict['freq']).copy(), 'uvw': {}, 'vis': {}, 'wgt': {}, 'msk': {}, 'bls': {}, 'jd': {}}
    if dataDict['isMasked']:
        sclUVData['isMasked'] = True
    else:
        sclUVData['isMasked'] = False
    fq = dataDict['freq'] / 1e9
    
    cGains = []
    for i in xrange(len(amps)):
        cGains.append( amps[i]*numpy.exp(2j*numpy.pi*fq*delays[i] + 1j*phaseOffsets[i]) )

    # Apply the scales and delays for all polarization pairs found in the original data
    for pol in dataDict['vis'].keys():
        sclUVData['bls'][pol] = []
        sclUVData['uvw'][pol] = []
        sclUVData['vis'][pol] = []
        sclUVData['wgt'][pol] = copy.copy(dataDict['wgt'][pol])
        sclUVData['msk'][pol] = copy.copy(dataDict['msk'][pol])
        sclUVData['jd'][pol] = copy.copy(dataDict['jd'][pol])

        for (i,j),uvw,vis in zip(dataDict['bls'][pol], dataDict['uvw'][pol], dataDict['vis'][pol]):
            sclUVData['bls'][pol].append( (i,j) )
            sclUVData['uvw'][pol].append( uvw )
            sclUVData['vis'][pol].append( vis*cGains[j].conj()*cGains[i] )

    return sclUVData


def _buildAmplitudeA(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix A for amplitude correction.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    A = numpy.zeros((nBLs, nStands))
    for i,(l,m) in enumerate(dataDict['bls'][pol]):
        A[i,l] = 1.0
        A[i,m] = 1.0
        
    A = numpy.matrix(A)
    return A


def _buildAmplitudeC(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix C for the amplitude correction.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    obsVis = []
    for vis in dataDict['vis'][pol]:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in simDict['vis'][pol]:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.log(numpy.abs(simVis / obsVis))
    
    Cp = numpy.zeros(nBLs)
    for i in xrange(C.shape[0]):
        Cp[i] = robust.mean(C[i,:])
    
    return Cp


def _buildPhaseOnlyA(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix A for phase correction with a frequency independent
    phase.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    A = numpy.zeros((nBLs*fq.size, nStands-1))
    for i in xrange(fq.size):
        for j,(l,m) in enumerate(dataDict['bls'][pol]):
            if l < refAnt:
                A[j+i*nBLs,l]   =  1.0
            elif l > refAnt:
                A[j+i*nBLs,l-1] =  1.0
            else:
                pass
                
            if m < refAnt:
                A[j+i*nBLs,m]   = -1.0
            elif m > refAnt:
                A[j+i*nBLs,m-1] = -1.0
            else:
                pass
        
    A = numpy.matrix(A)
    return A


def _buildPhaseOnlyC(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix C for phase correction with a frequency independent
    phase.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    obsVis = []
    for vis in dataDict['vis'][pol]:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in simDict['vis'][pol]:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.angle(simVis / obsVis)
    
    Cp = numpy.zeros(nBLs*fq.size)
    for i in xrange(fq.size):
        for j in xrange(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _buildDelayOnlyA(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix A for phase correction with a delay.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    A = numpy.zeros((nBLs*fq.size, nStands-1))
    for i in xrange(fq.size):
        for j,(l,m) in enumerate(dataDict['bls'][pol]):
            if l < refAnt:
                A[j+i*nBLs,l]   =  2*numpy.pi*fq[i]
            elif l > refAnt:
                A[j+i*nBLs,l-1] =  2*numpy.pi*fq[i]
            else:
                pass
                
            if m < refAnt:
                A[j+i*nBLs,m]   = -2*numpy.pi*fq[i]
            elif m > refAnt:
                A[j+i*nBLs,m-1] = -2*numpy.pi*fq[i]
            else:
                pass
        
    A = numpy.matrix(A)
    return A


def _buildDelayOnlyC(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix C for phase correction with a delay.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    obsVis = []
    for vis in dataDict['vis'][pol]:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in simDict['vis'][pol]:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.angle(simVis / obsVis)
    
    Cp = numpy.zeros(nBLs*fq.size)
    for i in xrange(fq.size):
        for j in xrange(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _buildDelayAndPhaseA(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix A for phase correction with a delay and a phase offset.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    A = numpy.zeros((nBLs*fq.size, 2*(nStands-1)))
    for i in xrange(fq.size):
        for j,(l,m) in enumerate(dataDict['bls'][pol]):
            if l < refAnt:
                A[j+i*nBLs,l]   =  2*numpy.pi*fq[i]
                A[j+i*nBLs,l+(nStands-1)] = 1.0
            elif l > refAnt:
                A[j+i*nBLs,l-1] =  2*numpy.pi*fq[i]
                A[j+i*nBLs,l-1+(nStands-1)] = 1.0
                
            else:
                pass
                
            if m < refAnt:
                A[j+i*nBLs,m]   = -2*numpy.pi*fq[i]
                A[j+i*nBLs,m+(nStands-1)] = -1.0
            elif m > refAnt:
                A[j+i*nBLs,m-1] = -2*numpy.pi*fq[i]
                A[j+i*nBLs,m-1+(nStands-1)] = -1.0
            else:
                pass
        
    A = numpy.matrix(A)
    return A


def _buildDelayAndPhaseC(aa, dataDict, simDict, chan, pol, refAnt=0):
    """
    Build the matrix C for phase correction with a delay and a phase offset.
    """
    
    # Get the baseline and stand counts
    nBLs = len(dataDict['bls'][pol])
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    obsVis = []
    for vis in dataDict['vis'][pol]:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in simDict['vis'][pol]:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.angle(simVis / obsVis)
    
    Cp = numpy.zeros(nBLs*fq.size)
    for i in xrange(fq.size):
        for j in xrange(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _selfCal(aa, dataDict, simDict, chan, pol, refAnt=0, nIter=30, amplitude=False, phaseOnly=False, delayOnly=False, delayAndPhase=False, amplitudeCutoff=1.001, phaseCutoff=0.01, delayCutoff=0.2, verbose=True):
    """
    Function used to perform a variety of self-calibration strategies on 
    data stored in a readUVData dictionary and a model sky stored in a 
    lsl.sim.vis.buildSimSky dictionary for a given polarization and 
    channel(s).
    
    The supported self-cal. schemes are:
    * phase-only
    * amplitude and phase
    * delay-only
    * amplitude and delay
    * delay/phase offset
    * amplitude and delay/phase offset
    
    The function exits when either the maximum number of iterations is 
    reached (nIter) or the maximum absolute value of the relevant quantity
    drops below the quantity's cutoff value.  The cutoff keywords and 
    their default values are:
    * amplitudeCutoff - 1.001, 
    * phaseCutoff - 0.01 radians, 
    * delayCutoff - 0.2 ns (0.01 radians over 10 MHz),
    """

    # Make sure we have the right polarization
    if pol not in dataDict['bls'].keys() and pol.lower() not in dataDict['bls'].keys():
        raise RuntimeError("Data dictionary does not have data for polarization '%s'" % pol)
    if pol not in simDict['bls'].keys() and pol.lower() not in simDict['bls'].keys():
        raise RuntimeError("Simulation dictionary does not have data for polarization '%s'" % pol)

    # Make sure that `chan' is an array by trying to find its length
    try:
        junk = len(chan)
    except TypeError:
        chan = [chan]
        
    N = len(aa.ants)
    found = False
    tauNames = []
    phsNames = []
    if refAnt != 0:
        origRefAnt = refAnt
        for i,ant in enumerate(aa.ants):
            if origRefAnt == ant.stand:
                refAnt = i
                found = True
            else:
                tauNames.append('tau%i' % ant.stand)
                phsNames.append('phs%i' % ant.stand)
    else:
        found = True
    if not found:
        raise RuntimeError("Stand #%i not found in the array provided" % refAnt)
    if verbose:
        print "Using antenna #%i as a reference (Stand #%i)" % (refAnt, aa.ants[refAnt].stand)
        
    # Frequency in GHz so that the delays can be in ns
    fq = dataDict['freq'][chan] / 1e9
    
    tempGains = numpy.ones(N)
    tempDelays = numpy.zeros(N)
    tempPhaseOffsets = numpy.zeros(N)
    for i in xrange(nIter):
        #
        # Amplitude
        #
        if amplitude:
            if verbose:
                print '  %iA' % (i+1,)
                
            A = _buildAmplitudeA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            C = _buildAmplitudeC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestGains, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestGains)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            bestGains = numpy.exp(bestGains)
            tempGains *= bestGains
            
            valid = numpy.where( numpy.abs(bestGains) < 1e6 )[0]
            metric = (numpy.abs(bestGains[valid])).max()
            if verbose:
                print '    ', metric
            if metric < amplitudeCutoff:
                amplitude = False
                
            dataDict = _scaleData(dataDict, bestGains, numpy.zeros_like(bestGains), numpy.zeros_like(bestGains))
        
        #
        # Delay and/or phase
        #
        if phaseOnly:
            if verbose:
                print '  %iP' % (i+1,)
                
            A = _buildPhaseOnlyA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            C = _buildPhaseOnlyC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestPhaseOffsets, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestPhaseOffsets)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestPhaseOffsets = list(bestPhaseOffsets)
            bestPhaseOffsets.insert(refAnt, 0.0)
            bestPhaseOffsets = numpy.array(bestPhaseOffsets)
            tempPhaseOffsets += bestPhaseOffsets
            
            valid = valid = numpy.where( numpy.abs(bestPhaseOffsets) < 1e6 )[0]
            metric = (numpy.abs(bestPhaseOffsets[valid])).max()
            if verbose:
                print '    ', metric
            if metric < phaseCutoff:
                phaseOnly = False
                
            dataDict = _scaleData(dataDict, numpy.ones_like(bestPhaseOffsets), numpy.zeros_like(bestPhaseOffsets), bestPhaseOffsets)
            
        elif delayOnly:
            if verbose:
                print '  %iD' % (i+1,)
                
            A = _buildDelayOnlyA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            C = _buildDelayOnlyC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestDelays, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestDelays)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestDelays = list(bestDelays)
            bestDelays.insert(refAnt, 0.0)
            bestDelays = numpy.array(bestDelays)
            tempDelays += bestDelays
            
            valid = numpy.where( numpy.abs(bestDelays) < 1e6 )[0]
            metric = (numpy.abs(bestDelays[valid])).max()
            if verbose:
                print '    ', metric
            if metric < delayCutoff:
                delayOnly = False
                
            dataDict = _scaleData(dataDict, numpy.ones_like(bestDelays), bestDelays, numpy.zeros_like(bestDelays))
            
        elif delayAndPhase:
            if verbose:
                print '  %iD+P' % (i+1,)
                
            A = _buildDelayAndPhaseA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            C = _buildDelayAndPhaseC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestDelaysAndPhaseOffsets, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestDelaysAndPhaseOffsets)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestDelays = list(bestDelaysAndPhaseOffsets[:(N-1)])
            bestDelays.insert(refAnt, 0.0)
            bestDelays = numpy.array(bestDelays)
            tempDelays += bestDelays
            
            bestPhaseOffsets = list(bestDelaysAndPhaseOffsets[(N-1):])
            bestPhaseOffsets.insert(refAnt, 0.0)
            bestPhaseOffsets = numpy.array(bestPhaseOffsets)
            tempPhaseOffsets += bestPhaseOffsets
            
            valid = numpy.where( numpy.abs(bestDelays) < 1e6 )[0]
            metric1 = (numpy.abs(bestDelays[valid])).max()
            metric2 = (numpy.abs(bestPhaseOffsets[valid])).max()
            if verbose:
                print '    ', metric1, metric2
            if metric1 < delayCutoff and metric2 < phaseCutoff:
                delayAndPhase = False
                
            dataDict = _scaleData(dataDict, numpy.ones_like(bestDelays), bestDelays, bestPhaseOffsets)
            
        else:
            pass
            
    # Make sure the phase is (-pi, pi]
    tempPhaseOffsets %= 2*numpy.pi
    tempPhaseOffsets[numpy.where( tempPhaseOffsets >  numpy.pi )] -= 2*numpy.pi
    
    if verbose:
        print 'Best Gains: ', tempGains
    bestGains = tempGains
    
    if verbose:
        print 'Best Delays: ', tempDelays
    bestDelays = tempDelays
    
    if verbose:
        print 'Best Phase Offsets: ', tempPhaseOffsets
    bestPhaseOffsets = tempPhaseOffsets
    
    return dataDict, bestGains, bestDelays, bestPhaseOffsets


def phaseOnly(aa, dataDict, simDict, chan, pol, refAnt=0, nIter=30, amplitude=False, amplitudeCutoff=1.001, phaseCutoff=0.01, verbose=True):
    """
    Function to apply a phase-only (and, optionally, a amplitude) self-
    calibration to data stored in a readUVData dictionary and a model sky 
    stored in a lsl.sim.vis.buildSimSky dictionary for the given 
    polarization and channel(s).
    
    .. note::
        If the "amplitude" keyword is set to True, a three-element tuple of 
        self-cal'd data, gain coefficients, and phase offsets is returned 
        rather than the standard two-element tuple.
    """
    
    caldDict, gains, delays, phaseOffsets = _selfCal(aa, dataDict, simDict, chan, pol, refAnt=refAnt, nIter=nIter, 
                                            amplitude=amplitude, phaseOnly=True, 
                                            amplitudeCutoff=amplitudeCutoff, phaseCutoff=phaseCutoff, 
                                            verbose=verbose)
                                            
    if amplitude:
        return caldDict, gains, phaseOffsets
    else:
        return caldDict, phaseOffsets


def delayOnly(aa, dataDict, simDict, chan, pol, refAnt=0, nIter=30, amplitude=False, amplitudeCutoff=1.001, delayCutoff=0.2, verbose=True):
    """
    Function to apply a delay-only (and, optionally, a amplitude) self-
    calibration to data stored in a readUVData dictionary and a model sky 
    stored in a lsl.sim.vis.buildSimSky dictionary for the given 
    polarization and channel(s).
    
    .. note::
        If the "amplitude" keyword is set to True, a three-element tuple of 
        self-cal'd data, gain coefficients, and delays in ns is returned rather 
        than the standard two-element tuple.
    """
    
    caldDict, gains, delays, phaseOffsets = _selfCal(aa, dataDict, simDict, chan, pol, refAnt=refAnt, nIter=nIter, 
                                            amplitude=amplitude, delayOnly=True, 
                                            amplitudeCutoff=amplitudeCutoff, delayCutoff=delayCutoff,
                                            verbose=verbose)
                                            
    if amplitude:
        return caldDict, gains, delays
    else:
        return caldDict, delays


def delayAndPhase(aa, dataDict, simDict, chan, pol, refAnt=0, nIter=30, amplitude=False, amplitudeCutoff=1.001, delayCutoff=0.2, phaseCutoff=0.01, verbose=True):
    """
    Function to apply a delay and phase offset (and, optionally, a amplitude)
    self-calibration to data stored in a readUVData dictionary and a model 
    sky stored in a lsl.sim.vis.buildSimSky dictionary for the given 
    polarization and channel(s).
    
    .. note::
        If the "amplitude" keyword is set to True, a four-element tuple of 
        self-cal'd data, gain coefficients, delays in ns, and phase offsets 
        is returned rather than the standard three-element tuple.
    """
    
    caldDict, gains, delays, phaseOffsets = _selfCal(aa, dataDict, simDict, chan, pol, refAnt=refAnt, nIter=nIter, 
                                            amplitude=amplitude, delayAndPhase=True, 
                                            amplitudeCutoff=amplitudeCutoff, delayCutoff=delayCutoff, phaseCutoff=phaseCutoff, 
                                            verbose=verbose)
                                            
    if amplitude:
        return caldDict, gains, delays, phaseOffsets
    else:
        return caldDict, delays, phaseOffsets
