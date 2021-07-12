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

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import numpy

from lsl.statistics import robust

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['phase_only', 'delay_only', 'delay_and_phase']


def _scale_data(dataSet, amps, delays, phase_offsets):
    """
    Apply a set of antenna-based real gain value, phase delays in ns, phase 
    offsets in radians to a data dictionary.  Returned the new scaled and 
    delayed dictionary.
    """

    import copy

    # Build the data dictionary to hold the scaled and delayed data
    sclData = dataSet.copy(include_pols=True)
    fq = dataSet.freq / 1e9
    
    cGains = []
    for i in range(len(amps)):
        cGains.append( amps[i]*numpy.exp(2j*numpy.pi*fq*delays[i] + 1j*phase_offsets[i]) )

    # Apply the scales and delays for all polarization pairs found in the original data
    for pds in sclData:
        for b,(i,j) in enumerate(sclData.baselines):
            pds.data[b,:] *= cGains[j].conj()*cGains[i]
            
    return sclData


def _build_amplitude_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix A for amplitude correction.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    A = numpy.zeros((nBLs, nStands))
    for i,(l,m) in enumerate(dataSet.baselines):
        A[i,l] = 1.0
        A[i,m] = 1.0
        
    A = numpy.matrix(A)
    return A


def _build_amplitude_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix C for the amplitude correction.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    obsVis = []
    for vis in getattr(dataSet, pol).data:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.log(numpy.abs(simVis / obsVis))
    
    Cp = numpy.zeros(nBLs)
    for i in range(C.shape[0]):
        try:
            Cp[i] = robust.mean(C[i,:])
        except ValueError:
            Cp[i] = numpy.mean(C[i,:])
            
    return Cp


def _build_phaseonly_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix A for phase correction with a frequency independent
    phase.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    A = numpy.zeros((nBLs*fq.size, nStands-1))
    for i in range(fq.size):
        for j,(l,m) in enumerate(dataSet.baselines):
            if l < ref_ant:
                A[j+i*nBLs,l]   =  1.0
            elif l > ref_ant:
                A[j+i*nBLs,l-1] =  1.0
            else:
                pass
                
            if m < ref_ant:
                A[j+i*nBLs,m]   = -1.0
            elif m > ref_ant:
                A[j+i*nBLs,m-1] = -1.0
            else:
                pass
        
    A = numpy.matrix(A)
    return A


def _build_phaseonly_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix C for phase correction with a frequency independent
    phase.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    obsVis = []
    for vis in getattr(dataSet, pol).data:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.angle(simVis / obsVis)
    
    Cp = numpy.zeros(nBLs*fq.size)
    for i in range(fq.size):
        for j in range(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _build_delayonly_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix A for phase correction with a delay.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    A = numpy.zeros((nBLs*fq.size, nStands-1))
    for i in range(fq.size):
        for j,(l,m) in enumerate(dataSet.baselines):
            if l < ref_ant:
                A[j+i*nBLs,l]   =  2*numpy.pi*fq[i]
            elif l > ref_ant:
                A[j+i*nBLs,l-1] =  2*numpy.pi*fq[i]
            else:
                pass
                
            if m < ref_ant:
                A[j+i*nBLs,m]   = -2*numpy.pi*fq[i]
            elif m > ref_ant:
                A[j+i*nBLs,m-1] = -2*numpy.pi*fq[i]
            else:
                pass
        
    A = numpy.matrix(A)
    return A


def _build_delayonly_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix C for phase correction with a delay.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    obsVis = []
    for vis in getattr(dataSet, pol).data:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.angle(simVis / obsVis)
    
    Cp = numpy.zeros(nBLs*fq.size)
    for i in range(fq.size):
        for j in range(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _build_delayandphase_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix A for phase correction with a delay and a phase offset.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    A = numpy.zeros((nBLs*fq.size, 2*(nStands-1)))
    for i in range(fq.size):
        for j,(l,m) in enumerate(dataSet.baselines):
            if l < ref_ant:
                A[j+i*nBLs,l]   =  2*numpy.pi*fq[i]
                A[j+i*nBLs,l+(nStands-1)] = 1.0
            elif l > ref_ant:
                A[j+i*nBLs,l-1] =  2*numpy.pi*fq[i]
                A[j+i*nBLs,l-1+(nStands-1)] = 1.0
                
            else:
                pass
                
            if m < ref_ant:
                A[j+i*nBLs,m]   = -2*numpy.pi*fq[i]
                A[j+i*nBLs,m+(nStands-1)] = -1.0
            elif m > ref_ant:
                A[j+i*nBLs,m-1] = -2*numpy.pi*fq[i]
                A[j+i*nBLs,m-1+(nStands-1)] = -1.0
            else:
                pass
        
    A = numpy.matrix(A)
    return A


def _build_delayandphase_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix C for phase correction with a delay and a phase offset.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    obsVis = []
    for vis in getattr(dataSet, pol).data:
        obsVis.append( numpy.array(vis[chan]) )
    obsVis = numpy.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( numpy.array(vis[chan]) )
    simVis = numpy.array(simVis)
    
    C = numpy.angle(simVis / obsVis)
    
    Cp = numpy.zeros(nBLs*fq.size)
    for i in range(fq.size):
        for j in range(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=0, max_iter=30, amplitude=False, phase_only=False, delay_only=False, delay_and_phase=False, amplitude_cutoff=1.001, phase_cutoff=0.01, delay_cutoff=0.2, verbose=True):
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
    reached (max_iter) or the maximum absolute value of the relevant quantity
    drops below the quantity's cutoff value.  The cutoff keywords and 
    their default values are:
     * amplitude_cutoff - 1.001, 
     * phase_cutoff - 0.01 radians, 
     * delay_cutoff - 0.2 ns (0.01 radians over 10 MHz),
    """

    # Make sure we have the right polarization
    if pol not in dataSet.pols:
        raise RuntimeError("Data set does not have data for polarization '%s'" % pol)
    if pol not in simSet.pols:
        raise RuntimeError("Simulation set does not have data for polarization '%s'" % pol)

    # Make sure that `chan' is an array by trying to find its length
    try:
        junk = len(chan)
    except TypeError:
        chan = [chan]
        
    N = len(aa.ants)
    found = False
    tauNames = []
    phsNames = []
    if ref_ant != 0:
        origRefAnt = ref_ant
        for i,ant in enumerate(aa.ants):
            if origRefAnt == ant.stand:
                ref_ant = i
                found = True
            else:
                tauNames.append('tau%i' % ant.stand)
                phsNames.append('phs%i' % ant.stand)
    else:
        found = True
    if not found:
        raise RuntimeError("Stand #%i not found in the array provided" % ref_ant)
    if verbose:
        print("Using antenna #%i as a reference (Stand #%i)" % (ref_ant, aa.ants[ref_ant].stand))
        
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    tempGains = numpy.ones(N)
    tempDelays = numpy.zeros(N)
    tempPhaseOffsets = numpy.zeros(N)
    for i in range(max_iter):
        #
        # Amplitude
        #
        if amplitude:
            if verbose:
                print('  %iA' % (i+1,))
                
            A = _build_amplitude_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_amplitude_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
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
                print('    ', metric)
            if metric < amplitude_cutoff:
                amplitude = False
                
            dataSet = _scale_data(dataSet, bestGains, numpy.zeros_like(bestGains), numpy.zeros_like(bestGains))
        
        #
        # Delay and/or phase
        #
        if phase_only:
            if verbose:
                print('  %iP' % (i+1,))
                
            A = _build_phaseonly_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_phaseonly_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestPhaseOffsets, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestPhaseOffsets)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestPhaseOffsets = list(bestPhaseOffsets)
            bestPhaseOffsets.insert(ref_ant, 0.0)
            bestPhaseOffsets = numpy.array(bestPhaseOffsets)
            tempPhaseOffsets += bestPhaseOffsets
            
            valid = valid = numpy.where( numpy.abs(bestPhaseOffsets) < 1e6 )[0]
            metric = (numpy.abs(bestPhaseOffsets[valid])).max()
            if verbose:
                print('    ', metric)
            if metric < phase_cutoff:
                phase_only = False
                
            dataSet = _scale_data(dataSet, numpy.ones_like(bestPhaseOffsets), numpy.zeros_like(bestPhaseOffsets), bestPhaseOffsets)
            
        elif delay_only:
            if verbose:
                print('  %iD' % (i+1,))
                
            A = _build_delayonly_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_delayonly_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestDelays, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestDelays)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestDelays = list(bestDelays)
            bestDelays.insert(ref_ant, 0.0)
            bestDelays = numpy.array(bestDelays)
            tempDelays += bestDelays
            
            valid = numpy.where( numpy.abs(bestDelays) < 1e6 )[0]
            metric = (numpy.abs(bestDelays[valid])).max()
            if verbose:
                print('    ', metric)
            if metric < delay_cutoff:
                delay_only = False
                
            dataSet = _scale_data(dataSet, numpy.ones_like(bestDelays), bestDelays, numpy.zeros_like(bestDelays))
            
        elif delay_and_phase:
            if verbose:
                print('  %iD+P' % (i+1,))
                
            A = _build_delayandphase_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_delayandphase_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = numpy.where( numpy.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            bestDelaysAndPhaseOffsets, resid, rank, s = numpy.linalg.lstsq(A, C)
            resid = numpy.array(C - numpy.dot(A, bestDelaysAndPhaseOffsets)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestDelays = list(bestDelaysAndPhaseOffsets[:(N-1)])
            bestDelays.insert(ref_ant, 0.0)
            bestDelays = numpy.array(bestDelays)
            tempDelays += bestDelays
            
            bestPhaseOffsets = list(bestDelaysAndPhaseOffsets[(N-1):])
            bestPhaseOffsets.insert(ref_ant, 0.0)
            bestPhaseOffsets = numpy.array(bestPhaseOffsets)
            tempPhaseOffsets += bestPhaseOffsets
            
            valid = numpy.where( numpy.abs(bestDelays) < 1e6 )[0]
            metric1 = (numpy.abs(bestDelays[valid])).max()
            metric2 = (numpy.abs(bestPhaseOffsets[valid])).max()
            if verbose:
                print('    ', metric1, metric2)
            if metric1 < delay_cutoff and metric2 < phase_cutoff:
                delay_and_phase = False
                
            dataSet = _scale_data(dataSet, numpy.ones_like(bestDelays), bestDelays, bestPhaseOffsets)
            
        else:
            pass
            
    # Make sure the phase is (-pi, pi]
    tempPhaseOffsets %= 2*numpy.pi
    tempPhaseOffsets[numpy.where( tempPhaseOffsets >  numpy.pi )] -= 2*numpy.pi
    
    if verbose:
        print('Best Gains: ', tempGains)
    bestGains = tempGains
    
    if verbose:
        print('Best Delays: ', tempDelays)
    bestDelays = tempDelays
    
    if verbose:
        print('Best Phase Offsets: ', tempPhaseOffsets)
    bestPhaseOffsets = tempPhaseOffsets
    
    return dataSet, bestGains, bestDelays, bestPhaseOffsets


def phase_only(aa, dataSet, simSet, chan, pol, ref_ant=0, max_iter=30, amplitude=False, amplitude_cutoff=1.001, phase_cutoff=0.01, verbose=True):
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
    
    caldDict, gains, delays, phaseOffsets = _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, max_iter=max_iter, 
                                            amplitude=amplitude, phase_only=True, 
                                            amplitude_cutoff=amplitude_cutoff, phase_cutoff=phase_cutoff, 
                                            verbose=verbose)
                                            
    if amplitude:
        return caldDict, gains, phaseOffsets
    else:
        return caldDict, phaseOffsets


def delay_only(aa, dataSet, simSet, chan, pol, ref_ant=0, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, verbose=True):
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
    
    caldDict, gains, delays, phaseOffsets = _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, max_iter=max_iter, 
                                            amplitude=amplitude, delay_only=True, 
                                            amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff,
                                            verbose=verbose)
                                            
    if amplitude:
        return caldDict, gains, delays
    else:
        return caldDict, delays


def delay_and_phase(aa, dataSet, simSet, chan, pol, ref_ant=0, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, phase_cutoff=0.01, verbose=True):
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
    
    caldDict, gains, delays, phaseOffsets = _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, max_iter=max_iter, 
                                            amplitude=amplitude, delay_and_phase=True, 
                                            amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff, phase_cutoff=phase_cutoff, 
                                            verbose=verbose)
                                            
    if amplitude:
        return caldDict, gains, delays, phaseOffsets
    else:
        return caldDict, delays, phaseOffsets
