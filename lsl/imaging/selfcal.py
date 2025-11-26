"""
Simple self-calibration module for correlated TBW and TBN data.  The 
supported self-calibration methods are:
 * phase-only
 * amplitude and phase
 * delay-only
 * amplitude and delay
 * delay/phase offset
 * amplitude and delay/phase offset

.. versionchanged:: 3.0.8
   Added support for Tikhonov regularization and exposed the loop gain

..versionchanged:: 0.6.3
    Reworked the module to make it more versatile
    
..versionadded:: 0.5.5
"""

import numpy as np

from scipy.linalg import lstsq, decomp, decomp_svd, pinv

from lsl.statistics import robust

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
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
        cGains.append( amps[i]*np.exp(2j*np.pi*fq*delays[i] + 1j*phase_offsets[i]) )

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
    
    A = np.zeros((nBLs, nStands))
    for i,(l,m) in enumerate(dataSet.baselines):
        A[i,l] = 1.0
        A[i,m] = 1.0
        
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
        obsVis.append( np.array(vis[chan]) )
    obsVis = np.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( np.array(vis[chan]) )
    simVis = np.array(simVis)
    
    C = np.log(np.abs(simVis / obsVis))
    
    Cp = np.zeros(nBLs)
    for i in range(C.shape[0]):
        try:
            Cp[i] = robust.mean(C[i,:])
        except ValueError:
            Cp[i] = np.mean(C[i,:])
            
    return Cp


def _build_phase_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix A for phase correction with a frequency independent
    phase.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    A = np.zeros((nBLs*fq.size, nStands-1))
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
        
    return A


def _build_phase_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
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
        obsVis.append( np.array(vis[chan]) )
    obsVis = np.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( np.array(vis[chan]) )
    simVis = np.array(simVis)
    
    C = np.angle(simVis / obsVis)
    
    Cp = np.zeros(nBLs*fq.size)
    for i in range(fq.size):
        for j in range(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i]
    
    return Cp


def _build_delay_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix A for phase correction with a delay.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    A = np.zeros((nBLs*fq.size, nStands-1))
    for i in range(fq.size):
        for j,(l,m) in enumerate(dataSet.baselines):
            if l < ref_ant:
                A[j+i*nBLs,l]   =  1
            elif l > ref_ant:
                A[j+i*nBLs,l-1] =  1
            else:
                pass
                
            if m < ref_ant:
                A[j+i*nBLs,m]   = -1
            elif m > ref_ant:
                A[j+i*nBLs,m-1] = -1
            else:
                pass
        
    return A


def _build_delay_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
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
        obsVis.append( np.array(vis[chan]) )
    obsVis = np.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( np.array(vis[chan]) )
    simVis = np.array(simVis)
    
    C = np.angle(simVis / obsVis)
    
    Cp = np.zeros(nBLs*fq.size)
    for i in range(fq.size):
        for j in range(C.shape[0]):
            Cp[j+i*nBLs] = C[j,i] / (2*np.pi*fq[i])
    
    return Cp


def _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, phase_only=False, delay_only=False, delay_and_phase=False, amplitude_cutoff=1.001, phase_cutoff=0.01, delay_cutoff=0.2, inv_epsilon=0.0, verbose=True):
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
     * phase_cutoff - 0.01 radians, and
     * delay_cutoff - 0.2 ns (0.01 radians over 10 MHz)
     
    Tikhonov regularization is also supported with the `inv_epsilon` keyword
    that controls the degree regularization.  For details see:
      https://en.wikipedia.org/wiki/Ridge_regpression
    Setting `inv_epsilon` to zero (the default) disables this feature.
    """

    # Make sure we have the right polarization
    if pol not in dataSet.pols:
        raise RuntimeError(f"Data set does not have data for polarization '{pol}'")
    if pol not in simSet.pols:
        raise RuntimeError(f"Simulation set does not have data for polarization '{pol}'")

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
        raise RuntimeError(f"Stand #{ref_ant} not found in the array provided")
    if verbose:
        print(f"Using antenna #{ref_ant} as a reference (Stand #{aa.ants[ref_ant].stand})")
        
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    origSet = dataSet
    
    converged = False
    tempGains = np.ones(N)
    tempDelays = np.zeros(N)
    tempPhaseOffsets = np.zeros(N)
    for i in range(max_iter):
        if converged:
            break
            
        #
        # Amplitude
        #
        if amplitude:
            if verbose:
                print('  %iA' % (i+1,))
                
            A = _build_amplitude_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_amplitude_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = np.where( np.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            if inv_epsilon > 0:
                R = 1/inv_epsilon * np.eye(*A.shape)
                trc, rank = pinv(np.dot(A.T,A) + np.dot(R.T,R), return_rank=True)
                bestGains = np.dot(trc, np.dot(A.T,C))
                
            else:
                bestGains, resid, rank, s = np.linalg.lstsq(A, C, rcond=None)
            resid = np.array(C - np.dot(A, bestGains)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            bestGains = np.exp(bestGains)
            tempGains *= bestGains
            
            valid = np.where( np.abs(bestGains) < 1e6 )[0]
            metric = (np.abs(bestGains[valid])).max()
            if verbose:
                print('    ', metric)
            if metric < amplitude_cutoff:
                amplitude = False
                converged = True
                
            dataSet = _scale_data(origSet, tempGains, np.zeros_like(tempGains), np.zeros_like(tempGains))
        
        #
        # Delay and/or phase
        #
        if phase_only:
            if verbose:
                print('  %iP' % (i+1,))
                
            A = _build_phase_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_phase_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = np.where( np.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            if inv_epsilon > 0:
                R = 1/inv_epsilon * np.eye(*A.shape)
                trc, rank = pinv(np.dot(A.T,A) + np.dot(R.T,R), return_rank=True)
                bestPhaseOffsets = np.dot(trc, np.dot(A.T,C))
                
            else:
                bestPhaseOffsets, resid, rank, s = np.linalg.lstsq(A, C, rcond=None)
            resid = np.array(C - np.dot(A, bestPhaseOffsets)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestPhaseOffsets = list(bestPhaseOffsets)
            bestPhaseOffsets.insert(ref_ant, 0.0)
            bestPhaseOffsets = np.array(bestPhaseOffsets)
            tempPhaseOffsets += loop_gain*bestPhaseOffsets
            
            valid = np.where( np.abs(bestPhaseOffsets) < 1e6 )[0]
            metric = (np.abs(bestPhaseOffsets[valid])).max()
            if verbose:
                print('    ', metric)
            if metric < phase_cutoff:
                phase_only = False
                converged = True
                
            dataSet = _scale_data(origSet, np.ones_like(tempPhaseOffsets), np.zeros_like(tempPhaseOffsets), tempPhaseOffsets)
            
        elif delay_only:
            if verbose:
                print('  %iD' % (i+1,))
                
            A = _build_delay_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_delay_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = np.where( np.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            if inv_epsilon > 0:
                R = 1/inv_epsilon * np.eye(*A.shape)
                trc, rank = pinv(np.dot(A.T,A) + np.dot(R.T,R), return_rank=True)
                bestDelays = np.dot(trc, np.dot(A.T,C))
                
            else:
                bestDelays, resid, rank, s = np.linalg.lstsq(A, C, rcond=None)
            resid = np.array(C - np.dot(A, bestDelays)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestDelays = list(bestDelays)
            bestDelays.insert(ref_ant, 0.0)
            bestDelays = np.array(bestDelays)
            tempDelays += loop_gain*bestDelays
            
            valid = np.where( np.abs(bestDelays) < 1e6 )[0]
            metric = (np.abs(bestDelays[valid])).max()
            if verbose:
                print('    ', metric)
            if metric < delay_cutoff:
                delay_only = False
                converged = True
                
            dataSet = _scale_data(origSet, np.ones_like(tempDelays), tempDelays, np.zeros_like(tempDelays))
            
        elif delay_and_phase:
            if verbose:
                print('  %iD+P' % (i+1,))
                
            ## Delay first
            A = _build_delay_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_delay_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = np.where( np.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            if inv_epsilon > 0:
                R = 1/inv_epsilon * np.eye(*A.shape)
                trc, rank = pinv(np.dot(A.T,A) + np.dot(R.T,R), return_rank=True)
                bestDelays = np.dot(trc, np.dot(A.T,C))
                
            else:
                bestDelays, resid, rank, s = np.linalg.lstsq(A, C, rcond=None)
            resid = np.array(C - np.dot(A, bestDelays)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestDelays = list(bestDelays)
            bestDelays.insert(ref_ant, 0.0)
            bestDelays = np.array(bestDelays)
            tempDelays += loop_gain*bestDelays
            
            dataSet = _scale_data(origSet, np.ones_like(tempDelays), tempDelays, tempPhaseOffsets)
            
            ## Then phase
            A = _build_phase_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_phase_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            good = np.where( np.isfinite(C) == 1 )[0]
            A = A[good,:]
            C = C[good]
            
            if inv_epsilon > 0:
                R = 1/inv_epsilon * np.eye(*A.shape)
                trc, rank = pinv(np.dot(A.T,A) + np.dot(R.T,R), return_rank=True)
                bestPhaseOffsets = np.dot(trc, np.dot(A.T,C))
                
            else:
                bestPhaseOffsets, resid, rank, s = np.linalg.lstsq(A, C, rcond=None)
            resid = np.array(C - np.dot(A, bestPhaseOffsets)).ravel()
            resid = (C**2).sum(), (resid**2).sum()
            
            bestPhaseOffsets = list(bestPhaseOffsets)
            bestPhaseOffsets.insert(ref_ant, 0.0)
            bestPhaseOffsets = np.array(bestPhaseOffsets)
            tempPhaseOffsets += 0.5*bestPhaseOffsets
            
            valid = np.where( (np.abs(bestDelays) < 1e6) & (np.abs(bestPhaseOffsets) < 1e6) )[0]
            metric1 = (np.abs(bestDelays[valid])).max()
            metric2 = (np.abs(bestPhaseOffsets[valid])).max()
            if verbose:
                print('    ', metric1, metric2)
            if metric1 < delay_cutoff and metric2 < phase_cutoff:
                delay_and_phase = False
                converged = True
                
            dataSet = _scale_data(origSet, np.ones_like(tempDelays), tempDelays, tempPhaseOffsets)
            
        else:
            pass
            
    # Make sure the phase is (-pi, pi]
    tempPhaseOffsets %= 2*np.pi
    tempPhaseOffsets[np.where( tempPhaseOffsets >  np.pi )] -= 2*np.pi
    
    if verbose:
        print('Best Gains: ', tempGains)
    bestGains = tempGains
    
    if verbose:
        print('Best Delays: ', tempDelays)
    bestDelays = tempDelays
    
    if verbose:
        print('Best Phase Offsets: ', tempPhaseOffsets)
    bestPhaseOffsets = tempPhaseOffsets
    
    return dataSet, bestGains, bestDelays, bestPhaseOffsets, converged


def phase_only(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, phase_cutoff=0.01, inv_epsilon=0.0, return_convergence=False, verbose=True):
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
    
    caldDict, gains, delays, phaseOffsets, converged = \
      _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant,
                loop_gain=loop_gain, max_iter=max_iter, 
                amplitude=amplitude, phase_only=True, 
                amplitude_cutoff=amplitude_cutoff, phase_cutoff=phase_cutoff, 
                inv_epsilon=inv_epsilon, verbose=verbose)
    
                           
    if amplitude:
        output = (caldDict, gains, phaseOffsets)
    else:
        output = (caldDict, phaseOffsets)
    if return_convergence:
        output = output+(converged,)
    return output


def delay_only(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, inv_epsilon=0.0, return_convergence=False, verbose=True):
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
    
    caldDict, gains, delays, phaseOffsets, converged = \
      _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, 
                loop_gain=loop_gain, max_iter=max_iter, 
                amplitude=amplitude, delay_only=True, 
                amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff,
                inv_epsilon=inv_epsilon, verbose=verbose)
                                            
    if amplitude:
        output = (caldDict, gains, delays)
    else:
        output = (caldDict, delays)
    if return_convergence:
        output = output+(converged,)
    return output


def delay_and_phase(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, phase_cutoff=0.01, inv_epsilon=0.0,  return_convergence=False, verbose=True):
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
    
    caldDict, gains, delays, phaseOffsets, converged = \
      _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, 
                loop_gain=loop_gain, max_iter=max_iter, 
                amplitude=amplitude, delay_and_phase=True, 
                amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff, phase_cutoff=phase_cutoff, 
                inv_epsilon=inv_epsilon, verbose=verbose)
                                            
    if amplitude:
        output = (caldDict, gains, delays, phaseOffsets)
    else:
        output = (caldDict, delays, phaseOffsets)
    if return_convergence:
        output = output+(converged,)
    return output
