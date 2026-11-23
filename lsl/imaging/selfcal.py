"""
Simple self-calibration module for correlated LWA data.  The 
supported self-calibration methods are:
 * phase-only
 * amplitude and phase
 * delay-only
 * amplitude and delay
 * delay/phase offset
 * amplitude and delay/phase offset

.. versionchanged:: 4.0.2
   Fixed the amplitude solutions when combined with phase or delay
   Solve for the delays and phase offsets jointly in delay_and_phase
   Switched to sparse matrices and scipy.sparse.linalg.lsqr for the solutions

.. versionchanged:: 3.0.8
   Added support for Tikhonov regularization and exposed the loop gain

..versionchanged:: 0.6.3
    Reworked the module to make it more versatile
    
..versionadded:: 0.5.5
"""

import logging
import numpy as np

from scipy import sparse
from scipy.sparse.linalg import lsqr

from lsl.statistics import robust
from lsl.logger import LSL_LOGGER

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
    Build the sparse matrix A for amplitude correction.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    rows, cols, vals = [], [], []
    for i,(l,m) in enumerate(dataSet.baselines):
        rows.extend([i, i])
        cols.extend([l, m])
        vals.extend([1.0, 1.0])
        
    return sparse.csr_matrix((vals, (rows, cols)), shape=(nBLs, nStands))


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
    Build the sparse matrix A for phase correction with a frequency 
    independent phase.
    """
    
    # Get the baseline and stand counts
    nBLs = dataSet.nbaseline
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    # Build the entries for a single channel...
    rows, cols, vals = [], [], []
    for j,(l,m) in enumerate(dataSet.baselines):
        if l < ref_ant:
            rows.append(j)
            cols.append(l)
            vals.append( 1.0)
        elif l > ref_ant:
            rows.append(j)
            cols.append(l-1)
            vals.append( 1.0)
        else:
            pass
            
        if m < ref_ant:
            rows.append(j)
            cols.append(m)
            vals.append(-1.0)
        elif m > ref_ant:
            rows.append(j)
            cols.append(m-1)
            vals.append(-1.0)
        else:
            pass
            
    # ... and then repeat them for each channel
    rows = np.concatenate([np.array(rows) + i*nBLs for i in range(fq.size)])
    cols = np.tile(cols, fq.size)
    vals = np.tile(vals, fq.size)
    
    return sparse.csr_matrix((vals, (rows, cols)), shape=(nBLs*fq.size, nStands-1))


def _build_phase_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix C for phase correction with a frequency independent
    phase.
    """
    
    obsVis = []
    for vis in getattr(dataSet, pol).data:
        obsVis.append( np.array(vis[chan]) )
    obsVis = np.array(obsVis)
    
    simVis = []
    for vis in getattr(simSet, pol).data:
        simVis.append( np.array(vis[chan]) )
    simVis = np.array(simVis)
    
    C = np.angle(simVis / obsVis)
    
    # Order by channel, then baseline, to match A
    return C.T.ravel()


def _build_delay_a(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the sparse matrix A for phase correction with a delay.  This is 
    the same as the matrix for the phase correction since the frequency 
    dependence is included in C.
    """
    
    return _build_phase_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)


def _build_delay_c(aa, dataSet, simSet, chan, pol, ref_ant=0):
    """
    Build the matrix C for phase correction with a delay.
    """
    
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
    
    C = np.angle(simVis / obsVis) / (2*np.pi*fq)
    
    # Order by channel, then baseline, to match A
    return C.T.ravel()


def _solve(A, C, inv_epsilon=0.0):
    """
    Solve A x = C in a least squares sense for the finite values of C.  If
    inv_epsilon is greater than zero, Tikhonov regularization is applied
    with a damping factor of 1/inv_epsilon.
    """
    
    good = np.where( np.isfinite(C) )[0]
    damp = 1/inv_epsilon if inv_epsilon > 0 else 0.0
    
    return lsqr(A[good,:], C[good], damp=damp, atol=1e-10, btol=1e-10)[0]


def _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, phase_only=False, delay_only=False, delay_and_phase=False, amplitude_cutoff=1.001, phase_cutoff=0.01, delay_cutoff=0.2, inv_epsilon=0.0):
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
    LSL_LOGGER.info(f"Using antenna #{ref_ant} as a reference (Stand #{aa.ants[ref_ant].stand})")
        
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
            LSL_LOGGER.debug(f'  {i+1}A')
            
            A = _build_amplitude_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_amplitude_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            bestGains = _solve(A, C, inv_epsilon=inv_epsilon)
            bestGains = np.exp(bestGains)
            tempGains *= bestGains
            
            valid = np.where( np.abs(bestGains) < 1e6 )[0]
            metric = (np.abs(bestGains[valid])).max()
            LSL_LOGGER.debug(f'    {metric}')
            if metric < amplitude_cutoff:
                amplitude = False
                converged = not (phase_only or delay_only or delay_and_phase)
                
            dataSet = _scale_data(origSet, tempGains, tempDelays, tempPhaseOffsets)
            
        #
        # Delay and/or phase
        #
        if phase_only:
            LSL_LOGGER.debug(f'  {i+1}P')
            
            A = _build_phase_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_phase_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            bestPhaseOffsets = _solve(A, C, inv_epsilon=inv_epsilon)
            
            bestPhaseOffsets = list(bestPhaseOffsets)
            bestPhaseOffsets.insert(ref_ant, 0.0)
            bestPhaseOffsets = np.array(bestPhaseOffsets)
            tempPhaseOffsets += loop_gain*bestPhaseOffsets
            
            valid = np.where( np.abs(bestPhaseOffsets) < 1e6 )[0]
            metric = (np.abs(bestPhaseOffsets[valid])).max()
            LSL_LOGGER.debug(f'    {metric}')
            if metric < phase_cutoff:
                phase_only = False
                converged = not amplitude
                
            dataSet = _scale_data(origSet, tempGains, tempDelays, tempPhaseOffsets)
            
        elif delay_only:
            LSL_LOGGER.debug(f'  {i+1}D')
            
            A = _build_delay_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            C = _build_delay_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            bestDelays = _solve(A, C, inv_epsilon=inv_epsilon)
            
            bestDelays = list(bestDelays)
            bestDelays.insert(ref_ant, 0.0)
            bestDelays = np.array(bestDelays)
            tempDelays += loop_gain*bestDelays
            
            valid = np.where( np.abs(bestDelays) < 1e6 )[0]
            metric = (np.abs(bestDelays[valid])).max()
            LSL_LOGGER.debug(f'    {metric}')
            if metric < delay_cutoff:
                delay_only = False
                converged = not amplitude
                
            dataSet = _scale_data(origSet, tempGains, tempDelays, tempPhaseOffsets)
            
        elif delay_and_phase:
            LSL_LOGGER.debug(f'  {i+1}D+P')
            
            ## Solve for the delays and phase offsets jointly.  The delay block
            ## is scaled by 2*pi*f since the fit is done in phase, not delay.
            Ad = _build_delay_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            Ap = _build_phase_a(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            A = sparse.hstack([sparse.diags(np.repeat(2*np.pi*fq, dataSet.nbaseline)) @ Ad, Ap], format='csr')
            C = _build_phase_c(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant)
            
            bestDelaysAndPhaseOffsets = _solve(A, C, inv_epsilon=inv_epsilon)
            
            bestDelays = list(bestDelaysAndPhaseOffsets[:N-1])
            bestDelays.insert(ref_ant, 0.0)
            bestDelays = np.array(bestDelays)
            tempDelays += loop_gain*bestDelays
            
            bestPhaseOffsets = list(bestDelaysAndPhaseOffsets[N-1:])
            bestPhaseOffsets.insert(ref_ant, 0.0)
            bestPhaseOffsets = np.array(bestPhaseOffsets)
            tempPhaseOffsets += loop_gain*bestPhaseOffsets
            
            valid = np.where( (np.abs(bestDelays) < 1e6) & (np.abs(bestPhaseOffsets) < 1e6) )[0]
            metric1 = (np.abs(bestDelays[valid])).max()
            metric2 = (np.abs(bestPhaseOffsets[valid])).max()
            LSL_LOGGER.debug(f'    {metric1} {metric2}')
            if metric1 < delay_cutoff and metric2 < phase_cutoff:
                delay_and_phase = False
                converged = not amplitude
                
            dataSet = _scale_data(origSet, tempGains, tempDelays, tempPhaseOffsets)
            
        else:
            pass
            
    # Did we converge?
    if converged:
        LSL_LOGGER.info('Fit converged')
    else:
        LSL_LOGGER.warning('Fit failed to convege')
        
    # Make sure the phase is (-pi, pi]
    tempPhaseOffsets %= 2*np.pi
    tempPhaseOffsets[np.where( tempPhaseOffsets >  np.pi )] -= 2*np.pi

    if amplitude:
        LSL_LOGGER.info(f'Best Gains: {tempGains}')
    bestGains = tempGains

    if delay_only or delay_and_phase:
        LSL_LOGGER.info(f'Best Delays: {tempDelays}')
    bestDelays = tempDelays

    if phase_only or delay_and_phase:
        LSL_LOGGER.info(f'Best Phase Offsets: {tempPhaseOffsets}')
    bestPhaseOffsets = tempPhaseOffsets
    
    return dataSet, bestGains, bestDelays, bestPhaseOffsets, converged


def phase_only(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, phase_cutoff=0.01, inv_epsilon=0.0, return_convergence=False):
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
                inv_epsilon=inv_epsilon)
    
                           
    if amplitude:
        output = (caldDict, gains, phaseOffsets)
    else:
        output = (caldDict, phaseOffsets)
    if return_convergence:
        output = output+(converged,)
    return output


def delay_only(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, inv_epsilon=0.0, return_convergence=False):
    """
    Function to apply a delay-only (and, optionally, a amplitude) self-
    calibration to data stored in a readUVData dictionary and a model sky 
    stored in a lsl.sim.vis.buildSimSky dictionary for the given 
    polarization and channel(s).
    
    .. note::
        If the "amplitude" keyword is set to True, a three-element tuple of 
        self-cal'd data, gain coefficients, and delays in ns is returned rather 
        than the standard two-element tuple.
        
    The delays returned are positive for an antenna whose signal arrives 
    late, and can be passed to lsl.sim.vis.scale_data to correct the data.
    """
    
    caldDict, gains, delays, phaseOffsets, converged = \
      _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, 
                loop_gain=loop_gain, max_iter=max_iter, 
                amplitude=amplitude, delay_only=True, 
                amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff,
                inv_epsilon=inv_epsilon)
                                            
    if amplitude:
        output = (caldDict, gains, delays)
    else:
        output = (caldDict, delays)
    if return_convergence:
        output = output+(converged,)
    return output


def delay_and_phase(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, phase_cutoff=0.01, inv_epsilon=0.0,  return_convergence=False):
    """
    Function to apply a delay and phase offset (and, optionally, a amplitude)
    self-calibration to data stored in a readUVData dictionary and a model 
    sky stored in a lsl.sim.vis.buildSimSky dictionary for the given 
    polarization and channel(s).
    
    .. note::
        If the "amplitude" keyword is set to True, a four-element tuple of 
        self-cal'd data, gain coefficients, delays in ns, and phase offsets 
        is returned rather than the standard three-element tuple.
        
    The delays returned are positive for an antenna whose signal arrives 
    late, and can be passed to lsl.sim.vis.scale_data to correct the data.
    """
    
    caldDict, gains, delays, phaseOffsets, converged = \
      _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, 
                loop_gain=loop_gain, max_iter=max_iter, 
                amplitude=amplitude, delay_and_phase=True, 
                amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff, phase_cutoff=phase_cutoff,
                inv_epsilon=inv_epsilon)
                                            
    if amplitude:
        output = (caldDict, gains, delays, phaseOffsets)
    else:
        output = (caldDict, delays, phaseOffsets)
    if return_convergence:
        output = output+(converged,)
    return output
