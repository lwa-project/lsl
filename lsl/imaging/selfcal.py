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
   Added `check_delays` to check `delay_only` solutions for cycle ambiguities

.. versionchanged:: 3.0.8
   Added support for Tikhonov regularization and exposed the loop gain

..versionchanged:: 0.6.3
    Reworked the module to make it more versatile
    
..versionadded:: 0.5.5
"""

import enum
import logging
import numpy as np

from scipy import sparse
from scipy.sparse.linalg import lsqr

from lsl.statistics import robust
from lsl.logger import LSL_LOGGER

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
__all__ = ['phase_only', 'delay_only', 'delay_and_phase', 'DelayStatus', 'check_delays']


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


class DelayStatus(enum.IntEnum):
    """
    Per-antenna results from `check_delays`.
    """
    
    OK        = 0    # delay is on the best cycle
    MOVED     = 1    # delay was shifted to a better cycle
    UNCERTAIN = 2    # delay could not be checked or the best cycle is not clear


def _build_cross_spectra(aa, dataSet, simSet, chan, pol, delays):
    """
    Build the per-antenna cross spectra between the data, corrected by
    `delays` in ns, and the model.  Each baseline/channel is reduced to a unit
    phasor.  Returns the cross spectra and the number of valid terms in each.
    """
    
    # Get the stand count
    nStands = len(aa.ants)
    
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    # Unit phasors of the delay corrected data relative to the model, using
    # the same sign convention as _scale_data
    bl = np.array(dataSet.baselines)
    blDelays = delays[bl[:,0]] - delays[bl[:,1]]
    obsVis = getattr(dataSet, pol).data[:,chan] * np.exp(2j*np.pi*fq*blDelays[:,None])
    simVis = getattr(simSet, pol).data[:,chan]
    
    with np.errstate(invalid='ignore', divide='ignore'):
        r = obsVis*simVis.conj()
        r /= np.abs(r)
    valid = np.isfinite(r)
    r[~valid] = 0.0
    
    # A delay shift on antenna k multiplies baseline (k,j) by exp(2 pi i f d)
    # and baseline (i,k) by exp(-2 pi i f d) so conjugate the second set
    S = np.zeros((nStands, fq.size), dtype=np.complex128)
    M = np.zeros((nStands, fq.size))
    for k in range(nStands):
        first, second = (bl[:,0] == k), (bl[:,1] == k)
        S[k,:] = r[first,:].sum(axis=0) + r[second,:].conj().sum(axis=0)
        M[k,:] = valid[first,:].sum(axis=0) + valid[second,:].sum(axis=0)
        
    return S, M


def check_delays(aa, dataSet, simSet, chan, pol, delays, search_range=None, min_coherence=0.5, min_improvement=0.3, ambiguity_ratio=0.9, return_details=False):
    """
    Check a set of delays in ns from `delay_only` for antennas that are off
    by one or more cycles.  Each antenna is tested on its own over trial
    shifts of +/- `search_range` ns (default is two cycles at the lowest
    frequency) using the full-band real coherence with the model.
    
    An antenna is flagged as DelayStatus.UNCERTAIN and left unchanged if its
    peak coherence is below `min_coherence` of the array median or if a
    second peak, at least half a cycle at the highest frequency away, is
    above `ambiguity_ratio` of its best peak.  Otherwise, it is shifted and
    flagged as DelayStatus.MOVED if the best shift improves the coherence by
    more than `min_improvement` of the array median.
    
    Returns the updated delays and a list of DelayStatus values.  If
    `return_details` is True a dictionary of the trial shifts, coherence
    curves, and thresholds is also returned.
    
    Since each antenna is checked independently, any moves should be
    followed by another solver pass, e.g.:
    
        >>> cal, delays = selfcal.delay_only(aa, data, sim, chan, pol)
        >>> delays, status = selfcal.check_delays(aa, data, sim, chan, pol, delays)
        >>> if selfcal.DelayStatus.MOVED in status:
        ...     data = vis.scale_data(data, np.ones(len(delays)), delays)
        ...     cal, resid = selfcal.delay_only(aa, data, sim, chan, pol)
        ...     delays += resid
        
    .. note::
        The check needs the same (u,v) cut as the solver since short baselines
        that are dominated by unmodeled diffuse emission produce false moves.
        
    .. versionadded:: 4.0.2
    """
    
    # Make sure that `chan' is an array by trying to find its length
    try:
        junk = len(chan)
    except TypeError:
        chan = [chan]
        
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    delays = np.array(delays, dtype=np.float64)
    S, M = _build_cross_spectra(aa, dataSet, simSet, chan, pol, delays)
    nTerms = M.sum(axis=1)
    
    # Coherence for all trial shifts.  A pure delay model has no free phase so
    # the real part only peaks on the correct cycle.
    if search_range is None:
        search_range = 2 / fq.min()
    search_step = 1 / (32*fq.max())
    trials = np.arange(-search_range, search_range+search_step/2, search_step)
    trials -= trials[np.argmin(np.abs(trials))]
    phasors = np.exp(2j*np.pi*fq[None,:]*trials[:,None])
    coh = np.real(np.dot(S, phasors.T)) / np.maximum(nTerms, 1)[:,None]
    
    best = np.argmax(coh, axis=1)
    cohBest = coh[np.arange(coh.shape[0]),best]
    cohZero = coh[:,np.argmin(np.abs(trials))]
    
    # Runner-up peak, excluding the main lobe of the best peak
    far = np.abs(trials[None,:] - trials[best][:,None]) >= 1 / (2*fq.max())
    cohNext = np.where(far, coh, -np.inf).max(axis=1)
    
    # Thresholds relative to the typical coherence of the array
    valid = nTerms > 0
    typical = np.median(cohBest[valid]) if valid.any() else 0.0
    cohCut = min_coherence*typical
    improvementCut = min_improvement*typical
    
    newDelays = delays.copy()
    status = [DelayStatus.OK for ant in aa.ants]
    for k,ant in enumerate(aa.ants):
        if not valid[k] or cohBest[k] < cohCut:
            status[k] = DelayStatus.UNCERTAIN
            LSL_LOGGER.info(f"Delay check cannot verify stand {ant.stand} (peak coherence {cohBest[k]:.3f} < {cohCut:.3f})")
        elif cohNext[k] > ambiguity_ratio*cohBest[k]:
            status[k] = DelayStatus.UNCERTAIN
            LSL_LOGGER.info(f"Delay check found an ambiguous delay for stand {ant.stand} (peaks {cohBest[k]:.3f} and {cohNext[k]:.3f})")
        elif cohBest[k] - cohZero[k] > improvementCut:
            newDelays[k] += trials[best[k]]
            status[k] = DelayStatus.MOVED
            LSL_LOGGER.info(f"Delay check moved stand {ant.stand} by {trials[best[k]]:+.2f} ns (coherence {cohZero[k]:.3f} -> {cohBest[k]:.3f})")
            
    if return_details:
        details = {'trials': trials, 'coherence': coh, 'n_terms': nTerms,
                   'best_shift': trials[best], 'coherence_zero': cohZero,
                   'coherence_best': cohBest, 'coherence_next': cohNext,
                   'typical': typical, 'coherence_cut': cohCut,
                   'improvement_cut': improvementCut}
        return newDelays, status, details
    return newDelays, status
