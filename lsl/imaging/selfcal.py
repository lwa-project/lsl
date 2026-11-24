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
   Added a post-convergence delay ambiguity check to `delay_only`

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
__all__ = ['DelayStatus', 'phase_only', 'delay_only', 'delay_and_phase', 'check_delays']


class DelayStatus(enum.IntFlag):
    """
    Per-antenna flags reported by the delay ambiguity check.
    """
    
    OK           = 0
    MOVED        = 1    # delay was shifted by the check
    UNVERIFIABLE = 2    # too little coherence with the model to check
    UNRESOLVED   = 4    # still wanted to move after the last check


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


def _get_reference_index(aa, ref_ant):
    """
    Convert a reference stand number into an index in aa.ants.  A value of
    zero selects the first antenna.
    """
    
    if ref_ant == 0:
        return 0
    for i,ant in enumerate(aa.ants):
        if ant.stand == ref_ant:
            return i
    raise RuntimeError(f"Stand #{ref_ant} not found in the array provided")


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


def check_delays(aa, dataSet, simSet, chan, pol, delays, search_range=None, search_step=0.25, min_gain=0.3, min_coherence=0.5, return_details=False):
    """
    Check a set of converged delays in ns from `delay_only` for antennas that
    are off by one or more cycles.  Each antenna is tested on its own over
    trial shifts of +/- `search_range` ns (default is two cycles at the lowest
    frequency) using the full-band real coherence with the model.
    
    The thresholds are relative to the median peak coherence of the array:
     * antennas whose peak is below `min_coherence` of this are flagged as
       DelayStatus.UNVERIFIABLE and left unchanged;
     * antennas whose best shift improves the coherence by more than
       `min_gain` of this are shifted and flagged as DelayStatus.MOVED.
    
    Returns the updated delays and a list of DelayStatus flags.  If
    `return_details` is True a dictionary of the trial shifts, coherence
    curves, and thresholds is also returned.
    
    .. note::
        All antennas are tested against the same input delays so any shifts
        are independent and should be followed by another solver pass.
        
    .. note::
        The check needs the same (u,v) cut as the solver since short baselines
        that are dominated by unmodeled diffuse emission produce false moves.
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
    trials = np.arange(-search_range, search_range+search_step/2, search_step)
    phasors = np.exp(2j*np.pi*fq[None,:]*trials[:,None])
    coh = np.real(np.dot(S, phasors.T)) / np.maximum(nTerms, 1)[:,None]
    
    best = np.argmax(coh, axis=1)
    cohBest = coh[np.arange(coh.shape[0]),best]
    cohZero = coh[:,np.argmin(np.abs(trials))]
    
    # Thresholds relative to the typical coherence of the array
    valid = nTerms > 0
    typical = np.median(cohBest[valid]) if valid.any() else 0.0
    gainCut = min_gain*typical
    cohCut = min_coherence*typical
    
    newDelays = delays.copy()
    status = [DelayStatus.OK for ant in aa.ants]
    for k,ant in enumerate(aa.ants):
        if not valid[k] or cohBest[k] < cohCut:
            status[k] = DelayStatus.UNVERIFIABLE
            LSL_LOGGER.info(f"Delay check cannot verify stand {ant.stand} (peak coherence {cohBest[k]:.3f} < {cohCut:.3f})")
        elif cohBest[k] - cohZero[k] > gainCut:
            newDelays[k] += trials[best[k]]
            status[k] = DelayStatus.MOVED
            LSL_LOGGER.info(f"Delay check moved stand {ant.stand} by {trials[best[k]]:+.2f} ns (coherence {cohZero[k]:.3f} -> {cohBest[k]:.3f})")
            
    if return_details:
        details = {'trials': trials, 'coherence': coh, 'n_terms': nTerms,
                   'best_shift': trials[best], 'coherence_zero': cohZero,
                   'coherence_best': cohBest, 'typical': typical,
                   'gain_cut': gainCut, 'coherence_cut': cohCut}
        return newDelays, status, details
    return newDelays, status


def _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, phase_only=False, delay_only=False, delay_and_phase=False, amplitude_cutoff=1.001, phase_cutoff=0.01, delay_cutoff=0.2, inv_epsilon=0.0, initial_delays=None):
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
    
    The delay solutions start from `initial_delays` in ns if it is provided,
    otherwise from zero.
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
    ref_ant = _get_reference_index(aa, ref_ant)
    LSL_LOGGER.info(f"Using antenna #{ref_ant} as a reference (Stand #{aa.ants[ref_ant].stand})")
        
    # Frequency in GHz so that the delays can be in ns
    fq = dataSet.freq[chan] / 1e9
    
    origSet = dataSet
    
    converged = False
    tempGains = np.ones(N)
    tempDelays = np.zeros(N)
    tempPhaseOffsets = np.zeros(N)
    if initial_delays is not None:
        tempDelays = np.array(initial_delays, dtype=np.float64)
        dataSet = _scale_data(origSet, tempGains, tempDelays, tempPhaseOffsets)
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


def delay_only(aa, dataSet, simSet, chan, pol, ref_ant=0, loop_gain=0.75, max_iter=30, amplitude=False, amplitude_cutoff=1.001, delay_cutoff=0.2, inv_epsilon=0.0, return_convergence=False, check_ambiguity=False, search_range=None, search_step=0.25, min_gain=0.3, min_coherence=0.5, max_checks=5, return_status=False):
    """
    Function to apply a delay-only (and, optionally, a amplitude) self-
    calibration to data stored in a readUVData dictionary and a model sky 
    stored in a lsl.sim.vis.buildSimSky dictionary for the given 
    polarization and channel(s).
    
    If `check_ambiguity` is True the converged delays are checked for
    antennas that are off by one or more cycles with `check_delays` and the
    solver is re-run from the corrected delays, up to `max_checks` times.
    The delays are then re-referenced to the reference antenna.  See
    `check_delays` for the `search_range`, `search_step`, `min_gain`, and
    `min_coherence` keywords.
    
    .. note::
        If the "amplitude" keyword is set to True, a three-element tuple of 
        self-cal'd data, gain coefficients, and delays in ns is returned rather 
        than the standard two-element tuple.
        
    The delays returned are positive for an antenna whose signal arrives 
    late, and can be passed to lsl.sim.vis.scale_data to correct the data.
    
    .. note::
        If the "return_status" keyword is set to True a list of per-antenna
        DelayStatus flags from the ambiguity check is appended to the output.
        
    .. versionchanged:: 4.0.2
        Added the `check_ambiguity` and `return_status` keywords
    """
    
    caldDict, gains, delays, phaseOffsets, converged = \
      _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, 
                loop_gain=loop_gain, max_iter=max_iter, 
                amplitude=amplitude, delay_only=True, 
                amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff,
                inv_epsilon=inv_epsilon)
                
    status = [DelayStatus.OK for ant in aa.ants]
    if check_ambiguity:
        refIndex = _get_reference_index(aa, ref_ant)
        moved = np.zeros(len(aa.ants), dtype=bool)
        for c in range(max_checks+1):
            newDelays, status = check_delays(aa, dataSet, simSet, chan, pol, delays,
                                             search_range=search_range, search_step=search_step,
                                             min_gain=min_gain, min_coherence=min_coherence)
            pending = np.array([bool(s & DelayStatus.MOVED) for s in status])
            if not pending.any() or c == max_checks:
                break
                
            # Re-reference and re-solve starting from the shifted delays
            moved |= pending
            newDelays -= newDelays[refIndex]
            caldDict, gains, delays, phaseOffsets, converged = \
              _self_cal(aa, dataSet, simSet, chan, pol, ref_ant=ref_ant, 
                        loop_gain=loop_gain, max_iter=max_iter, 
                        amplitude=amplitude, delay_only=True, 
                        amplitude_cutoff=amplitude_cutoff, delay_cutoff=delay_cutoff,
                        inv_epsilon=inv_epsilon, initial_delays=newDelays)
                        
        # Shifts from the last check have not been applied
        if pending.any():
            LSL_LOGGER.warning(f"Delay check still moving {pending.sum()} antennas after {max_checks} passes")
        for k in range(len(aa.ants)):
            status[k] &= DelayStatus.UNVERIFIABLE
            if moved[k]:
                status[k] |= DelayStatus.MOVED
            if pending[k]:
                status[k] |= DelayStatus.UNRESOLVED
                
    if amplitude:
        output = (caldDict, gains, delays)
    else:
        output = (caldDict, delays)
    if return_convergence:
        output = output+(converged,)
    if return_status:
        output = output+(status,)
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
    
    .. note::
        The delay ambiguity check in `delay_only` is not available here.
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
