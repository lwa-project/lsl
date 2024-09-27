import os
import sys
import h5py
import ephem
import numpy as np
from functools import lru_cache
from collections import OrderedDict
from scipy.interpolate import interp1d, RegularGridInterpolator

from astropy.coordinates import Angle as AstroAngle

from lsl.common.data_access import DataAccess

from lsl.misc import telemetry
telemetry.track_module()

__version__ = '0.1'
__all__ = ['mueller_matrix', 'beam_response', 'get_avaliable_models']


_MODELS = OrderedDict()
_MODELS['empirical'] = 'antenna/lwa1-dipole-emp.npz'
_MODELS['nec_002']   = 'antenna/NEC_er-3p0_sig-0p002.hdf5'
_MODELS['nec_004']   = 'antenna/NEC_er-3p0_sig-0p004.hdf5'
_MODELS['nec_008']   = 'antenna/NEC_er-3p0_sig-0p008.hdf5'
_MODELS['nec_015']   = 'antenna/NEC_er-3p0_sig-0p015.hdf5'
_MODELS['nec_030']   = 'antenna/NEC_er-3p0_sig-0p030.hdf5'


@lru_cache(maxsize=10)
def _load_response_fitted(frequency, corrected=False):
    """
    Given an observing frequency in Hz, load in the empirical model (see LWA
    memos #175 and #178) for an isolated LWA dipole and compute its response at
    a collection of azimuths and altiudes.  The `corrected` keyword enables an 
    altitude-based correction to the model based on work by Dowell et al. (2017).
    Returns a seven-element list of:
     * frequencies represented (1D; Hz)
     * azimuths computed (1D; degrees)
     * altitudes computed (1D; degrees)
     * X pol. E plane response (3D; freq x alt x az; normalized)
     * X pol. H plane response (3D; freq x alt x az; normalized)
     * Y pol. E plane response (3D; freq x alt x az; normalized)
     * Y pol. H plane response (3D; freq x alt x az; normalized)
     
    .. note:: This empirical model is only valid for XX and YY and the returned
              values the E and H plane responses are sqrt(XX/2) or sqrt(YY/2).
    """
    
    for pol in ('X', 'Y'):
        # Get the empirical model of the beam and compute it for the correct frequencies
        with DataAccess.open(_MODELS['empirical'], 'rb') as fh:
            beamDict = np.load(fh)
            beamCoeff = beamDict[f"fit{pol.upper()}"]
            alphaE = np.polyval(beamCoeff[0,0,:], frequency)
            betaE =  np.polyval(beamCoeff[0,1,:], frequency)
            gammaE = np.polyval(beamCoeff[0,2,:], frequency)
            deltaE = np.polyval(beamCoeff[0,3,:], frequency)
            alphaH = np.polyval(beamCoeff[1,0,:], frequency)
            betaH =  np.polyval(beamCoeff[1,1,:], frequency)
            gammaH = np.polyval(beamCoeff[1,2,:], frequency)
            deltaH = np.polyval(beamCoeff[1,3,:], frequency)
            
        if corrected:
            with DataAccess.open('antenna/lwa1-dipole-cor.npz', 'rb') as fh:
                corrDict = np.load(fh)
                cFreqs = corrDict['freqs'][...]
                cAlts  = corrDict['alts'][...]
                if corrDict['degrees'].item():
                    cAlts *= np.pi / 180.0
                cCorrs = corrDict['corrs'][...]
                
            if frequency/1e6 < cFreqs.min() or frequency/1e6 > cFreqs.max():
                print("WARNING: Input frequency of %.3f MHz is out of range, skipping correction" % (frequency/1e6,))
                corrFnc = None
            else:
                fCors = cAlts*0.0
                for j in range(fCors.size):
                    ffnc = interp1d(cFreqs, cCorrs[:,j], bounds_error=False)
                    fCors[j] = ffnc(frequency/1e6)
                corrFnc = interp1d(cAlts, fCors, bounds_error=False)
                
        else:
            corrFnc = None
            
        def compute_beam_pattern(az, alt, corr=corrFnc):
            zaR = np.pi/2 - alt*np.pi / 180.0 
            azR = az*np.pi / 180.0
            
            c = 1.0
            if corrFnc is not None:
                c = corrFnc(alt*np.pi / 180.0)
                c = np.where(np.isfinite(c), c, 1.0)
                
            pE = (1-(2*zaR/np.pi)**alphaE)*np.cos(zaR)**betaE + gammaE*(2*zaR/np.pi)*np.cos(zaR)**deltaE
            pH = (1-(2*zaR/np.pi)**alphaH)*np.cos(zaR)**betaH + gammaH*(2*zaR/np.pi)*np.cos(zaR)**deltaH

            return c*np.sqrt((pE*np.cos(azR))**2 + (pH*np.sin(azR))**2)

        # Calculate the beam
        az = np.arange(361)
        alt = np.arange(91)
        maz, malt = np.meshgrid(az, alt)
        resp = np.sqrt(compute_beam_pattern(maz, malt)/2)
        resp.shape = (1,)+resp.shape
        
        try:
            output.extend([resp, resp])
        except NameError:
            output = [np.array([frequency]), az, alt, resp, resp]
            
    return output


@lru_cache(maxsize=10)
def _load_response_full(frequency, model='nec_004'):
    """
    Given an observing frequency in Hz, load in data from an electromagnetic
    model for an isolated LWA dipole as a function of frequency, azimuth, and
    altitude.  Returns a seven-element list of:
     * frequencies represented (1D; Hz)
     * azimuths computed (1D; degrees)
     * altitudes computed (1D; degrees)
     * X pol. E plane response (3D; freq x alt x az; normalized)
     * X pol. H plane response (3D; freq x alt x az; normalized)
     * Y pol. E plane response (3D; freq x alt x az; normalized)
     * Y pol. H plane response (3D; freq x alt x az; normalized)
    """
    
    model = model.lower()
    try:
        filename = _MODELS[model]
    except KeyError:
        raise ValueError(f"Unknown dipole model source '{model}'")
        
    for pol in ('X', 'Y'):
        # Open the file and select the mode relevant frequencies
        with DataAccess.open(filename, 'rb') as fh:
            h = h5py.File(fh)
            mfreq = h['freq'][...]*1e6
            best = np.where(np.abs(mfreq - frequency) < 10e6)[0]
            maz = h['phi'][...]
            malt = 90 - h['theta'][...] # zenith angle -> altitude
            mflds = h[f"{pol.upper()}-pol_Efields"]
            E = mflds['Etheta(Mag)'][best,:,:maz.size]*np.exp(1j*mflds['Etheta(Phase)'][best,:,:maz.size])
            H = mflds['Ephi(Mag)'][best,:,:maz.size]*np.exp(1j*mflds['Ephi(Phase)'][best,:,:maz.size])
            s = max([E.max(), H.max()])
            E /= s
            H /= s
            
        if maz.max() < 360:
            ## Catch to make sure we can wrap around the north
            maz = np.append(maz, 360)
            E = np.append(E, E[...,[0,]], axis=2)
            H = np.append(H, H[...,[0,]], axis=2)
            
        try:
            output.extend([E,H])
        except NameError:
            output = [mfreq[best], maz, malt, E, H]
            
    return output


def _jones_matrix(model, frequency=74e6):
    """
    Given a EM model name, load in the model and return a four element tuple of:
     * frequencyy (1D; Hz)
     * azimuth (1D; degrees)
     * altitude (1D; degrees)
     * Jones matrix (5D; 2 x 2 x freq x alt x az)
    """
    
    # Load
    model = model.lower()
    if model == 'empirical':
        mfreq, maz, malt, XE, XH, YE, YH = _load_response_fitted(frequency, corrected=False)
    elif model == 'llfss':
        mfreq, maz, malt, XE, XH, YE, YH = _load_response_fitted(frequency, corrected=True)
    else:
        mfreq, maz, malt, XE, XH, YE, YH = _load_response_full(frequency, model=model)
        
    # Create the Jones matrix
    J = np.array([[XE, XH],
                  [YE, YH]])
                  
    # Return
    return mfreq, maz, malt, J


def mueller_matrix(model, az, alt, frequency=74e6, degrees=True):
    """
    Given a model source, an array of azimuth values and an array of altitude
    values, compute the Mueller matrix for an isolated LWA dipole.  Returns the
    Mueller matrix with shape (4,4)+az.shape.
    
    .. note:: If az/alt are scalar then a (4,4) matrix is returned.
    """
    
    # Convert
    if isinstance(az, AstroAngle):
        az = az.deg
    elif isinstance(az, ephem.Angle) or not degrees:
        az = az * 180/np.pi
    
    if isinstance(alt, AstroAngle):
        alt = alt.deg
    elif isinstance(az, ephem.Angle) or not degrees:
        alt = alt * 180/np.pi
        
    is_1d = False
    if not isinstance(az, np.ndarray) or not isinstance(alt, np.ndarray):
        is_1d = True
        az, alt = np.array(az), np.array(alt)
        
    # Load in correct the Jones matrix
    mfreq, maz, malt, J = _jones_matrix(model, frequency=frequency)
    
    # Set the "A" matrix
    A = np.array([[1, 0,   0,  1],
                  [1, 0,   0, -1],
                  [0, 1,   1,  0],
                  [0, 1j, -1j, 0]])
                  
    # Compute the Kronecker product of J and J*
    M = np.zeros((4,4,J.shape[2],J.shape[3],J.shape[4]), dtype=J.dtype)
    M[:2,:2,...] = J[0,0,...]*J.conj()
    M[:2,2:,...] = J[0,1,...]*J.conj()
    M[2:,:2,...] = J[1,0,...]*J.conj()
    M[2:,2:,...] = J[1,1,...]*J.conj()
    
    # Swap the axis order to make the matrix multiplication easier
    M = M.transpose(2,3,4,0,1)
    
    # Compute the Mueller matrix
    M = np.matmul(A, np.matmul(M, np.linalg.inv(A)))
    
    # Another axis swap to put the polarization axes first
    M = M.transpose(3,4,0,1,2)
    
    # Evaluate at each azimuth/altitude using a 3D interpolation function
    ffreq = np.ones(az.size)*frequency
    faz = az.ravel()
    falt = alt.ravel()
    resp = np.zeros((4,4)+faz.shape, dtype=np.float64)
    for i in range(4):
        for j in range(4):
            P = M[i,j,...]
            if mfreq.size > 1:
                pfunc = RegularGridInterpolator((mfreq, malt, maz), P.real, bounds_error=False, fill_value=np.nan)
                resp[i,j,:] = pfunc((ffreq, falt, faz))
            else:
                pfunc = RegularGridInterpolator((malt, maz), P[0,:,:].real, bounds_error=False, fill_value=np.nan)
                resp[i,j,:] = pfunc((falt, faz))
    resp.shape = (4,4)+az.shape
    
    return resp


def _compute_from_feeds(pol, XE, XH, YE, YH):
    """
    Given a desired polarization product and a collection of feed response in
    the E and H planes, return the combined response for that product.  The
    returned array as the same dimensionallity as the input feed responses.
    """
    
    pol = pol.upper()
    if pol == 'XX':
        return np.abs(XE)**2 + np.abs(XH)**2
    elif pol == 'YY':
        return np.abs(YE)**2 + np.abs(YH)**2
    elif pol == 'XY':
        return XE*YE.conj() + XH*YH.conj()
    elif pol == 'YX':
        return YE*XE.conj() + YH*XH.conj()
    elif pol == 'I':
        return _compute_from_feeds('XX', XE, XH, YE, YH) \
               + _compute_from_feeds('YY', XE, XH, YE, YH)
    elif pol == 'Q':
        return _compute_from_feeds('XX', XE, XH, YE, YH) \
               - _compute_from_feeds('YY', XE, XH, YE, YH)
    elif pol == 'U':
        return _compute_from_feeds('XY', XE, XH, YE, YH) \
               + _compute_from_feeds('YX', XE, XH, YE, YH)
    elif pol == 'V':
        return _compute_from_feeds('XY', XE, XH, YE, YH) \
               - _compute_from_feeds('YX', XE, XH, YE, YH)
    else:
        raise ValueError(f"Unknown polarization produce '{pol}'")


def beam_response(model, pol, az, alt, frequency=74e6, degrees=True):
    """
    Return the response of an isolated LWA dipole in a particular azimuth/
    altitude, for the specified module source, polarization product, and
    frequency.
    """
    
    # Convert
    if isinstance(az, AstroAngle):
        az = az.deg
    elif isinstance(az, ephem.Angle) or not degrees:
        az = az * 180/np.pi
    
    if isinstance(alt, AstroAngle):
        alt = alt.deg
    elif isinstance(az, ephem.Angle) or not degrees:
        alt = alt * 180/np.pi
        
    is_1d = False
    if not isinstance(az, np.ndarray) or not isinstance(alt, np.ndarray):
        is_1d = True
        az, alt = np.array(az), np.array(alt)
        
    # Load
    model = model.lower()
    if model == 'empirical':
        mfreq, maz, malt, XE, XH, YE, YH = _load_response_fitted(frequency, corrected=False)
    elif model == 'llfss':
        mfreq, maz, malt, XE, XH, YE, YH = _load_response_fitted(frequency, corrected=True)
    else:
        mfreq, maz, malt, XE, XH, YE, YH = _load_response_full(frequency, model=model)
        
    # Build up the correct polarization product
    P = _compute_from_feeds(pol, XE, XH, YE, YH)
    if mfreq.size > 1:
        is_3d = True
        pfunc = RegularGridInterpolator((mfreq, malt, maz), P.real, bounds_error=False, fill_value=np.nan)
    else:
        is_3d = False
        pfunc = RegularGridInterpolator((malt, maz), P[0,:,:].real, bounds_error=False, fill_value=np.nan)
    
    # Evaluate
    ffreq = np.ones(az.size)*frequency
    faz = az.ravel()
    falt = alt.ravel()
    if is_3d:
        resp = pfunc((ffreq, falt, faz))
    else:
        resp = pfunc((falt, faz))
    resp.shape = az.shape
    
    if is_1d:
        resp = resp.item()
        
    return resp


def get_avaliable_models():
    return list(_MODELS.keys()) + ['llfss']
