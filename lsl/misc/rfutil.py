# -*- coding: utf-8 -*-

"""
A collection of utilities to help convert RF engineering/communications terms
into radio astronomy terms.  This module is geared toward taking antenna 
parameters in dBi or dBd and converting them into meaningful RA gains in
K/Jy and then using those to get a system equivalent flux density.

.. versionadded:: 1.0.3
"""

import math

from lsl.common.constants import c as speedOfLight, kB

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['dBd_to_dBi', 'dBd_to_dBi', 'dBi_to_gain', 'dBd_to_gain', 'gain_to_dBi', 'gain_to_dBd', 'calculate_sefd', 
        'calculate_effective_area', 'Jy_to_dBm', 'dBm_to_Jy']


def dBd_to_dBi(dBd):
    """
    Convert an antenna gain in dBb (gain relative to the maximum gain of a
    half-wave dipole) into a gain in dBi (gain relative to an isotropic 
    antenna).
    """
    
    return dBd + 2.14844


def dBi_to_dBd(dBi):
    """
    Convert and antenna gain dBi (gain relative to an isotropic antenna)
    into a gain in dBd (gain relative to the maximum gain of a half-wave 
    dipole).
    """
    
    return dBi - 2.14844


def dBi_to_gain(dBi, frequency):
    """
    Convert an antenna gain in dBi (gain relative to an isotropic antenna) 
    at a particular frequency in Hz into a gain in K/Jy.
    """
    
    # 10^(dBi/10) is the ratio of the effective area of the antenna relative 
    # to an isotropic antenna where the area of the isotropic antenna is 
    # given by:  wavelength^2 / 4 / pi
    return speedOfLight**2 / (8.0*math.pi*kB*frequency**2) * 10**(dBi/10.0)


def dBd_to_gain(dBd, frequency):
    """
    Convert an antenna gain in dBd (gain relative to the maximum gain of a
    half-wave dipole) at a particular frequency in Hz into a gain in K/Jy
    """
    
    dBi = dBd_to_dBi(dBd)
    return dBi_to_gain(dBi, frequency)


def gain_to_dBi(gain, frequency):
    """
    Convert an antenna gain in K/Jy at a particular frequency in Hz into a
    gain in dBi (gain relative to an isotropic antenna).
    """
    
    # Calculate the area ratio
    areaRatio = 8.0*math.pi*kB*frequency**2 / speedOfLight**2 * gain
    
    # Convert to dBi
    return 10.0*math.log10(areaRatio)


def gain_to_dBd(gain, frequency):
    """
    Convert an antenna gain in K/Jy at a particular frequency in Hz into a
    gain in dBd (gain relative to the maximum gain of a half-wave dipole).
    """
    
    dBi = gain_to_dBi(gain, frequency)
    return dBi_to_dBd(dBi)


def calculate_sefd(Tsys, gain=None, area=None, efficiency=None):
    """
    Given a variety of parameters, calculate the system equivalent flux 
    density in Jy for an antenna.  The parameters are:
     * Tsys - system temperature in K - required
     * gain - antenna gain in K/Jy - optional
     * area - antenna collecting area in m^2 - optional
     * efficiency - aperture efficiency - optional
    
    Of the optional parameters either 'gain' needs to be supplied or
    both 'area' and 'efficiency'.
    """
    
    # Have we been supplied a gain or do we need to calculate it?
    if gain is None:
        ## Do we have everything we need to calculate the gain?
        if area is None or efficiency is None:
            raise RuntimeError("Too few parameters supplied to calculate the SEFD")
            
        ## Looks like it
        gain = efficiency*area / 2.0 / kB
        
    # Calculate the SEFD
    SEFD = Tsys / gain
    
    return SEFD


def calculate_effective_area(gain):
    """
    Given the gain of an antenna in K/Jy, calculate the effective collecting
    area in square meters.
    """
    
    # Calculate the area
    area = 2.0*kB*gain
    
    return area


def Jy_to_dBm(flux, bandwidth, gain):
    """
    Convert a flux density in Jy into a received signal strength in dBm 
    under the assumptions of:
     * signal bandwidth in Hz
     * antenna gain in K/Jy
    """
    
    # Antenna parameters
    area = calculate_effective_area(gain)
    
    # Power in mW
    P =  flux * 10**-26 * bandwidth * area * 1000.0
    
    # To dBm
    return 10.0*math.log10(P)


def dBm_to_Jy(dBm, bandwidth, gain):
    """
    Convert a received signal strength in dBm into a flux density in Jy under
    the assumptions of:
     * signal bandwidth in Hz
     * antenna gain in K/Jy
    """
    
    # Antenna parameters
    area = calculate_effective_area(gain)
    
    # Power in mW
    P = 10**(dBm/10.0) / 1000.0
    
    # To Jy
    return P / (10**-26 * bandwidth * area)
    
