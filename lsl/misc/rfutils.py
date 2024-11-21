"""
A collection of utilities to help convert RF engineering/communications terms
into radio astronomy terms.  This module is geared toward taking antenna 
parameters in dBi or dBd and converting them into meaningful RA gains in
K/Jy and then using those to get a system equivalent flux density.

.. versionadded:: 1.0.3
"""

import math
from astropy import units
from astropy.constants import c as speedOfLight, k_B as kB

from lsl.misc import telemetry
telemetry.track_module()

from typing import Optional, Union


__version__ = '0.1'
__all__ = ['dBd_to_dBi', 'dBd_to_dBi', 'dBi_to_gain', 'dBd_to_gain', 'gain_to_dBi', 'gain_to_dBd', 'calculate_sefd', 
        'calculate_effective_area', 'Jy_to_dBm', 'dBm_to_Jy']


def _make_quantity(value: Union[str,float,units.quantity.Quantity], default_units: units.Unit) -> Union[units.quantity.Quantity,None]:
    """
    Helper function for taking values and giving them units.
    """
    
    if value is not None:
        if isinstance(value, str):
            value = units.quantity.Quantity(value)
        elif not isinstance(value, units.quantity.Quantity):
            value = value*default_units
    return value


def dBd_to_dBi(dBd: float) -> float:
    """
    Convert an antenna gain in dBb (gain relative to the maximum gain of a
    half-wave dipole) into a gain in dBi (gain relative to an isotropic 
    antenna).
    """
    
    return dBd + 2.14844


def dBi_to_dBd(dBi: float) -> float:
    """
    Convert and antenna gain dBi (gain relative to an isotropic antenna)
    into a gain in dBd (gain relative to the maximum gain of a half-wave 
    dipole).
    """
    
    return dBi - 2.14844


def dBi_to_gain(dBi: float, frequency: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Convert an antenna gain in dBi (gain relative to an isotropic antenna) 
    at a particular frequency in Hz into a gain in K/Jy.
    """
    
    # Deal with units
    frequency = _make_quantity(frequency, units.Hz)
     
    # 10^(dBi/10) is the ratio of the effective area of the antenna relative 
    # to an isotropic antenna where the area of the isotropic antenna is 
    # given by:  wavelength^2 / 4 / pi
    gain = speedOfLight**2 / (8.0*math.pi*kB*frequency**2) * 10**(dBi/10.0) # type: ignore
    gain = gain.to('K/Jy').value
    
    return gain


def dBd_to_gain(dBd: float, frequency: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Convert an antenna gain in dBd (gain relative to the maximum gain of a
    half-wave dipole) at a particular frequency in Hz into a gain in K/Jy
    """
    
    dBi = dBd_to_dBi(dBd)
    return dBi_to_gain(dBi, frequency)


def gain_to_dBi(gain:  Union[str,float,units.quantity.Quantity], frequency: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Convert an antenna gain in K/Jy at a particular frequency in Hz into a
    gain in dBi (gain relative to an isotropic antenna).
    """
    
    # Deal with units
    gain = _make_quantity(gain, units.K/units.Jy)
    frequency = _make_quantity(frequency, units.Hz)
        
    # Calculate the area ratio
    areaRatio = 8.0*math.pi*kB*frequency**2 / speedOfLight**2 * gain    # type: ignore
    
    # Convert to dBi
    return 10.0*math.log10(areaRatio)


def gain_to_dBd(gain:  Union[str,float,units.quantity.Quantity], frequency: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Convert an antenna gain in K/Jy at a particular frequency in Hz into a
    gain in dBd (gain relative to the maximum gain of a half-wave dipole).
    """
    
    dBi = gain_to_dBi(gain, frequency)
    return dBi_to_dBd(dBi)


def calculate_sefd(Tsys: Union[str,float,units.quantity.Quantity], gain: Optional[Union[str,float,units.quantity.Quantity]]=None, area: Optional[Union[str,float,units.quantity.Quantity]]=None, efficiency: Optional[float]=None) -> float:
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
    
    # Deal with units
    Tsysq = _make_quantity(Tsys, units.K)
    gainq = _make_quantity(gain, units.K/units.Jy)
    areaq = _make_quantity(area, units.meter**2)
    
    # Have we been supplied a gain or do we need to calculate it?
    if gainq is None:
        ## Do we have everything we need to calculate the gain?
        if area is not None and efficiency is not None:
            ## Looks like it
            gainq = efficiency*areaq / 2.0 / kB # type: ignore
        else:
            ## Can't
            raise RuntimeError("Too few parameters supplied to calculate the SEFD")
            
    # Calculate the SEFD
    SEFD = Tsysq / gainq
    SEFD = SEFD.to('Jy').value
    
    return SEFD


def calculate_effective_area(gain: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Given the gain of an antenna in K/Jy, calculate the effective collecting
    area in square meters.
    """
    
    # Deal with units
    gain = _make_quantity(gain, units.K/units.Jy)
    
    # Calculate the area
    area = 2.0*kB*gain
    area = area.to('m^2').value
    
    return area


def Jy_to_dBm(flux: Union[str,float,units.quantity.Quantity], bandwidth: Union[str,float,units.quantity.Quantity], gain: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Convert a flux density in Jy into a received signal strength in dBm 
    under the assumptions of:
     * signal bandwidth in Hz
     * antenna gain in K/Jy
    """
    
    # Deak with units
    flux = _make_quantity(flux, units.Jy)
    bandwidth = _make_quantity(bandwidth, units.Hz)
    gain = _make_quantity(gain, units.K/units.Jy)
        
    # Antenna parameters
    area = calculate_effective_area(gain)
    area = _make_quantity(area, units.meter**2) # type: ignore
    
    # Power in mW
    P = flux*bandwidth*area # type: ignore
    P = P.to('mW').value
    
    # To dBm
    return 10.0*math.log10(P)


def dBm_to_Jy(dBm: float, bandwidth: Union[str,float,units.quantity.Quantity], gain: Union[str,float,units.quantity.Quantity]) -> float:
    """
    Convert a received signal strength in dBm into a flux density in Jy under
    the assumptions of:
     * signal bandwidth in Hz
     * antenna gain in K/Jy
    """
    
    # Deak with units
    bandwidth = _make_quantity(bandwidth, units.Hz)
    gain = _make_quantity(gain, units.K/units.Jy)
    
    # Antenna parameters
    area = calculate_effective_area(gain)
    area = _make_quantity(area, units.meter**2) # type: ignore
    
    # Power in mW
    P = 10**(dBm/10.0)*(units.W/1000)
    
    # To Jy
    flux = P / (bandwidth*area)     # type: ignore
    flux = flux.to('Jy').value
    
    return flux
