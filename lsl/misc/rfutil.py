# -*- coding: utf-8 -*-
"""
A collection of utilities to help convert RF engineering/communications terms
into radio astronomy terms.  This module is geared toward taking antenna 
parameters in dBi or dBd and converting them into meaningful RA gains in
K/Jy and then using those to get a system equivalent flux density.

.. versionadded:: 1.0.3
"""

import math

from lsl.common.constants import c as speedOfLight

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['dBd2dBi', 'dBi2gain', 'dBd2gain', 'calculateSEFD', 'jy2dBm', '__version__', '__revision__', '__all__']


def dBd2dBi(dBd):
	"""
	Convert an antenna gain in dBb (gain relative to the maximum gain of a
	half-wave dipole) into a gain in dBi (gain relative to an isotropic 
	antenna).
	"""
	
	return dBd + 2.15


def dBi2gain(dBi, frequency):
	"""
	Convert an antenna gain in dBi (gain relative to an isotropic antenna) 
	at a particular frequency in Hz into a gain in K/Jy.
	"""
	
	return 2.56 / (frequency/1e6)**2 * 10**(dBi/10.0)


def dBd2gain(dBd, frequency):
	"""
	Convert an antenna gain in dBd (gain relative to the maximum gain of a
	half-wave dipole) at a particular frequency in Hz into a gain in K/Jy
	"""
	
	dBi = dBd2dBi(dBd)
	return dBi2gain(dBi, frequency)


def calculateSEFD(Tsys, gain=None, area=None, efficiency=None):
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
	
	# Boltzmann's constant
	kB = 1380.0	# Jy m^2 / K
	
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


def jy2dBm(flux, frequency, bandwidth):
	"""
	Convert a flux density in Jy into a received signal strength in dBm 
	under the assumptions of:
	  * center frequency in Hz
	  * signal bandwidth in Hz
	  * received using a isotropic antenna with area wavelength^2 / (4 pi)
	"""
	
	# Antenna parameters
	wavelength = speedOfLight / frequency	# m
	areaIsotropic = wavelength**2 / 4 / math.pi	# m^2
	
	# Power in mW
	P = 10**-26 * bandwidth * areaIsotropic * 1000.0
	
	# To dBm
	return 10.0*math.log10(P)
