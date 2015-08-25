# -*- coding: utf-8 -*-

"""
Python module to handle the channelization and cross-correlation of TBW and
TBN data.  The main python functions in this module are:
  * calcSpectra - calculate power spectra for a collection of signals
  * FXCorrelator - calculate cross power spectra for a collection of signals
  * FXStokes - calculate Stokes cross power spectra for a collection of signals
both of which have been deprecated in favor of the new C extension based 
routines listed below.

The main python/C extension functions in this module are:
  * SpecMaster - similar to calcSpectra but uses the _spec module for all 
    computations and does not support automatic sub-integration
  * SpecMasterP - SpecMaster with a 64-tap uniform DFT filter bank
  * StokesMaster - similar to SpecMaster but computes all four Stokes parameters
  * FXMaster - calculate cross power spectra for a collection of signals

Each function is set up to process the signals in parallel using the 
multiprocessing module and accepts a variety of options controlling the processing
of the data, including various window functions and time averaging.

.. versionchanged:: 1.0.1
	Removed SpecMasterP.

.. versionchanged:: 1.0.0
	All of the functions here now return all 'LFFT' channels.
"""

import os
import sys
import ephem
import numpy

from lsl.common.constants import c as vLight
from lsl.common import dp as dp_common
from lsl.common.constants import *
from lsl.correlator import uvUtils

import _spec
import _stokes
import _core

__version__ = '1.0'
__revision__ = '$Rev$'
__all__ = ['pol2pol', 'noWindow', 'SpecMaster', 'StokesMaster', 'FXMaster', 'FXStokes', '__version__', '__revision__', '__all__']


def pol2pol(pol):
	"""
	Convert a polarization string, e.g., XX/XY or RR/LL, to a numeric :class:`lsl.common.stations.Antenna`
	instance polarization.
	"""
	
	pol = pol.upper()
	out = []
	for p in pol:
		if p in ('X', 'R'):
			out.append(0)
		elif p in ('Y', 'L'):
			out.append(1)
		else:
			raise RuntimeError("Unknown polarization code '%s'" % pol)
		
	return out


def noWindow(L):
	"""
	Default "empty" windowing function for use with the various routines.  This
	function returned a numpy array of '1's of the specified length.
	"""

	return numpy.ones(L)


def SpecMaster(signals, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0, ClipLevel=0):
	"""
	A more advanced version of calcSpectra that uses the _spec C extension 
	to handle all of the P.S.D. calculations in parallel.  Returns a two-
	element tuple of the frequencies (in Hz) and PSDs in linear power/RBW.
	
	.. note::
		SpecMaster currently average all data given and does not support the
		SampleAverage keyword that calcSpectra does.
	"""
	
	# Figure out if we are working with complex (I/Q) data or only real.  This
	# will determine how the FFTs are done since the real data mirrors the pos-
	# itive and negative Fourier frequencies.
	if signals.dtype.kind == 'c':
		lFactor = 1
		doFFTShift = True
		CentralFreq = float(CentralFreq)
	else:
		lFactor = 2
		doFFTShift = False

	# Calculate the frequencies of the FFTs.  We do this for twice the FFT length
	# because the real-valued signal only occupies the positive part of the 
	# frequency space.
	if SampleRate is None:
		SampleRate = dp_common.fS
	freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/SampleRate)
	# Deal with TBW and TBN data in the correct way
	if doFFTShift:
		freq += CentralFreq
		freq = numpy.fft.fftshift(freq)
	freq = freq[:LFFT]
		
	if window is noWindow:
		# Data without a window function provided
		if signals.dtype.kind == 'c':
			output = _spec.FPSDC2(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel)
		else:
			output = _spec.FPSDR2(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			output = _spec.FPSDC3(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel, window=window)
		else:
			output = _spec.FPSDR3(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel, window=window)
	
	return (freq, output)


def StokesMaster(signals, antennas, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0, ClipLevel=0):
	"""
	Similar to SpecMaster, but accepts an array of signals and a list of 
	antennas in order to compute the PSDs for the four Stokes parameters: 
	I, Q, U, and V.  Returns a two-element tuple of the frequencies (in Hz) 
	and PSDs in linear power/RBW.  The PSD are three dimensional with 
	dimensions Stokes parameter (0=I, 1=Q, 2=U, 3=V) by stand by channel).
	"""
	
	# Figure out if we are working with complex (I/Q) data or only real.  This
	# will determine how the FFTs are done since the real data mirrors the pos-
	# itive and negative Fourier frequencies.
	if signals.dtype.kind == 'c':
		lFactor = 1
		doFFTShift = True
		CentralFreq = float(CentralFreq)
	else:
		lFactor = 2
		doFFTShift = False
		
	# Match the X and Y stand data
	antennas1 = [a for a in antennas if a.pol == 0]
	signalsIndex1 = [i for (i, a) in enumerate(antennas) if a.pol == 0]
	antennas2 = [a for a in antennas if a.pol == 1]
	signalsIndex2 = [i for (i, a) in enumerate(antennas) if a.pol == 1]
	if len(signalsIndex1) != len(signalsIndex2):
		raise RuntimeError("Supplied data does not contain an equal number of X and Y signals.")

	# Calculate the frequencies of the FFTs.  We do this for twice the FFT length
	# because the real-valued signal only occupies the positive part of the 
	# frequency space.
	if SampleRate is None:
		SampleRate = dp_common.fS
	freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/SampleRate)
	# Deal with TBW and TBN data in the correct way
	if doFFTShift:
		freq += CentralFreq
		freq = numpy.fft.fftshift(freq)
	freq = freq[:LFFT]
	
	if window is noWindow:
		# Data without a window function provided
		if signals.dtype.kind == 'c':
			output = _stokes.FPSDC2(signals[signalsIndex1], signals[signalsIndex2], LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel)
		else:
			output = _stokes.FPSDR2(signals[signalsIndex1], signals[signalsIndex2], LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			output = _stokes.FPSDC3(signals[signalsIndex1], signals[signalsIndex2], LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel, window=window)
		else:
			output = _stokes.FPSDR3(signals[signalsIndex1], signals[signalsIndex2], LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel, window=window)
	
	return (freq, output)


def FXMaster(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, verbose=False, window=noWindow, SampleRate=None, CentralFreq=0.0, Pol='XX', GainCorrect=False, ReturnBaselines=False, ClipLevel=0, phaseCenter='z'):
	"""
	A more advanced version of FXCorrelator for TBW and TBN data.  Given an 
	2-D array of signals (stands, time-series) and an array of stands, compute 
	the cross-correlation of the data for all baselines.  Return the frequencies 
	and visibilities as a two-elements tuple.
	
	.. versionchanged:: 1.1.0
		Made the 'phaseCenter' keyword more flexible.  It can now be either:
		 * 'z' to denote the zenith,
		 * a ephem.Body instances which has been computed for the observer, or
		 * a two-element tuple of azimuth, elevation in degrees.
		
	.. versionchanged:: 1.0.0
		Added a phase-center keyword that accept a two-element tuple of 
		azimuth and elelvation (in degrees) to change where the 
		correlations are phased to
		
	.. versionchanged:: 0.4.0
		Switched over to passing in Antenna instances generated by the
		:mod:`lsl.common.stations` module instead of a list of stand ID
		numbers.
	"""
	
	# Decode the polarization product into something that we can use to figure 
	# out which antennas to use for the cross-correlation
	pol1, pol2 = pol2pol(Pol)
	
	antennas1 = [a for a in antennas if a.pol == pol1]
	signalsIndex1 = [i for (i, a) in enumerate(antennas) if a.pol == pol1]
	antennas2 = [a for a in antennas if a.pol == pol2]
	signalsIndex2 = [i for (i, a) in enumerate(antennas) if a.pol == pol2]

	nStands = len(antennas1)
	baselines = uvUtils.getBaselines(antennas1, antennas2=antennas2, IncludeAuto=IncludeAuto, Indicies=True)
	nBL = len(baselines)

	# Figure out if we are working with complex (I/Q) data or only real.  This
	# will determine how the FFTs are done since the real data mirrors the pos-
	# itive and negative Fourier frequencies.
	if signals.dtype.kind == 'c':
		lFactor = 1
		doFFTShift = True
		CentralFreq = float(CentralFreq)
	else:
		lFactor = 2
		doFFTShift = False
		
	if SampleRate is None:
		SampleRate = dp_common.fS
	freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/SampleRate)
	if doFFTShift:
		freq += CentralFreq
		freq = numpy.fft.fftshift(freq)
	freq = freq[:LFFT]
	
	# Get the location of the phase center in radians and create a 
	# pointing vector
	if phaseCenter == 'z':
		azPC = 0.0
		elPC = numpy.pi/2.0
	else:
		if isinstance(phaseCenter, ephem.Body):
			azPC = phaseCenter.az * 1.0
			elPC = phaseCenter.alt * 1.0
		else:
			azPC = phaseCenter[0]*numpy.pi/180.0
			elPC = phaseCenter[1]*numpy.pi/180.0
			
	source = numpy.array([numpy.cos(elPC)*numpy.sin(azPC), 
					  numpy.cos(elPC)*numpy.cos(azPC), 
					  numpy.sin(elPC)])
					  
	# Define the cable/signal delay caches to help correlate along and compute 
	# the delays that we need to apply to align the signals
	dlyRef = len(freq)/2
	delays1 = numpy.zeros((nStands,LFFT))
	delays2 = numpy.zeros((nStands,LFFT))
	for i in list(range(nStands)):
		xyz1 = numpy.array([antennas1[i].stand.x, antennas1[i].stand.y, antennas1[i].stand.z])
		xyz2 = numpy.array([antennas2[i].stand.x, antennas2[i].stand.y, antennas2[i].stand.z])
		
		delays1[i,:] = antennas1[i].cable.delay(freq) - numpy.dot(source, xyz1) / vLight
		delays2[i,:] = antennas2[i].cable.delay(freq) - numpy.dot(source, xyz2) / vLight
	if delays1[:,dlyRef].min() < delays2[:,dlyRef].min():
		minDelay = delays1[:,dlyRef].min()
	else:
		minDelay = delays2[:,dlyRef].min()
	delays1 -= minDelay
	delays2 -= minDelay
	
	# F - defaults to running parallel in C via OpenMP
	if window is noWindow:
		# Data without a window function provided
		if signals.dtype.kind == 'c':
			FEngine = _core.FEngineC2
		else:
			FEngine = _core.FEngineR2
		signalsF1, validF1 = FEngine(signals[signalsIndex1,:], freq, delays1, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			FEngine = _core.FEngineC3
		else:
			FEngine = _core.FEngineR3
		signalsF1, validF1 = FEngine(signals[signalsIndex1,:], freq, delays1, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel, window=window)
	
	if pol2 == pol1:
		signalsF2 = signalsF1
		validF2 = validF1
	else:
		if window is noWindow:
			signalsF2, validF2 = FEngine(signals[signalsIndex2,:], freq, delays2, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel)
		else:
			signalsF2, validF2 = FEngine(signals[signalsIndex2,:], freq, delays2, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel, window=window)

	# X
	output = _core.XEngine2(signalsF1, signalsF2, validF1, validF2)
	if not IncludeAuto:
		# Remove auto-correlations from the output of the X engine if we don't 
		# need them.  To do this we need to first build the full list of baselines
		# (including auto-correlations) and then prune that.
		baselinesFull = uvUtils.getBaselines(antennas1, antennas2=antennas2, IncludeAuto=True, Indicies=True)
		fom = numpy.array([a1-a2 for (a1,a2) in baselinesFull])
		nonAuto = numpy.where( fom != 0 )[0]
		output = output[nonAuto,:]

	# Divide the cross-multiplied data by the number of channels used
	output /= LFFT
	if signals.dtype.kind != 'c':
		output /= 2.0
		
	# Apply cable gain corrections (if needed)
	if GainCorrect:
		for bl in xrange(output.shape[0]):
			cableGain1 = antennas1[baselines[bl][0]].cable.gain(freq)
			cableGain2 = antennas2[baselines[bl][1]].cable.gain(freq)
			
			output[bl,:] /= numpy.sqrt(cableGain1*cableGain2)
			
	# Create antenna baseline list (if needed)
	if ReturnBaselines:
		antennaBaselines = []
		for bl in xrange(output.shape[0]):
			antennaBaselines.append( (antennas1[baselines[bl][0]], antennas2[baselines[bl][1]]) )
		returnValues = (antennaBaselines, freq, output)
	else:
		returnValues = (freq, output)

	return returnValues


def FXStokes(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, verbose=False, window=noWindow, SampleRate=None, CentralFreq=0.0,  GainCorrect=False, ReturnBaselines=False, ClipLevel=0, phaseCenter='z'):
	"""
	A more advanced version of FXCorrelator for TBW and TBN data.  Given an 
	2-D array of signals (stands, time-series) and an array of stands, compute 
	the cross-correlation of the data for all baselines.  Return the frequencies 
	and visibilities as a two-elements tuple.
	
	.. versionchanged:: 1.1.0
		Made the 'phaseCenter' keyword more flexible.  It can now be either:
		 * 'z' to denote the zenith,
		 * a ephem.Body instances which has been computed for the observer, or
		 * a two-element tuple of azimuth, elevation in degrees.
		 
	.. versionchanged:: 1.0.0
		Added a phase-center keyword that accept a two-element tuple of 
		azimuth and elelvation (in degrees) to change where the 
		correlations are phased to
	"""
	
	# Since we want to compute Stokes parameters, we need both pols
	pol1 = 0
	pol2 = 1
	
	antennas1 = [a for a in antennas if a.pol == pol1]
	signalsIndex1 = [i for (i, a) in enumerate(antennas) if a.pol == pol1]
	antennas2 = [a for a in antennas if a.pol == pol2]
	signalsIndex2 = [i for (i, a) in enumerate(antennas) if a.pol == pol2]

	nStands = len(antennas1)
	baselines = uvUtils.getBaselines(antennas1, antennas2=antennas2, IncludeAuto=IncludeAuto, Indicies=True)
	nBL = len(baselines)

	# Figure out if we are working with complex (I/Q) data or only real.  This
	# will determine how the FFTs are done since the real data mirrors the pos-
	# itive and negative Fourier frequencies.
	if signals.dtype.kind == 'c':
		lFactor = 1
		doFFTShift = True
		CentralFreq = float(CentralFreq)
	else:
		lFactor = 2
		doFFTShift = False

	if SampleRate is None:
		SampleRate = dp_common.fS
	freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/SampleRate)
	if doFFTShift:
		freq += CentralFreq
		freq = numpy.fft.fftshift(freq)
	freq = freq[:LFFT]
	
	# Get the location of the phase center in radians and create a 
	# pointing vector
	if phaseCenter == 'z':
		azPC = 0.0
		elPC = numpy.pi/2.0
	else:
		if isinstance(phaseCenter, ephem.Body):
			azPC = phaseCenter.az * 1.0
			elPC = phaseCenter.alt * 1.0
		else:
			azPC = phaseCenter[0]*numpy.pi/180.0
			elPC = phaseCenter[1]*numpy.pi/180.0
	source = numpy.array([numpy.cos(elPC)*numpy.sin(azPC), 
					  numpy.cos(elPC)*numpy.cos(azPC), 
					  numpy.sin(elPC)])
					  
	# Define the cable/signal delay caches to help correlate along and compute 
	# the delays that we need to apply to align the signals
	dlyRef = len(freq)/2
	delays1 = numpy.zeros((nStands,LFFT))
	delays2 = numpy.zeros((nStands,LFFT))
	for i in list(range(nStands)):
		xyz1 = numpy.array([antennas1[i].stand.x, antennas1[i].stand.y, antennas1[i].stand.z])
		xyz2 = numpy.array([antennas2[i].stand.x, antennas2[i].stand.y, antennas2[i].stand.z])
		
		delays1[i,:] = antennas1[i].cable.delay(freq) - numpy.dot(source, xyz1) / vLight
		delays2[i,:] = antennas2[i].cable.delay(freq) - numpy.dot(source, xyz2) / vLight
	if delays1[:,dlyRef].min() < delays2[:,dlyRef].min():
		minDelay = delays1[:,dlyRef].min()
	else:
		minDelay = delays2[:,dlyRef].min()
	delays1 -= minDelay
	delays2 -= minDelay

	# F - defaults to running parallel in C via OpenMP
	if window is noWindow:
		# Data without a window function provided
		if signals.dtype.kind == 'c':
			FEngine = _core.FEngineC2
		else:
			FEngine = _core.FEngineR2
		signalsF1, validF1 = FEngine(signals[signalsIndex1,:], freq, delays1, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			FEngine = _core.FEngineC3
		else:
			FEngine = _core.FEngineR3
		signalsF1, validF1 = FEngine(signals[signalsIndex1,:], freq, delays1, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel, window=window)
	
	if window is noWindow:
		signalsF2, validF2 = FEngine(signals[signalsIndex2,:], freq, delays2, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel)
	else:
		signalsF2, validF2 = FEngine(signals[signalsIndex2,:], freq, delays2, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate, ClipLevel=ClipLevel, window=window)

	# X
	output = _stokes.XEngine2(signalsF1, signalsF2, validF1, validF2)
	if not IncludeAuto:
		# Remove auto-correlations from the output of the X engine if we don't 
		# need them.  To do this we need to first build the full list of baselines
		# (including auto-correlations) and then prune that.
		baselinesFull = uvUtils.getBaselines(antennas1, antennas2=antennas2, IncludeAuto=True, Indicies=True)
		fom = numpy.array([a1-a2 for (a1,a2) in baselinesFull])
		nonAuto = numpy.where( fom != 0 )[0]
		output = output[:,nonAuto,:]

	# Divide the cross-multiplied data by the number of channels used
	output /= LFFT
	if signals.dtype.kind != 'c':
		output /= 2.0
		
	# Apply cable gain corrections (if needed)
	if GainCorrect:
		for bl in xrange(output.shape[0]):
			cableGain1 = antennas1[baselines[bl][0]].cable.gain(freq)
			cableGain2 = antennas2[baselines[bl][1]].cable.gain(freq)
			
			output[:,bl,:] /= numpy.sqrt(cableGain1*cableGain2)
			
	# Create antenna baseline list (if needed)
	if ReturnBaselines:
		antennaBaselines = []
		for bl in xrange(output.shape[1]):
			antennaBaselines.append( (antennas1[baselines[bl][0]], antennas2[baselines[bl][1]]) )
		returnValues = (antennaBaselines, freq, output)
	else:
		returnValues = (freq, output)

	return returnValues	
