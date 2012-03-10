# -*- coding: utf-8 -*-

"""
Python module to handle the channelization and cross-correlation of TBW and
TBN data.  The main python functions in this module are:
  * calcSpectra - calculate power spectra for a collection of signals
  * FXCorrelator - calculate cross power spectra for a collection of signals
both of which have been deprecated in favor of the new C extension based 
routines listed below.

The main python/C extension functions in this module are:
  * SpecMaster - similar to calcSpectra but uses the _spec module for all 
    computations and does not support automatic sub-integration
  * SpecMasterP - SpecMaster with a 64-tap uniform DFT filter bank
  * FXMaster - calculate cross power spectra for a collection of signals

Each function is set up to process the signals in parallel using the 
multiprocessing module and accepts a variety of options controlling the processing
of the data, including various window functions and time averaging.
"""

import os
import sys
import numpy

from lsl.common.constants import c as vLight
from lsl.common import dp as dp_common
from lsl.common.constants import *
from lsl.common.warns import warnDeprecated
from lsl.correlator import uvUtils

import _spec
import _stokes
import _core

__version__ = '0.6'
__revision__ = '$Rev$'
__all__ = ['pol2pol', 'noWindow', 'calcSpectrum', 'calcSpectra', 'SpecMaster', 'SpecMasterP', 'StokesMaster', 'FXMaster', 'FXStokes', '__version__', '__revision__', '__all__']


def pol2pol(pol):
	"""
	Convert a polarization string, e.g., XX or XY, to a numeric :mod:`lsl.common.stations.Antena`
	instance polarization.
	"""
	
	pol = pol.upper()
	out = []
	for p in pol:
		if p == 'X':
			out.append(0)
		elif p == 'Y':
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


def calcSpectrum(signal, LFFT=64, window=noWindow, verbose=False):
	"""
	Worker function for calcSpectra.
	"""

	# Figure out if we are working with complex (I/Q) data or only real.  This
	# will determine how the FFTs are done since the real data mirrors the pos-
	# itive and negative Fourier frequencies.
	if signal.dtype.kind == 'c':
		lFactor = 1
	else:
		lFactor = 2

	out = numpy.zeros(LFFT-1)

	# Remove the mean
	intSignal = 1.0*signal
	intSignal -= intSignal.mean()

	nChunks = intSignal.shape[0]/lFactor/LFFT
	if verbose:
		print "  %i samples -> %i chunks" % (intSignal.shape[0], nChunks)

	for i in range(nChunks):
		cs = intSignal[i*lFactor*LFFT:(i+1)*lFactor*LFFT]
		cs *= window(lFactor*LFFT)
		cf = numpy.fft.fft(cs, lFactor*LFFT)
		cp = numpy.abs(cf)**2.0 / lFactor / LFFT
			
		out += cp[1:LFFT]
	out /= float(nChunks)

	return out


def calcSpectra(signals, LFFT=64, SampleAverage=None, window=noWindow, DisablePool=False, verbose=False, SampleRate=None, CentralFreq=0.0):
	"""
	Given a collection of time series data with inputs on the first dimension
	and data on the second, compute the spectra for all inputs.  By default, 
	all data in the time series are average into a single spectrum.  However, this 
	behavior can be modified if the SampleAverage keyword is set.  SampleAverage 
	specifies how many individual spectra of length LFFT are to be averaged.

	.. versionchanged:: 0.3.4
		Prior to LSL version 0.3.4, the window functions available for calcSpectra 
		were limited to Blackman and an (untest) Polyphase filter.  With version
		0.4.0, the window to be used is passed to the function call via the 'window'
		keyword and an "empty" window is provided by the module.  This allows for
		the various window function defined in numpy (i.e., bartlett, blackman, 
		hamming, etc.) to be used.  It also makes it easier to filter the data using
		a custom window.  For example, a Kaiser window with a shape factor of 5 
		could be made with:
			
			>>> import numpy
			>>> def newWindow(L):
			...      return numpy.kaiser(L, 5)
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
	freq = freq[1:LFFT]

	# Calculate the number of output bins with respect to time.  If SampleAverage 
	# is None, average up everything.
	if SampleAverage is None:
		SampleAverage = signals.shape[1] / lFactor / LFFT
	nSamples = signals.shape[1] / lFactor / LFFT / int(SampleAverage)

	# The multiprocessing module allows for the creation of worker pools to help speed
	# things along.  If the processing module is found, use it.  Otherwise, set
	# the 'usePool' variable to false and run single threaded.
	try:
		from multiprocessing import Pool, cpu_count
		
		# To get results pack from the pool, you need to keep up with the workers.  
		# In addition, we need to keep up with which workers goes with which 
		# signal since the workers are called asynchronously.  Thus, we need a 
		# taskList array to hold tuples of signal ('count') and workers.
		taskPool = Pool(processes=int(numpy.ceil(cpu_count()*0.70)))
		taskList = []

		usePool = True
	except ImportError:
		usePool = False
	
	# Turn off the thread pool if we are explicitly told not to use it.
	if DisablePool:
		usePool = False

	output = numpy.zeros( (signals.shape[0], nSamples, freq.shape[0]) )
	for i in range(signals.shape[0]):
		for j in range(nSamples):
			# Figure out which section of each signal to average
			bchan = j*lFactor*LFFT*SampleAverage
			echan = (j+1)*lFactor*LFFT*SampleAverage
			if echan >= signals.shape[1]:
				echan = -1

			if verbose:
				print "Working on signal %i of %i, section %i to %i" % (i+1, signals.shape[0], bchan, echan)

			# If pool, pool...  Otherise don't
			if usePool:
				task = taskPool.apply_async(calcSpectrum, args=(signals[i,bchan:echan], ), 
									kwds={'LFFT': LFFT, 'window': window, 'verbose': verbose})
				taskList.append((i,j,task))
			else:
				tempPS = calcSpectrum(signals[i,bchan:echan], LFFT=LFFT, window=window, verbose=verbose)

				if doFFTShift:
					output[i,j,:] = numpy.fft.fftshift(tempPS)
				else:
					output[i,j,:] = tempPS

	# If pooling... Close the pool so that it knows that no ones else is joining.  
	# Then, join the workers together and wait on the last one to finish before 
	# saving the results.
	if usePool:
		taskPool.close()
		taskPool.join()

		# This is where he taskList list comes in handy.  We now know who did what
		# when we unpack the various results
		for i,j,task in taskList:
			if doFFTShift:
				output[i,j,:] = numpy.fft.fftshift(task.get())
			else:
				output[i,j,:] = task.get()

		# Destroy the taskPool
		del(taskPool)

	# Trim single dimensions from the output if they are found.  This is the 
	# default behavior
	output = numpy.squeeze(output)

	return (freq, output)


def SpecMaster(signals, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0, ClipLevel=0):
	"""
	A more advanced version of calcSpectra that uses the _spec C extension 
	to handle all of the P.S.D. calculations in parallel.  Returns a two-
	element tuple of the frequencies (in Hz) and PSDs in dB/RBW.
	
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
	freq = freq[1:LFFT]
	
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


def SpecMasterP(signals, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0, ClipLevel=0):
	"""
	Similar to SpecMaster but uses a 4-tap polyphase filter bank instead
	of a FFT.  Returns a two-element tuple of the frequencies (in Hz) and 
	PSDs in dB/RBW.
	
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
	freq = freq[1:LFFT]
	
	if window is noWindow:
		# Data without a window function provided
		if signals.dtype.kind == 'c':
			output = _spec.PPSDC2(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel)
		else:
			output = _spec.PPSDR2(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			output = _spec.PPSDC3(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel, window=window)
		else:
			output = _spec.PPSDR3(signals, LFFT=LFFT, Overlap=1, ClipLevel=ClipLevel, window=window)
	
	return (freq, output)


def StokesMaster(signals, antennas, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0, ClipLevel=0):
	"""
	Similar to SpecMaster, but accepts an array of signals and a list of 
	antennas in order to compute the PSDs for the four Stokes parameters: 
	I, Q, U, and V.  Returns a two-element tuple of the frequencies (in Hz) 
	and PSDs in dB/RBW.  The PSD are three dimensional with dimensions 
	Stokes parameter (0=I, 1=Q, 2=U, 3=V) by stand by channel).
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
	freq = freq[1:LFFT]
	
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


def FXMaster(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, verbose=False, window=noWindow, SampleRate=None, CentralFreq=0.0, Pol='XX', GainCorrect=False, ReturnBaselines=False, ClipLevel=0):
	"""
	A more advanced version of FXCorrelator for TBW and TBN data.  Given an 
	2-D array of signals (stands, time-series) and an array of stands, compute 
	the cross-correlation of the data for all baselines.  Return the frequencies 
	and visibilities as a two-elements tuple.
	
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
	freq = freq[1:LFFT]

	# Define the cable/signal delay caches to help correlate along and compute 
	# the delays that we need to apply to align the signals
	dlyRef = len(freq)/2
	delays1 = numpy.zeros((nStands,LFFT-1))
	delays2 = numpy.zeros((nStands,LFFT-1))
	for i in list(range(nStands)):
		delays1[i,:] = antennas1[i].cable.delay(freq) - antennas1[i].stand.z / vLight
		delays2[i,:] = antennas2[i].cable.delay(freq) - antennas2[i].stand.z / vLight
	if delays1[:,dlyRef].max() > delays2[:,dlyRef].max():
		maxDelay = delays1[:,dlyRef].max()
	else:
		maxDelay = delays2[:,dlyRef].max()
	delays1 = maxDelay - delays1
	delays2 = maxDelay - delays2

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


def FXStokes(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, verbose=False, window=noWindow, SampleRate=None, CentralFreq=0.0,  GainCorrect=False, ReturnBaselines=False, ClipLevel=0):
	"""
	A more advanced version of FXCorrelator for TBW and TBN data.  Given an 
	2-D array of signals (stands, time-series) and an array of stands, compute 
	the cross-correlation of the data for all baselines.  Return the frequencies 
	and visibilities as a two-elements tuple.
	
	.. versionchanged:: 0.4.0
		Switched over to passing in Antenna instances generated by the
		:mod:`lsl.common.stations` module instead of a list of stand ID
		numbers.
	"""
	
	# Since we want to compute Stokes parameters, we need both pols
	pol1, pol2 = 0, 1
	
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
	freq = freq[1:LFFT]

	# Define the cable/signal delay caches to help correlate along and compute 
	# the delays that we need to apply to align the signals
	dlyRef = len(freq)/2
	delays1 = numpy.zeros((nStands,LFFT-1))
	delays2 = numpy.zeros((nStands,LFFT-1))
	for i in list(range(nStands)):
		delays1[i,:] = antennas1[i].cable.delay(freq) - antennas1[i].stand.z / vLight
		delays2[i,:] = antennas2[i].cable.delay(freq) - antennas2[i].stand.z / vLight
	if delays1[:,dlyRef].max() > delays2[:,dlyRef].max():
		maxDelay = delays1[:,dlyRef].max()
	else:
		maxDelay = delays2[:,dlyRef].max()
	delays1 = maxDelay - delays1
	delays2 = maxDelay - delays2

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
