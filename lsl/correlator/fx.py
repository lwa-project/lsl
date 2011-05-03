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

from lsl.common import dp as dp_common
from lsl.common.constants import *
from lsl.common.warns import warnDeprecated
from lsl.correlator import uvUtils

import _spec
import _core

__version__ = '0.5'
__revision__ = '$ Revision: 25 $'
__all__ = ['noWindow', 'calcSpectrum', 'calcSpectra', 'SpecMaster', 'SpecMasterP', 'correlate', 'FXCorrelator', 'FXMaster', '__version__', '__revision__', '__all__']


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


def SpecMaster(signals, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0):
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
			output = _spec.FPSDC2(signals, LFFT=LFFT, Overlap=1)
		else:
			output = _spec.FPSDR2(signals, LFFT=LFFT, Overlap=1)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			output = _spec.FPSDC3(signals, LFFT=LFFT, Overlap=1, window=window)
		else:
			output = _spec.FPSDR3(signals, LFFT=LFFT, Overlap=1, window=window)
	
	return (freq, output)


def SpecMasterP(signals, LFFT=64, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0):
	"""
	Similar to SpecMaster but uses a 4-tap polyphase filter bank instead
	of a FFT.  Returns a two-element tuple of the frequencies (in Hz) and 
	PSDs in dB/RBW.
	
	.. note::
		SpecMaster currently average all data given and does not support the
		SampleAverage keyword that calcSpectra does.
		
	.. seealso::
		:mod:`lsl.correlator.filterbank`
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
			output = _spec.PPSDC2(signals, LFFT=LFFT, Overlap=1)
		else:
			output = _spec.PPSDR2(signals, LFFT=LFFT, Overlap=1)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			output = _spec.PPSDC3(signals, LFFT=LFFT, Overlap=1, window=window)
		else:
			output = _spec.PPSDR3(signals, LFFT=LFFT, Overlap=1, window=window)
	
	return (freq, output)


def correlate(signal1, signal2, antenna1, antenna2, LFFT=64, Overlap=1, window=noWindow, verbose=False, SampleRate=None, CentralFreq=0.0):
	"""
	Channalize and cross-correlate signals from two antennae.  Both the signals 
	and the stand numbers are needed to implement time delay and phase corrections.
	The resulting visibilities from the cross-correlation are time average and a 
	LFFT-1 length numpy array is returned.  This array does not contain the DC 
	component of the signal.
	"""

	import aipy

	# Figure out if we are working with complex (I/Q) data or only real.  This
	# will determine how the FFTs are done since the real data mirrors the pos-
	# itive and negative Fourier frequencies.
	if signal1.dtype.kind == 'c':
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
	delayRef = len(freq)/2

	delay1 = antenna1.cable.delay(freq)
	delay2 = antenna2.cable.delay(freq)

	if delay2[delayRef] > delay1[delayRef]:
		delay1 = delay2 - delay1
		delay2 = delay2 - delay2
	else:
		delay2 = delay1 - delay2 
		delay1 = delay1 - delay1

	start1 = int(round(delay1[delayRef]*SampleRate))
	start2 = int(round(delay2[delayRef]*SampleRate))

	if verbose:
		print "  Delay for stands #%i, pol. %i: %.2f - %.2f ns (= %i samples)" % (antenna1.stand.id, antenna1.pol, delay1[-1]*1e9, delay1[0]*1e9, start1)
		print "  Delay for stands #%i, pol. %i: %.2f - %.2f ns (= %i samples)" % (antenna2.stand.id, antenna2.pol, delay2[-1]*1e9, delay2[0]*1e9, start2)

	# Remove the mean
	is1 = 1.0*signal1
	is1 -= is1.mean()
	is2 = 1.0*signal2
	is2 -= is2.mean()

	# Compute the length of each time series and find the shortest one.
	length1 = (is1.shape[0] - start1) / LFFT / lFactor
	length2 = (is2.shape[0] - start2) / LFFT / lFactor
	# Factor in  the overlap amount for all but the last full sample
	length = numpy.array([length1, length2]).min() * Overlap - Overlap + 1

	# Begin computing the visibilities and loop over the signal chunks in order 
	# to average in time.
	visibility = numpy.zeros(LFFT-1, dtype=numpy.complex_)
	for i in range(length):
		# Current chunks
		s1 = is1[(start1+int(lFactor*LFFT*(1.0*i)/Overlap)):(start1+int(lFactor*LFFT*((1.0*i)/Overlap+1)))]
		s2 = is2[(start2+int(lFactor*LFFT*(1.0*i)/Overlap)):(start2+int(lFactor*LFFT*((1.0*i)/Overlap+1)))]

		s1 *= window(lFactor*LFFT)
		s2 *= window(lFactor*LFFT)
		
		# The FFTs
		fft1 = numpy.fft.fft(s1, lFactor*LFFT)
		fft2 = numpy.fft.fft(s2, lFactor*LFFT)
		
		if doFFTShift:
			fft1 = numpy.fft.fftshift(fft1)
			fft2 = numpy.fft.fftshift(fft2)

		# Calculate the visibility of the desired part of the frequency space
		tempVis = fft1[1:LFFT]*(fft2[1:LFFT].conj()) / lFactor / LFFT

		# Add it on to the master output visibility for averaging
		visibility += tempVis
	
	# Apply the phase rotator to the visibility
	visibility *= numpy.exp(-2j*numpy.pi*freq*((delay1-start1/SampleRate)-(delay2-start2/SampleRate)))

	# Average and return
	visibility /= float(length)

	return visibility


def FXCorrelator(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, window=noWindow, CrossPol=None, DisablePool=False, verbose=False, SampleRate=None, CentralFreq=0.0):
	"""
	A basic FX correlators for the TBW data.  Given an 2-D array of signals
	(stands, time-series) and an array of stands, compute the cross-correlation of
	the data for all baselines.  If cross-polarizations need to be calculated, the
	CrossPol keyword allows for the other polarization data to be entered into the
	correlator.  Return the frequencies and visibilities as a two-elements tuple.

	.. deprecated:: 0.3.4
		FXCorrelator has been deprecated as of LSL version 0.4.0 in favor of the new
		C-based correlator "FXMaster".  The new correlator, however, does not 
		currently support windowing the data.

	.. versionchanged:: 0.3.4
		Prior to LSL version 0.4.0, the window functions available for FXCorrelator
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
			
	.. versionchanged:: 0.4.0
		Switched over to passing in Antenna instances generated by the
		:mod:`lsl.common.stations` module instead of a list of stand ID
		numbers.
	"""

	warnDeprecated("FXCorrelator", "Please use the C-based FX correlator 'FXMaster'")

	N = len(antennas)
	baselines = uvUtils.getBaselines(antennas, IncludeAuto=IncludeAuto, Indicies=True)
	Nbase = len(baselines)
	output = numpy.zeros( (Nbase, LFFT-1), dtype=numpy.complex64)

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

	# The multiprocessing module allows for the creation of worker pools to help speed
	# things along.  If the processing module is found, use it.  Otherwise, set
	# the 'usePool' variable to false and run single threaded.
	try:
		from multiprocessing import Pool, cpu_count
		
		# To get results pack from the pool, you need to keep up with the workers.  
		# In addition, we need to keep up with which workers goes with which 
		# baseline since the workers are called asynchronously.  Thus, we need a 
		# taskList array to hold tuples of baseline ('count') and workers.
		taskPool = Pool(processes=int(numpy.ceil(cpu_count()*0.70)))
		taskList = []

		usePool = True
	except ImportError:
		usePool = False

	# Turn off the thread pool if we are explicitly told not to use it.
	if DisablePool:
		usePool = False

	# Loop over baselines
	count = 0
	for i,j in baselines:
		if verbose:
			if i == j:
				print "Auto-correlating stand %i" % stands[i]
			else:
				print "Cross-correlating stands %i and %i" % (stands[i], stands[j])
			
		# If pool, pool...  Otherwise don't
		if usePool:
			if CrossPol is not None:
				task = taskPool.apply_async(correlate, args=(signals[i,:], CrossPol[j,:], antennas[i], antennas[j]), 
								kwds={'LFFT': LFFT, 'Overlap': Overlap, 'window': window, 'verbose': verbose, 'SampleRate': SampleRate, 'CentralFreq': CentralFreq})
			else:
				task = taskPool.apply_async(correlate, args=(signals[i,:], signals[j,:], antennas[i], antennas[j]), 
								kwds={'LFFT': LFFT, 'Overlap': Overlap, 'window': window, 'verbose': verbose, 'SampleRate': SampleRate, 'CentralFreq': CentralFreq})
			taskList.append((count,task))
		else:
			if CrossPol is not None:
				tempCPS = correlate(signals[i,:], CrossPol[j,:], antennas[i], antennas[j], LFFT=LFFT, Overlap=Overlap, window=window, verbose=verbose, SampleRate=SampleRate, CentralFreq=CentralFreq)
			else:
				tempCPS = correlate(signals[i,:], signals[j,:], antennas[i], antennas[j], LFFT=LFFT, Overlap=Overlap, window=window, verbose=verbose, SampleRate=SampleRate, CentralFreq=CentralFreq)

			if doFFTShift:
				output[count,:] = tempCPS
			else:
				output[count,:] = tempCPS
		count = count + 1

	# If pooling... Close the pool so that it knows that no ones else is joining.  
	# Then, join the workers together and wait on the last one to finish before 
	# saving the results.
	if usePool:
		taskPool.close()
		taskPool.join()

		# This is where he taskList list comes in handy.  We now know who did what
		# when we unpack the various results
		for count,task in taskList:
			if doFFTShift:
				output[count,:] = task.get()
			else:
				output[count,:] = task.get()

		# Destroy the taskPool
		del(taskPool)

	return (freq, output)


def FXMaster(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, verbose=False, window=noWindow, SampleRate=None, CentralFreq=0.0, GainCorrect=False):
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

	nStands = len(antennas)
	baselines = uvUtils.getBaselines(antennas, IncludeAuto=IncludeAuto, Indicies=True)
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
	delays = numpy.zeros((nStands,LFFT-1))
	for i in list(range(nStands)):
		delays[i,:] = antennas[i].cable.delay(freq)
	delays = delays[:,dlyRef].max() - delays

	# F - defaults to running parallel in C via OpenMP
	if window is noWindow:
		# Data without a window function provided
		if signals.dtype.kind == 'c':
			signalsF = _core.FEngineC2(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate)
		else:
			signalsF = _core.FEngineR2(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate)
	else:
		# Data with a window function provided
		if signals.dtype.kind == 'c':
			signalsF = _core.FEngineC3(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, 
									SampleRate=SampleRate, window=window)
		else:
			signalsF = _core.FEngineR3(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, 
									SampleRate=SampleRate, window=window)
	signalsFC = signalsF.conj()

	# X
	output = _core.XEngine2(signalsF, signalsFC)

	# Divide the cross-multiplied data by the number of channels used
	output /= LFFT
	if signals.dtype.kind != 'c':
		output /= 2.0
		
	# Apply cable gain corrections (if needed)
	if GainCorrect:
		for bl in xrange(output.shape[0]):
			cableGain1 = antennas[baselines[bl][0]].cable.gain(freq)
			cableGain2 = antennas[baselines[bl][1]].cable.gain(freq)
			
			output[bl,:] /= numpy.sqrt(cableGain1*cableGain2)

	return (freq, output)


#def PXMaster(signals, antennas, LFFT=64, Overlap=1, IncludeAuto=False, verbose=False, window=noWindow, SampleRate=None, CentralFreq=0.0):
	#"""
	#A version of FXMaster that uses a 64-tap polyphase filter for the FFT step
	#rather than a normal FFT.  Returns the frequencies and visibilities as a 
	#two-elements tuple.
	#"""

	#nStands = stands.shape[0]
	#baselines = uvUtils.getBaselines(stands, IncludeAuto=IncludeAuto, Indicies=True)
	#nBL = len(baselines)

	## Figure out if we are working with complex (I/Q) data or only real.  This
	## will determine how the FFTs are done since the real data mirrors the pos-
	## itive and negative Fourier frequencies.
	#if signals.dtype.kind == 'c':
		#lFactor = 1
		#doFFTShift = True
		#CentralFreq = float(CentralFreq)
	#else:
		#lFactor = 2
		#doFFTShift = False

	#if SampleRate is None:
		#SampleRate = dp_common.fS
	#freq = numpy.fft.fftfreq(lFactor*LFFT, d=1.0/SampleRate)
	#if doFFTShift:
		#freq += CentralFreq
		#freq = numpy.fft.fftshift(freq)
	#freq = freq[1:LFFT]

	## Define the cable/signal delay caches to help correlate along and compute 
	## the delays that we need to apply to align the signals
	#dlyRef = len(freq)/2
	#delays = numpy.zeros((nStands,LFFT-1))
	#for i in list(range(nStands)):
		#delays[i,:] = antennas.cable.delay(freq)
	#delays = delays[:,dlyRef].max() - delays

	## F - defaults to running parallel in C via OpenMP
	#if window is noWindow:
		#print "No Window"
		## Data without a window function provided
		#if signals.dtype.kind == 'c':
			#signalsF = _core.PEngineC2(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate)
		#else:
			#signalsF = _core.PEngineR2(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, SampleRate=SampleRate)
	#else:
		## Data with a window function provided
		#if signals.dtype.kind == 'c':
			#signalsF = _core.PEngineC3(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, 
									#SampleRate=SampleRate, window=window)
		#else:
			#signalsF = _core.PEngineR3(signals, freq, delays, LFFT=LFFT, Overlap=Overlap, 
									#SampleRate=SampleRate, window=window)
	#signalsFC = signalsF.conj()

	## X
	#output = _core.XEngine2(signalsF, signalsFC)

	## Divide the cross-multiplied data by the number of channels used
	#output /= LFFT
	#if signals.dtype.kind != 'c':
		#output /= 2.0

	#return (freq, output)
