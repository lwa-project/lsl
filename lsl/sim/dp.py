# -*- coding: utf-8 -*-

"""
Module to simulate observations made with the DP system.
"""

import time
import numpy
from aipy import coord as aipycoord

from lsl import astro
from lsl.common import dp as dp_common
from lsl.common import stations as lwa_common
from lsl.sim import tbw
from lsl.sim import tbn
from lsl.sim import drx
from lsl.sim import vis
from lsl.correlator import uvUtils
from lsl.reader.tbn import filterCodes as TBNFilters
from lsl.reader.drx import filterCodes as DRXFilters

__version__ = '0.3'
__revision__ = '$Rev$'
__all__ = ['basicSignal', 'pointSource', '__version__', '__revision__', '__all__']


def __basicTBW(fh, stands, nFrames, **kwargs):
	"""
	Private function for generating a basic TBW signal.
	"""

	tStart = kwargs['tStart']
	bits = kwargs['bits']
	verbose = kwargs['verbose']
	noiseStrength = kwarg['noiseStrength']
	
	if bits == 12:
		maxValue = 2047
		samplesPerFrame = 400
	else:
		maxValue =  7
		samplesPerFrame = 1200

	if verbose:
		print "Simulating %i captures of %i-bit TBW data for %i stands:" % (int(numpy.ceil(nFrames / 30000.0)), bits, len(stands))

	nCaptures = int(numpy.ceil(nFrames / 30000.0))
	for capture in range(nCaptures):
		for stand1, stand2 in zip(stands[0::2], stands[1::2]):
			FramesThisBatch = nFrames - capture*30000
			if FramesThisBatch > 30000:
				FramesThisBatch = 30000
			if verbose:
				print " capture %i, stands %i and %i" % (capture+1, stand1, stand2)

			for i in range(FramesThisBatch):
				t = long(tStart*dp_common.fS) + i*samplesPerFrame
				t += long(60*dp_common.fS*capture)
				tFrame = t/dp_common.fS - tStart + numpy.arange(samplesPerFrame, dtype=numpy.float32) / dp_common.fS
				
				cFrame = tbw.SimFrame(stand=stand1, frameCount=i+1, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn(2, samplesPerFrame)
				cFrame.xy[0,:] *= maxValue*noiseStrength
				cFrame.xy[0,:] += maxValue*numpy.cos(2*numpy.pi*40.0e6*tFrame)
				cFrame.xy[1,:] *= maxValue*noiseStrength
				cFrame.xy[1,:] += maxValue*numpy.cos(2*numpy.pi*60.0e6*tFrame)
				cFrame.writeRawFrame(fh)

				cFrame = tbw.SimFrame(stand=stand2, frameCount=i+1, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn(2, samplesPerFrame)
				cFrame.xy[0,:] *= maxValue*noiseStrength
				cFrame.xy[0,:] += maxValue*numpy.cos(2*numpy.pi*30.0e6*tFrame)
				cFrame.xy[1,:] *= maxValue*noiseStrength
				cFrame.xy[1,:] += maxValue*numpy.cos(2*numpy.pi*50.0e6*tFrame)
				cFrame.writeRawFrame(fh)


def __basicTBN(fh, stands, nFrames, **kwargs):
	"""
	Private function for generating a basic TBN signal.
	"""

	tStart = kwargs['tStart']
	filter = kwargs['filter']
	verbose = kwargs['verbose']
	noiseStrength = kwarg['noiseStrength']
	sampleRate = TBNFilters[filter]
	
	maxValue = 127
	samplesPerFrame = 512
	upperSpike = sampleRate / 4.0
	lowerSpike = -sampleRate / 4.0
	
	if verbose:
		print "Simulating %i frames of TBN Data @ %.2f kHz for %i stands:" % (nFrames, sampleRate/1e3, len(stands))
	
	for i in range(nFrames):
		if i % 1000 == 0 and verbose:
			print " frame %i" % (i+1)
		t = long(tStart*dp_common.fS) + long(i*dp_common.fS*samplesPerFrame/sampleRate)
		tFrame = t/dp_common.fS - tStart + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sampleRate
		for stand in stands:
			cFrame = tbn.SimFrame(stand=stand, pol=0, freq=40e6, gain=20, frameCount=i+1, obsTime=t)
			cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
			cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
			cFrame.iq *= maxValue*noiseStrength
			cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*upperSpike*tFrame)
			cFrame.writeRawFrame(fh)

			cFrame = tbn.SimFrame(stand=stand, pol=1, freq=40e6, gain=20, frameCount=i+1, obsTime=t)
			cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
			cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
			cFrame.iq *= maxValue*noiseStrength
			cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*lowerSpike*tFrame)
			cFrame.writeRawFrame(fh)

def __basicDRX(fh, stands, nFrames, **kwargs):
	"""
	Private function for generating a basic TBN signal.
	"""

	tStart = kwargs['tStart']
	filter = kwargs['filter']
	nTuning = kwargs['nTuning']
	verbose = kwargs['verbose']
	noiseStrength = kwarg['noiseStrength']
	sampleRate = DRXFilters[filter]
	
	maxValue = 7
	samplesPerFrame = 4096
	upperSpike1 = sampleRate / 4.0
	lowerSpike1 = -sampleRate / 4.0
	upperSpike2 = sampleRate / 3.0
	lowerSpike2 = -sampleRate / 3.0

	if verbose:
		print "Simulating %i frames of DRX Data @ %.2f MHz for %i beams, %i tunings each:" % (nFrames, sampleRate/1e6, len(stands), nTuning)

	beams = stands
	for i in range(nFrames):
		if i % 1000 == 0 and verbose:
			print " frame %i" % i
		t = long(tStart*dp_common.fS) + long(i*dp_common.fS*samplesPerFrame/sampleRate)
		tFrame = t/dp_common.fS - tStart + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sampleRate
		for beam in beams:
			for tune in range(1, nTuning+1):
				if tune == 1:
					# Tuning 1:
					cFrame = drx.SimFrame(beam=beam, tune=1, pol=0, frameCount=i+1, filterCode=filter, timeOffset=0, obsTime=t, flags=0)
					cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
					cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
					cFrame.iq *= maxValue*noiseStrength
					cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*upperSpike1*tFrame)
					cFrame.writeRawFrame(fh)
			
					cFrame = drx.SimFrame(beam=beam, tune=1, pol=1, frameCount=i+1, filterCode=filter, timeOffset=0, obsTime=t, flags=0)
					cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
					cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
					cFrame.iq *= maxValue*noiseStrength
					cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*lowerSpike1*tFrame)
					cFrame.writeRawFrame(fh)
				else:
					# Tuning 2:
					cFrame = drx.SimFrame(beam=beam, tune=2, pol=0, frameCount=i+1, filterCode=filter, timeOffset=0, obsTime=t, flags=0)
					cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
					cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
					cFrame.iq *= maxValue*noiseStrength
					cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*lowerSpike2*tFrame)
					cFrame.writeRawFrame(fh)
			
					cFrame = drx.SimFrame(beam=beam, tune=2, pol=1, frameCount=i+1, filterCode=filter, timeOffset=0, obsTime=t, flags=0)
					cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
					cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
					cFrame.iq *= maxValue*noiseStrength
					cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*upperSpike2*tFrame)
					cFrame.writeRawFrame(fh)


def basicSignal(fh, stands, nFrames, mode='DRX', filter=6, nTuning=2, bits=12, tStart=0, noiseStength=0.1, verbose=False):
	"""
	Generate a collection of frames with a basic test signal for TBW, TBN, 
	and DRX.  The signals for the three modes are:
	
	TBW
	  * noise + 40 MHz signal for x-pol.
	  * noise + 60 MHz signal for y-pol.
	    -> odd stands
	  * noise + 30 MHz signal for x-pol.
	  * noise + 50 MHz signal for ypol.
	    -> even stands
	
	TBN
	  * noise + (sampleRate/4) kHz signal for x-pol. and noise + 
	    (-sampleRate/4) for y-pol.

	DRX
	  * noise + (sampleRate/4) kHz signal for x-pol. and noise + 
	    (-sampleRate/4) for y-pol. -> tuning 1
	  * noise + (-sampleRate/3) kHz signal for x-pol. and noise + 
	    (sampleRate/3) for y-pol. -> tuning 2
	    
	All modes need to have stands (beams in the case of DRX) and number of
	frames to generate.  TBW also needs to 'bits' keyword set to generate 
	either 12-bit or 4-bit data.  The TBN and DRX frames need the 'filter'
	keyword set to specify the filter width.  In addition, the 'stands' 
	argument is interpreted as beam numbers for DRX.
	
	.. versionchanged:: 0.4.4
		Added the `noiseStrength` keyword to control how much noise is added to 
		the data.
	"""

	if tStart == 0:
		tStart = time.time()

	if mode == 'TBW':
		__basicTBW(fh, stands, nFrames, bits=bits, tStart=tStart, noiseStrengh=noiseStrength, verbose=verbose)
	elif mode == 'TBN':
		__basicTBN(fh, stands, nFrames, filter=filter, tStart=tStart, noiseStrengh=noiseStrength, verbose=verbose)
	elif mode == 'DRX':
		__basicDRX(fh, stands, nFrames, filter=filter, nTuning=nTuning, tStart=tStart, noiseStrengh=noiseStrength, verbose=verbose)
	else:
		raise RuntimeError("Unknown observations mode: %s" % mode)


def __getAntennaArray(station, stands, time, freqs):
	"""
	Given a LWA station object, a list of stands, an observation time, and
	a list of frequencies in Hz, build an aipy AntennaArray object.
	"""

	return vis.buildSimArray(station, stands, freqs/1e9, jd=astro.unix_to_utcjd(time))


def __getSourceParameters(aa, timestamp, srcs):
	"""
	Given an aipy AntennaArray object, an observation time, and aipy.src 
	object, return all of the parameters needed for a simulation.
	"""
	
	# Set the time for the array
	aa.set_unixtime(timestamp)

	# Compute the source parameters
	srcs_tp = []
	srcs_mp = []
	srcs_jy = []
	srcs_fq = []
	for name,src in srcs.iteritems():
		## Update the source's coordinates
		src.compute(aa)

		## Get parameters
		top = src.get_crds(crdsys='top', ncrd=3)	# topo. coords.
		map = src.map							# equitorial -> topo. rotation matrix
		jys = src.get_jys()						# F_nu
		frq = aa.get_afreqs()					# nu

		## Fix the lowest frequencies to avoid problems with the flux blowing up
		## at nu = 0 Hz by replacing flux values below 1 MHz with the flux at 
		## 1 MHz
		Jyat1MHz = jys[ numpy.where( numpy.abs(frq-0.001) == numpy.abs(frq-0.001).min() ) ]
		jys = numpy.where( frq >= 0.001, jys, Jyat1MHz )

		## Filter out sources that are below the horizon or have no flux
		srcAzAlt = aipycoord.top2azalt(top)
		if srcAzAlt[1] <= 0 or jys.sum() <= 0:
			continue

		## Save values into the source arrays
		srcs_tp.append( top )
		srcs_mp.append( map )
		srcs_jy.append( jys )
		srcs_fq.append( frq )

	# Return the values as a dictionary
	return {'topo': srcs_tp, 'trans': srcs_mp, 'flux': srcs_jy, 'freq': srcs_fq}


def __buildSignals(aa, stands, srcParams, times, pol='x', phaseCenter='z'):
	"""
	Given an aipy AntennaArray, a list of stand numbers, a dictionary of source 
	parameters, and an array of times in ns, return a numpy array of the simulated 
	signals that is Nstands x Ntimes in shape
	."""

	# Find out how many stands, srcs, and samples (times) we are working with
	Nstand = len(aa.ants)
	Nsrc = len(srcParams['topo'])
	Ntime = len(times)
	
	# Get the topocentric coorindates for the zenith
	zen = aipycoord.azalt2top(numpy.array([[numpy.pi/4],[numpy.pi/2]]))
	zen = numpy.squeeze(zen)
	
	# Update the phase center if necessary
	if phaseCenter == 'z':
		phaseCenterMap = aipycoord.eq2top_m(0.0, aa.lat)
	else:
		phaseCenter.compute(aa)
		phaseCenterMap = phaseCenter.map

	# Setup a temporary array to hold the signals per source, stand, and time.  
	# This array is complex so that it can accomidate both TBW and TBN data at
	# the same time
	temp = numpy.zeros((Nsrc, Nstand, Ntime), dtype=numpy.complex64)

	# Loop over sources and stands to build up the signals
	srcCount = 0
	for topo,trans,flux,freq in zip(srcParams['topo'], srcParams['trans'], srcParams['flux'], srcParams['freq']):
		antCount = 0
		for ant,std in zip(aa.ants, stands):
			# Zeroth, get the beam response in the direction of the current source for all frequencies
			antResponse = numpy.squeeze( ant.bm_response(topo, pol=pol) )

			# First, do the geometric delay
			geoDelay = ( numpy.dot(trans, ant.pos).transpose() )[2] 
			geoDelayPC = ( numpy.dot(phaseCenterMap, ant.pos).transpose() )[2]

			# Second, do the cable delay
			Delayat1MHz = std.cable.delay(frequency=1.0e6) * 1e9 # s -> ns
			cblDelay = std.cable.delay(frequency=aa.get_afreqs()*1e9) * 1e9 # s -> ns
			# NB: Replace the cable delays below 1 MHz with the 1 MHz value to keep the 
			# delays from blowing up for small f
			cblDelay = numpy.where( freq >= 0.001, cblDelay, Delayat1MHz )

			for a,j,f,d in zip(antResponse, flux, freq, cblDelay):
				factor = a * numpy.sqrt(j)
				angle = 2*numpy.pi*f*(times + (d - (geoDelay-geoDelayPC)))
				temp[srcCount,antCount,:] += factor*(numpy.cos(angle) + 1j*numpy.sin(angle))
			antCount = antCount + 1
		srcCount = srcCount + 1

	# Sum over sources and done
	tdSignals = temp.sum(axis=0)
	return tdSignals


def __pointSourceTBW(fh, stands, src, nFrames, **kwargs):
	"""
	Private function to build TBW point sources.
	"""
	
	bits = kwargs['bits']
	tStart = kwargs['tStart']
	phaseCenter = kwargs['phaseCenter']
	verbose = kwargs['verbose']
	noiseStrength = kwargs['noiseStrength']
	
	sampleRate = dp_common.fS
	freqs = (numpy.fft.fftfreq(1024, d=1.0/sampleRate))[1:512]
	aa = __getAntennaArray(lwa_common.lwa1, stands, tStart, freqs)
	
	if bits == 12:
		maxValue = 2047
		samplesPerFrame = 400
	else:
		maxValue =  7
		samplesPerFrame = 1200

	if verbose:
		print "Simulating %i captures of %-bit TBW data for %i stands:" % (int(numpy.ceil(nFrames / 30000.0)), bits, len(stands))

	nCaptures = int(numpy.ceil(nFrames / 30000.0))
	for capture in range(nCaptures):
		j = 0
		k = 1
		for stand1, stand2 in zip(stands[0::2], stands[1::2]):
			FramesThisBatch = nFrames - capture*30000
			if FramesThisBatch > 30000:
				FramesThisBatch = 30000
			if verbose:
				print " capture %i, stands %i and %i" % (capture+1, stand1, stand2)

			for i in range(FramesThisBatch):
				t = long(tStart*dp_common.fS) + i*samplesPerFrame
				t += long(60*dp_common.fS*capture)
				tFrame = t/dp_common.fS - tStart + numpy.arange(samplesPerFrame, dtype=numpy.float32) / dp_common.fS

				# Get the source parameters
				srcParams = __getSourceParameters(aa, tFrame[0], src)

				# Generate the time series response of each signal at each frequency
				tdSignalsX = __buildSignals(aa, stands, srcParams, tFrame*1e9, pol='x', phaseCenter=phaseCenter)
				tdSignalsY = __buildSignals(aa, stands, srcParams, tFrame*1e9, pol='y', phaseCenter=phaseCenter)

				cFrame = tbw.SimFrame(stand=stand1.stand.id, frameCount=i+1, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn(2, samplesPerFrame)
				cFrame.xy *= maxValue*noiseStrength
				cFrame.xy[0,:] += maxValue*tdSignalsX.real[j,:]
				cFrame.xy[1,:] += maxvalue*tdSignalsY.real[j,:]
				
				cFrame.writeRawFrame(fh)

				cFrame = tbw.SimFrame(stand=stand2.stand.id, frameCount=i+1, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn(2, samplesPerFrame)
				cFrame.xy *= maxValue*noiseStrength
				cFrame.xy[0,:] += maxValue*tdSignalsX.real[k,:]
				cFrame.xy[1,:] += maxValue*tdSignalsY.real[k,:]
				
				cFrame.writeRawFrame(fh)


def __pointSourceTBN(fh, stands, src, nFrames, **kwargs):
	"""
	Private function to build TBN point sources.
	"""
	
	CentralFreq = kwargs['CentralFreq']
	filter = kwargs['filter']
	tStart = kwargs['tStart']
	phaseCenter = kwargs['phaseCenter']
	verbose = kwargs['verbose']
	noiseStrength = kwargs['noiseStrength']
	
	sampleRate = TBNFilters[filter]
	maxValue = 127
	samplesPerFrame = 512
	freqs = (numpy.fft.fftfreq(samplesPerFrame, d=1.0/sampleRate)) + CentralFreq
	freqs = numpy.fft.fftshift(freqs)
	aa = __getAntennaArray(lwa_common.lwa1, stands, tStart, freqs)
	
	if verbose:
		print "Simulating %i frames of TBN Data @ %.2f kHz for %i stands:" % (nFrames, sampleRate/1e3, len(stands))
	
	for i in range(nFrames):
		if i % 1000 == 0 and verbose:
			print " frame %i" % (i+1)
		t = long(tStart*dp_common.fS) + long(i*dp_common.fS*samplesPerFrame/sampleRate)
		tFrame = t/dp_common.fS - tStart + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sampleRate
		
		# Get the source parameters
		srcParams = __getSourceParameters(aa, tFrame[0], src)
		
		# Generate the time series response of each signal at each frequency
		tdSignalsX = __buildSignals(aa, stands, srcParams, tFrame*1e9, pol='x', phaseCenter=phaseCenter)
		tdSignalsY = __buildSignals(aa, stands, srcParams, tFrame*1e9, pol='y', phaseCenter=phaseCenter)
		
		j = 0
		for stand in stands:
			cFrame = tbn.SimFrame(stand=stand.stand.id, pol=0, freq=CentralFreq, gain=19, frameCount=i+1, obsTime=t)
			cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
			cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
			cFrame.iq *= maxValue*noiseStrength
			cFrame.iq += mavValue*tdSignalsX[j,:].astype(numpy.singlecomplex)
			
			cFrame.writeRawFrame(fh)

			cFrame = tbn.SimFrame(stand=stand.stand.id, pol=1, freq=CentralFreq, gain=19, frameCount=i+1, obsTime=t)
			cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
			cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
			cFrame.iq *= maxValue*noiseStrength
			cFrame.iq += maxValue*tdSignalsY[j,:].astype(numpy.singlecomplex)
			
			cFrame.writeRawFrame(fh)
			
			j += 1


def pointSource(fh, stands, src, nFrames, mode='TBN', CentralFreq=49.0e6, filter=7, bits=12, tStart=0, phaseCenter='z', noiseStrength=0.1, verbose=False):
	"""
	Generate a collection of frames with a point source signal for TBW
	and TBN.  The point source is specified as a aipy.src object.
	    
	All modes need to have stands (beams in the case of DRX) and number of
	frames to generate.  TBW also needs to `bits' keyword set to generate 
	either 12-bit or 4-bit data.  The TBN frames need the `filter' keyword 
	set to specify the filter width.
	
	.. versionchanged:: 0.4.4
		Added the `noiseStrength` keyword to control how much noise is added to 
		the data.
	"""

	if tStart == 0:
		tStart = time.time()

	if mode == 'TBW':
		__pointSourceTBW(fh, stands, src, nFrames, bits=bits, tStart=tStart, phaseCenter=phaseCenter, noiseStrength=noiseStrength, verbose=verbose)
	elif mode == 'TBN':
		__pointSourceTBN(fh, stands, src, nFrames, CentralFreq=CentralFreq, filter=filter, tStart=tStart, phaseCenter=phaseCenter, noiseStrength=noiseStrength, verbose=verbose)
	else:
		raise RuntimeError("Unknown observations mode: %s" % mode)
