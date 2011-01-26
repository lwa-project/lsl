# -*- coding: utf-8 -*-

"""Module to simulate observations made with the DP system."""

import functools
import itertools
import time
import numpy
import scipy.stats
from aipy import coord as aipycoord

from lsl import astro
from lsl.common import dp as dp_common
from lsl.common import stations as lwa_common
from lsl.sim import tbw
from lsl.sim import tbn
from lsl.sim import drx
from lsl.sim import vis
from lsl.reader.tbn import filterCodes as TBNFilters
from lsl.reader.drx import filterCodes as DRXFilters

__version__ = '0.1'
__revision__ = '$ Revision: 1 $'
__all__ = ['basicSignal', '__version__', '__revision__', '__all__']


def __basicTBW(fh, stands, nFrames, **kwargs):
	"""Private function for generating a basic TBW signal."""

	tStart = kwargs['tStart']
	bits = kwargs['bits']
	if bits == 12:
		maxValue = 2047
		samplesPerFrame = 400
	else:
		maxValue =  7
		samplesPerFrame = 1200

	nCaptures = int(numpy.ceil(nFrames / 30000.0))
	for capture in range(nCaptures):
		for stand1, stand2 in zip(stands[0::2], stands[1::2]):
			print "Simulating TBW capture %i, stands %i and %i" % (capture+1, stand1, stand2)
			FramesThisBatch = nFrames - capture*30000
			if FramesThisBatch > 30000:
				FramesThisBatch = 30000
			print "-> Frames %i" % FramesThisBatch

			for i in range(FramesThisBatch):
				t = long(tStart*dp_common.fS) + i*samplesPerFrame
				t += long(60*dp_common.fS*capture)
				tFrame = numpy.arange(samplesPerFrame, dtype=numpy.float32)
				
				cFrame = tbw.SimFrame(stand=stand1, frameCount=i+1, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn(2, samplesPerFrame)
				cFrame.xy[0,:] *= maxValue/15.0
				cFrame.xy[0,:] += maxValue*numpy.cos(2*numpy.pi*0.2041*tFrame)
				cFrame.xy[1,:] *= maxValue/15.0
				cFrame.xy[1,:] += maxValue*numpy.cos(2*numpy.pi*0.3061*tFrame)
				
				cFrame.writeRawFrame(fh)

				cFrame = tbw.SimFrame(stand=stand2, frameCount=i+1, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn(2, samplesPerFrame)
				cFrame.xy[0,:] *= maxValue/15.0
				cFrame.xy[0,:] += maxValue*numpy.cos(2*numpy.pi*0.1531*tFrame)
				cFrame.xy[1,:] *= maxValue/15.0
				cFrame.xy[1,:] += maxValue*numpy.cos(2*numpy.pi*0.2551*tFrame)
				
				cFrame.writeRawFrame(fh)


def __basicTBN(fh, stands, nFrames, **kwargs):
	"""Private function for generating a basic TBN signal."""

	tStart = kwargs['tStart']
	filter = kwargs['filter']
	sampleRate = TBNFilters[filter]
	maxValue = 127
	samplesPerFrame = 512
	upperSpike = sampleRate / 4.0
	lowerSpike = -sampleRate / 4.0
	
	for i in range(nFrames):
		if i % 1000 == 0:
			print "Simulating TBN frame %i" % (i+1)
		t = long(tStart*dp_common.fS) + long(i*dp_common.fS*samplesPerFrame/sampleRate)
		tFrame = t/dp_common.fS + numpy.arange(samplesPerFrame, dtype=numpy.float32) / sampleRate
		for stand in stands:
			cFrame = tbn.SimFrame(stand=stand, pol=0, frameCount=i+1, obsTime=t)
			cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
			cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
			cFrame.iq *= maxValue/15.0
			cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*upperSpike*tFrame)
			
			cFrame.writeRawFrame(fh)

			cFrame = tbn.SimFrame(stand=stand, pol=1, frameCount=i+1, obsTime=t)
			cFrame.iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
			cFrame.iq += numpy.random.randn(samplesPerFrame) + 1j*numpy.random.randn(samplesPerFrame)
			cFrame.iq *= maxValue/15.0
			cFrame.iq += maxValue*numpy.exp(2j*numpy.pi*lowerSpike*tFrame)
			
			cFrame.writeRawFrame(fh)

def __basicDRX(fh, stands, nFrames, **kwargs):
	"""Private function for generating a basic TBN signal."""

	tStart = kwargs['tStart']
	filter = kwargs['filter']
	sampleRate = DRXFilters[filter]
	maxValue = 7
	samplesPerFrame = 4096
	beams = stands
	
	part1 = scipy.stats.norm(loc=1000, scale=100, size=samplesPerFrame)
	part2 = scipy.stats.norm(loc=3000, scale=100, size=samplesPerFrame)
	part3 = scipy.stats.norm(loc=1500, scale=100, size=samplesPerFrame)
	part4 = scipy.stats.norm(loc=2500, scale=100, size=samplesPerFrame)

	signal1 = part1.pdf(numpy.arange(samplesPerFrame)) + 0.3*part2.pdf(numpy.arange(samplesPerFrame))
	signal1 *= 0.2 / signal1.max() 
	signal2 = 0.3*part1.pdf(numpy.arange(samplesPerFrame)) + part2.pdf(numpy.arange(samplesPerFrame))
	signal2 *= 0.2 / signal2.max() 
	signal3 = part3.pdf(numpy.arange(samplesPerFrame)) + 0.3*part4.pdf(numpy.arange(samplesPerFrame))
	signal3 *= 0.2 / signal3.max() 
	signal4 = 0.3*part3.pdf(numpy.arange(samplesPerFrame)) + part4.pdf(numpy.arange(samplesPerFrame))
	signal4 *= 0.2 / signal4.max() 
	signal = [signal1, signal2, signal3, signal4]

	norm = functools.partial(numpy.random.normal, loc=8)

	for i in range(nFrames):
		if i % 1000 == 0:
			print "Simulating DRX frame %i" % i
		t = long(tStart*dp_common.fS) + long(i*dp_common.fS*samplesPerFrame/sampleRate)
		for beam in beams:
			for tune in [1, 2]:
				for pol in [0, 1]:
					cFrame = drx.SimFrame(beam=beam, tune=tune, pol=pol, frameCount=i+1, filterCode=filter, timeOffset=0, obsTime=t, flags=0)

					iq = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
					for chan, scale in itertools.izip(numpy.arange(samplesPerFrame), signal[2*(tune-1)+pol]):
						iq.real[chan] = norm(scale=scale)
						iq.imag[chan] = norm(scale=scale)
					cFrame.iq = iq

					cFrame.writeRawFrame(fh)

def basicSignal(fh, stands, nFrames, mode='DRX', filter=6, bits=12, tStart=0):
	"""Generate a collection of frames with a basic test signal for TBW, TBN, 
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
	    * same test signal used in the original lwa_dp_sim
	    
	All modes need to have stands (beams in the case of DRX) and number of
	frames to generate.  TBW also needs to 'bits' keyword set to generate 
	either 12-bit or 4-bit data.  The TBN and DRX frames need the 'filter'
	keyword set to specify the filter width.  In addition, the 'stands' 
	argument is interpreted as beam numbers for DRX."""

	if tStart == 0:
		tStart = time.time()

	if mode == 'TBW':
		__basicTBW(fh, stands, nFrames, mode=mode, filter=filter, bits=bits, tStart=tStart)
	elif mode == 'TBN':
		__basicTBN(fh, stands, nFrames, mode=mode, filter=filter, bits=bits, tStart=tStart)
	elif mode == 'DRX':
		__basicDRX(fh, stands, nFrames, mode=mode, filter=filter, bits=bits, tStart=tStart)
	else:
		raise RuntimeError("Unknown observations mode: %s" % mode)


def __getAntennaArray(station, stands, time, freqs):
	"""Given a LWA station object, a list of stands, an observation time, and
	a list of frequencies in Hz, build an aipy AntennaArray object."""

	return vis.buildSimArray(station, stands, freqs/1e9, jd=astro.unix_to_utcjd(time))


def __getSourceParameters(aa, time, srcs):
	"""Given an aipy AntennaArray object, an observation time, and aipy.src 
	object, return all of the parameters needed for a simulation."""
	
	# Set the time for the array
	aa.set_unixtime(time)

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
		srcAzAlt = aipycoord.top2azalt(top) * 180/numpy.pi
		if srcAzAlt[1] <= 0 or jys.sum() <= 0:
			continue

		## Save values into the source arrays
		srcs_tp.append( top )
		srcs_mp.append( map )
		srcs_jy.append( jys )
		srcs_fq.append( frq )

	# Return the values as a dictionary
	return {'topo': srcs_tp, 'trans': srcs_mp, 'flux': srcs_jy, 'freq': srcs_fq}


def __buildSignals(aa, srcParams, times):
	"""Given an aipy AntennaArray, a list of stand numbers, a dictionary of source parameters, and an 
	array of times in ns, return a numpy array of the simulated signals that is 
	Nstands x Ntimes in shape."""

	# Find out how many stands, srcs, and samples (times) we are working with
	Nstand = len(aa.ants)
	Nsrc = len(srcParams['topo'])
	Ntime = len(times)

	# Define the stand position and cable/signal delay caches to the simulator 
	# move along faster
	dlyCache = uvUtils.SignalCache(aa.get_afreqs()*1e9, applyDispersion=True)

	# Setup a temporary array to hold the signals per source, stand, and time.  
	# This array is complex so that it can accomidate both TBW and TBN data at
	# the same time
	temp = n.zeros((Nsrc, Nstand, Ntime), dtype=n.complex64)

	# Loop over sources and stands to build up the signals
	srcCount = 0
	for topo,trans,flux,freq in zip(srcParams['topo'], srcParams['trans'], srcParams['flux'], srcParams['freq']):
		antCount = 0
		for ant in aa.ants:
			# Zeroth, get the beam response in the direction of the current source for all frequencies
			antResponse = numpy.squeeze( ant.bm_response(tp, pol=pol) )

			# First, do the geometric delay
			geoDelay = ( numpy.dot(trans, ant.pos).transpose() )[2] 

			# Second, do the cable delay
			Delayat1MHz = dlyCache.cableDelay(ant.stand, freq=1.0e6)[0] * 1e9 # s -> ns
			cblDelay = dlyCache.cableDelay(ant.stand) * 1e9 # s -> ns
			# NB: Replace the cable delays below 1 MHz with the 1 MHz value to keep the 
			# delays from blowing up for small f
			cblDelay = numpy.where( freq >= 0.001, cblDelay, Delayat1MHz )

			for a,j,f,d in zip(antResponse, flux, freq, cblDelay):
				factor = a * n.sqrt(j)
				angle = 2*n.pi*f*(t + (d - geoDelay))
				temp[srcCount,antCount,:] += factor*(n.cos(angle) + 1j*n.sin(angle))
			antCount = antCount + 1
		srcCount = srcCount + 1

	# Sum over sources and done
	tdSignals = temp.sum(axis=0)
	return tdSignals


def pointSource(fh, stands, src, nFrames, mode='TBN', filter=7, bits=12, tStart=0):
	"""Generate a collection of frames with a point source signal for TBW
	and TBN.  The point source is specified as a aipy.src object.
	    
	All modes need to have stands (beams in the case of DRX) and number of
	frames to generate.  TBW also needs to `bits' keyword set to generate 
	either 12-bit or 4-bit data.  The TBN frames need the `filter' keyword 
	set to specify the filter width."""

	if tStart == 0:
		tStart = time.time()
	else:
		tStart = float(tStart)

	if mode == 'TBW':
		sampleRate = dp_common.fS
		freqs = (numpy.fft.fftfreq(1024, d=1.0/sampleRate))[1:512]
		aa = __buildAntennaArray(lwa_common.lwa1(), stands, tStart, freqs)
		
		if bits == 12:
			maxValue = 2047
			samplesPerFrame = 400
		else:
			maxValue =  7
			samplesPerFrame = 1200

		nCaptures = int(numpy.ceil(nFrames / 30000.0))
		for capture in range(nCaptures):
			j = 0
			k = 1
			for stand1, stand2 in zip(stands[0::2], stands[1::2]):
				print "Simulating TBW capture %i, stands %i and %i" % (capture+1, stand1, stand2)
				FramesThisBatch = nFrames - capture*30000
				if FramesThisBatch > 30000:
					FramesThisBatch = 30000
				print "-> Frames %i" % FramesThisBatch

				for i in range(FramesThisBatch):
					frameT = tStart + i*samplesPerFrame/sampleRate + 60.0*capture

					# Get the source parameters
					srcParams = __getSourceParameters(aa, frameT, src)

					# Generate the time series response of each signal at each frequency
					t = frameT + n.arange(samplesPerFrame)/sampleRate * 1e9 # s -> ns
					tdSignals = __buildSignals(aa, srcParams, t)

					cFrame = tbw.SimFrame(stand=stand1, frameCount=i+1, dataBits=bits, obsTime=frameT)
					cFrame.xy = numpy.random.randn(2, samplesPerFrame)
					cFrame.xy[0,:] = tdSignals.real[j,:]
					cFrame.xy[1,:] = tdSignals.real[j,:]
					
					cFrame.writeRawFrame(fh)

					cFrame = tbw.SimFrame(stand=stand2, frameCount=i+1, dataBits=bits, obsTime=frameT)
					cFrame.xy = numpy.random.randn(2, samplesPerFrame)
					cFrame.xy[0,:] = tdSignals.real[k,:]
					cFrame.xy[1,:] = tdSignals.real[k,:]
					
					cFrame.writeRawFrame(fh)
			j = j + 2
			k = k + 2

	if mode == 'TBN':
		sampleRate = TBNFilters[filter]
		maxValue = 127
		samplesPerFrame = 512
		freqs = (numpy.fft.fftfreq(samplesPerFrame, d=1.0/sampleRate)) + CentralFreq
		freqs = numpy.fft.fftshift(freqs)
		aa = __buildAntennaArray(lwa_common.lwa1(), stands, tStart, freqs)

		for i in range(nFrames):
			if i % 1000 == 0:
				print "Simulating TBN frame %i" % (i+1)
			frameT = tStart + i*samplesPerFrame/sampleRate

			# Get the source parameters
			srcParams = __getSourceParameters(aa, frameT, src)

			# Generate the time series response of each signal at each frequency
			t = frameT + n.arange(samplesPerFrame)/sampleRate * 1e9 # s -> ns
			tdSignals = __buildSignals(aa, srcParams, t)

			j = 0
			for stand in stands:
				cFrame = tbn.SimFrame(stand=stand, pol=0, frameCount=i+1, obsTime=frameT)
				cFrame.iq = tdSignals[j,:].astype(numpy.singlecomplex)
				
				cFrame.writeRawFrame(fh)

				cFrame = tbn.SimFrame(stand=stand, pol=1, frameCount=i+1, obsTime=frameT)
				cFrame.iq = tdSignals[j,:].astype(numpy.singlecomplex)
				
				cFrame.writeRawFrame(fh)
				
				j = j + 1
