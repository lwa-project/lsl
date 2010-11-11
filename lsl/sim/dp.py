# -*- coding: utf-8 -*-

"""Module to simulate DP observations."""

import functools
import itertools
import math
import numpy
import scipy.stats

from lsl.common import dp as dp_common
from lsl.sim import tbw
from lsl.sim import tbn
from lsl.sim import drx
from lsl.reader.tbn import filterCodes as TBNFilters
from lsl.reader.drx import filterCodes as DRXFilters

__version__ = '0.1'
__revision__ = '$ Revision: 1 $'
__all__ = ['__version__', '__revision__', '__all__']


def basicSignal(fh, stands, nFrames, mode='DRX', filter=6, bits=12):
	"""Generate a collection of frames with a basic test signal for TBW, TBN, 
	and DRX.  The signls are:
	  TBW:
	    noise + 40 MHz signal for x-pol.; noise + 60 MHz signal for y-pol.-> odd stands
	    noise + 30 MHz signal for x-pol.; noise + 50 MHz signal for ypol. -> even stands
	  TBN:
	    noise + 
	  DRX:
	    same test signal used in the original lwa_dp_sim"""

	if mode == 'TBW':
		sampleRate = dp_common.fS
		
		if bits == 12:
			maxValue = 2047
			samplesPerFrame = 400
		else:
			maxValue =  7
			samplesPerFrame = 1200

		nCaptures = int(math.ceil(nFrames / 30000.0))
		for capture in range(nCaptures):
			for stand1, stand2 in zip(stands[0::2], stands[1::2]):
				print "Simulating capture %i, stands %i and %i" % (capture+1, stand1, stand2)
				FramesThisBatch = nFrames - capture*30000
				if FramesThisBatch > 30000:
					FramesThisBatch = 30000
				print "-> Frames %i" % FramesThisBatch
				for i in range(FramesThisBatch):
					t = i*samplesPerFrame/sampleRate + 60.0*capture
					
					cFrame = tbw.SimFrame(stand=stand1, frameCount=i+1, dataBits=bits, obsTime=t)
					cFrame.xy = numpy.random.randn(2, samplesPerFrame)
					cFrame.xy[0,:] *= maxValue/15.0
					cFrame.xy[0,:] += maxValue*numpy.cos(2*numpy.pi*0.2041*numpy.arange(samplesPerFrame))
					cFrame.xy[1,:] *= maxValue/15.0
					cFrame.xy[1,:] += maxValue*numpy.cos(2*numpy.pi*0.3061*numpy.arange(samplesPerFrame))
					
					cFrame.writeRawFrame(fh)

					cFrame = tbw.SimFrame(stand=stand2, frameCount=i+1, dataBits=bits, obsTime=t)
					cFrame.xy = numpy.random.randn(2, samplesPerFrame)
					cFrame.xy[0,:] *= maxValue/15.0
					cFrame.xy[0,:] += maxValue*numpy.cos(2*numpy.pi*0.1531*numpy.arange(samplesPerFrame))
					cFrame.xy[1,:] *= maxValue/15.0
					cFrame.xy[1,:] += maxValue*numpy.cos(2*numpy.pi*0.2551*numpy.arange(samplesPerFrame))
					
					cFrame.writeRawFrame(fh)
	
	if mode == 'TBN':
		sampleRate = TBNFilters[filter]
		maxValue = 127
		samplesPerFrame = 512
		
		for i in range(nFrames):
			if i % 1000 == 0:
				print "Simulating frame %i" % (i+1)
			t = 1.0*i*samplesPerFrame/sampleRate
			for stand in stands:
				cFrame = tbn.SimFrame(stand=stand, pol=0, frameCount=i+1, obsTime=t)
				temp = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
				temp.real = numpy.random.randn(samplesPerFrame)
				temp.real *= maxValue/15.0
				temp.real += maxValue*numpy.cos(2*numpy.pi/100.0*(i*samplesPerFrame+numpy.arange(samplesPerFrame)))
				temp.imag = numpy.random.randn(samplesPerFrame)
				temp.imag *= maxValue/15.0
				temp.imag += maxValue*numpy.sin(2*numpy.pi/100.0*(i*samplesPerFrame+numpy.arange(samplesPerFrame)))
				cFrame.iq = temp
				
				cFrame.writeRawFrame(fh)

				cFrame = tbn.SimFrame(stand=stand, pol=1, frameCount=i+1, obsTime=t)
				temp = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
				temp.real = numpy.random.randn(samplesPerFrame)
				temp.real *= maxValue/15.0
				temp.real += maxValue*numpy.cos(2*numpy.pi/200.0*(i*samplesPerFrame+numpy.arange(samplesPerFrame)))
				temp.imag = numpy.random.randn(samplesPerFrame)
				temp.imag *= maxValue/15.0
				temp.imag += maxValue*numpy.sin(2*numpy.pi/200.0*(i*samplesPerFrame+numpy.arange(samplesPerFrame)))
				cFrame.iq = temp
				
				cFrame.writeRawFrame(fh)
	
	if mode == 'DRX':
		sampleRate = DRXFilters[filter]
		maxValue = 7
		samplesPerFrame = 4096
		
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
				print "Simulating frame %i" % i
			t = 1.0*i*samplesPerFrame/sampleRate
			for beam in beams:
				for tune in [1, 2]:
					for pol in [0, 1]:
						cFrame = drx.SimFrame(beam=beam, tune=tune, pol=pol, frameCount=i+1, filterCode=filter, timeOffset=0, obsTime=t, flags=0)

						iq = numpy.zeros(samplePerFrame, dtype=numpy.singlecomplex)
						for chan, scale in itertools.izip(numpy.arange(samplesPerFrame), signal[2*(tune-1)+pol]):
							iq[chan].real = norm(scale=scale)
							iq[chan].imag = norm(scale=scale)
						cFrame.iq = iq

						cFrame.writeRawFrame(fh)
							