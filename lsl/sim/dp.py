# -*- coding: utf-8 -*-

"""Module to simulate DP observations."""

import numpy

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
	if mode == 'TBW':
		sampleRate = dp_common.fS
		
		if nFrames > 30000:
			nFrames = 30000
		if bits == 12:
			maxValue = 2047
			samplesPerFrame = 400
		else:
			maxValue =  7
			samplesPerFrame = 1200

		for stand1,stand2 in zip(stands[0::2], stands[1::2]):
			for i in range(nFrames):
				t = i*samplesPerFrame*sampleRate
				
				cFrame = tbw.SimFrame(stand=stand, frameCount=i, dataBits=bits, obsTime=t)
				cFrame.xy = numpy.random.randn((2,samplesPerFrame)) * maxValue/4.0
				
				cFrame.writeFrame(fh)
	
	if mode == 'TBN':
		sampleRate = TBNFilters[filter]
		maxValue = 127
		samplesPerFrame = 512
		
		for i in range(nFrames):
			t = i*samplesPerFrame*sampleRate
			for stand in stands:
				for pol in [0, 1]:
					cFrame = tbn.SimFrame(stand=stand, pol=pol, frameCount=i, obsTime=t)
					temp = numpy.zeros(samplesPerFrame, dtype=numpy.singlecomplex)
					temp.real = numpy.random.randn(samplesPerFrame) * maxValue/4.0
					temp.imag = numpy.random.randn(samplesPerFrame) * maxValue/3.0
					cFrame.data.iq = temp
				
					cFrame.writeFrame(fh)
	
	if mode == 'DRX':
		sampleRate = DRXFilters[filter]
		maxValue = 7
		samplesPerFrame = 4096
		
	