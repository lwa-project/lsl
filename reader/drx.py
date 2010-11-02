#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python module to read in DRX data."""

import os
import sys
import math
import time
import numpy
import array
import struct

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


DRXFrameSize = 4128


class DRXFrameHeader(object):
	"""Class that stores the information found in the header of a DRX 
	frame.  All six fields listed in the DP IDC version H are stored as 
	well as the original binary header data."""
	
	def __init__(self, frameCount=None, drxID=None, secondsCount=None, decimation=None, timeOffset=None, raw=None):
		self.frameCount = frameCount
		self.drxID = drxID
		self.secondsCount = secondsCount
		self.decimation = decimation
		self.timeOffset = timeOffset
		self.raw = raw
	
	def parseID(self):
		"""Parse the DRX ID into a tuple containing the beam (1 through
		4), tunning (1 and 2), and polarization (0 and 1)."""
		
		beam = self.drxID&7
		tune = (self.drxID>>3)&7
		pol  = (self.drxID>>7)&1

		return (beam, tune, pol)
	
	def getSampleRate(self):
		"""Return the sample rate of the data in samples/second."""
		
		sampleRate = 196.0 / self.decimation
		return sampleRate


class DRXFrameData(object):
	"""Class that stores the information found in the data section of a DRX
	frame.  All three fields listed in the DP IDC version H are stored."""
	def __init__(self, timeTag=None, flags=None, iq=None):
		self.timeTag = timeTag
		self.flags = flags
		self.iq = iq


class DRXFrame(object):
	"""Class that stores the information contained within a single DRX 
	frame.  It's properties are DRXFrameHeader and DRXFrameData objects."""
	def __init__(self, header=DRXFrameHeader(), data=DRXFrameData()):
		self.header = header
		self.data = data


class DRXObservingBlock(object):
	"""Class that stores all frames associates with a particular beam at a
	particular time."""
	def __init__(self, x1=DRXFrame(), y1=DRXFrame(), x2=DRXFrame(), y2=DRXFrame()):
		self.x1 = x1
		self.y1 = y1
		self.x2 = x2
		self.y2 = y2
		

def readDRXHeader(filehandle, Verbose=False):
	rawHeader = ''
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	sync4, sync3, sync2, sync1 = struct.unpack(">BBBB", s)
	#print sync1,sync2,sync3,sync4
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	drxID, frameCount = struct.unpack(">B3s", s)
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	secondsCount = struct.unpack(">L", s)
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	decimation, timeOffset = struct.unpack(">HH", s)

	newHeader = DRXFrameHeader()
	newHeader.frameCount = frameCount
	newHeader.drxID = drxID
	newHeader.secondsCount = secondsCount
	newHeader.decimation = decimation
	newHeader.timeOffset = timeOffset
	newHeader.raw = rawHeader

	if Verbose:
		beam, tune, pol = newHeader.parseID()
		print "Header: ", drxID, decimation, timeOffset
		print "  Beam: ", beam
		print "  Tuning: ", tune
		print "  Polarization: ", pol

	return newHeader


def readDRXData(filehandle):
	s = filehandle.read(8)
	timeTag = struct.unpack(">Q", s)
	s = filehandle.read(8)
	flags = struct.unpack(">Q", s)
		
	s = filehandle.read(4096)
	rawData = struct.unpack(">1024L", s)
	data = numpy.zeros(4096, dtype=numpy.complex64)
	
	data[0::4] = [complex((word>>28)&15, (word>>24)&15) for word in rawData]
	data[1::4] = [complex((word>>20)&15, (word>>16)&15) for word in rawData]
	data[2::4] = [complex((word>>12)&15, (word>>8)&15) for word in rawData]
	data[3::4] = [complex((word>>4)&15, word&15) for word in rawData]
	
	newData = DRXFrameData()
	newData.timeTag = timeTag[0]
	newData.flags = flags
	newData.iq = data - complex(8.0, 8.0)

	return newData


def readDRXFrame(filehandle):
	"""Function to read in a single DRX frame (header+data) and store the 
	contents as a DRXFrame object.  This function wraps readerHeader and 
	readData."""
	
	try:
		# Read header
		hdr = readDRXHeader(filehandle)
		# Read data
		dat = readDRXData(filehandle)
	except IOError (errorNumber, errorString):
		print 'file I/O error {0}: {1}' % (errorNumber, errorString) 
		raise
	
	# Create the new frame object and return
	newFrame = DRXFrame()
	newFrame.header = hdr
	newFrame.data = dat

	return newFrame


def readDRXBlock(filehandle):
	"""Function to read in a single DRX block (four frames) and store the 
	contents as a DRXObservingBlock object.  This function wraps 
	readDRXFrame."""
	
	try:
		x1 = readDRXFrame(filehandle)
		y1 = readDRXFrame(filehandle)
		x2 = readDRXFrame(filehandle)
		y2 = readDRXFrame(filehandle)
	except IOError (errorNumber, errorString):
		print 'file I/O error {0}: {1}' % (errorNumber, errorString) 
		raise

	block = DRXObservingBlock(x1=x1, y1=y1, x2=x2, y2=y2)
	return block


def getBeamCount(filehandle):
	"""Find out how many beams are present by examining the next 32 DRX
	records.  Return the number of beams found."""

	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Build up the list-of-lists that store ID codes and loop through 32
	# frames.  In each case, parse pull the DRX ID, extract the beam number, 
	# and append the DRX ID to the relevant beam array if it is not already 
	# there.
	beams = []
	for i in range(32):
		cFrame = readDRXFrame(filehandle)
		cID = cFrame.header.drxID
		beam = cID&7
		if beam not in beams:
			beams.append(beam)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)

	# Return the number of beams found
	return len(beams)


def getFramesPerObs(filehandle):
	"""Find out how many frames are present per beam by examining the next 
	32 DRX records.  32 elements are read so that this function can be 
	called at anytime during the reading and the correct values returned.  
	Return the number of frames per observations as a four-element tuple, 
	one for each beam."""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Build up the list-of-lists that store ID codes and loop through 32
	# frames.  In each case, parse pull the DRX ID, extract the beam number, 
	# and append the DRX ID to the relevant beam array if it is not already 
	# there.
	idCodes = [[], [], [], []]
	for i in range(32):
		cFrame = readDRXFrame(filehandle)
		cID = cFrame.header.drxID
		beam = cID&7
		if cID not in idCodes[beam-1]:
			idCodes[beam-1].append(cID)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Get the length of each beam list and return them as a tuple
	return (len(idCodes[0]), len(idCodes[1]), len(idCodes[2]), len(idCodes[3]))
	

def averageObservations(Observations):
	"""Given a list of ObservingBlock objects, average the observations 
	together on a per tuning, per polarization basis.  A new ObsevingBlock 
	object is returned that contains the averaged data."""
	
	newBlock = Observations[0]
	tempx1 = newBlock.x1.data.iq * complex(0.0, 0.0)
	tempy1 = newBlock.y1.data.iq * complex(0.0, 0.0)
	tempx2 = newBlock.x2.data.iq * complex(0.0, 0.0)
	tempy2 = newBlock.y2.data.iq * complex(0.0, 0.0)
	
	obsCount = 0
	for Observation in Observations:
		tempx1 = tempx1 + Observation.x1.data.iq
		tempy1 = tempy1 + Observation.y1.data.iq
		tempx2 = tempx2 + Observation.x2.data.iq
		tempy2 = tempy2 + Observation.y2.data.iq
		obsCount = obsCount + 1
	
	tempx1 = tempx1 / float(obsCount)
	tempy1 = tempy1 / float(obsCount)
	tempx2 = tempx2 / float(obsCount)
	tempy2 = tempy2 / float(obsCount)
			
	newBlock.x1.data.iq = tempx1
	newBlock.y1.data.iq = tempy1
	newBlock.x2.data.iq = tempx2
	newBlock.y2.data.iq = tempy2
			
	return newBlock


def averageObservations2(Observations, timeAvg=1, chanAvg=1):
	"""Given a list of ObservingBlock objects, average the observations 
	together on a per tuning, per polarization basis.  A new ObsevingBlock 
	object is returned that contains the averaged data."""
	
	caBlocks = []
	if chanAvg > 1:
		for Observation in Observations:
			currBlock = Observation
			for i in range(len(Observation.x1.data.iq)/chanAvg):
				currBlock.x1.data.iq[i] = numpy.average( Observation.x1.data.iq[i*chanAvg:(i+1)*chanAvg] )
				currBlock.y1.data.iq[i] = numpy.average( Observation.y1.data.iq[i*chanAvg:(i+1)*chanAvg] )
				currBlock.x2.data.iq[i] = numpy.average( Observation.x2.data.iq[i*chanAvg:(i+1)*chanAvg] )
				currBlock.y2.data.iq[i] = numpy.average( Observation.y2.data.iq[i*chanAvg:(i+1)*chanAvg] )
			currBlock.x1.data.iq = currBlock.x1.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			currBlock.y1.data.iq = currBlock.y1.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			currBlock.x2.data.iq = currBlock.x2.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			currBlock.y2.data.iq = currBlock.y2.data.iq[0:len(Observation.x1.data.iq)/chanAvg]
			caBlocks.append(currBlock)
	else:
		caBlocks = Observations
			
	taBlocks = []
	if timeAvg > 1:
		for i in range(len(Observations)/timeAvg):
			taBlocks.append( averageObservations(caBlocks[i*timeAvg:(i+1)*timeAvg]) )
	else:
		taBlocks = caBlocks
			
	return taBlocks


def writefits(data):
	import pyfits

	# Data
	nTime = data.shape[0]
	c1 = pyfits.Column(name='data', format='4096C', array=numpy.zeros((nTime, 4096), dtype=numpy.complex64))

	# Frequency
	c2 = pyfits.Column(name='ctype1', format='A8')
	c3 = pyfits.Column(name='crval1', format='1D')
	c4 = pyfits.Column(name='cdelt1', format='1D')
	c5 = pyfits.Column(name='crpix1', format='1D')

	# RA/Dec
	c6 = pyfits.Column(name='ctype2', format='A8')
	c7 = pyfits.Column(name='crval2', format='1D', unit='degrees')
	c8 = pyfits.Column(name='ctype3', format='A8')
	c9 = pyfits.Column(name='crval3', format='1D', unit='degrees')

	# Tunning/Polarization
	c10 = pyfits.Column(name='ctype4', format='A8')
	c11 = pyfits.Column(name='crval4', format='1D')
	c12 = pyfits.Column(name='ctype5', format='A8')
	c13 = pyfits.Column(name='crval5', format='1D')

	# Time
	c14 = pyfits.Column(name='stime', format='1D')

	# Obs. Mode/Channels
	c15 = pyfits.Column(name='obsmode', format='A3')
	c16 = pyfits.Column(name='nchan', format='1D')

	colDefs = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16])

	prim = pyfits.PrimaryHDU()
	sdfits11 = pyfits.new_table(colDefs)
	sdfits12 = pyfits.new_table(colDefs)
	sdfits21 = pyfits.new_table(colDefs)
	sdfits22 = pyfits.new_table(colDefs)

	for i in range(data.shape[0]):
		print "%i of %s" % (i, data.shape[0])
		print sdfits11.data.field('data').shape, numpy.squeeze(data[i,:]).shape, numpy.reshape(data[i,:].real, (1,4096)).shape
		sdfits11.data.field('data')[i,:] = numpy.reshape(data[i,:], (1,4096))
	sdfits11.data.field('ctype1')[:] = 'FREQ'
	sdfits11.data.field('crval1')[:] = 0
	sdfits11.data.field('cdelt1')[:] = 1
	sdfits11.data.field('crpix1')[:] = 1
	sdfits11.data.field('crval2')[:] = 180.0
	sdfits11.data.field('crval3')[:] = -35.0
	sdfits11.data.field('obsmode')[:] = 'DRX'
	sdfits11.data.field('nchan')[:] = 4096

	print sdfits11.header.ascardlist()

	masterList = pyfits.HDUList([prim, sdfits11, sdfits12, sdfits21, sdfits22])

	masterList.writeto('test-sdfits.fits')
	pyfits.info('test-sdfits.fits')
	

def main(args):
	nSamples = os.path.getsize(args[0]) / DRXFrameSize
	print "Samples in file: ", nSamples
	fh = open(args[0], "rb", buffering=DRXFrameSize)
	nFpO = getFramesPerObs(fh)
	nBeams = getBeamCount(fh)
	print "Beams: ", nBeams
	print "Frames per Observations: ", nFpO
	blockBuffer = []
	blocks = []

	tStart = time.time()

	nSamples = 256

	fig = plt.figure()

	for i in range(0,nSamples):
		currBlock = readDRXBlock(fh)
		blockBuffer.append(currBlock)

		if len(blockBuffer) == 16:
			avgBlock = averageObservations(blockBuffer)
			#avgBlock = averageObservations2(blockBuffer, timeAvg=16, chanAvg=2)
			blocks.append(avgBlock)
			blockBuffer = []
	
	nChan = blocks[0].x1.data.iq.shape[0]
	outSpec = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec[row,:] = block.x1.data.iq
	outSpec2 = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec2[row,:] = block.y1.data.iq
	outSpec3 = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec3[row,:] = block.x2.data.iq
	outSpec4 = numpy.zeros((nSamples/16, nChan), dtype=numpy.complex64)
	for row,block in zip(range(nSamples),blocks):
		outSpec4[row,:] = block.y2.data.iq

	tEnd = time.time()
	print 'Read %i frames in %0.3f s (%0.1f frames/s)' % (4*nSamples, (tEnd-tStart), 4*nSamples/(tEnd-tStart))
	
	writefits(outSpec)

	ax = fig.add_subplot(221)
	dB = outSpec - outSpec.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 1, Pol. 0')
	ax.set_ylabel('Channel')
	ax.axis('auto')
	
	ax = fig.add_subplot(222)
	dB = outSpec2 - outSpec2.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 1, Pol. 1')
	ax.axis('auto')

	ax = fig.add_subplot(223)
	dB = outSpec3 - outSpec3.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 2, Pol. 0')
	ax.set_xlabel('Time')
	ax.set_ylabel('Channel')
	ax.axis('auto')

	ax = fig.add_subplot(224)
	dB = outSpec4 - outSpec4.mean(axis=0)
	dB = numpy.log10( (dB*dB.conj()).real )*10.0
	
	ax.imshow(numpy.transpose(dB), origin='lower')
	ax.set_title('Tuning 2, Pol. 1')
	ax.set_xlabel('Time')
	ax.axis('auto')

	plt.show()
	fig.savefig("readDRX.png")

	fh.close()

if __name__ == "__main__":
	main(sys.argv[1:])

