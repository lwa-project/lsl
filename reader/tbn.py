#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy
import array
import struct


TBNFrameSize = 1024


class TBNFrameHeader(object):
	def __init__(self, frameCount=None, secondsCount=None, tbnID=None,  raw=None):
		self.frameCount = frameCount
		self.secondsCount = secondsCount
		self.tbnID = tbnID
		self.raw = raw
		
	def parseID(self):
		if self.tbnID % 2 == 0:
			stand = self.tbnID / 2
			pol = 1
		else:
			stand = (self.tbnID-1)/2 + 1
			pol = 0
			
		return (stand, pol)


class TBNFrameData(object):
	def __init__self(self, timeTag=None, iq=None):
		self.timeTag = timeTage
		self.iq = iq


class TBNFrame(object):
	def __init__(self, header=TBNFrameHeader(), data=TBNFrameData()):
		self.header = header
		self.data = data


class TBNObservingBlock(object):
	def __init__(self, x=[], y=[]):
		self.x = x
		self.y = y


def readTBNHeader(filehandle, Verbose=True):
	rawHeader = ''
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	sync4, sync3, sync2, sync1 = struct.unpack(">cccc", s)
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	m5cID, frameCount = struct.unpack(">B3s", s)
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	secondsCount = struct.unpack(">L", s)
	s = filehandle.read(4)
	rawHeader = rawHeader + s
	tbnID, junk = struct.unpack(">HH", s)

	newHeader = TBNFrameHeader()
	newHeader.frameCount = frameCount
	newHeader.secondsCount = secondsCount
	newHeader.tbnID = tbnID
	newHeader.raw = rawHeader

	return newHeader


def readTBNData(filehandle):
	s = filehandle.read(8)
	timeTag = struct.unpack(">Q", s)
	
	s = filehandle.read(1000)
	rawData = struct.unpack(">250L", s)
	data = numpy.zeros(500, dtype=numpy.complex64)
	
	data[0::2] = [complex((word>>24)&255, (word>>16)&255) for word in rawData]
	data[1::2] = [complex((word>>8)&255, word&255) for word in rawData]
	
	newData = TBNFrameData()
	newData.timeTag = timeTag[0]
	newData.iq = data - complex(128.0, 128.0)

	return newData


def readTBNFrame(filehandle):
	hdr = readTBNHeader(filehandle)
	dat = readTBNData(filehandle)
	
	newFrame = TBNFrame()
	newFrame.header = hdr
	newFrame.data = dat

	return newFrame


def getFramesPerObs(filehandle):
	"""Find out how many frames are present per observation by examining 
	the next 512 TBN frames.  512 elements are read so that this function 
	can be called at anytime during the reading and the correct values 
	returned.  Return the number of frames per observations as a two-
	element tuple, one for each polarization."""
	
	# Save the current position in the file so we can return to that point
	fhStart = filehandle.tell()
	
	# Build up the list-of-lists that store ID codes and loop through 512
	# frames.  In each case, parse pull the TBN ID, extract the stand 
	# number, and append the stand number to the relevant polarization array 
	# if it is not already there.
	idCodes = [[], []]
	for i in range(512):
		cFrame = readTBNFrame(filehandle)
		cID = cFrame.header.parseID()
		if cID not in idCodes[cPol]:
			idCodes[cPol].append(cID)
			
	# Return to the place in the file where we started
	filehandle.seek(fhStart)
	
	# Get the length of each beam list and return them as a tuple
	return (len(idCodes[0]), len(idCodes[1]))


def main(args):
	nSamples = os.path.getsize(args[0]) / TBNFrameSize
	print "Samples in file: ", nSamples
	fh = open(args[0], "rb", buffering=TBNFrameSize)
	print getFramesPerObs(fh)
	
	for i in range(0,2048):
		cFrame = readTBNFrame(fh)
		print cFrame.header.parseID(),cFrame.data.timeTag
	
	fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])

