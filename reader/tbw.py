#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python module to read in TBW data."""

import os
import sys
import numpy
import array
import struct


TBWFrameSize = 1224


class TBWFrameHeader(object):
	def __init__(self, frameCount=None, secondsCount=None, tbwID=None,  raw=None):
		self.frameCount = frameCount
		self.secondsCount = secondsCount
		self.tbwID = tbwID
		self.timeOffset = timeOffset
		self.raw = raw

	def getBits(self):
		# I don't think that this is correct as-is
		bits = (self.twbID>>1)&1
		if bits == 1:
			dataBits = 12
		else:
			dataBits = 4

		return dataBits


class TBWFrameData(object):
	def __init__self(self, timeTag=None, samples=400, xy=None):
		self.timeTag = timeTage
		self.xy = xy


class TBWFrame(object):
	def __init__(self, header=TBWFrameHeader(), data=TBWFrameData()):
		self.header = header
		self.data = data


def readTBWHeader(filehandle, Verbose=True):
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
	tbwID, junk = struct.unpack(">HH", s)

	newHeader = TBWFrameHeader()
	newHeader.frameCount = frameCount
	newHeader.secondsCount = secondsCount
	newHeader.tbwID = tbwID

	return newHeader


def readTBWData12(filehandle, dataBits=12):
	s = filehandle.read(8)
	timeTag = struct.unpack(">Q", s)
	
	s = filehandle.read(1200)
	unalignedData = struct.unpack(">1200B", s)
	rawData = [part[0]*65536+part[1]*256+part[2] for part in zip(unalignedData[0::3], unalignedData[1::3], unalignedData[2::3])]
	data = numpy.zeros((400, 2), dtype=numpy.float32) # ???
	
	data[0::2,0] = [float((word24>>18)&63) for word24 in rawData]
	data[0::2,1] = [float((word24>>12)&63) for word24 in rawData]
	data[1::2,0] = [float((word24>>6)&63) for word24 in rawData]
	data[1::2,1] = [float(word24&63) for word24 in rawData]

	newData = TBWFrameData()
	newData.timeTag = timeTag[0]
	newData.xy = data - 32.0

	return newData
	
	
def readTWBData4(filehandle):
	s = filehandle.read(8)
	timeTag = struct.unpack(">Q", s)
	
	s = filehandle.read(1200)
	rawData = struct.unpack(">300L", s)
	data = numpy.zeros((1200, 2), dtype=numpy.float32)
	
	data[0::4,0] = [float((word>>28)&15) for word in rawData]
	data[0::4,1] = [float((word>>24)&15) for word in rawData]
	data[1::4,0] = [float((word>>20)&15) for word in rawData]
	data[1::4,1] = [float((word>>16)&15) for word in rawData]
	data[2::4,0] = [float((word>>12)&15) for word in rawData]
	data[2::4,1] = [float((word>>8)&15) for word in rawData]
	data[3::4,0] = [float((word>>4)&15) for word in rawData]
	data[3::4,1] = [float(word&15) for word in rawData]
	
	newData = TBWFrameData()
	newData.timeTag = timeTag[0]
	newData.xy = data - 8.0
	
	return newData


def readTBWFrame(filehandle):
	hdr = readTBWHeader(filehandle)
	if hdr.getBits() == 12:
		dat = readTBWData12(filehandle)
	else:
		dat = readTBWData4(filehandle)
	
	newFrame = TBWFrame()
	newFrame.header = hdr
	newFrame.data = dat

	return newFrame
