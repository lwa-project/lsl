# -*- coding: utf-8 -*-

"""Python module for creating creating and writing simulated S60 frames to 
a file."""

import numpy

from lsl.common import dp as dp_common
from lsl.reader import s60
from errors import *

__version__ = '0.1'
__revision__ = '$ Revision: 2 $'
__all__ = ['frame2frame', 'chunk2frame', '__version__', '__revision__', '__all__']


def frame2frame(s60Data):
	"""Convert a complex numpy array with the same length as a S60 frame to a 
	raw S60 frame."""

	# The raw frame
	rawFrame = numpy.zeros(s60.FrameSize, dtype=numpy.uint8)

	# Part 1: Sanity checks
	## Compared to the DP data formats (TBW, TBN, and DRX), the S60 data is 
	## relatively barebones.  These checks replace the SimFrame.isValid()
	## function in the DP data formats.
	### Does the data type make sense?
	if s60Data.dtype.kind != 'c':
		raise invalidDataType()
	### Does the data length make sense?
	if s60Data.shape[0] != s60.FrameSize/2:
		raise invalidDataSize()
	
	# Part 2: The data
	i = s60Data.real + 128.0
	q = s60Data.imag + 128.0
	## Round, clip, and convert to unsigned numbers
	i = i.round()
	i = i.clip(0, 255)
	i = i.astype(numpy.uint8)
	q = q.round()
	q = q.clip(0, 255)
	q = q.astype(numpy.uint8)

	rawFrame[0::2] = i
	rawFrame[1::2] = q

	return rawFrame


def chunk2frame(s60Data):
	"""Convert a complex numpy array to a raw S60 frame."""

	# The raw frame
	rawFrame = numpy.zeros(2*len(s60Data), dtype=numpy.uint8)

	# Part 1: Sanity checks
	## Compared to the DP data formats (TBW, TBN, and DRX), the S60 data is 
	## relatively barebones.  These checks replace the SimFrame.isValid()
	## function in the DP data formats.
	### Does the data type make sense?
	if s60Data.dtype.kind != 'c':
		raise invalidDataType()
	### Does the data length make sense?
	if s60Data.shape[0] < 1:
		raise invalidDataSize()
	
	# Part 2: The data
	i = s60Data.real + 128.0
	q = s60Data.imag + 128.0
	## Round, clip, and convert to unsigned numbers
	i = i.round()
	i = i.clip(0, 255)
	i = i.astype(numpy.uint8)
	q = q.round()
	q = q.clip(0, 255)
	q = q.astype(numpy.uint8)

	rawFrame[0::2] = i
	rawFrame[1::2] = q

	return rawFrame
