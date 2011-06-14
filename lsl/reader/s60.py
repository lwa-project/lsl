# -*- coding: utf-8 -*-

"""
Python module for reading data in from S60 files.  Unlike TBW or TBN data, 
the S60 data consists of only packed binary values.  Thus, the functions for 
reading in the S60 data do not return objects describing each frame but rather
return numpy arrays of the data.  The two reading functions defined are:

readFrame
  return a complex array containing the data from one UDP packet (734 samples)

readChunk
  return a complex array containing a specified number of samples, i.e., the 
  'chunk' size.

In addition to the two readers, there is also a function, getBandpassModel, 
which reads in Joe Craig's model of the S60 bandpass. See `the S60 data page <http://www.phys.unm.edu/~lwa/dokuwiki/doku.php?id=s60_data#view_and_plot_with_octave>`_
for information about the model.
"""

import numpy

from lsl.common.paths import data as dataPath
from errors import numpyError

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['readFrame', 'readChunk', 'getBandpassModel', 'FrameSize', 'SampleRate', '__version__', '__revision__', '__all__']

# Packet length seems to describe how many bytes of data are in 
# each UDP packet.  1468 bytes == 734 samples.
FrameSize = 1468
# Packets per second seens to indicate how many UDP packets are 
# sent out by the S60 system to be recorded every second.  Thus in
# one second, it would seem that there are 5109*1468 = 7.50 M byes 
# or 3.75 M samples.  
__PacketsPerSecond = 5109
# Sample rate would then be the product of FrameSize and 
# PacketsPerSecond / 2 (two bytes used per I/Q sample)
SampleRate = FrameSize*__PacketsPerSecond / 2.0


def readFrame(filehandle):
	"""
	Function to read in a single S60 frame and return the data as a numpy
	single-precision complex array with 734 elements.
	"""

	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=FrameSize)
	data = numpy.zeros(FrameSize/2, dtype=numpy.csingle)
	if rawData.shape[0] < 2*data.shape[0]:
		raise numpyError()

	data.real = rawData[0::2] - 128.0
	data.imag = rawData[1::2] - 128.0

	return data


def readChunk(filehandle, Chunk=4096):
	"""
	Function to read a certain number of chunks (samples) from the input file-
	handle and return them as a single-precision complex numpy array.  The output 
	is in samples so twice as many bytes need to be read in as samples.
	"""

	rawData = numpy.fromfile(filehandle, dtype=numpy.uint8, count=2*Chunk)
	data = numpy.empty(Chunk, dtype=numpy.csingle)
	if rawData.shape[0] < 2*data.shape[0]:
		raise numpyError()
	
	data.real = rawData[0::2] - 128.0
	data.imag = rawData[1::2] - 128.0

	return data


def getBandpassModel():
	"""
	Read in Joe's model for the S60 badpass that is posted on the wiki.
	"""

	from scipy.io import loadmat

	modelFile = os.path.join(dataPath, 's60-bandpass-model.mat')
	model = loadmat(modelFile, struct_as_record=True, squeeze_me=True)
	return model['y']
