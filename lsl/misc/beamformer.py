# -*- coding: utf-8 -*-

"""Module to allow for post-acquisition delay-and-sum beamforming with TBW 
and TBN data using both integer sample delays for time series data 
(intDelayAndSum) and fractional delays for time series data transformed to 
the frequency domain (fftDelayAndSum).  fftDelayAndSum is still under 
development."""

import os
import sys
import aipy
import math
import numpy

from lsl.common.paths import data as dataPath
from lsl.common.constants import c
from lsl.common import dp as dp_common
from lsl.correlator import uvUtils

__version__ = '0.2'
__revision__ = '$ Revision: 5 $'
__all__ = ['BeamformingError', 'calcDelay', 'intDelayAndSum', 'intBeamShape', 'fftDelayAndSum', 'fftBeamShape', '__version__', '__revision__', '__all__']

class BeamformingError(Exception):
	"""Base class for all exceptions in this file."""

	def __init__(self, strerror, errno='-1'):
		self.errno = errno
		self.strerror = strerror
		self.filename = None
		self.args = (errno, strerror)

	def __str__(self):
		return "%s" % self.strerror


def __loadStandResponse(freq=49.0e6):
	"""Create an aipy.amp.beam object that holds the reponse for a single 
	isolated stand.  The stand response is based on NEC4 models at a variety
	of frequencies within the LWA frequency range."""

	# Read in the spherical harmonic representation of the beam distributed with
	# LSL
	dd = numpy.load(os.path.join(dataPath, 'beam-shape.npz'))
	coeffs = dd['coeffs']

	# Calculate how many harmonics are stored in the data set and reorder the data
	# to aipy's liking
	deg = coeffs.shape[0]-1
	lmax = int((math.sqrt(1+8*coeffs.shape[1])-3)/2)
	beamShapeDict = {}
	for i in range(deg+1):
		beamShapeDict[i] = numpy.squeeze(coeffs[-1-i,:])

	# Build the beam object and done
	return aipy.amp.BeamAlm(numpy.array([freq/1e9]), lmax=lmax, mmax=lmax, deg=deg, nside=128, coeffs=beamShapeDict)


def calcDelay(stands, freq=49.0e6, azimuth=0.0, elevation=90.0):
	"""Calculate the time delays for delay-and-sum beam forming a collection of 
	stands looking in at a particular azimuth and elevation (both in degrees).  
	A numpy array of the needed delays in seconds is returned."""

	# Make sure the pointing coordinates make sense
	if elevation < 0 or elevation > 90:
		raise BeamformingError("Pointing elevation (%.2f deg) is out of range [0, 90]" % elevation)
	if azimuth < 0 or azimuth > 360:
		raise BeamformingError("Pointing azimuth (%.2f deg) is out of range [0, 360]" % azimuth)

	# Get the positions of the stands and compute the mean center of the array
	xyz = uvUtils.getXYZ(stands)
	arrayX = xyz[:,0].mean()
	arrayY = xyz[:,1].mean()
	arrayZ = xyz[:,2].mean()

	# Build up a unit vector that points in the direction azimuth,elevation
	rAz = azimuth*numpy.pi/180.0
	rEl = elevation*numpy.pi/180.0
	source = numpy.array([numpy.cos(rEl)*numpy.sin(rAz), 
						numpy.cos(rEl)*numpy.cos(rAz), 
						numpy.sin(rEl)])

	# Compute the stand positions relative to the average and loop over stands
	# to compute the time delays in seconds
	arrayXYZ = xyz - numpy.array([arrayX, arrayY, arrayZ])
	delays = numpy.zeros((len(stands),))
	for i in list(range(len(stands))):
		delays[i] = numpy.dot(source, arrayXYZ[i,:]) / c

	# Get the cable delays for each stand and add that in as well
	dlyCache = uvUtils.CableCache(freq, applyDispersion=True)
	for i in list(range(len(stands))):
		delays[i] = dlyCache.cableDelay(stands[i]) - delays[i]

	# Done
	return delays


def intDelayAndSum(stands, data, sampleRate=dp_common.fS, azimuth=0.0, elevation=90.0):
	"""Given a list of stands and a data stream of the form stands x times, 
	delay and sum the data stream into one beam.  The delays applied are 
	integer sample delays.  Return a numpy array of the time series data 
	associated with the formed beam."""

	# Get the stand delays and convert the delay times from seconds to samples
	delays = calcDelay(stands, azimuth=azimuth, elevation=elevation)
	delays = numpy.round(delays*sampleRate).astype(numpy.int16)

	# Make the delays into something meaningful for the shifting of the data 
	# streams
	delays = delays.max() - delays

	# Delay and sum by looping over stands inside of looping over times
	output = numpy.zeros((data.shape[1]-delays.max()), dtype=data.dtype)
	for s in list(range(len(stands))):
		start = delays[s]
		stop = data.shape[1] - delays.max() + start
		output = output + data[s,start:stop]

	# Done
	return output


def intBeamShape(stands, sampleRate=dp_common.fS, azimuth=0.0, elevation=90.0, progress=False):
	"""Given a list of stands, compute the on-sky response of the delay-and-sum
	scheme implemented in intDelayAndSum.  A 360x90 numpy array spaning azimuth
	and elevation is returned."""

	# Get the stand delays and convert the delay times from seconds to samples
	delays = calcDelay(stands, freq=49.0e6, azimuth=azimuth, elevation=elevation)
	delays = numpy.round(delays*sampleRate).astype(numpy.int16)

	# Build up a base time array, load in the cable delays, and get the stand 
	# positions for geometric delay calculations.
	t = numpy.arange(0,1000)/sampleRate
	dlyCache = uvUtils.CableCache(49.0e6, applyDispersion=True)
	xyz = uvUtils.getXYZ(stands)
	arrayX = xyz[:,0].mean()
	arrayY = xyz[:,1].mean()
	arrayZ = xyz[:,2].mean()
	arrayXYZ = xyz - numpy.array([arrayX, arrayY, arrayZ])

	# Load in the respoonse of a single isolated stand
	standBeam = __loadStandResponse(freq=49.0e6)

	# Build the output array and loop over all azimuths and elevations
	output = numpy.zeros((360,90))
	for az in list(range(360)):
		rAz = az*numpy.pi/180.0
		for el in list(range(90)):
			rEl = el*numpy.pi/180.0

			# Display the progress meter if the `progress' keyword is set to True.  The
			# progress meter displays a `.' every 2% complete and the percentages every
			# 10%.  At 100%, `Done' is displayed.
			if progress:
				fracDone = (az*90+el) / 32400.0 * 100
				if fracDone % 10 == 0 and round(fracDone,2) != 100:
					sys.stdout.write("%i%%" % fracDone)
				elif round(fracDone,2) == 100:
					sys.stdout.write("Done\n")
				elif round(fracDone,3) % 2 == 0:
					sys.stdout.write(".")
				else:
					pass
				sys.stdout.flush()

			# Unit vector for the currect on-sky location
			currPos = numpy.array([numpy.cos(rEl)*numpy.sin(rAz), 
							numpy.cos(rEl)*numpy.cos(rAz), 
							numpy.sin(rEl)])
			# Stand response in this direction
			currResponse = standBeam.response(aipy.coord.azalt2top(numpy.concatenate([[rAz], [rEl]])))[0][0]

			# Loop over stands to build the simulated singnals
			signals = numpy.zeros((len(stands), 1000))
			for i in list(range(len(stands))):
				currDelay = dlyCache.cableDelay(stands[i]) - numpy.dot(currPos, arrayXYZ[i,:]) / c
				signals[i,:] = currResponse * numpy.cos(2*numpy.pi*49.0e6*(t + currDelay))

			# Beamform with delay-and-sum and store the RMS result
			beam = intDelayAndSum(stands, signals, sampleRate=sampleRate, azimuth=azimuth, elevation=elevation)
			output[az,el] = numpy.sqrt((beam**2).mean())

	# Done
	return output


def fftDelayAndSum(stands, data, sampleRate=dp_common.fS, LFFT=256, CentralFreq=49.0e6, azimuth=0.0, elevation=90.0):
	"""Given a list of stands and a data stream of the form stands x times, 
	delay and sum the data stream into one beam.  The delays first applied 
	as integer sample delays.  Then, the data are transformed to the 
	frequency domain and the remainder of the delay is applied as a phase 
	rotation.  A numpy array of the frequency-domain data over time is
	returned."""

	# Get the stand delays in seconds
	delays = calcDelay(stands, azimuth=azimuth, elevation=elevation)

	# Make the delays into something meaningful for the shifting of the data 
	# streams.  Then, get the integer delay and sub-sample delay for each stand
	delays = delays.max() - delays
	intDelays = numpy.round(delays*sampleRate).astype(numpy.int16)
	frcDelays = delays - intDelays/sampleRate

	# Compute the frequencies of the FFT
	freq = numpy.fft.fftfreq(LFFT, d=1.0/sampleRate)
	if data.dtype.kind == 'c':
		freq = freq + CentralFreq

	# Loop over stands and FFT sections to compute the formed beam
	output = numpy.zeros((LFFT, (data.shape[1]-intDelays.max())/LFFT), dtype=numpy.complex64)
	for s in list(range(len(stands))):
		for l in list(range(output.shape[1])):
			section = data[s,(l*LFFT+intDelays[s]):((l+1)*LFFT+intDelays[s])]
			sectionF = numpy.fft.fft(section)
			sectionF *= numpy.exp(2j*numpy.pi*freq*frcDelays[s])
			output[:,l] = output[:,l] + sectionF

	# Done
	return output


def fftBeamShape(stands, sampleRate=dp_common.fS, LFFT=256, CentralFreq=49.0e6, azimuth=0.0, elevation=90.0, progress=False):
	"""Given a list of stands, compute the on-sky response of the delay-and-sum
	scheme implemented in intDelayAndSum.  A 360x90 numpy array spaning azimuth
	and elevation is returned."""

	# Get the stand delays and convert the delay times from seconds to samples
	delays = calcDelay(stands, freq=CentralFreq, azimuth=azimuth, elevation=elevation)
	delays = numpy.round(delays*sampleRate).astype(numpy.int16)

	# Build up a base time array, load in the cable delays, and get the stand 
	# positions for geometric delay calculations.
	t = numpy.arange(0,1000)/sampleRate
	dlyCache = uvUtils.CableCache(CentralFreq, applyDispersion=True)
	xyz = uvUtils.getXYZ(stands)
	arrayX = xyz[:,0].mean()
	arrayY = xyz[:,1].mean()
	arrayZ = xyz[:,2].mean()
	arrayXYZ = xyz - numpy.array([arrayX, arrayY, arrayZ])

	# Load in the respoonse of a single isolated stand
	standBeam = __loadStandResponse(freq=CentralFreq)

	# Build the output array and loop over all azimuths and elevations
	output = numpy.zeros((360,90))
	for az in list(range(360)):
		rAz = az*numpy.pi/180.0
		for el in list(range(90)):
			rEl = el*numpy.pi/180.0

			# Display the progress meter if the `progress' keyword is set to True.  The
			# progress meter displays a `.' every 2% complete and the percentages every
			# 10%.  At 100%, `Done' is displayed.
			if progress:
				fracDone = (az*90+el) / 32400.0 * 100
				if fracDone % 10 == 0 and round(fracDone,2) != 100:
					sys.stdout.write("%i%%" % fracDone)
				elif round(fracDone,2) == 100:
					sys.stdout.write("Done\n")
				elif round(fracDone,3) % 2 == 0:
					sys.stdout.write(".")
				else:
					pass
				sys.stdout.flush()

			# Unit vector for the currect on-sky location
			currPos = numpy.array([numpy.cos(rEl)*numpy.sin(rAz), 
							numpy.cos(rEl)*numpy.cos(rAz), 
							numpy.sin(rEl)])
			# Stand response in this direction
			currResponse = standBeam.response(aipy.coord.azalt2top(numpy.concatenate([[rAz], [rEl]])))[0][0]

			# Loop over stands to build the simulated singnals
			signals = numpy.zeros((len(stands), 1000))
			for i in list(range(len(stands))):
				currDelay = dlyCache.cableDelay(stands[i]) - numpy.dot(currPos, arrayXYZ[i,:]) / c
				signals[i,:] = currResponse * numpy.cos(2*numpy.pi*CentralFreq*(t + currDelay))

			# Beamform with delay-and-sum and store the RMS result
			beam = fftDelayAndSum(stands, signals, sampleRate=sampleRate, LFFT=LFFT, CentralFreq=CentralFreq, 
								azimuth=azimuth, elevation=elevation)
			output[az,el] = (numpy.abs(beam)**2).max()

	# Done
	return output