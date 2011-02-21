# -*- coding: utf-8 -*-

"""Module for generating simulated arrays and visilibity data.  The chief 
functions of this module are:

buildSimArray
  given a station object, a list of stands, and a list of frequencies, build 
  a AIPY AntennaArray-like object.  This module can also generate AntennaArray 
  objects with positional errors by setting the 'PosError' keyword to a 
  positive value.

buildSimData
  given a SimArray and a list of aipy.src sources, build up a collection of 
  visibilities for a given set of Julian dates

scaleData
  given a dictionary of simulated visibilities from buildSimData, apply 
  antenna-based gains and delays to the visibilities

shiftData
  given a dictionary of simulated visibilities from buildSimData, shift the uvw 
  coordinates of the visibilities.
  .. note::

	This only changes the uvw values and does not phase-shift the data.

The format of the data dictionaries mentioned above is:

primary keys
  The primary keys store the major aspects of the visiblity data, e.g., 
  frequency coverage, baseline pairs, uvw coordinates, etc.  Valid keys are:
    * *freq* - list of frequencies used in Hz
    * *isMasked* - whether or not the visibility data have been masked 
      (numpy.compress'd)
    * *bls* - list of baselines in (stand 1, stand2) format
    * *uvw* - list of uvw coordinates as 3-element numpy arrays
    * *vis* - list of visibility numpy arrays
    * *wgt* - list of weight arrays with the same length as the visilitity arrays
    * *msk* - list of mask arrays used for the data.  1 = masked, 0 = valid
    * *jd*  - list of Julian dates associated with each list element

secondary keys
  The bls, uvw, vis, wgt, msk, and jd primary keys also have secondary keys that 
  indicate which polarizations are being stored.  Valid keys are:
    * *xx*
    * *yy*
    * *xy*
    * *yx*

In addition to simulation functions, this module includes buildGriddedImage
which takes a dictionary of visibilities and returns and aipy.im.ImgW object.

.. versionchanged:: 0.3.0
	This module was formerly called lsl.sim.sim
"""

import os
import sys
import aipy
import math
import numpy as n

from lsl import astro
from lsl.common import dp as dp_common
from lsl.common.paths import data as dataPath
from lsl.correlator import uvUtils

__version__ = '0.2'
__revision__ = '$ Revision: 13 $'
__all__ = ['srcs', 'Antenna', 'AntennaArray', 'buildSimArray', 'buildSimData', 'scaleData', 'shiftData', 'buildGriddedImage', '__version__', '__revision__', '__all__']


# A dictionary of bright sources in the sky to use for simulations
srcs = aipy.src.get_catalog(srcs=['Sun', 'Jupiter', 'cas', 'crab', 'cyg', 'her', 'sgr', 'vir'])


class Antenna(aipy.amp.Antenna):
	"""Modification to the aipy.amp.Antenna class to also store the stand ID 
	number in the Antenna.stand attribute.  This also add a getBeamShape 
	attribute that pulls in the old vis.getBeamShape function."""

	def __init__(self, x, y, z, beam, phsoff=[0.,0.], bp_r=n.array([1]), bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0), stand=0, **kwargs):
		"""New init function that include stand ID number support.  From aipy.amp.Antenna:

		  * x,y z = antenna coordinates in equatorial (ns) coordinates
		  * beam = Beam object (implements response() function)
		  * phsoff = polynomial phase vs. frequency.  Phs term that is linear
		    with freq is often called 'delay'.
		  * bp_r = polynomial (in freq) modeling real component of passband
		  * bp_i = polynomial (in freq) modeling imaginary component of passband
		  * amp = overall multiplicative scaling of gain
		  * pointing = antenna pointing (az,alt).  Default is zenith.
		"""

		aipy.phs.Antenna.__init__(self, x,y,z, beam=beam, phsoff=phsoff)
		self.set_pointing(*pointing)
		self.bp_r = bp_r
		self.bp_i = bp_i
		self.amp = amp
		self.stand = stand
		self._update_gain()

	def get_beam_shape(self, pol='x'):
		"""Return a 360 by 90 by nFreqs numpy array showning the beam pattern of a
		particular antenna in the array.  The first two dimensions of the output 
		array contain the azimuth (from 0 to 359 degrees in 1 degree steps) and 
		altitlude (from 0 to 89 degrees in 1 degree steps)."""

		# Build azimuth and altitude arrays.  Be sure to convert to radians
		az = n.zeros((360,90))
		for i in range(360):
			az[i,:] = i*n.pi/180.0
		alt = n.zeros((360,90))
		for i in range(90):
			alt[:,i] = i*n.pi/180.0

		# The beam model is computed in terms of topocentric coordinates, so make that
		# converseion right quick using the aipy.coord module.
		xyz = aipy.coord.azalt2top(n.concatenate([[az],[alt]]))
		
		# I cannot figure out how to do this all at once, so loop through azimuth/
		# altitude pairs
		resp = n.zeros((360,90,len(self.beam.freqs)))
		for i in range(360):
			for j in range(90):
				resp[i,j,:] = n.squeeze( self.bm_response(n.squeeze(xyz[:,i,j]), pol=pol) )

		return resp


class AntennaArray(aipy.amp.AntennaArray):
	"""Modification to the aipy.ant.AntennaArray class to add a fuction to 
	retrieve the stands stored in the AntennaArray.ants attribute.  Also add 
	a function to set the array time from a UNIX timestamp."""

	def get_stands(self):
		"""Return a numpy array listing the stands found in the AntennaArray 
		object."""

		stands = []
		for ant in self.ants:
			stands.append(ant.stand)
		
		return numpy.array(stands)

	def set_unixtime(self, timestamp):
		"""Set the array time using a UNIX timestamp (epoch 1970)."""

		self.set_jultime(astro.unix_to_utcjd(time))


def buildSimArray(station, stands, freq, jd=None, PosError=0.0, ForceFlat=False, verbose=False):
	"""Build a AIPY AntennaArray for simulation purposes.  Inputs are a station 
	object defined from the lwa_common module, a numpy array of stand numbers, 
	and a numpy array of frequencies in either Hz of GHz.  Optional inputs are
	a Julian Date to set the array to and a positional error terms that perturbs
	each of the stands in x, y, and z.  The output of this module is an AIPY
	AntennaArray object.

	The shape of the antenna response is either flat (gain of 1 in all 
	directions) or modeled by a collection of spherical harmonics that are poly-
	nomials in frequency.  The spherical harmonics are used if the file 
	'beam_shape.npz' is found in the current directory."""

	# If the frequencies are in Hz, we need to convert to GHz
	freqs = freq.copy()
	if freqs.min() > 1e6:
		freqs /= 1.0e9
	
	# If the beam Alm coefficient file is present, build a more relatistc beam 
	# response.  Otherwise, assume a flat beam
	if not ForceFlat and os.path.exists(os.path.join(dataPath, 'beam-shape.npz')):
		dd = n.load(os.path.join(dataPath, 'beam-shape.npz'))
		coeffs = dd['coeffs']

		deg = coeffs.shape[0]-1
		lmax = int((math.sqrt(1+8*coeffs.shape[1])-3)/2)
		beamShapeDict = {}
		for i in range(deg+1):
			beamShapeDict[i] = n.squeeze(coeffs[-1-i,:])

		if verbose:
			print "Using Alm beam model with %i-order freq. polynomial and %i-order sph. harmonics" % (deg, lmax)
		beam = aipy.amp.BeamAlm(freqs, lmax=lmax, mmax=lmax, deg=deg, nside=128, coeffs=beamShapeDict)
	else:
		if verbose:
			print "Using flat beam model"
		beam = aipy.amp.Beam(freqs)

	if PosError != 0:
		print "WARNING:  Creating array with positional errors between %.3f and %.3f m" % (-PosError, PosError)

	# Build an array of AIPY Antenna objects
	ants = []
	posCache = uvUtils.PositionCache()
	for stand in stands:
		top = posCache.getXYZ(stand)
		top += (2*PosError*n.random.rand(3)-PosError)	# apply a random positional error if needed
		top.shape = (3,)
		eq = n.dot( aipy.coord.top2eq_m(0.0, station.lat), top )
		eq *= 100	# m -> cm
		eq /= aipy.const.c	# cm -> s
		eq *= 1e9	# s -> ns

		delayCoeff = n.zeros(2)

		amp = uvUtils.cableGain(stand, freq=(freq*1e9))
		
		ants.append( Antenna(eq[0], eq[1], eq[2], beam, phsoff=delayCoeff, amp=amp, stand=stand) )

	# Combine the array of antennas with the array's location to generate an
	# AIPY AntennaArray object
	simAA = AntennaArray(ants=ants, location=station.getAIPYLocation())

	# Set the Julian Data for the AntennaArray object if it is provided.  The try...except
	# clause is used to deal with people who may want to pass an array of JDs in rather than
	# just one.
	if jd is not None:
		try:
			junk = len(jd)
		except:
			simAA.set_jultime(jd)
		else:
			simAA.set_jultime(jd[0])

	return simAA


def __buildSimData(aa, srcs, pols=['xx', 'yy', 'xy', 'yx'], jd=None, phaseCenter='z', baselines=None, mask=None, verbose=False, count=None, max=None):
	"""Helper function for buildSimData so that buildSimData can be called with 
	a list of Julian Dates and reconstruct the data appropriately."""

	nFreq = (n.squeeze((aa.get_afreqs()))).shape[0]

	# Update the JD if necessary
	if jd is not None:
		if verbose:
			if count is not None and max is not None:
				print "Setting Julian Date to %.5f (%i of %i)" % (jd, count, max)
			else:
				print "Setting Julian Date to %.5f" % jd
		aa.set_jultime(jd)
	else:
		jd = aa.get_jultime()

	# Compute the source parameters
	srcs_eq = []
	srcs_jy = []
	srcs_fq = []
	srcs_sh = []
	if verbose:
		print "Sources Used for Simulation:"
	for name,src in srcs.iteritems():
		## Update the source's coordinates
		src.compute(aa)

		## RA/dec -> equitorial 
		srcEQ = src.get_crds(crdsys='eq', ncrd=3)
		srcEQ.shape = (3,1)
		srcs_eq.append( srcEQ )

		## Source flux over the bandpass
		jys = src.get_jys()
		jys.shape = (1,nFreq)
		srcs_jy.append( jys )

		## Frequencies that the fluxes correspond to
		frq = aa.get_afreqs()
		frq.shape = (1,nFreq)
		srcs_fq.append( frq )

		## Source shape parameters
		shp = n.array(src.srcshape)
		shp.shape = (3,1)
		srcs_sh.append( shp )

		if verbose:
			srcAzAlt = aipy.coord.top2azalt(src.get_crds(crdsys='top', ncrd=3))
			if srcAzAlt[1] < 0:
				print "  %s: below horizon" % name
			else:
				print "  %s: %+.2f altitude; flux mean %.2f kJy" % (name, srcAzAlt[1]*180.0/math.pi, jys.mean()/1e3)

	# Build the simulation cache
	aa.sim_cache( n.concatenate(srcs_eq, axis=1), 
				jys=n.concatenate(srcs_jy, axis=0),
				mfreqs=n.concatenate(srcs_fq, axis=0),
				srcshapes=n.concatenate(srcs_sh, axis=1) )

	# Build the simulated data.  If no baseline list is provided, build all 
	# baselines avaliable
	if baselines is None:
		baselines = uvUtils.getBaselines(n.zeros(len(aa.ants)), Indicies=True)

	# Define output data structure
	UVData = {'freq': n.squeeze(aa.get_afreqs()*1e9), 'uvw': {}, 'vis': {}, 'wgt': {}, 'msk': {}, 'bls': {}, 'jd': {}}
	if mask is None:
		UVData['isMasked'] = False
	else:
		UVData['isMasked'] = True

	# Go!
	if phaseCenter is not 'z':
		phaseCenter.compute(aa)

	for pol in pols:
		UVData['bls'][pol] = []
		UVData['uvw'][pol] = []
		UVData['vis'][pol] = []
		UVData['wgt'][pol] = []
		UVData['msk'][pol] = []
		UVData['jd'][pol] = []

		if mask is None:
			for (i,j) in baselines:
				if verbose:
					print "Simulating data for baseline %i-%i, pol. %s" % (i,j,pol)

				try:
					crd = aa.gen_uvw(i, j, src=phaseCenter)
					d = aa.sim(i, j, pol=pol)
					d = aa.phs2src(d, phaseCenter, i, j)
				except aipy.phs.PointingError:
					continue

				UVData['bls'][pol].append( (i,j) )
				UVData['uvw'][pol].append( n.squeeze(crd) )
				UVData['vis'][pol].append( d )
				UVData['wgt'][pol].append( n.ones_like(d)*len(UVData['vis'][pol][-1]) )
				UVData['msk'][pol].append( n.zeros_like(UVData['vis'][pol][-1]) )
				UVData['jd'][pol].append( jd )
		else:
			for (i,j),msk in zip(baselines,mask):
				if verbose:
					print "Simulating data for baseline %i-%i, pol. %s" % (i,j,pol)
				
				try:
					crd = aa.gen_uvw(i, j, src=phaseCenter)
					d = aa.sim(i, j, pol=pol)
					d = aa.phs2src(d, phaseCenter, i, j)
				except aipy.phs.PointingError:
					continue

				uvw = aa.gen_uvw(i, j, src=phaseCenter)
				uvw = n.squeeze(crd.compress(n.logical_not(msk), axis=2))
				vis = d.compress(n.logical_not(msk))
				wgt = n.ones_like(vis) * len(vis)

				UVData['bls'][pol].append( (i,j) )
				UVData['uvw'][pol].append( uvw )
				UVData['vis'][pol].append( vis )
				UVData['wgt'][pol].append( wgt )
				UVData['msk'][pol].append( msk )
				UVData['jd'][pol].append( jd )

	return UVData


def buildSimData(aa, srcs, pols=['xx', 'yy', 'xy', 'yx'], jd=None, phaseCenter='z', baselines=None, mask=None, verbose=False):
	"""Given an AIPY AntennaArray object and a dictionary of sources from 
	aipy.src.get_catalog, returned a data dictionary of simulated data taken at 
	zenith.  Optinally, the data can be masked using some referenced (observed) 
	data set or only a specific sub-set of baselines.
	
	..versionchanged:: 0.4.0
		Added the 'pols' keyword to only compute certain polarization components
	"""

	nFreq = (aa.get_afreqs()).shape[0]

	# Update the JD if necessary
	if jd is not None:
		try:
			junk = len(jd)
		except TypeError:
			jd = [jd]
	else:
			jd = [aa.get_jultime()]

	# Build up output dictionary
	UVData = {'freq': n.squeeze(aa.get_afreqs()*1e9), 'uvw': {}, 'vis': {}, 'wgt': {}, 'msk': {}, 'bls': {}, 'jd': {}}
	if mask is None:
		UVData['isMasked'] = False
	else:
		UVData['isMasked'] = True
	
	# Loop over polarizations to populate the simulated data set
	for pol in pols:
		UVData['bls'][pol] = []
		UVData['uvw'][pol] = []
		UVData['vis'][pol] = []
		UVData['wgt'][pol] = []
		UVData['msk'][pol] = []
		UVData['jd'][pol] = []

	# Loop over Julian days to fill in the simulated data set
	jdCounter = 1
	for juldate in jd:
		oBlk = __buildSimData(aa, srcs, pols=pols, jd=juldate, phaseCenter=phaseCenter, baselines=baselines, mask=mask, verbose=verbose, count=jdCounter, max=len(jd))
		jdCounter = jdCounter + 1

		for pol in oBlk['bls']:
			for b,u,v,w,m,j in zip(oBlk['bls'][pol], oBlk['uvw'][pol], oBlk['vis'][pol], oBlk['wgt'][pol], oBlk['msk'][pol], oBlk['jd'][pol]):
				UVData['bls'][pol].append( b )
				UVData['uvw'][pol].append( u )
				UVData['vis'][pol].append( v )
				UVData['wgt'][pol].append( w )
				UVData['msk'][pol].append( m )
				UVData['jd'][pol].append( j )

	return UVData


def buildSimSignals(aa, stands, srcs, length=1024000, jd=None, phaseCenter='z', pol='x', SampleRate=dp_common.fS, IQ=False, verbose=True):
	"""Build up a collection of time series data that represents the output
	voltage of each stand as a function of time for a given collection of
	sources.  Note:  This only seems to work for single sources."""

	eta = 119.9169832*n.pi	# impedance of free space (Ohm)

	Nstand = len(aa.ants)
	freq = n.squeeze((aa.get_afreqs())) * 1e9 # GHz to Hz
	Nfreq = freq.shape[0]
	
	# Compute the source parameters
	srcs_tp = []
	srcs_mp = []
	srcs_jy = []
	srcs_fq = []
	if verbose:
		print "Sources Used for Simulation:"
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
		Jyat1MHz = jys[ n.where( n.abs(frq-0.001) == n.abs(frq-0.001).min() ) ]
		jys = n.where( frq >= 0.001, jys, Jyat1MHz )

		## Filter out sources that are below the horizon or have no flux
		srcAzAlt = aipy.coord.top2azalt(top) * 180/math.pi
		if verbose:
			if srcAzAlt[1] < 0:
				print "  %s: below horizon" % name
			else:
				print "  %s: %+.2f altitude; flux mean %.2f kJy" % (name, srcAzAlt[1], jys.mean()/1e3)
		if srcAzAlt[1] <= 0 or jys.sum() <= 0:
			continue

		## Save values into the source arrays
		srcs_tp.append( top )
		srcs_mp.append( map )
		srcs_jy.append( jys )
		srcs_fq.append( frq )

	# Get the topocentric coorindates for the zenith
	zen = aipy.coord.azalt2top(n.array([[n.pi/4],[n.pi/2]]))
	zen = n.squeeze(zen)

	# Define the stand position and cable/signal delay caches to the simulator 
	# move along faster
	dlyCache = uvUtils.SignalCache(freq, applyDispersion=True)

	# Update the phase center if necessary
	if phaseCenter == 'z':
		phaseCenterMap = aipy.coord.eq2top_m(0.0, aa.lat)
	else:
		phaseCenter.compute(aa)
		phaseCenterMap = phaseCenter.map

	# Generate the time series response of each signal at each frequency
	t = n.arange(length)/SampleRate * 1e9 # s -> ns
	tdSignals = n.zeros((Nstand, length), dtype=n.complex64)
	temp = n.zeros((len(srcs_tp), Nstand, length), dtype=n.complex64)
	for sc,tp,mp,jy,fq in zip(range(temp.shape[0]), srcs_tp, srcs_mp, srcs_jy, srcs_fq):
		for i,ant in zip(range(Nstand), aa.ants):
			# Zeroth, get the beam response in the direction of the current source for all frequencies
			antResponse = n.squeeze( ant.bm_response(tp, pol=pol) )

			# First, do the geometric delay
			geoDelay = ( n.dot(mp, ant.pos).transpose() )[2] 
			geoDelayPC = ( n.dot(phaseCenterMap, ant.pos).transpose() )[2]
			print geoDelay, geoDelayPC

			# Second, do the cable delay
			Delayat1MHz = dlyCache.cableDelay(stands[i], freq=1.0e6)[0] * 1e9 # s -> ns
			cblDelay = dlyCache.cableDelay(stands[i]) * 1e9 # s -> ns
			# NB: Replace the cable delays below 1 MHz with the 1 MHz value to keep the 
			# delays from blowing up for small f
			cblDelay = n.where( fq >= 0.001, cblDelay, Delayat1MHz )

			for a,j,f,d in zip(antResponse, jy, fq, cblDelay):
				factor = a * n.sqrt(j)
				angle = 2*n.pi*f*(t + (d - (geoDelay-geoDelayPC)))
				temp[sc,i,:] += factor*(n.cos(angle) + 1j*n.sin(angle))

	print temp.shape
	tdSignals = temp.sum(axis=0)
	if not IQ:
		return tdSignals.real
	else:
		return tdSginals


def scaleData(dataDict, amps, delays):
	"""Apply a set of antenna-based real gain values and phase delays in ns to a 
	data dictionary.  Returned the new scaled and delayed dictionary.
	
	..versionchanged:: 0.4.0
		The delays are now expected to be in nanoseconds rather than radians.
	"""

	import copy

	# Build the data dictionary to hold the scaled and delayed data
	sclUVData = {'freq': (dataDict['freq']).copy(), 'uvw': {}, 'vis': {}, 'wgt': {}, 'msk': {}, 'bls': {}, 'jd': {}}
	if dataDict['isMasked']:
		sclUVData['isMasked'] = True
	else:
		sclUVData['isMasked'] = False
	fq = dataDict['freq'] / 1e9

	# Apply the scales and delays for all polarization pairs found in the original data
	for pol in dataDict['vis'].keys():
		sclUVData['bls'][pol] = copy.copy(dataDict['bls'][pol])
		sclUVData['uvw'][pol] = copy.copy(dataDict['uvw'][pol])
		sclUVData['vis'][pol] = []
		sclUVData['wgt'][pol] = copy.copy(dataDict['wgt'][pol])
		sclUVData['msk'][pol] = copy.copy(dataDict['msk'][pol])
		sclUVData['jd'][pol] = copy.copy(dataDict['jd'][pol])

		for (i,j),vis in zip(dataDict['bls'][pol], dataDict['vis'][pol]):
			sclUVData['vis'][pol].append( vis*amps[i]*amps[j]*n.exp(2j*math.pi*fq*(delays[j]-delays[i])) )

	return sclUVData
	

def shiftData(dataDict, aa):
	"""Shift the uvw coordinates in one data dictionary to a new set of uvw 
	coordinates that correspond to a new AntennaArray object.  This is useful
	for looking at how positional errors in the array affect the data."""

	import copy
	
	# Build the data dictionary to hold the shifted data
	sftUVData = {'freq': (dataDict['freq']).copy(), 'uvw': {}, 'vis': {}, 'wgt': {}, 'msk': {}, 'bls': {}, 'jd': {}}
	if dataDict['isMasked']:
		sftUVData['isMasked'] = True
	else:
		sftUVData['isMasked'] = False

	# Apply the coordinate shift all polarization pairs found in the original data
	for pol in dataDict['vis'].keys():
		sftUVData['bls'][pol] = copy.copy(dataDict['bls'][pol])
		sftUVData['uvw'][pol] = []
		sftUVData['vis'][pol] = copy.copy(dataDict['vis'][pol])
		sftUVData['wgt'][pol] = copy.copy(dataDict['wgt'][pol])
		sftUVData['msk'][pol] = copy.copy(dataDict['msk'][pol])
		sftUVData['jd'][pol] = copy.copy(dataDict['jd'][pol])

		for (i,j),m in zip(dataDict['bls'][pol], dataDict['msk'][pol]):
			crds = aa.gen_uvw(i, j, src='z')
			if dataDict['isMasked']:
				crds = crds.compress(n.logical_not(n.logical_not(m), axis=2))
			crds = n.squeeze(crds)
			sftUVData['uvw'][pol].append( crds )

	return sftUVData


def buildGriddedImage(dataDict, MapSize=30, MapRes=0.50, MapWRes=0.10, pol='xx', chan=None):
	"""Given a data dictionary, build an aipy.img.ImgW object of gridded uv data 
	which can be used for imaging.  The ImgW object itself is returned by this 
	function to make it more versatile."""

	im = aipy.img.ImgW(size=MapSize, res=MapRes, wres=MapWRes)

	if chan is not None:
		# Make sure that `chan' is an array by trying to find its length
		try:
			junk = len(chan)
		except TypeError:
			chan = [chan]

		# Build up the data using only the specified channels
		uvw = []
		vis = []
		wgt = []
		for (i,j),d,u,w,m in zip(dataDict['bls'][pol], dataDict['vis'][pol], dataDict['uvw'][pol], dataDict['wgt'][pol], dataDict['msk'][pol]):
			u = u[:,chan]
			u.shape = (3, len(chan))

			uvw.append(u)
			vis.append(n.array([d[chan]]))
			wgt.append(n.array([w[chan]]))
	else:
		uvw = dataDict['uvw'][pol]
		vis = dataDict['vis'][pol]
		wgt = dataDict['wgt'][pol]

	uvw = n.concatenate(uvw, axis=1)
	vis = n.concatenate(vis)
	wgt = n.concatenate(wgt)

	uvw, vis, wgt = im.append_hermitian(uvw, vis, wgts=wgt)
	im.put(uvw, vis, wgts=wgt)

	return im
