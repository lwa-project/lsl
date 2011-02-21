# -*- coding: utf-8 -*-

import math
import numpy as n

try:
	from sim import scaleData
	from visUtils import argument
except ImportError, err:
	moduleName = (err.args[0]).split()[-1]
	print "The '%s' module is needed by this file." % moduleName
	sys.exit(-1)

def residualsSCP(delays, amps, dataDict, simDict, chan, pol='yy', MapSize=30, MapRes=0.50):
	"""Function used to perform a phase self-calibration of data stored in a 
	readUVData dictionary and a model sky stored in a sim.buildSimSky dictionary.
	This function returns the residuals between a scaled and delayed data set 
	and the model sky.  Only the delays are allowed to be updated by 
	scipy.optimize.leastsq in this function."""

	# Scale and delay the observed visibility data according to the 'delays' and 
	# 'amps' inputs
	scaledVis = []
	for obsVis,(i,j) in zip(dataDict['vis'][pol],dataDict['bls'][pol]):
		scaledVis.append( n.array([obsVis[chan]*amps[i]*amps[j]*n.exp(-1j*2.0*math.pi*(delays[j]-delays[i]))]) )
	obsVis = n.concatenate(scaledVis)

	# Mask the simulated sky visibility data with the actual data mask if needed.
	if not simDict['isMasked']:
		maskedVis = []
		for simVis,mask in zip(simDict['vis'][pol], dataDict['msk'][pol]):
			maskedVis.append( simVis.compress(n.logical_not(mask)) )
		simVis = n.concatenate(maskedVis)
	else:
		maskedVis = []
		for simVis in simDict['vis'][pol]:
			maskedVis.append( n.array([simVis[chan]]) )
		simVis = n.concatenate(maskedVis)

	# Compute the difference of the two sets of visibilities and return the 
	err = argument( simVis / obsVis )
	errOut = n.where( n.abs(err) < n.pi, err, n.pi-err )
	#print err.max(), errOut.max(), errOut.mean()
	
	return (errOut**2).sum()


def residualsSCA(amps, delays, dataDict, simDict, chan, pol='yy', MapSize=30, MapRes=0.50):
	"""Function used to perform a phase self-calibration of data stored in a 
	readUVData dictionary and a model sky stored in a sim.buildSimSky dictionary.
	This function returns the residuals between a scaled and delayed data set 
	and the model sky.  This function is the counterpart ot residualsSCP and only 
	the gains are allowed to be updated by scipy.optimize.leastsq ."""

	# Scale and delay the observed visibility data according to the 'delays' and 
	# 'amps' inputs
	scaledVis = []
	for obsVis,(i,j) in zip(dataDict['vis'][pol],dataDict['bls'][pol]):
		scaledVis.append( n.array([obsVis[chan]*amps[i]*amps[j]*n.exp(-1j*2.0*math.pi*(delays[j]-delays[i]))]) )
	obsVis = n.concatenate(scaledVis)

	# Mask the simulated sky visibility data with the actual data mask if needed.
	if not simDict['isMasked']:
		maskedVis = []
		for simVis,mask in zip(simDict['vis'][pol], dataDict['msk'][pol]):
			maskedVis.append( simVis.compress(n.logical_not(mask)) )
		simVis = n.concatenate(maskedVis)
	else:
		maskedVis = []
		for simVis in simDict['vis'][pol]:
			maskedVis.append( n.array([simVis[chan]]) )
		simVis = n.concatenate(maskedVis)

	# Compute the difference of the two sets of visibilities and return the 
	# absolute value since scipy.optimize.leastsq is looking for real values.
	err = n.abs(obsVis) - n.abs(simVis)
	
	return n.abs(err)


def selfCal(dataDict, simDict, MapSize=30, MapRes=0.5, pol='xx', chan=22, doGain=False, returnDelays=False):
	"""Function used to perform a phase self-calibration of data stored in a 
	readUVData dictionary and a model sky stored in a sim.buildSimSky dictionary 
	for a given polarization and channel."""

	from scipy.optimize import leastsq, fmin, fmin_cg

	# Initial guesses are no phase delays and blaket gains of 0.184
	initialDelays = n.zeros(20)
	initialAmps = n.ones(20)*0.388
	bestFit = fmin_cg(residualsSCP, initialDelays, args=(initialAmps, dataDict, simDict, chan), full_output=True)#, xtol=1e-15, ftol=1e-15, maxiter=1000000, maxfun=1000000)#, maxfev=5000, factor=100000)

	# Break the results array into its various components
	bestDelays = bestFit[0][0:20]
	bad = n.where( bestDelays < 0 )
	while len(bad[0]) > 0:
		bestDelays[bad[0]] += 1.0
		bad = n.where( bestDelays < 0 )
	bad = n.where( bestDelays > 1 )
	while len(bad[0]) > 0:
		bestDelays[bad[0]] -= 1.0
		bad = n.where( bestDelays > 1 )

	# If we want to run a gain calibration as well, run the minimization using 
	# residualsSCA.  Otherwise, just use the gain guesses from above to provide
	# the "best gains".
	if doGain:
		bestFit = leastsq(residualsSCA, initialAmps, args=(bestDelays, dataDict, simDict, chan), full_output=True, xtol=1e-7, ftol=1e-7, maxfev=5000, factor=1)
		bestAmps =  bestFit[0][0:20]

		bestFit = leastsq(residualsSCP, bestDelays, args=(bestAmps, dataDict, simDict, chan), full_output=True, xtol=1e-7, ftol=1e-7, maxfev=5000, factor=100)

		# Break the results array into its various components
		bestDelays = bestFit[0][0:20]
		bad = n.where( bestDelays < 0 )
		while len(bad[0]) > 0:
			bestDelays[bad[0]] += 1.0
			bad = n.where( bestDelays < 0 )
		bad = n.where( bestDelays > 1 )
		while len(bad[0]) > 0:
			bestDelays[bad[0]] -= 1.0
			bad = n.where( bestDelays > 1 )
	else:
		bestAmps = initialAmps

	# Report on the fits
	print 'Best Delays: ', bestDelays
	print 'Best Gains:  ', bestAmps
	print 'Message: ', bestFit[3], bestFit[4]

	if returnDelays:
		return (scaleData(dataDict, bestAmps, bestDelays), bestDelays)
	else:
		return scaleData(dataDict, bestAmps, bestDelays)