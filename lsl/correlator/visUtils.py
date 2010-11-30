# -*- coding: utf-8 -*-

"""Module of methods for visualizing visibility data.  Routines include:
  plots of visibility amplitude as a function of uv radius
  plots of phase for specific baselines as a function of time
  plots of closure phase for three antennae as a function of time
In addition, most of the functions have the ability to compare the observed 
data with simulations."""

import math
import numpy as n

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

__version__ = '0.1'
__revision__ = '$ Revision: 15 $'
__all__ = ['argument', 'unwrap', 'unmaskCompressedData', 'plotVisibilities', 'plotPhases', 'fitPhases', 'plotClosure', '__version__', '__revision__', '__all__']


def argument(data):
	"""Return the argument (phase) of a complex number.  This function is 
	a 'short cut' for the numpy.angle function."""
	
	return n.angle(data)


def unwrap(theta):
	"""Given a numpy array of phases (complex number arguments), unwrap the 
	phases by looking for +/- pi jumps in the data.  Return an unwrapped 
	numpy array of the same length."""

	output = 1.0*theta
	for i in range(1, len(theta), 1):
		if output[i] - output[i-1] >= n.pi:
			output[i:] = output[i:] - 2*n.pi
		elif output[i] - output[i-1] <= -n.pi:
			output[i:] = output[i:] + 2*n.pi
	return output


def unmaskCompressedData(data, mask):
	"""Given a numpy array that has been compressed and the mask used to 
	compress it, rebuild the original full-sized array.  Values that were 
	previously masked are replaced with NaNs."""

	origData = n.zeros(len(mask), dtype=data.dtype)
	
	maskCount = 0
	dataCount = 0
	for maskValue in mask:
		if maskValue == 0:
			origData[maskCount] = data[dataCount]
			dataCount = dataCount + 1
		else:
			origData[maskCount] = n.nan
		maskCount = maskCount + 1

	return origData


def plotVisibilities(dataDict,  pol='yy', chan=None, jd=None, simDict=None, SaveFig=False):
	"""Given a data dictionary from readUVData, create a plot of visibility
	amplitude as a function of uv radius.  Alternately, a Julian Date can be 
	specified and only data matched that date are plotted.  If a dictionary of 
	simulated data are also supplied via the simDict keyword, a two-panel plot
	is created comparing the real and simulated data."""
	
	fig = plt.figure()
	if simDict is not None:
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 1, 2)
	else:
		ax1 = fig.add_subplot(1, 1, 1)

	# Calculate the uv radius, normalize the amplitudes, and plot.  Each
	# baseline is done all at once so they all have the same colors
	for t, uvw, data in zip(dataDict['jd'][pol], dataDict['uvw'][pol], dataDict['vis'][pol]):
		if jd is not None:
			if t != jd:
				continue
		if chan is None:
			uvDist = n.sqrt( uvw[0,:]**2 + uvw[1,:]**2)
			visAmp = n.abs(data)
		else:
			uvDist = n.array([n.sqrt( uvw[0,chan]**2 + uvw[1,chan]**2)])
			visAmp = n.array(n.abs(data[chan]))
		
		ax1.plot(uvDist, visAmp, marker='o', linestyle=' ', alpha=0.50)
		ax1.semilogy()

	if simDict is not None:
		for t, uvw, data in zip(simDict['jd'][pol], simDict['uvw'][pol], simDict['vis'][pol]):
			if jd is not None:
				if t != jd:
					continue
			if chan is None:
				uvDist = n.sqrt( uvw[0,:]**2 + uvw[1,:]**2)
				visAmp = n.abs(data)
			else:
				uvDist = n.array([n.sqrt( uvw[0,chan]**2 + uvw[1,chan]**2)])
				visAmp = n.array(n.abs(data[chan]))
			
			ax2.plot(uvDist, visAmp, marker='o', linestyle=' ', alpha=0.50)
			ax2.semilogy()

	# Set axis title and labels
	ax1.set_title('Pol.: %s' % pol)
	ax1.set_xlabel('$uv$ Radius [$\lambda$]')
	ax1.set_ylabel('Amplitude')

	if simDict is not None:
		ax2.set_title('Pol.: %s (Sim.)' % pol)
		ax2.set_xlabel('$uv$ Radius [$\lambda$]')
		ax2.set_ylabel('Amplitude')

	plt.show()

	if SaveFig:
		fig.savefig('visibilities.png')


def plotPhases(dataDict, baselines=[0], pol='yy', jd=None, stands=None, simDict=None, Unwrap=False, SaveFig=False):
	"""Given a data dictionary from readUVData, create a plot of visibility
	phase as a function of frequency.  If a dictionary of simulated data are 
	also supplied via the simDict keyword, a two-panel plot is created comparing 
	the real and simulated data."""

	fig = plt.figure()
	if simDict is not None:
		ax1 = fig.add_subplot(3, 1, 1)
		ax2 = fig.add_subplot(3, 1, 2)
		ax3 = fig.add_subplot(3, 1, 3)
	else:
		ax1 = fig.add_subplot(1, 1, 1)

	# Loop over baselines since we care about only certain baselines
	for baseline in baselines:
		# Load in the frequency and the data
		freq = n.ma.array(dataDict['freq'], mask=dataDict['msk'][pol][baseline])
		data = dataDict['vis'][pol][baseline]

		# Calculate the phase from the visibility data and generate a label.  If
		# the stands keyword has been used in the call, the label contains actual 
		# stand numbers.
		phs = argument(dataDict['vis'][pol][baseline])
		if Unwrap:
			phs = unwrap(phs)
		if stands is not None:
			i,j = dataDict['bls'][pol][baseline]
			lbl = '%i-%i' % (stands[i], stands[j])
		else:
			lbl = '%i-%i' % dataDict['bls'][pol][baseline]
		
		# Plot depending on if the data are masked or not
		if phs.shape[0] != dataDict['msk'][pol][baseline].shape[0]:
			obsLine = ax1.plot(freq.compressed()/1e6, phs, marker='o', alpha=0.50, label=lbl)
		else:
			obsLine = ax1.plot(freq/1e6, phs, marker='o', alpha=0.50, label=lbl)

		plotColor = obsLine[-1].get_color()
		if simDict is not None:
			# Load in the frequency and the data
			freq = n.ma.array(simDict['freq'], mask=simDict['msk'][pol][baseline])
			sdata = simDict['vis'][pol][baseline]

			# Calculate the phase from the visibility data and generate a label.  If
			# the stands keyword has been used in the call, the label contains actual 
			# stand numbers.
			phs = argument(sdata)
			if Unwrap:
				phs = unwrap(phs)
			if stands is not None:
				i,j = dataDict['bls'][pol][baseline]
				lbl = '%i-%i (S)' % (stands[i], stands[j])
				lbl2 = '%i-%i (D)' % (stands[i], stands[j])
			else:
				lbl = '%i-%i (S)' % dataDict['bls'][pol][baseline]
				lbl2 = '%i-%i (D)' % dataDict['bls'][pol][baseline]
			
			rdata = data / sdata
			rphs = argument(rdata)
			if Unwrap:
				rphs = unwrap(rphs)

			# Plot depending on if the data are masked or not
			if phs.shape[0] != dataDict['msk'][pol][baseline].shape[0]:
				ax2.plot(freq.compressed()/1e6, phs, marker='^', alpha=0.50, label=lbl, color=plotColor)
				ax3.plot(freq.compressed()/1e6, rphs, marker='x', alpha=0.50, label=lbl2, color=plotColor)
			else:
				ax2.plot(freq/1e6, phs, marker='^', alpha=0.50, label=lbl, color=plotColor)
				ax3.plot(freq/1e6, rphs, marker='x', alpha=0.50, label=lbl2, color=plotColor)

			if Unwrap:
				delayCoeffs = n.polyfit(freq, rphs, 1)
				delayFit = n.polyval(delayCoeffs, freq)
				delay = delayCoeffs[0] / 2.0 / n.pi * 1e9
				ax3.plot(freq/1e6, delayFit, linestyle='-', label="Delay: %.2f ns" % delay, color=plotColor)

	if not Unwrap:
		ax1.set_yticks([-math.pi, -math.pi/2, 0, math.pi/2, math.pi])
		ax1.set_yticklabels(['-$\pi$', '-$\pi$/2', '0', '$\pi$/2', '$\pi$'])

	ax1.set_title('Pol.: %s' % pol)
	ax1.set_xlabel('Frequency [MHz]')
	ax1.set_ylabel('Phase [rad]')
		
	ax1.legend(loc=0)

	if simDict is not None:
		if not Unwrap:
			ax2.set_ylim([-math.pi, math.pi])
			ax2.set_yticks([-math.pi, -math.pi/2, 0, math.pi/2, math.pi])
			ax2.set_yticklabels(['-$\pi$', '-$\pi$/2', '0', '$\pi$/2', '$\pi$'])

		ax2.set_title('Pol.: %s (Sim.)' % pol)
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('Phase [rad]')
			
		if not Unwrap:
			ax3.set_ylim([-math.pi, math.pi])
			ax3.set_yticks([-math.pi, -math.pi/2, 0, math.pi/2, math.pi])
			ax3.set_yticklabels(['-$\pi$', '-$\pi$/2', '0', '$\pi$/2', '$\pi$'])

		ax3.set_title('Pol.: %s (Sim.)' % pol)
		ax3.set_xlabel('Frequency [MHz]')
		ax3.set_ylabel('Phase Diff. [rad]')

		ax2.legend(loc=0)
		ax3.legend(loc=0)

	plt.show()

	if SaveFig:
		fig.savefig('phases.png')


def fitPhases(dataDict, simDict, baseline=0, pol='yy', minCorr=0.90, returnDispersion=False, stands=None, NoPlot=False):
	"""Given a data dictionary from readUVData, fit delays by examining the
	change of the visiblility phases with frequency.  This routine examines
	each baseline of interest and looks for flat regions that are 3 MHz wide 
	and have R correlation values greater than 'minCorr'.  In theses flat 
	regions, a line is fit to the data and average slope found from the fits
	yields a delay.  The delays are returned as a numpy array with values in
	ns.  If the returnDispersion keyword is set, a numpy array with the 
	standard deviation of the fits is also returned."""

	import pylab

	# Build a list of unique Julian dates present in the data
	junk = {}
	for dt in dataDict['jd'][pol]:
		junk[dt] = 1
	jds = sorted(junk.keys())
	nJD = len(jds)
	del(junk)

	# Number of baselines in the dataset
	nBL = len(dataDict['vis'][pol]) / nJD
	refBaseline = baseline % nBL

	# Setup the output figure
	if not NoPlot:
		fig = pylab.figure()
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 1, 2)

	# Loop over baselines
	allDelay = []
	baseDelay = []
	for baseline in range(refBaseline, len(dataDict['vis'][pol]), nBL):
		# Load the frequency and removed flagged values
		freq = n.ma.array(dataDict['freq'], mask=dataDict['msk'][pol][baseline])
		freq = freq.compressed()

		# Load the visibility data and calculate the phase
		phs = unwrap(argument(dataDict['vis'][pol][baseline] / simDict['vis'][pol][baseline]))

		# Plot
		if not NoPlot:
			obsLine = ax1.plot(freq/1e6, phs, marker='x')
			plotColor = obsLine[-1].get_color()

		# Fit the delay of the unwrapped phase difference
		delayCoeffs = n.polyfit(freq, phs, 1)
		delayFit = n.polyval(delayCoeffs, freq)
		delay = delayCoeffs[0] / 2.0 / n.pi * 1e9
		R = n.corrcoef(freq, phs)[0,1]

		allDelay.append( delay )
		if n.abs(R) >= minCorr:
			baseDelay.append( delay )
			if not NoPlot:
				ax1.plot(freq/1e6, delayFit, color=plotColor)

	# Convert baseDelay to a numpy array and remove bad fits characterized by a NaN
	allDelay = n.array(allDelay)
	allDelay = allDelay.compress(n.isfinite(allDelay))
	baseDelay = n.array(baseDelay)
	baseDelay = baseDelay.compress(n.isfinite(baseDelay))
	print refBaseline, len(allDelay), allDelay.mean(), allDelay.std()
	print refBaseline, len(baseDelay), baseDelay.mean(), baseDelay.std()

	delays = baseDelay.mean()
	disprs = baseDelay.std()
	
	if not NoPlot:
		dataRange = baseDelay.max() - baseDelay.min()
		nBins = int(round(n.abs(dataRange)/15.0))
		ax2.hist(allDelay, bins=nBins, range=[allDelay.min(), allDelay.max()], label='All')
		ax2.hist(baseDelay, bins=nBins, range=[allDelay.min(), allDelay.max()], label='R$\geq$%.3f' % minCorr)

		title = 'Delays for Baseline #%i' % refBaseline
		if stands is not None:
			i, j = dataDict['bls'][pol][refBaseline]
			title = "%s (%i-%i)" % (title, stands[i], stands[j])
		ax1.set_title(title)
		ax1.set_xlabel('Frequency [MHz]')
		ax1.set_ylabel('$\Delta$ Phase [rad]')
		ax2.set_xlabel('Delay [ns]')
		ax2.set_ylabel('Number of Samples')
		ax2.legend(loc=0)

		pylab.show()

	if returnDispersion:
		return (delays, disprs)
	else:
		return delays


def plotClosure(dataDict, pol='yy', antennas=[1, 2, 3], individualPhases=False, stands=None, SaveFig=False):
	"""Given a data dictionary from readUVData, plot the closure phase of three 
	antennae as a function of time.  If the individualPhases keyword is set, also 
	plot the phases over time of the three antennae that comprise the closure 
	triplet."""

	fig = plt.figure()
	if not individualPhases:
		ax1 = fig.add_subplot(1, 1, 1)
	else:
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 3, 4)
		ax3 = fig.add_subplot(2, 3, 5)
		ax4 = fig.add_subplot(2, 3, 6)

	# Build a list of unique Julian dates present in the data
	junk = {}
	for dt in dataDict['jd'][pol]:
		junk[dt] = 1
	jds = sorted(junk.keys())
	print jds
	del(junk)

	# Loop over the JDs to build up the closure phase as a function of time
	for jd in jds:
		print "Starting on JD %.5f" % jd
		# Initialize closerPhase to 0 for all frequencies
		trispectrum = n.ones_like(dataDict['freq'])
		individPhase = n.zeros((3,len(dataDict['freq'])))
		closureCount = 0
		
		# Loop over the data
		for (i,j),vis,msk,t in zip(dataDict['bls'][pol], dataDict['vis'][pol], dataDict['msk'][pol], dataDict['jd'][pol]):
			if t < jd:
				continue
			if t > jd:
				break

			# If we come across antennae that are part of the closure triplet, get their
			# phases.  Note:  one pair in the triplet needs it Hermitian conjugate to 
			# to properly close the phase.  This pair is taken to be the baseline formed
			# by the first and last antennae specified.
			if i in antennas and j in antennas:
				if i == antennas[0] and j == antennas[2]:
					temp = j
					j = i
					i = temp
					vis = vis.conj()
				vis = unmaskCompressedData(vis, msk)

				if closureCount == 0:
					trispectrum = vis
				else:
					trispectrum *= vis
				if i == antennas[0]:
					individPhase[0,:] = argument(vis)
				elif i == antennas[1]:
					individPhase[1,:] = argument(vis)
				else:
					individPhase[2,:] = argument(vis)
				closureCount += 1

			# If clouseCount indicates that we have added all 3 phases together, plot the
			# result
			if closureCount == 3:
				toUse = [len(dataDict['freq'])/4, len(dataDict['freq'])/2, 3*len(dataDict['freq'])/4]

				closurePhase = argument(trispectrum)

				time = (n.ones_like(dataDict['freq'])*t - jds[0])
				c = ax1.scatter(time[toUse]*24.0*60.0, closurePhase[toUse]*180.0/n.pi, c=dataDict['freq'][toUse])
				
				if individualPhases:
					bad = n.where( individPhase > n.pi )
					if len(bad[0]) > 0:
						individPhase -= (2*n.pi)

					ax2.scatter(time[toUse]*24.0*60.0, individPhase[0,toUse]*180.0/n.pi, c=dataDict['freq'][toUse])
					ax3.scatter(time[toUse]*24.0*60.0, individPhase[1,toUse]*180.0/n.pi, c=dataDict['freq'][toUse])
					ax4.scatter(time[toUse]*24.0*60.0, individPhase[2,toUse]*180.0/n.pi, c=dataDict['freq'][toUse])
				break

	if stands is None:
		title = 'Closure Phase for %i-%i-%i' % (antennas[0],antennas[1],antennas[2])
	else:
		title = 'Closure Phase for %i-%i-%i' % (stands[antennas[0]],stands[antennas[1]],stands[antennas[2]])
	ax1.set_title(title)
	ax1.set_xlabel('Time from First Capture [min]')
	ax1.set_ylabel('Closure Phase [deg]')
	step = len(dataDict['freq'])/10
	b = fig.colorbar(c, ax=ax1, ticks=dataDict['freq'][0::step])
	b.ax.set_yticklabels(["%.2f MHz" % (1e-6*f) for f in (dataDict['freq'][0::step])])

	ax1.set_ylim([-180, 180])
	if individualPhases:
		ax2.set_ylabel('Phase [deg]')
		ax2.set_ylim([-180, 180])
		ax3.set_ylim([-180, 180])
		ax3.yaxis.set_major_formatter( NullFormatter() )
		ax4.set_ylim([-180, 180])
		ax4.yaxis.set_major_formatter( NullFormatter() )

	plt.show()

	if SaveFig:
		fig.savefig('closure.png')


def plotWaterfall(dataDict, baseline=21, pol='yy', stands=None, simDict=None, SaveFig=False):
	"""Create a waterfall diagram (visibility amplitude as a function of time and 
	frequency) for a particular baseline.  If simulated data are also provided, 
	plot those on a seperate axis."""
	
	fig = plt.figure()
	if simDict is None:
		ax1 = fig.add_subplot(1, 1, 1)
	else:
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 1, 2)

	# Number of channels in the data set
	nChan = len(dataDict['freq'])

	# Build a list of unique Julian dates present in the data
	junk = {}
	for dt in dataDict['jd'][pol]:
		junk[dt] = 1
	jds = sorted(junk.keys())
	nJD = len(jds)
	del(junk)

	# Number of baselines in the dataset
	nBL = len(dataDict['vis'][pol]) / nJD
	baseline = baseline % nBL

	# Create the waterfall diagram for the data
	wfd = n.zeros((nJD,nChan))
	for i in range(nJD):
		if dataDict['isMasked']:
			wfd[i,:] = unmaskCompressedData(n.abs( dataDict['vis'][pol][baseline+nBL*i] ), dataDict['msk'][pol][baseline+nBL*i])
		else:
			wfd[i,:] = n.abs( dataDict['vis'][pol][baseline+nBL*i] )
	
	# Do the same for the simulated data if it is provided
	if simDict is not None:
		swfd = n.zeros((nJD,nChan))
		for i in range(nJD):
			if simDict['isMasked']:
				swfd[i,:] = unmaskCompressedData(n.abs( simDict['vis'][pol][baseline+nBL*i] ), simDict['msk'][pol][baseline+nBL*i])
			else:
				swfd[i,:] = n.abs( simDict['vis'][pol][baseline+nBL*i] )

	# Determine intelligent clipping levels
	hist, bins = n.histogram(wfd.compress(n.isfinite(wfd.ravel())), bins=1000)
	cdf = 1.0*hist.cumsum() / hist.sum()
	lowerClip = (n.where( cdf >= 0.05, bins[1:], n.zeros_like(bins[1:])+1e12)).min()
	upperClip = (n.where( cdf <= 0.90, bins[1:], n.zeros_like(bins[1:])-1e12)).max()
	print lowerClip, upperClip

	# Plot it all up
	ax1.imshow(wfd, origin='lower', vmin=lowerClip, vmax=upperClip, 
				extent=(dataDict['freq'][0]/1e6, dataDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))
	if simDict is not None:
		# Determine intelligent clipping levels
		hist, bins = n.histogram(swfd.compress(n.isfinite(swfd.ravel())), bins=1000)
		cdf = 1.0*hist.cumsum() / hist.sum()
		lowerClip = (n.where( cdf >= 0.05, bins[1:], n.zeros_like(bins[1:])+1e12)).min()
		upperClip = (n.where( cdf <= 0.95, bins[1:], n.zeros_like(bins[1:])-1e12)).max()
		print lowerClip, upperClip

		ax2.imshow(swfd, origin='lower', vmin=lowerClip, vmax=upperClip, 
					extent=(simDict['freq'][0]/1e6, simDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))

	if stands is None:
		title = 'Visibility Amplitudes for %i-%i' % dataDict['bls'][pol][baseline]
	else:
		i,j = dataDict['bls'][pol][baseline]
		title = 'Visibility Amplitudes for %i-%i' % (stands[i], stands[j])
	ax1.set_title(title)	
	ax1.set_xlabel('Frequency [MHz]')
	ax1.set_ylabel('Time [JD-%.4f]' % jds[0])
	ax1.axis('auto')

	if simDict is not None:
		if stands is None:
			title = 'Visibility Amplitudes for %i-%i (Sim.)' % simDict['bls'][pol][baseline]
		else:
			i,j = simDict['bls'][pol][baseline]
			title = 'Visibility Amplitudes for %i-%i (Sim.)' % (stands[i], stands[j])
		ax2.set_title(title)	
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('Time [JD-%.4f]' % jds[0])
		ax2.axis('auto')

	plt.show()

	if SaveFig:
		fig.savefig('waterfall.png')


def matchAmps(dataDict, baseline=21, pol='yy', stands=None, simDict=None, SaveFig=False):
	"""Create a waterfall diagram (visibility amplitude as a function of time and 
	frequency) for a particular baseline.  If simulated data are also provided, 
	plot those on a seperate axis."""
	
	fig = plt.figure()
	if simDict is None:
		ax1 = fig.add_subplot(1, 1, 1)
	else:
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 1, 2)

	# Number of channels in the data set
	nChan = len(dataDict['freq'])

	# Build a list of unique Julian dates present in the data
	junk = {}
	for dt in dataDict['jd'][pol]:
		junk[dt] = 1
	jds = sorted(junk.keys())
	nJD = len(jds)
	del(junk)

	# Number of baselines in the dataset
	nBL = len(dataDict['vis'][pol]) / nJD
	baseline = baseline % nBL

	# Create the waterfall diagram for the data
	wfd = n.zeros((nJD,nChan))
	for i in range(nJD):
		if dataDict['isMasked']:
			wfd[i,:] = unmaskCompressedData(n.abs( dataDict['vis'][pol][baseline+nBL*i] ), dataDict['msk'][pol][baseline+nBL*i])
		else:
			wfd[i,:] = n.abs( dataDict['vis'][pol][baseline+nBL*i] )
	

	bestMatch = -1
	bestChi2 = 1e90

	for b in range(nBL):
		# Do the same for the simulated data if it is provided
		if simDict is not None:
			swfd = n.zeros((nJD,nChan))
			for i in range(nJD):
				if simDict['isMasked']:
					swfd[i,:] = unmaskCompressedData(n.abs( simDict['vis'][pol][b+nBL*i] ), simDict['msk'][pol][b+nBL*i])
				else:
					swfd[i,:] = n.abs( simDict['vis'][pol][b+nBL*i] )

		Chi2 = ((wfd - swfd*wfd.std()/swfd.std())**2).sum()
		if Chi2 < bestChi2:
			print 'New best match is %i @ %f' % (b, Chi2)
			bestMatch = b
			bestChi2 = Chi2

	# Do the same for the simulated data if it is provided
	if simDict is not None:
		swfd = n.zeros((nJD,nChan))
		for i in range(nJD):
			if simDict['isMasked']:
				swfd[i,:] = unmaskCompressedData(n.abs( simDict['vis'][pol][b+nBL*i] ), simDict['msk'][pol][b+nBL*i])
			else:
				swfd[i,:] = n.abs( simDict['vis'][pol][b+nBL*i] )

	# Determine intelligent clipping levels
	hist, bins = n.histogram(wfd.compress(n.isfinite(wfd.ravel())), bins=1000)
	cdf = 1.0*hist.cumsum() / hist.sum()
	lowerClip = (n.where( cdf >= 0.05, bins[1:], n.zeros_like(bins[1:])+1e12)).min()
	upperClip = (n.where( cdf <= 0.90, bins[1:], n.zeros_like(bins[1:])-1e12)).max()
	print lowerClip, upperClip

	# Plot it all up
	ax1.imshow(wfd, origin='lower', vmin=lowerClip, vmax=upperClip, 
				extent=(dataDict['freq'][0]/1e6, dataDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))
	if simDict is not None:
		# Determine intelligent clipping levels
		hist, bins = n.histogram(swfd.compress(n.isfinite(swfd.ravel())), bins=1000)
		cdf = 1.0*hist.cumsum() / hist.sum()
		lowerClip = (n.where( cdf >= 0.05, bins[1:], n.zeros_like(bins[1:])+1e12)).min()
		upperClip = (n.where( cdf <= 0.95, bins[1:], n.zeros_like(bins[1:])-1e12)).max()
		print lowerClip, upperClip

		ax2.imshow(swfd, origin='lower', vmin=lowerClip, vmax=upperClip, 
					extent=(simDict['freq'][0]/1e6, simDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))

	if stands is None:
		title = 'Visibility Amplitudes for %i-%i' % dataDict['bls'][pol][baseline]
	else:
		i,j = dataDict['bls'][pol][baseline]
		title = 'Visibility Amplitudes for %i-%i' % (stands[i], stands[j])
	ax1.set_title(title)	
	ax1.set_xlabel('Frequency [MHz]')
	ax1.set_ylabel('Time [JD-%.4f]' % jds[0])
	ax1.axis('auto')

	if simDict is not None:
		if stands is None:
			title = 'Visibility Amplitudes for %i-%i (Sim.)' % simDict['bls'][pol][baseline]
		else:
			i,j = simDict['bls'][pol][baseline]
			title = 'Visibility Amplitudes for %i-%i (Sim.)' % (stands[i], stands[j])
		ax2.set_title(title)	
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('Time [JD-%.4f]' % jds[0])
		ax2.axis('auto')

	plt.show()

	if SaveFig:
		fig.savefig('match-waterfall.png')


def plotDiffWaterfall(dataDict, simDict, baseline=21, pol='yy', stands=None, SaveFig=False):
	"""Create a waterfall diagram (visibility amplitude as a function of time and 
	frequency) for a particular baseline.  If simulated data are also provided, 
	plot those on a seperate axis."""
	
	fig = plt.figure()
	ax1 = fig.add_subplot(3, 1, 1)
	ax2 = fig.add_subplot(3, 1, 2)
	ax3 = fig.add_subplot(3, 1, 3)

	# Number of channels in the data set
	nChan = len(dataDict['freq'])

	# Build a list of unique Julian dates present in the data
	junk = {}
	for dt in dataDict['jd'][pol]:
		junk[dt] = 1
	jds = sorted(junk.keys())
	nJD = len(jds)
	del(junk)

	# Number of baselines in the dataset
	nBL = len(dataDict['vis'][pol]) / nJD
	baseline = baseline % nBL

	# Create the waterfall diagram for the data
	wfd = n.zeros((nJD,nChan))
	for i in range(nJD):
		if dataDict['isMasked']:
			wfd[i,:] = unmaskCompressedData(n.abs( dataDict['vis'][pol][baseline+nBL*i] ), dataDict['msk'][pol][baseline+nBL*i])
		else:
			wfd[i,:] = n.abs( dataDict['vis'][pol][baseline+nBL*i] )
	
	# Do the same for the simulated data if it is provided
	swfd = n.zeros((nJD,nChan))
	for i in range(nJD):
		if simDict['isMasked']:
			swfd[i,:] = unmaskCompressedData(n.abs( simDict['vis'][pol][baseline+nBL*i] ), simDict['msk'][pol][baseline+nBL*i])
		else:
			swfd[i,:] = n.abs( simDict['vis'][pol][baseline+nBL*i] )

	# Scale and difference
	diff = wfd - swfd*n.median(wfd)/n.median(swfd)
	print diff.min(), diff.mean(), diff.max()

	# Determine intelligent clipping levels
	hist, bins = n.histogram(wfd.compress(n.isfinite(wfd.ravel())), bins=1000)
	cdf = 1.0*hist.cumsum() / hist.sum()
	lowerClip = (n.where( cdf >= 0.05, bins[1:], n.zeros_like(bins[1:])+1e12)).min()
	upperClip = (n.where( cdf <= 0.95, bins[1:], n.zeros_like(bins[1:])-1e12)).max()
	print lowerClip, upperClip

	# Plot it all up
	ax1.imshow(wfd, origin='lower', vmin=lowerClip, vmax=upperClip, 
				extent=(dataDict['freq'][0]/1e6, dataDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))
	
	# Determine intelligent clipping levels
	hist, bins = n.histogram(swfd.compress(n.isfinite(swfd.ravel())), bins=1000)
	cdf = 1.0*hist.cumsum() / hist.sum()
	lowerClip = (n.where( cdf >= 0.05, bins[1:], n.zeros_like(bins[1:])+1e12)).min()
	upperClip = (n.where( cdf <= 0.95, bins[1:], n.zeros_like(bins[1:])-1e12)).max()
	print lowerClip, upperClip

	ax2.imshow(swfd, origin='lower', vmin=lowerClip, vmax=upperClip, 
				extent=(simDict['freq'][0]/1e6, simDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))

	ax3.imshow(diff, origin='lower', 
				extent=(simDict['freq'][0]/1e6, simDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))

	if stands is None:
		title = 'Visibility Amplitudes for %i-%i' % dataDict['bls'][pol][baseline]
	else:
		i,j = dataDict['bls'][pol][baseline]
		title = 'Visibility Amplitudes for %i-%i' % (stands[i], stands[j])
	ax1.set_title(title)	
	ax1.set_xlabel('Frequency [MHz]')
	ax1.set_ylabel('Time [JD-%.4f]' % jds[0])
	ax1.axis('auto')

	if stands is None:
		title = 'Visibility Amplitudes for %i-%i (Sim.)' % simDict['bls'][pol][baseline]
	else:
		i,j = simDict['bls'][pol][baseline]
		title = 'Visibility Amplitudes for %i-%i (Sim.)' % (stands[i], stands[j])
	ax2.set_title(title)	
	ax2.set_xlabel('Frequency [MHz]')
	ax2.set_ylabel('Time [JD-%.4f]' % jds[0])
	ax2.axis('auto')

	if stands is None:
		title = 'Visibility Amplitudes for %i-%i (Diff.)' % simDict['bls'][pol][baseline]
	else:
		i,j = simDict['bls'][pol][baseline]
		title = 'Visibility Amplitudes for %i-%i (Diff.)' % (stands[i], stands[j])
	ax3.set_title(title)	
	ax3.set_xlabel('Frequency [MHz]')
	ax3.set_ylabel('Time [JD-%.4f]' % jds[0])
	ax3.axis('auto')

	plt.show()

	if SaveFig:
		fig.savefig('diff-waterfall.png')


def plotPhaseWaterfall(dataDict, baseline=21, pol='yy', stands=None, simDict=None, SaveFig=False):
	"""Create a waterfall diagram (visibility phase as a function of time and 
	frequency) for a particular baseline.  If simulated data are also provided, 
	plot those on a seperate axis."""
	
	fig = plt.figure()
	if simDict is None:
		ax1 = fig.add_subplot(1, 1, 1)
	else:
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 1, 2)

	# Number of channels in the data set
	nChan = len(dataDict['freq'])

	# Build a list of unique Julian dates present in the data
	junk = {}
	for dt in dataDict['jd'][pol]:
		junk[dt] = 1
	jds = sorted(junk.keys())
	nJD = len(jds)
	del(junk)

	# Number of baselines in the dataset
	nBL = len(dataDict['vis'][pol]) / nJD
	baseline = baseline % nBL

	# Create the waterfall diagram for the data
	wfd = n.zeros((nJD,nChan))
	for i in range(nJD):
		phs = argument(dataDict['vis'][pol][baseline+nBL*i])
		if dataDict['isMasked']:
			wfd[i,:] = unmaskCompressedData(phs, dataDict['msk'][pol][baseline+nBL*i])
		else:
			wfd[i,:] = phs
	
	# Do the same for the simulated data if it is provided
	if simDict is not None:
		swfd = n.zeros((nJD,nChan))
		for i in range(nJD):
			phs = argument(simDict['vis'][pol][baseline+nBL*i])
			if simDict['isMasked']:
				swfd[i,:] = unmaskCompressedData(phs, simDict['msk'][pol][baseline+nBL*i])
			else:
				swfd[i,:] = phs

	# Plot it all up
	ax1.imshow(wfd, origin='lower', vmin=-n.pi, vmax=n.pi, 
				extent=(dataDict['freq'][0]/1e6, dataDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))
	if simDict is not None:
		ax2.imshow(swfd, origin='lower', vmin=-n.pi, vmax=n.pi, 
					extent=(simDict['freq'][0]/1e6, simDict['freq'][-1]/1e6, 0, jds[-1]-jds[0]))

	if stands is None:
		title = 'Visibility Phases for %i-%i' % dataDict['bls'][pol][baseline]
	else:
		i,j = dataDict['bls'][pol][baseline]
		title = 'Visibility Phases for %i-%i' % (stands[i], stands[j])
	ax1.set_title(title)	
	ax1.set_xlabel('Frequency [MHz]')
	ax1.set_ylabel('Time [JD-%.4f]' % jds[0])
	ax1.axis('auto')

	if simDict is not None:
		if stands is None:
			title = 'Visibility Phases for %i-%i (Sim.)' % simDict['bls'][pol][baseline]
		else:
			i,j = simDict['bls'][pol][baseline]
			title = 'Visibility Phases for %i-%i (Sim.)' % (stands[i], stands[j])
		ax2.set_title(title)	
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('Time [JD-%.4f]' % jds[0])
		ax2.axis('auto')

	plt.show()

	if SaveFig:
		fig.savefig('phase-waterfall.png')