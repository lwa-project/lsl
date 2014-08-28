#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pytz
import math
import ephem
import getopt

from datetime import datetime, date, time, timedelta

from lsl.common.stations import lwa1
from lsl.transform import Time
from lsl.astro import utcjd_to_unix, MJD_OFFSET
from lsl.common import sdf, metabundle


__version__ = "0.2"
__revision__ = "$Rev: 941 $"

# Date/time manipulation
_UTC = pytz.utc
formatString = '%Y/%m/%d %H:%M:%S.%f %Z'


def getObsStartStop(obs):
	"""
	Given an observation, get the start and stop times (returned as a two-
	element tuple).
	"""
	
	# UNIX timestamp for the start
	tStart = utcjd_to_unix(obs.mjd + MJD_OFFSET)
	tStart += obs.mpm / 1000.0
	
	# UNIX timestamp for the stop
	tStop = tStart +  obs.dur / 1000.0
	
	# Conversion to a timezone-aware datetime instance
	tStart = _UTC.localize( datetime.utcfromtimestamp(tStart) )
	tStop  = _UTC.localize( datetime.utcfromtimestamp(tStop ) )
	
	# Return
	return tStart, tStop


def main(args):
	# Get the observer
	observer = lwa1.getObserver()
	
	# Filenames in an easier format
	inputTGZ  = args[0]
	
	# Parse the input file and get the dates of the observations
	project = metabundle.getSessionDefinition(inputTGZ)
	obsImpl = metabundle.getObservationSpec(inputTGZ)
	
	nObs = len(project.sessions[0].observations)
	tStart = [None,]*nObs
	for i in xrange(nObs):
		tStart[i]  = utcjd_to_unix(project.sessions[0].observations[i].mjd + MJD_OFFSET)
		tStart[i] += project.sessions[0].observations[i].mpm / 1000.0
		tStart[i]  = datetime.utcfromtimestamp(tStart[i])
		tStart[i]  = _UTC.localize(tStart[i])
	
	# Get the LST at the start
	observer.date = (min(tStart)).strftime('%Y/%m/%d %H:%M:%S')
	lst = observer.sidereal_time()
	
	# Report on the file
	print "Filename: %s" % inputTGZ
	print " Project ID: %s" % project.id
	print " Session ID: %i" % project.sessions[0].id
	print " Observations appear to start at %s" % (min(tStart)).strftime(formatString)
	print " -> LST at %s for this date/time is %s" % (lwa1.name, lst)
	
	lastDur = project.sessions[0].observations[nObs-1].dur
	lastDur = timedelta(seconds=int(lastDur/1000), microseconds=(lastDur*1000) % 1000000)
	sessionDur = max(tStart) - min(tStart) + lastDur
	
	print " "
	print " Total Session Duration: %s" % sessionDur
	print " -> First observation starts at %s" % min(tStart).strftime(formatString)
	print " -> Last observation ends at %s" % (max(tStart) + lastDur).strftime(formatString)
	if project.sessions[0].observations[0].mode not in ('TBW', 'TBN'):
		drspec = 'No'
		if project.sessions[0].spcSetup[0] != 0 and project.sessions[0].spcSetup[1] != 0:
			drspec = 'Yes'
		drxBeam = project.sessions[0].drxBeam
		if drxBeam < 1:
			drxBeam = "MCS decides"
		else:
			drxBeam = "%i" % drxBeam
		print " DRX Beam: %s" % drxBeam
		print " DR Spectrometer used? %s" % drspec
		if drspec == 'Yes':
			print " -> %i channels, %i windows/integration" % tuple(project.sessions[0].spcSetup)
	else:
		tbnCount = 0
		tbwCount = 0
		for obs in project.sessions[0].observations:
			if obs.mode == 'TBW':
				tbwCount += 1
			else:
				tbnCount += 1
		if tbwCount > 0 and tbnCount == 0:
			print " Transient Buffer Mode: TBW"
		elif tbwCount == 0 and tbnCount > 0:
			print " Transient Buffer Mode: TBN"
		else:
			print " Transient Buffer Mode: both TBW and TBN"
	print " "
	print "File Information:"
	fileInfo = metabundle.getSessionMetaData(inputTGZ)
	for obsID in fileInfo.keys():
		print " Obs. #%i: %s" % (obsID, fileInfo[obsID]['tag'])
	
	print " "
	print "ASP Configuration:"
	for which in ('Beginning', 'End'):
		aspConfig = metabundle.getASPConfigurationSummary(inputTGZ, which=which)
		
		print '  %s' % which
		for k in aspConfig.keys():
			print '    %s: %i' % (k, aspConfig[k])
			
	print " "
	print " Number of observations: %i" % nObs
	print " Observation Detail:"
	for i in xrange(nObs):
		currDur = project.sessions[0].observations[i].dur
		currDur = timedelta(seconds=int(currDur/1000), microseconds=(currDur*1000) % 1000000)
		
		print "  Observation #%i" % (i+1,)
		currObs = None
		for j in xrange(len(obsImpl)):
			if obsImpl[j]['obsID'] == i+1:
				currObs = obsImpl[j]
				break
				
		## Basic setup
		print "   Target: %s" % project.sessions[0].observations[i].target
		print "   Mode: %s" % project.sessions[0].observations[i].mode
		print "   Start:"
		print "    MJD: %i" % project.sessions[0].observations[i].mjd
		print "    MPM: %i" % project.sessions[0].observations[i].mpm
		print "    -> %s" % getObsStartStop(project.sessions[0].observations[i])[0].strftime(formatString)
		print "   Duration: %s" % currDur
		
		## DP setup
		if project.sessions[0].observations[i].mode not in ('TBW',):
			print "   Tuning 1: %.3f MHz" % (project.sessions[0].observations[i].frequency1/1e6,)
		if project.sessions[0].observations[i].mode not in ('TBW', 'TBN'):
			print "   Tuning 2: %.3f MHz" % (project.sessions[0].observations[i].frequency2/1e6,)
		if project.sessions[0].observations[i].mode not in ('TBW',):
			print "   Filter code: %i" % project.sessions[0].observations[i].filter
		if currObs is not None:
			if project.sessions[0].observations[i].mode not in ('TBW',):
				if project.sessions[0].observations[i].mode == 'TBN':
					print "   Gain setting: %i" % currObs['tbnGain']
				else:
					print "   Gain setting: %i" % currObs['drxGain']
		else:
			print "   WARNING: observation specification not found for this observation"
			
		## Comments/notes
		print "   Observer Comments: %s" % project.sessions[0].observations[i].comments


if __name__ == "__main__":
	main(sys.argv[1:])
	
