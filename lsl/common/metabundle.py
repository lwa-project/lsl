# -*- coding: utf-8 -*-

"""
Module for working with a MCS meta-data tarball and extracting the useful bits out 
it and putting those bits into Python objects, e.g, :class:`lsl.common.stations.LWAStation` 
and :class:`lsl.common.sdm.SDM`.
"""

import os
import re
import glob
import shutil
import tarfile
import tempfile

from lsl.common.dp import fS
from lsl.common import stations, sdm

__version__ = "0.1"
__revision__ = "$ Revision: 1 $"
__all__ = ['readSESFile', 'readOBSFile', 'getSDM', 'getStation', 'getSessionMetaData', 'getSessionSpec', 'getObservationSpec', '__version__', '__revision__', '__all__']

# Regular expression for figuring out filenames
filenameRE = re.compile(r'(?P<projectID>[a-zA-Z0-9]{1,8})_(?P<sessionID>\d+)(_(?P<obsID>\d+)(_(?P<obsOutcome>\d+))?)?.*\..*')


def __guidedBinaryRead(fh, fmt):
	"""
	Function to wrap reading in packed binary data directrly from an open file
	handle.  This function calls struct.unpack() and struct.calcsize() to figure 
	out what to read and how.
	
	Return either a single item if a single item is requested or a list of items.
	
	Used by SDM.binaryFill()
	"""
	
	
	data = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))
	if len(data) == 1:
		return data[0]
	else:
		return list(data)


def __word2freq(word):
	"""
	Convert a LWA DP tuning word to a frequency in Hz.
	"""
	
	return word*fS / 2**32


def readSESFile(filename):
	"""
	Read in a session specification file (MCS0030, Section 5) and return the data
	as a dictionary.
	"""
	
	# Read the SES
	fh = open(filename, 'rb')

	version = __guidedBinaryRead(fh, "<H")
	projectID, sessionID = __guidedBinaryRead(fh, "<9sI")
	cra = __guidedBinaryRead(fh, "<h")
	drxBeam = __guidedBinaryRead(fh, "<h")
	sessionMJD, sessionMPM = __guidedBinaryRead(fh, "<QQ")
	sessionDur = __guidedBinaryRead(fh, "<Q")
	sessionObs = __guidedBinaryRead(fh, "<I")
	
	record = {}
	for k in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
		record[k] = __guidedBinaryRead(fh, "<h")
	update = {}
	for k in ['ASP', 'DP_', 'DR1', 'DR2', 'DR3', 'DR4', 'DR5', 'SHL', 'MCS']:
		update[k] = __guidedBinaryRead(fh, "<h")
	
	logSch, logExe = __guidedBinaryRead(fh, "<bb")
	incSMIF, incDesi = __guidedBinaryRead(fh, "<bb")

	fh.close()
	
	return {'version': version, 'projectID': projectID, 'sessionID': sessionID, 'CRA': cra, 
		   'MJD': sessionmjd, 'MPM': sessionMPM, 'Dur': sessionDur, 'nObs': sessionObs, 
		   'record': record, 'update': update, 'logSch': logSch, 'logExe': logExe, 
		   'incSMIF': incSMIF, 'incDesi': incDesi}


def readOBSFile(filename):
	"""
	Read in a observation specification file (MCS0030, Section 6) and return the
	data as a dictionary.
	"""
	
	# Read the OBS
	fh = open(filename, 'rb')
	
	version = __guidedBinaryRead(fh, "<H")
	projectID, sessionID = __guidedBinaryRead(fh, "<9sI")
	obsID = __guidedBinaryRead(fh, "<I")
	obsMJD, obsMPM = __guidedBinaryRead(fh, "<QQ")
	obsDur = __guidedBinaryRead(fh, "<Q")
	obsMode = __guidedBinaryRead(fh, "<H")
	obsRA, obsDec = __guidedBinaryRead(fh, "<ff")
	obsB = __guidedBinaryRead(fh, "<H")
	obsFreq1, obsFreq2, obsBW = __guidedBinaryRead(fh, "<IIH")
	nSteps, stepRADec = __guidedBinaryRead(fh, "<IH")

	steps = []
	for n in xrange(nSteps):
		c1, c2 = __guidedBinaryRead(fh, "<ff")
		t = __guidedBinaryRead(fh, "<I")
		f1, f2 = __guidedBinaryRead(fh, "<II")
		b = __guidedBinaryRead(fh, "<H")
		if b == 3:
			delay = __guidedBinaryRead(fh, "<520H")
			gain = __guidedBinaryRead(fh, "<1040h")
		else:
			delay = []
			gain = []
		
		steps.append([c1, c2, t, f1, f2, b, delay, gain])
		
		alignment = __guidedBinayRead(fh, "<I")
		if alignment != (2**32 - 2):
			raise IOError("Bytes alignment lost at bytes %i" % fh.tell())

	fee = __guidedBinaryRead(fh, "<520h")
	flt = __guidedBinaryRead(fh, "<260h")
	at1 = __guidedBinaryRead(fh, "<260h")
	at2 = __guidedBinaryRead(fh, "<260h")
	ats = __guidedBinaryRead(fh, "<260h")

	tbwBits, tbwSamps = __guidedBinaryRead(fh, "<HI")
	tbnGain, drxGain = __guidedBinaryRead(fh, "<hh")
	alignment = __guidedBinaryRead(fh, "<I")
	if aligment != (2**32 - 1):
		raise IOError("Bytes alignment lost at bytes %i" % fh.tell())

	fh.close()
	
	return {'version': version, 'projectID': projectID, 'sessionID': sessionID, 'obsID': obsID, 'MJD': obsMJD, 'MPM': obsMPM, 
		   'Dur': obsDur, 'Mode': obsMode, 'RA': obsRA, 'Dec': obsDec, 'Beam': obsB, 'Freq1': __word2freq(obsFreq1), 
		   'Freq2': __word2freq(obsFreq2), 'BW': obsBW, 'nSteps': nSteps, 'StepRADec': stepRADec,  'steps': steps, 
		   'fee': fee, 'flt': flt, 'at1': at1, 'at2': at2, 'ats': ats}


def getSDM(tarname):
	"""
	Given a MCS meta-data tarball, extract the information stored in the 
	dynamic/sdm.dat file and return a :class:`lsl.common.sdm.SDM` instance
	describing the dynamic condition of the station.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	
	# Extract the SDM file
	tf = tarfile.open(tarname, mode='r:gz')
	ti = fh.getmember('dynamic/sdm.dat')
	tf.extractall(path=tempDir, members=[ti,])
	
	# Parse the SDM file and build the SDM instance
	dynamic = sdm.SDM(filename=os.path.join(tempDir, 'dynamic', 'sdm.dat'))
	
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	return dynamic


def getStation(tarname, ApplySDM=True):
	"""
	Given a MCS meta-data tarball, extract the information stored in the ssmif.dat 
	file and return a :class:`lsl.common.stations.LWAStation` object.  Optionaly, 
	update the :class:`lsl.common.stations.Antenna` instances associated whith the
	LWAStation object using the included SDM file.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	
	# Extract the SSMIF and SDM files
	tf = tarfile.open(tarname, mode='r:gz')
	ti = fh.getmember('ssmif.dat')
	tf.extractall(path=tempDir, members=[ti,])
	
	# Read in the SSMIF
	station = stations.parseSSMIF(os.path.join(tempDir, 'ssmif.dat'))
	
	# Get the SDM (if we need to)
	if ApplySDM:
		try:
			dynamic = getSDM(tarname)
		except:
			dynamic = None
	else:
		dynamic = None
	
	# Update the SSMIF entries
	if dynamic is not None:
		newAnts = dynamic.updateAntennas(station.getAntennas())
		station.antennas = newAnts
	
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return station


def getSessionMetaData(tarname):
	"""
	Given a MCS meta-data tarball, extract the session meta-data file (MCS0030, 
	Section 7) and return a dictionary of observations that contain dictionaries 
	of the OP_TAG (tag), OBS_OUTCOME (outcome), and the MSG (msg).
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Extract the session meta-data file
	tf = tarfile.open(tarname, mode='r:gz')
	ti = fh.getmember('%s_metadata.txt' % basename)
	tf.extractall(path=tempDir, members=[ti,])
	
	# Read in the SMF
	fh = open(filename, 'r')

	result = {}
	for line in fh:
		line = line.replace('\n', '')
		if len(line) == 0:
			continue

		## I don't really know how the messages will look so we use this try...except
		## block should take care of the various situations.
		try:
			obsID, opTag, obsOutcome, msg = line.split(None, 3)
		except ValueError:
			obsID, opTag, obsOutcome = line.split(None, 2)
			msg = ''

		obsID = int(obsID)
		obsOutcome = int(obsOutcome)
		result[obsID] = {'tag': opTag, 'outcome': obsOutcome, 'msg': msg}
		
	fh.close()
	
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return result


def getSessionSpec(tarname):
	"""
	Given a MCS meta-data tarball, extract the session specification file (MCS0030, 
	Section 5) and return a dictionary of parameters.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Extract the session specification file
	tf = tarfile.open(tarname, mode='r:gz')
	ti = fh.getmember('%s.ses' % basename)
	tf.extractall(path=tempDir, members=[ti,])
	
	# Read in the SES
	ses = readSESFile(os.path.join(tempDir, '%s.ses' % basename))
	
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return ses


def getObservationSpec(tarname, selectObs=None):
	"""
	Given a MCS meta-data tarball, extract one or more observation specification 
	file (MCS0030, Section 6) and return a list of dictionaries corresponding to
	each OBS file.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Find all of the obs files and extract them
	tf = tarfile.open(tarname, mode='r:gz')
	tis = []
	for ti in tf.getmembers():
		if ti.name.find('.obs') != -1:
			tis.append(ti)
	tf.extractall(path=tempDir, members=tis)
	
	# Read in the OBS files
	obs = []
	for of in glob.glob(tempDir, '*.obs'):
		obs.append( readOBSFile(of) )
		
	# Cull the list based on the observation ID selection
	if obs is not None:
		outObs = []
		for o in obs:
			if o['obsID'] in selectObs:
				outObs.append(o)
	else:
		outObs = obs
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return outObs
