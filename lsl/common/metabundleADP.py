# -*- coding: utf-8 -*-

# Python3 compatiability
from __future__ import print_function
import sys
if sys.version_info > (3,):
	xrange = range
	long = int

"""
Module for working with an MCS meta-data tarball and extracting the useful bits out 
it and putting those bits into Python objects, e.g, :class:`lsl.common.stations.LWAStation` 
and :class:`lsl.common.sdm.SDM`.
"""

import os
import re
import copy
import glob
import shutil
import tarfile
import tempfile
from datetime import datetime, timedelta

from lsl.common import stations, sdmADP, sdfADP
from lsl.common.mcsADP import *
from lsl.common.adp import word2freq
from lsl.transform import Time
from lsl.misc.lru_cache import lru_cache

__version__ = '1.0'
__revision__ = '$Rev$'
__all__ = ['readSESFile', 'readOBSFile', 'readCSFile', 'getSDM', 'getStation', 'getSessionMetaData', 
		 'getSessionSpec', 'getObservationSpec', 'getSessionDefinition', 'getCommandScript', 
		 'getASPConfiguration', 'getASPConfigurationSummary', 'isValid', 
		 '__version__', '__revision__', '__all__']

# Regular expression for figuring out filenames
filenameRE = re.compile(r'(?P<projectID>[a-zA-Z0-9]{1,8})_(?P<sessionID>\d+)(_(?P<obsID>\d+)(_(?P<obsOutcome>\d+))?)?.*\..*')


@lru_cache(maxsize=1)
def _openTarball(tarname):
	return tarfile.open(tarname, mode='r:*')


@lru_cache(maxsize=1)
def _getMembers(tarname):
	tf = _openTarball(tarname)
	return tf.getmembers()


def readSESFile(filename):
	"""
	Read in a session specification file (MCS0030, Section 5) and return the data
	as a dictionary.
	"""
	
	# Read the SES
	with open(filename, 'rb') as fh:
		bses = parseCStruct(SSF_STRUCT, endianness='little')
		fh.readinto(bses)
		
		## LWA-1 check
		#if bses.FORMAT_VERSION not in (6,):
		#	fh.close()
		#	raise RuntimeError("Version mis-match: File appears to be from LWA-1")
			
		if bses.SESSION_NOBS > 150:
			## Pre SESSION_SPC
			fh.seek(0)
			
			newStruct = []
			for line in SSF_STRUCT.split('\n'):
				if line.find('SESSION_SPC') != -1:
					continue
				newStruct.append(line)
			newStruct = '\n'.join(newStruct)
			
			bses = parseCStruct(newStruct, endianness='little')
			fh.readinto(bses)
			bses.SESSION_SPC = ''
			
	record = {'ASP': bses.SESSION_MRP_ASP, 'ADP': bses.SESSION_MRP_DP_, 'SHL': bses.SESSION_MRP_SHL, 
			'MCS': bses.SESSION_MRP_MCS, 'DR1': bses.SESSION_MRP_DR1, 'DR2': bses.SESSION_MRP_DR2}
	
	update = {'ASP': bses.SESSION_MUP_ASP, 'ADP': bses.SESSION_MUP_DP_, 'SHL': bses.SESSION_MUP_SHL, 
			'MCS': bses.SESSION_MUP_MCS, 'DR1': bses.SESSION_MUP_DR1, 'DR2': bses.SESSION_MUP_DR2}
	
	return {'version': bses.FORMAT_VERSION, 'projectID': bses.PROJECT_ID.lstrip().rstrip(), 
		   'sessionID': bses.SESSION_ID,  'CRA': bses.SESSION_CRA,  'drxBeam': bses.SESSION_DRX_BEAM,
		   'spcSetup': bses.SESSION_SPC, 'MJD': bses.SESSION_START_MJD, 'MPM': bses.SESSION_START_MPM, 
		   'Dur': bses.SESSION_DUR, 'nObs': bses.SESSION_NOBS, 'record': record, 'update': update, 
		   'logSch': bses.SESSION_LOG_SCH, 'logExe': bses.SESSION_LOG_EXE, 'incSMIF': bses.SESSION_INC_SMIB,
		   'incDesi': bses.SESSION_INC_DES}


def readOBSFile(filename):
	"""
	Read in a observation specification file (MCS0030, Section 6) and return the
	data as a dictionary.
	"""
	
	# Read the OBS
	with open(filename, 'rb') as fh:
		bheader = parseCStruct(OSF_STRUCT, endianness='little')
		bstep   = parseCStruct(OSFS_STRUCT, endianness='little')
		bbeam   = parseCStruct(BEAM_STRUCT, endianness='little')
		bfooter = parseCStruct(OSF2_STRUCT, endianness='little')
		fh.readinto(bheader)
		
		## LWA-1 check
		#if bheader.FORMAT_VERSION not in (6,):
		#	fh.close()
		#	raise RuntimeError("Version mis-match: File appears to be from LWA-1")
			
		if bheader.OBS_ID > 150:
			## Pre SESSION_SPC and OBS_BDM
			fh.seek(0)
			
			newStruct = []
			for line in OSF_STRUCT.split('\n'):
				if line.find('OBS_BDM') != -1:
					continue
				if line.find('SESSION_SPC') != -1:
					continue
				newStruct.append(line)
			newStruct = '\n'.join(newStruct)
			
			bheader = parseCStruct(newStruct, endianness='little')
			fh.readinto(bheader)
			bheader.SESSION_SPC = ''
			bheader.OBS_BDM = ''
			
		elif bheader.OBS_B > 2:
			## Pre OBS_BDM
			fh.seek(0)
			
			newStruct = []
			for line in OSF_STRUCT.split('\n'):
				if line.find('OBS_BDM') != -1:
					continue
				newStruct.append(line)
			newStruct = '\n'.join(newStruct)
			
			bheader = parseCStruct(newStruct, endianness='little')
			fh.readinto(bheader)
			bheader.OBS_BDM = ''
			
		if IS_32BIT_PYTHON:
			skip = parseCStruct("""
			int junk;
			""", endianness='little')
			fh.readinto(skip)

		steps = []
		for n in xrange(bheader.OBS_STP_N):
			fh.readinto(bstep)
			if bstep.OBS_STP_B == 3:
				fh.readinto(bbeam)
				bstep.delay = copy.deepcopy(bbeam.OBS_BEAM_DELAY)
				bstep.gain  = copy.deepcopy(single2multi(bbeam.OBS_BEAM_GAIN, *bbeam.dims['OBS_BEAM_GAIN']))
			else:
				bstep.delay = []
				bstep.gain  = []
			
			steps.append(copy.deepcopy(bstep))
			
			alignment = parseCStruct("""
			unsigned int block;
			""", endianness='little')
			fh.readinto(alignment)
			
			if alignment.block != (2**32 - 2):
				raise IOError("Byte alignment lost at byte %i" % fh.tell())
				
		fh.readinto(bfooter)
		
		if bfooter.alignment != (2**32 - 1):
			raise IOError("Byte alignment lost at byte %i" % fh.tell())
			
	output = {'version': bheader.FORMAT_VERSION, 'projectID': bheader.PROJECT_ID.lstrip().rstrip(), 
		     'sessionID': bheader.SESSION_ID, 'drxBeam': bheader.SESSION_DRX_BEAM, 
		     'spcSetup': bheader.SESSION_SPC, 'obsID': bheader.OBS_ID,
		     'MJD': bheader.OBS_START_MJD, 'MPM': bheader.OBS_START_MPM, 'Dur': bheader.OBS_DUR, 
		     'Mode': bheader.OBS_MODE, 'beamDipole': bheader.OBS_BDM, 
		     'RA': bheader.OBS_RA, 'Dec': bheader.OBS_DEC, 'Beam': bheader.OBS_B, 
		     'Freq1': word2freq(bheader.OBS_FREQ1), 'Freq2': word2freq(bheader.OBS_FREQ2), 'BW': bheader.OBS_BW, 'nSteps': bheader.OBS_STP_N, 'StepRADec': bheader.OBS_STP_RADEC,  'steps': steps, 
		     'fee': single2multi(bfooter.OBS_FEE, *bfooter.dims['OBS_FEE']), 
		     'flt': list(bfooter.OBS_ASP_FLT), 'at1': list(bfooter.OBS_ASP_AT1), 
		     'at2': list(bfooter.OBS_ASP_AT2), 'ats': list(bfooter.OBS_ASP_ATS)}
	output['tbfSamples'] = bfooter.OBS_TBF_SAMPLES
	output['tbfGain'] = bfooter.OBS_TBF_GAIN
	output['tbnGain'] = bfooter.OBS_TBN_GAIN
	output['drxGain'] = bfooter.OBS_DRX_GAIN
	
	return output


def readCSFile(filename):
	"""
	Read in a command script file (MCS0030, currently undocumented) and return the
	data as a list of dictionaries.
	"""
	
	# Read the CS file
	with open(filename, 'rb') as fh:
		commands = []
		while True:
			action = parseCStruct("""
			long int tv[2];
			int bASAP;
			int sid;
			int cid;
			int len;
			""", endianness='little')
			
			try:
				fh.readinto(action)
				
				if action.tv[0] == 0:
					break
					
				if action.len > 0:
					data = parseCStruct("""
					char data[%i];
					""" % action.len, endianness='little')
					
					fh.readinto(data)
					data = data.data
				else:
					data = None
				
				actionPrime = {'time': action.tv[0] + action.tv[1]/1.0e6, 
							'ignoreTime': True if action.bASAP else False, 
							'subsystemID': sid2string(action.sid), 'commandID': cid2string(action.cid), 
							'commandLength': action.len, 'data': data}
							
				commands.append( actionPrime )
			except IOError:
				break
				
	return commands


def getSDM(tarname):
	"""
	Given an MCS meta-data tarball, extract the information stored in the 
	dynamic/sdm.dat file and return a :class:`lsl.common.sdm.SDM` instance
	describing the dynamic condition of the station.
	
	If a sdm.dat file cannot be found in the tarball, None is returned.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	
	# Extract the SDM file.  If the dynamic/sdm.dat file cannot be found, None
	# is returned via the try...except block.
	tf = _openTarball(tarname)
	try:
		ti = tf.getmember('dynamic/sdm.dat')
	except KeyError:
		return None
	tf.extractall(path=tempDir, members=[ti,])
	
	try:
		# Parse the SDM file and build the SDM instance
		dynamic = sdmADP.parseSDM(os.path.join(tempDir, 'dynamic', 'sdm.dat'))
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	return dynamic


def getStation(tarname, ApplySDM=True):
	"""
	Given an MCS meta-data tarball, extract the information stored in the ssmif.dat 
	file and return a :class:`lsl.common.stations.LWAStation` object.  Optionally, 
	update the :class:`lsl.common.stations.Antenna` instances associated whith the
	LWAStation object using the included SDM file.
	
	If a ssmif.dat file cannot be found in the tarball, None is returned.  
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	
	# Extract the SSMIF and SDM files.  If the ssmif.dat file cannot be found, None
	# is returned via the try...except block
	tf = _openTarball(tarname)
	try:
		ti = tf.getmember('ssmif.dat')
	except KeyError:
		return None
	tf.extractall(path=tempDir, members=[ti,])
	
	try:
		# Read in the SSMIF
		station = stations.parseSSMIF(os.path.join(tempDir, 'ssmif.dat'))
		
		# Get the SDM (if we need to)
		if ApplySDM:
			dynamic = getSDM(tarname)
		else:
			dynamic = None
		
		# Update the SSMIF entries
		if dynamic is not None:
			newAnts = dynamic.updateAntennas(station.getAntennas())
			station.antennas = newAnts
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return station


def getSessionMetaData(tarname):
	"""
	Given an MCS meta-data tarball, extract the session meta-data file (MCS0030, 
	Section 7) and return a dictionary of observations that contain dictionaries 
	of the OP_TAG (tag), DRSU Barcode (drsu), OBS_OUTCOME (outcome), and the 
	MSG (msg).
	
	.. versionchanged:: 0.6.5
		Update to the new _metadata.txt format
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Extract the session meta-data file (_metadata.txt)
	tf = _openTarball(tarname)
	try:
		ti = tf.getmember('%s_metadata.txt' % basename)
	except KeyError:
		for ti in _getMembers(tarname):
			if ti.name[-13:] == '_metadata.txt':
				break
	tf.extractall(path=tempDir, members=[ti,])
	
	try:
		# Read in the SMF
		filename = os.path.join(tempDir, ti.name)
		with open(filename, 'r') as fh:
			# Define a regular expresion to match the latest format
			lineRE = re.compile(r"\s*(?P<id>\d{1,}?)\s+\[(?P<tag>[\d_]+?)\]\s+\['?(?P<barcode>.+?)'?\]\s+(?P<outcome>\d)\s+\[(?P<msg>.*?)\]")
			
			result = {}
			for line in fh:
				line = line.replace('\n', '')
				if len(line) == 0:
					continue
					
				mtch = lineRE.search(line)
				if mtch is not None:
					## If it matches the new format
					obsID = mtch.group('id')
					opTag = mtch.group('tag')
					drsuBarcode = mtch.group('barcode')
					if drsuBarcode[:3] == 'Err':
						try:
							drsuBarcode = result[int(obsID)-1]['barcode']
						except KeyError:
							drsuBarcode = 'UNK'
					obsOutcome = mtch.group('outcome')
					msg = mtch.group('msg')
					
				else:
					## Otherwise, I don't really know how the messages will look so we use this try...except
					## block should take care of the various situations.
					try:
						obsID, opTag, drsuBarcode, obsOutcome, msg = line.split(None, 4)
						opTag = opTag.replace('[', '')
						opTag = opTag.replace(']', '')
						drsuBarcode = drsuBarcode.replace('[', '')
						drsuBarcode = drsuBarcode.replace(']', '')
						drsuBarcode = drsuBarcode.replace("'", '')
					except ValueError:
						try:
							obsID, opTag, drsuBarcode, obsOutcome = line.split(None, 3)
							msg = 'UNK'
						except ValueError:
							try:
								obsID, opTag, obsOutcome = line.split(None, 2)
								drsuBarcode = 'UNK'
								obsOutcome = '-1'
								msg = 'UNK'
							except ValueError:
								obsID, obsOutcome = line.split(None, 1)
								drsuBarcode = 'UNK'
								opTag = 'UNK'
								msg = 'UNK'
								
				obsID = int(obsID)
				obsOutcome = int(obsOutcome) if obsOutcome != 'Failed' else 1
				result[obsID] = {'tag': opTag, 'barcode': drsuBarcode, 'outcome': obsOutcome, 'msg': msg}
				
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return result


def getSessionSpec(tarname):
	"""
	Given an MCS meta-data tarball, extract the session specification file (MCS0030, 
	Section 5) and return a dictionary of parameters.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Extract the session specification file (.ses)
	tf = _openTarball(tarname)
	try:
		ti = tf.getmember('%s.ses' % basename)
	except KeyError:
		for ti in _getMembers(tarname):
			if ti.name[-4:] == '.ses':
				break
	tf.extractall(path=tempDir, members=[ti,])
	
	try:
		# Read in the SES
		ses = readSESFile(os.path.join(tempDir, ti.name))
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return ses


def getObservationSpec(tarname, selectObs=None):
	"""
	Given an MCS meta-data tarball, extract one or more observation specification 
	file (MCS0030, Section 6) and return a list of dictionaries corresponding to
	each OBS file.  If the `selectObs` keyword is set to a list of observation
	numbers, only observations matching the numbers in `selectObs` are returned.
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Find all of the .obs files and extract them
	tf = _openTarball(tarname)
	tis = []
	for ti in _getMembers(tarname):
		if ti.name[-4:] == '.obs':
			tis.append(ti)
	tf.extractall(path=tempDir, members=tis)
	
	try:
		# Read in the OBS files
		obsList = []
		for of in sorted(glob.glob(os.path.join(tempDir, '*.obs'))):
			obsList.append( readOBSFile(of) )
			
		# Cull the list based on the observation ID selection
		if selectObs is not None:
			outObs = []
			for o in obsList:
				try:
					if o['obsID'] in selectObs:
						outObs.append(o)
				except TypeError:
					if o['obsID'] == selectObs:
						outObs.append(o)
						
			if len(outObs) == 1:
				outObs = outObs[0]
		else:
			outObs = obsList
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return outObs


def getSessionDefinition(tarname):
	"""
	Given an MCS meta-data tarball, extract the session specification file, the 
	session meta-data file, and all observation specification files to build up
	a SDF-representation of the session.
	
	.. note::
		This function returns a full :class:`lsl.common.sdfADP.Project` instance 
		with the session in question stored under `project.sessions[0]` and the 
		observations under `project.sessions[0].observations`.
	"""
	
	# Find the SDF file contained in the tarball
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Find the right .txt file (not the metadata one) and extract it
	tf = _openTarball(tarname)
	for ti in _getMembers(tarname):
		if ti.name[-4:] == '.txt' and ti.name.find('metadata') == -1:
			break
	tf.extractall(path=tempDir, members=[ti,])
	
	try:
		# Parse it
		project = sdfADP.parseSDF(os.path.join(tempDir, ti.name))
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Clean up
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return the filled-in SDF instance
	return project


def getCommandScript(tarname):
	"""
	Given an MCS meta-data tarball, extract the command script and parse it.  The
	commands are returned as a list of dictionaries (one dictionary per command).
	"""
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Find the .cs file and extract it
	tf = _openTarball(tarname)
	for ti in _getMembers(tarname):
		if ti.name[-3:] == '.cs':
			break
	tf.extractall(path=tempDir, members=[ti,])
	
	try:
		# Read in the CS
		cs = readCSFile(os.path.join(tempDir, ti.name))
	except Exception as e:
		shutil.rmtree(tempDir, ignore_errors=True)
		raise e
		
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	# Return
	return cs


def getASPConfiguration(tarname, which='beginning'):
	"""
	Given an MCS meta-data tarball, extract the ASP MIB contained in it and return 
	a dictionary of values for the filter, AT1, AT2, and ATSplit.  The 'which'
	keyword is used to specify whether or not the configuration returned is at the
	beginning (default) or end of the session.
	
	.. versionadded:: 0.6.5
	"""
	
	which = which.lower()
	if which not in ('beginning', 'begin', 'end'):
		raise ValueError("Unknown configuration time '%s'" % which)
		
	# Stub ASP configuration
	aspConfig = {'filter': [-1 for i in xrange(264)],
			   'at1': [-1 for i in xrange(264)],
			   'at2': [-1 for i in xrange(264)],
			   'atsplit': [-1 for i in xrange(264)]}
	
	tempDir = tempfile.mkdtemp(prefix='metadata-bundle-')
	path, basename = os.path.split(tarname)
	basename, ext = os.path.splitext(basename)
	
	# Find the .dir/.pag file and extract it
	tf = _openTarball(tarname)
	mibs = []
	for ti in _getMembers(tarname):
		if ti.name.find('_ASP_%s' % which[:5]) != -1:
			if ti.name[-4:] == '.pag' or ti.name[-4:] == '.dir':
				mibs.append(ti)
				
	if len(mibs) > 0:
		tf.extractall(path=tempDir, members=mibs)
		
		try:
			# Read in the right MIB
			aspMIB = {}
			for mib in mibs:
				if mib.name[-4:] == '.dir':
					continue
				if which[:5] == 'begin' and mib.name.find('_ASP_begin') == -1:
					continue
				if which == 'end' and mib.name.find('_ASP_end') == -1:
					continue
					
				aspMIB = MIB()
				aspMIB.fromFile(os.path.join(tempDir, mib.name))
				break
				
			# Extract the configuration
			for key in aspMIB.keys():
				values = [int(v) for v in key.split('.')]
				
				if values[0] == 3:
					## Filter
					aspConfig['filter'][values[1]-1] = int(aspMIB[key].value)
					
				elif values[0] == 4:
					## Attenuators
					if values[1] == 1:
						### AT1
						aspConfig['at1'][values[2]-1] = int(aspMIB[key].value)
					elif values[1] == 2:
						### AT2
						aspConfig['at2'][values[2]-1] = int(aspMIB[key].value)
					elif values[1] == 3:
						### ATSPLIT
						aspConfig['atsplit'][values[2]-1] = int(aspMIB[key].value)
					else:
						pass
						
				else:
					pass
		except Exception as e:
			shutil.rmtree(tempDir, ignore_errors=True)
			raise e
			
	# Cleanup
	shutil.rmtree(tempDir, ignore_errors=True)
	
	return aspConfig


def getASPConfigurationSummary(tarname, which='beginning'):
	"""
	Similar to getASPConfiguration, but returns only a single value for each
	of the four ASP paramters:  filter, AT, AT2, and ATSplit.  The values
	are based off the mode of the parameter.
	
	.. versionadded:: 0.6.5
	"""
	
	# Get the full configuration
	aspConfig = getASPConfiguration(tarname, which=which)
	
	# Count
	count = {}
	for param in aspConfig.keys():
		count[param] = {}
		for ant in xrange(len(aspConfig[param])):
			value = aspConfig[param][ant]
			try:
				count[param][value] += 1
			except KeyError:
				count[param][value] = 1
				
	# Modes
	mode = {}
	for param in count.keys():
		best = 0
		mode[param] = 0
		
		for value in count[param].keys():
			num = count[param][value]
			if num > best:
				best = num
				mode[param] = value
				
	# Done
	return mode


def isValid(tarname, verbose=False):
	"""
	Given a filename, see if it is valid metadata tarball or not.
	
	.. versionadded:: 1.2.0
	"""
	
	passes = 0
	failures = 0
	try:
		getSessionSpec(tarname)
		passes += 1
		if verbose:
			print("Session specification - OK")
	except IOError as e:
		raise e
	except:
		failures += 1
		if verbose:
			print("Session specification - FAILED")
		
	try:
		getObservationSpec(tarname)
		passes += 1
		if verbose:
			print("Observation specification(s) - OK")
	except:
		failures += 1
		if verbose:
			print("Observation specification(s) - FAILED")
			
	try:
		getCommandScript(tarname)
		passes += 1
		if verbose:
			print("Command script - OK")
	except:
		failures += 1
		if verbose:
			print("Command script - FAILED")
			
	if verbose:
		print("---")
		print("%i passed / %i failed" % (passes, failures))
		
	return False if failures else True