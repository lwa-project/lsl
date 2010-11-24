# -*- coding: utf-8 -*-

import sys
import numpy
try:
	import cStringIO as StringIO
except ImportError:
	import StringIO

from lsl.misc import geodesy
from lsl import astro as astro
from lsl.common import stations as lwa_common
from lsl.correlator import uvUtils as uvUtils
from lsl.reader import tbw as tbw
from lsl.reader import tbn as tbn
from lsl.reader import errors as errors
from lsl.writer import vdif


class ConfigError(Exception):
	"""Base exception class for dealing with errors in creating the DiFX 
	configuration files."""

	def __init__(self,  strerror, errno='-1'):
		self.errno = errno
		self.strerror = strerror
		self.args = (errno, strerror)

	def __str__(self):
		return "%s" % self.strerror


class DiFXConfig(object):
	"""Object for creating DiFX input, calc, flag, and joblist configuration 
	files.  The class is initialized with a set of common observation parameters
	and the output files are written using the class's 'write' function.  By 
	default, the outfile files are:
	  lwa_<mode>_<mjd>.joblist
	  lwa_<mode>_<mjd>_1.input
	  lwa_<mode>_<mjd>_1.calc
	  lwa_<mode>_<mjd>_1.flag"""

	def __init__(self, site=lwa_common.lwa1(), stands=[], mode='TBW', refTime=0.0, intTime=0, nInts=1, bits=12, filter=7, centralFreq=49.0e6):
		"""Initialize all of the values needed to create a complete set of DiFX 
		configuration files.  Required parameters are:
		  + LWA station object, 
		  + list of stands, 
		  + observations mode (TBW or TBN), 
		  + time of first observation, 
		  + integration time in seconds, and
		  + number of integration times to compute.
		For TBW data, the additional information is required:
		  + number of bits in the data (12 or 4).
		For TBN data, the additional information is required:
		  + TBN filter code and 
		  + central frequency in Hz.
		"""

		# Basic info
		self.site = site
		self.stands = stands
		self.mode = mode.upper()
		self.refTime = self.parseRefTime(refTime)
		self.intTime = intTime
		self.nInts = 1

		# Make sure we have a valid mode
		if self.mode not in ['TBW', 'TBN']:
			raise ConfigError("Unknown LWA observing mode '%s' " % self.mode)

		# Default file basename creation
		refMJD = refJD - astro.MJD_OFFSET
		self.basename = "lwa_%s_%i" % (self.mode.lower(), int(refMJD))

		# Mode-specific frequency, bandwidth, and data setup
		if self.mode == 'TBW':
			self.centralFreq = 49.0e6
			self.bandwidth = 98e6
			self.bits = bits
		else:
			self.centralFreq = centralFreq
			self.bandwidth = tbn.filterCodes[filter]
			self.bits = 8

	def setBasename(self, basename):
		"""Override the default file basename of "lwa_<mode>_<MJD>" with a user-
		specified string."""

		self.basename = basename

	def parseRefTime(self, refTime):
		"""Given a time as either a integer, float, string, or datetime object, 
		convert it to a string in the formation 'YYYY-MM-DDTHH:MM:SS'."""

		# Valid time string (modulo the 'T')
		timeRE = re.compile(r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?')

		if type(refTime).__name__ in ['int', 'long', 'float']:
			refDateTime = datetime.utcfromtimestamp(refTime)
			refTime = refDateTime.strftime("%Y-%m-%dT%H:%M:%S")
		elif type(refTime).__name__ in ['datetime']:
			refTime = refTime.strftime("%Y-%m-%dT%H:%M:%S")
		elif type(refTime).__name__ in ['str']:
			# Make sure that the string times are of the correct format
			if re.match(timeRE, refTime) is None:
				raise RuntimeError("Malformed date/time provided: %s" % refTime)
			else:
				refTime = refTime.replace(' ', 'T', 1)
		else:
			raise RuntimeError("Unknown time format provided.")

		return refTime

	def refTime2AstroDate(self):
		"""Convert a reference time string to an astro.date object."""

		dateStr = self.refTime.replace('T', '-').replace(':', '-').split('-')
		return astro.date(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), int(dateStr[3]), int(dateStr[4]), float(dateStr[5]))

	def getJoblist(self):
		"""Create a string that contains the list of DiFX jobs that needed to be 
		run."""

		refDate = self.refTime2AstroDate()
		refJD = refDate.to_jd()
		refMJD = refJD - astro.MJD_OFFSET
		stpMJD = refMJD + (self.intTime*self.nInts/astro.SECS_IN_DAY)

		fh = StringIO.StringIO()
		fh.write('exper=%s  v2d=lwa_%s.v2d  pass=%s  mjd=%.7f DiFX=DiFX-1.5.4  vex2difx=1.0.4\n' % (self.mode, self.mode.lower(), self.mode, refMJD))
		fh.write('%s_1 %.7f %.7f 3 6.053 10.91  # %s\n' % (self.basename, refMJD, stpMJD, ' '.join(['LWA%03i' % i for i in self.stands])))
	
		joblistString = fh.getvalue()
		fh.close()

		return joblistString

	def getInput(self):
		"""Create a string that contains a DiFX input file."""

		refDate = self.refTime2AstroDate()
		refJD = refDate.to_jd()
		refMJD = refJD - astro.MJD_OFFSET
		stpMJD = refMJD + (self.intTime*self.nInts/astro.SECS_IN_DAY)

		fh = StringIO.String()
		fh.write('# COMMON SETTINGS ##!\n')
		fh.write('DELAY FILENAME:     %s_1.delay\n' % self.basename)
		fh.write('UVW FILENAME:       %s_1.uvw\n' % self.basename)
		fh.write('CORE CONF FILENAME: %s_1.threads\n' % self.basename)
		fh.write('EXECUTE TIME (SEC): 60\n')
		fh.write('START MJD:          %i\n' % int(refMJD))
		fh.write('START SECONDS:      %i\n' % int((refMJD-int(refMJD))*astro.SECS_IN_DAY))
		fh.write('ACTIVE DATASTREAMS: %i\n' % len(self.stands))
		fh.write('ACTIVE BASELINES:   %i\n' % len(self.stands))
		fh.write('VIS BUFFER LENGTH:  %i\n' % 32)
		fh.write('OUTPUT FORMAT:      SWIN\n')
		fh.write('OUTPUT FILENAME:    %s.difx\n' % self.basename)
		fh.write(' \n')
	
		fh.write('# CONFIGURATIONS ###!\n')
		fh.write('NUM CONFIGURATIONS: 1\n')
		fh.write('CONFIG SOURCE:      %s_default\n' % self.mode)
		fh.write('INT TIME (SEC):     %.4f\n' % self.intTime)
		fh.write('NUM CHANNELS:       512\n')
		fh.write('CHANNELS TO AVERAGE:1\n')
		fh.write('OVERSAMPLE FACTOR:  1\n')
		fh.write('DECIMATION FACTOR:  1\n')
		fh.write('BLOCKS PER SEND:    3125\n')
		fh.write('GUARD BLOCKS:       1\n')
		fh.write('POST-F FRINGE ROT:  FALSE\n')
		fh.write('QUAD DELAY INTERP:  TRUE\n')
		fh.write('WRITE AUTOCORRS:    TRUE\n')
		fh.write('PULSAR BINNING:     FALSE\n')

		for i in list(range(len(self.stands))):
			fh.write('DATASTREAM %i INDEX: %i\n' % (i, i))
		for i in list(range((len(self.stands)*(len(self.stands)+1)/2))):
			fh.write('BASELINE %i INDEX:   %i\n' % (i, i))
		fh.write(' \n')
	
		fh.write('# FREQ TABLE #######!\n')
		fh.write('FREQ ENTRIES:       1\n')
		fh.write('FREQ (MHZ) 0:       %.5f\n' % (self.centralFreq-self.bandwidth/2.0))
		fh.write('BW (MHZ) 0:         %.5f\n' % self.bandwidth)
		fh.write('SIDEBAND 0:         U\n')
		fh.write(' \n')
	
		fh.write('# TELESCOPE TABLE ##!\n')
		fh.write('TELESCOPE ENTRIES:  %i\n' % len(stands))
		i = 0
		for stand in self.stands:
			fh.write('TELESCOPE NAME %i:   LWA%03i\n' % (i, stand))
			fh.write('CLOCK DELAY (us) %i: %.9f\n' % (i, uvUtils.cableDelay(stand, freq=self.centralFreq)*1e6))
			fh.write('CLOCK RATE(us/s) %i: %.9f\n' % (i, 0.0*1e6))
			i = i + 1
		fh.write(' \n')
	
		fh.write('# DATASTREAM TABLE #!\n')
		fh.write('DATASTREAM ENTRIES: %i\n' % len(stands))
		fh.write('DATA BUFFER FACTOR: 32\n')
		fh.write('NUM DATA SEGMENTS:  %i\n' % len(stands))
		i = 0
		for stand in self.stands:
			fh.write('TELESCOPE INDEX:    %i\n' % i)
			fh.write('TSYS:               %.4f\n' % 0.0)
			fh.write('DATA FORMAT:        VDIF\n')
			fh.write('QUANTISATION BITS:  %i\n' % self.bits)
			fh.write('DATA FRAME SIZE:    40004096\n')
			fh.write('DATA SOURCE:        FILE\n')
			fh.write('FILTERBANK USED:    FALSE\n')
			fh.write('NUM FREQS:          1\n')
			fh.write('FREQ TABLE INDEX 0: 0\n')
			fh.write('CLK OFFSET 0 (us):  0.000000\n')
			fh.write('NUM POLS 0:         1\n')
			fh.write('INPUT BAND 0 POL:   X\n')
			fh.write('INPUT BAND 0 INDEX: 0\n')
			i = i + 1
		fh.write(' \n')
	
		fh.write('# BASELINE TABLE ###!\n')
		fh.write('BASELINE ENTRIES:   %i\n' % (len(stands)*(len(stands)+1)/2))
		c = 0
		for i in list(range((len(stands)-1))):
			for j in list(range(1,len(stands))):
				if i == j:
					continue
				fh.write('D/STREAM A INDEX %i: %i\n' % (c, i))
				fh.write('D/STREAM B INDEX %i: %i\n' % (c, j))
				fh.write('NUM FREQS %i:        1\n' % c)
				fh.write('POL PRODUCTS %i/0:   1\n' % c)
				for k in list(range(len(stands))):
					fh.write('D/STREAM A BAND %i: 0\n' % k)
					fh.write('D/STREAM B BAND %i: 0\n' % k)
				c = c + 1
		fh.write(' \n')
	
		fh.write('# DATA TABLE #######!\n')
		i = 0
		for stand in self.stands:
			fh.write('D/STREAM %i FILES:   1\n' % i)
			fh.write('FILE %i/0:           %s-%i.vdif\n' % (i, self.mode, stand))
			i = i + 1
		fh.write(' \n')

		inputString = fh.getvalue()
		fh.close()

		return inputString

	def getCalc(self):
		"""Create a string that contains the bulk of the DiFX calc file."""

		refDate = self.refTime2AstroDate()
		refJD = refDate.to_jd()
		refMJD = refJD - astro.MJD_OFFSET
		stpMJD = refMJD + (self.intTime*self.nInts/astro.SECS_IN_DAY)
		xyzs = uvUtils.getXYZ(self.stands)
		
		fh = StringIO.StringIO()
		fh.write('JOB ID:             %i\n' % 1)
		fh.write('JOB START TIME:     %.7f\n' % refMJD)
		fh.write('JOB STOP TIME:      %.7f\n' % stpMJD)
		fh.write('DUTY CYCLE:         1.000\n')
		fh.write('OBSCODE:            %s\n' self.mode)
		fh.write('DIFX VERSION:       DiFX-1.5.4\n')
		fh.write('SUBJOB ID:          0\n')
		fh.write('SUBARRAY ID:        0\n')
		fh.write('START MJD:          %.7f\n' % refMJD)
		fh.write('START YEAR:         %i\n' % refDate.years)
		fh.write('START MONTH:        %i\n' % refDate.months)
		fh.write('START DAY:          %i\n' % refDat.days)
		fh.write('START HOUR:         %i\n' % refDat.hours)
		fh.write('START MINUTE:       %i\n' % refDat.minutes)
		fh.write('START SECOND:       %i\n' % int(refDat.seconds))
		fh.write('INCREMENT (SECS):   %.f\n' % self.intTIme)
		fh.write('SPECTRAL AVG:       1\n')
		fh.write('TAPER FUNCTION:     UNIFORM\n')
		fh.write('NUM TELESCOPES:     %i\n' % len(self.stands))
		
		i = 0
		for stand in self.stands:
			siteECI = site.getGeocentricLocation()
			standECI = siteECI + numpy.dot(site.getECEFTransform(), xyzs[i,:])
			
			fh.write('TELESCOPE %i NAME:   LWA%03i\n' % (i, stand))
			fh.write('TELESCOPE %i MOUNT:  azel\n' % i)
			fh.write('TELESCOPE %i OFFSET (m):0.000000\n' % i)
			fh.write('TELESCOPE %i X (m):  %.6f\n' % (i, standECI[0]))
			fh.write('TELESCOPE %i Y (m):  %.6f\n' % (i, standECI[1]))
			fh.write('TELESCOPE %i Z (m):  %.6f\n' % (i, standECI[2]))
			fh.write('TELESCOPE %i SHELF:  NONE\n' % i)
			i = i + 1

		fh.write('NUM SCANS:          1\n')
		fh.write('SCAN 0 POINTS:      10\n')
		fh.write('SCAN 0 START PT:    0\n')
		fh.write('SCAN 0 SRC NAME:    %s_default\n' % self.mode)
		fh.write('SCAN 0 REAL NAME:   ZENITH\n')
		fh.write('SCAN 0 SRC RA:      %.9f\n' % 0.0)
		fh.write('SCAN 0 SRC DEC:     %.9f\n' % 0.0)
		fh.write('SCAN 0 CALCODE:     \n')
		fh.write('SCAN 0 QUAL:        0\n')

		refDate = self.refTime2AstroDate()
		refMJD = refDate.to_jd() - astro.MJD_OFFSET
		eops = geodesy.getEOPRange(refMJD, stpMJD)
		if eops is None:
			eops = [geodesy.EOP(mjd=refMJD)]

		fh.write('NUM EOP:            %i\n' % len(eops))
		i = 0
		for eop in eops:
			fh.write('EOP %i TIME (mjd):   %i\n' % (i, int(eop.mjd)))
			fh.write('EOP %i TAI_UTC (sec):%i\n' % (i, astro.leap_secs(refJD)))
			fh.write('EOP %i UT1_UTC (sec):%.6f\n' % (i, eop.utDiff))
			fh.write('EOP %i XPOLE (arcsec):%.6f\n' % (i, eop.x))
			fh.write('EOP %i YPOLE (arcsec):%.6f\n' % (i, eop.y))
			i = i + 1
		fh.write('NUM SPACECRAFT:     0\n')

		fh.write('DELAY FILENAME:     %s_1.delay\n' % self.basename)
		fh.write('UVW FILENAME:       %s_1.uvw\n' % self.basename)
		fh.write('RATE FILENAME:      %s_1.rate\n' % self.basename)
		fh.write('IM_FILENAME:        %s_1.im\n' % self.basename)
		
		calcString = fh.getvalue()
		fh.close()

		return calcString

	def getFlag(self):
		"""Create a string that contains the flag file associated with this DiFX 
		configuration."""

		fh = StringIO.StringIO()
		fh.write('0')

		flagString = fh.getvalue()
		fh.close()
	
		return flagString

	def write(self):
		inputFile = "%s_1.input" % self.basename
		inputCont = self.getInput()
		
		calcFile = "%s_1.calc" % self.basename
		calcCont = self.getCalc()
		
		flagFile = "%s_1.flag" % self.basename
		flagCont = self.getFlag()

		jobsFile = "%s.joblist" % self.basename
		jobsCont = self.getJoblist()

		names = [inputFile, calcFile, flagFile, jobsFile]
		conts = [inputCont, calcCont, flagCont, jobsCont]
		for name,cont = zip(names, conts):
			print "Writing %i lines to file '%s'" % (len(cont.split('\n')), name)
			fh = open(name, "w")
			fh.write(cont)
			fh.close()

	def convert(self, inputFile):
		"""Given a TBW/TBN file, build a list of VDIF files that are ready to go 
		into DiFX.  This function splits the input DP file into Nstands VDIF files 
		that are named in the format speficied by the configuration writer, i.e.:
		  <mode>-<stand>.vdif
		"""

		if self.mode == 'TBW':
			print "Converting TBW file '%s' to %i VDIF files" % (inputFile, len(self.stands))
			self.__convertTBW(inputFile)
		else:
			print "Converting TBN file '%s' to %i VDIF files" % (inputFile, len(self.stands))
			self.__convertTBN(inputFile)
		
	def __convertTBW(self, inputFile):
		"""Private function call by convert for working with TBW data."""

		fh = open(inputFile, 'rb')
	
		# Make sure that the data is TBW and determine the data length
		test = tbw.readFrame(fh)
		if not test.header.isTBW():
			raise errors.notTBWError()
		fh.seek(0)

		# Create the output files
		ofh = []
		outputFiles = []
		for stand in self.stands:
			outputFiles.append('%s-%i.vdif' % (self.mode, stand))
			ofh.append(open(outputFiles[-1], 'wb'))

		# Read in the data and add it to the VDIF files created above
		for i in range(300000):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = tbw.readFrame(fh, Verbose=False)
			except errors.eofError:
				break
			except errors.syncError:
				syncCount = syncCount + 1
				continue
			except errors.numpyError:
				break

			stand = cFrame.header.parseID()
			if cFrame.header.frameCount % 10000 == 0:
				print("%2i  %14i  %6.3f  %5i  %5i" % (stand, cFrame.data.timeTag, cFrame.getTime(), cFrame.header.frameCount, cFrame.header.secondsCount))
			
			frame1 = Frame(stand=2*(stand-1)+1, time=cFrame.getTime(), bits=cFrame.getDataBits(), 
							sampleRate = dp_common.fS, data=cFrame.data.xy[0,:])
			frame2 = Frame(stand=2*(stand-1)+2, time=cFrame.getTime(), bits=cFrame.getDataBits(), 
							sampleRate = dp_common.fS, data=cFrame.data.xy[1,:])
			
			frame1.writeRawFrame(ofh[2*(stand-1)+0])
			frame2.writeRawFrame(ofh[2*(stand-1)+1])
		fh.close()

		# Close all of the constituent VDIF files.
		for fh in ofh:
			fh.close()

	def __convertTBN(self, inputFile):
		"""Private function call by convert for working with TBN data."""
	
		fh = open(inputFile, "rb")
		
		test = tbn.readFrame(fh)
		if not test.header.isTBN():
			raise errors.notTBNError()
		fh.seek(0)

		nFpO = tbn.getFramesPerObs(fh)
		nFpO = nFpO[0] + nFpO[1]

		sampleRate = tbn.getSampleRate(fh, nFrames=nFpO)
		nSamples = 340000

		blocks = []
		count = {}
		syncCount = 0
		masterCount = 0
		for i in range(nSamples):
			try:
				cFrame = tbn.readBlock(fh, nFrames=nFpO, SampleRate=sampleRate)
			except errors.eofError:
				break
			except errors.syncError:
				syncCount = syncCount + 1
				continue
			except errors.numpyError:
				break
			blocks.append(cFrame)

			for frame in blocks[-1].x:
				if frame is None:
					print "sync error"
					continue
				stand, pol = frame.parseID()
				if stand not in count.keys():
					count[stand] = 0
				count[stand] = count[stand] + 1

	def make(self, inputFile):
		"""Convience wrapper for write and convert so that one call will create
		everything needed to run DiFX on LWA data."""

		self.write()
		self.convert(inputFile)