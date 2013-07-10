# -*- coding: utf-8 -*-

"""
LWA astronomical source catalogs.
"""


import os
import math
import abc
import collections
import logging

from lsl import astro
from lsl import transform
from lsl.common.paths import data as dataPath

__version__   = '0.1'
__revision__ = '$Rev$'
__author__    = 'D. L .Wood'
__maintainer__ = 'Jayce Dowell'


class CatalogEntry(object):
	"""
	Represents one source entry in a catalogue.
	
	Contains members:
	* name        - The source name.
	* position    - The source equatorial J2000 position as object
		of type transform.CelestialPosition.
	* alias_list  - A list of strings providing alternate names for
		the source.
	"""
	
	# limit the class attributes
	__slots__ = ('name', 'position', 'alias_list')
	
	def __init__(self, name, position):
		"""
		Create a catalog entry.
		"""
		
		self.name = name
		self.position = position
		self.alias_list = []
		
	def __repr__(self):
		"""
		Low level string representation.
		"""
		
		return "%s.%s(%s,%s)" % (type(self).__module__, type(self).__name__, repr(self.name), repr(self.position))


class Catalog(collections.Mapping):
	"""
	Class representing astronomical source catalog information.
	This is an abstract class; derived classes must provide a
	parse_file() method which populates the catalog object with
	information from file or other source.
	
	Catalog instances support the read-only collections.Mapping
	interface.  That is, they support the read-only methods of
	the dict built-in type.
	"""
	
	__metaclass__ = abc.ABCMeta
	
	# a logging instance for catalogs
	log = logging.getLogger('catalog')
	
	def __init__(self, name):
		"""
		Create a source catalog.
		"""
		
		# initialize catalog data structures
		self.name = name
		self.source_map = {}
		self.alias_map = {}
		
		# parse_file() is an abstract method which must be defined in
		# a concrete implementation for a particular catalog
		self.log.debug("parsing data file for catalog %s", name)    
		self.parse_file()
		
	@abc.abstractmethod
	def parse_file(self):
		"""
		Read catalog information from file into internal data
		structures.
		"""
		
		pass
		
	@staticmethod
	def get_directory():
		"""
		Returns the path to the catalog data file directory.
		"""
			
		return os.path.join(dataPath, 'catalog')
		
	def __repr__(self):   
		"""
		Low level string representation.
		"""
			
		return "%s.Catalog(%s)" % (type(self).__module__, repr(self.name))
		
	def __len__(self):
		"""
		Return the number of sources in the catalog.
		"""
		
		return len(self.source_map)
		
	def __getitem__(self, key):
		"""
		Access source by subscript name.  Raises KeyError if the source
		is not in the catalog.
		"""
		
		entry = self.lookup(key)
		if entry is None:
			raise KeyError("name %s not in catalog" % repr(key))
		return entry
		
	def __iter__(self):
		"""
		Return an iterator to the primary names in the catalog.
		"""
		
		return iter(self.source_map.keys())
		
	def lookup(self, name):
		"""
		Lookup a source in the catalog.
		
		Param: name - The primary name or alias of the source.
		
		Returns: An object of type CatalogEntry giving the source information,
				or None if the name is not found in the catalog.
		"""
		
		try:
			entry = self.source_map[name]
		except KeyError:
			try:
				entry = self.alias_map[name]
			except KeyError:
				entry = None
				
		return entry


class LWA_Catalog(Catalog):
	"""
	Specific definition for LWA observation source catalogue data file.
	"""         
	
	def __init__(self):
		"""
		Create a LWA catalog instance.
		"""
		
		Catalog.__init__(self, 'LWA')
		
	def parse_file(self):
		"""
		Read a source catalog data file.
		"""
		
		# open data file 
		fileName = os.path.join(self.get_directory(), 'lwa_catalog.dat')
		catFile = open(fileName, 'r')
		
		# read source info
		lineNum = 0
		for line in catFile:
			lineNum += 1
			if line.startswith('#') or line.isspace():
				continue
				
			try:
				name = line[0:8]
				raHours = int(line[9:11], 10)
				raMinutes = int(line[12:14], 10)
				raSeconds = float(line[15:20])
				decSign = line[21]
				decDegrees = int(line[22:24], 10)
				decMinutes = int(line[25:27], 10)
				decSeconds = float(line[28:32])
			except ValueError as err:
				raise RuntimeError("file %s, line %d incorectly formated: \'%s\' : %s]" % (fileName, lineNum, line, err))
				
			name = name.rstrip()           
			
			ra = astro.hms(raHours, raMinutes, raSeconds)
			if decSign == '-':
				sign = True
			else:
				sign = False
			dec = astro.dms(sign, decDegrees, decMinutes, decSeconds)
			
			entry = CatalogEntry(name, transform.CelestialPosition((ra, dec), name=name))
			self.source_map[name] = entry
			
			aliasList = line[34:].split()
			if len(aliasList) > 0:
				for alias in aliasList:
					entry.alias_list.append(alias)
					self.alias_map[alias] = entry
					
		catFile.close()


class PSR_Catalog(Catalog):
	"""
	Specific definition for ATNF Pulsar (PSRCAT) catalog.
	Data file is psrcat.db which can be retreived from:
	<http://www.atnf.csiro.au/research/pulsar/psrcat/download.html>
	""" 
	
	def __init__(self):
		"""
		Create an instance of the PSR catalog.
		"""
		
		Catalog.__init__(self, 'PSR')
		
	def parse_file(self):
		"""
		Read a source catalog data file.
		"""
		
		debug = False
		# open data file 
		fileName = os.path.join(self.get_directory(), 'psrcat.db')
		catFile = open(fileName, 'r')
		
		# read source info
		psrb = None
		psrj = None
		ra = None
		dec = None
		bad = False
		for line in catFile:
			if line.startswith('PSRB'):
				psrb = line.split()[1]
			if line.startswith('PSRJ'):
				psrj = line.split()[1]
			if line.startswith('RAJ'):
				rastr = line.split()[1]
				sp = rastr.split(':')
				if len(sp) == 3:
					(raHours, raMinutes, raSeconds) = sp
				elif len(sp) == 2:
					(raHours, raMinutes) = sp
					raSeconds = 0.0
				else:
					if debug: print('Bad format for RAJ line : '+line)
					bad = True
				raHours = int(raHours)
				raMinutes = int(raMinutes)
				raSeconds = float(raSeconds)
				try:
					ra = astro.hms(raHours, raMinutes, raSeconds)
				except:
					if debug: print('PSRCAT: Bad RA for ',psrj," : ",rastr)
					bad = True
			if line.startswith('DECJ'):
				decstr = line.split()[1]
				sp = decstr.split(':')
				if len(sp) == 3:
					(decDegrees, decMinutes, decSeconds) = sp
				elif len(sp) == 2:
					(decDegrees, decMinutes) = sp
					decSeconds = 0.0
				else:
					if debug: print('PSRCAT: Bad format for DECJ line : '+line)
					bad = True
					continue
				if decDegrees.startswith('-'):
					decDegrees = decDegrees[1:]
					sign = True
				else:
					sign = False
				decDegrees = int(decDegrees)
				decMinutes = int(decMinutes)
				decSeconds = float(decSeconds)
				try:
					dec = astro.dms(sign, decDegrees, decMinutes, decSeconds)
				except:
					if debug: print('PSRCAT: Bad DEC for ',psrj," : ",decstr)
					bad = True
					
			if line.startswith('@-'):
				# New source, save current source info
				if psrb is not None:
					name = psrb
					alias = psrj
				else:
					name = psrj
					alias = None
					
				if ra is None or dec is None:
					# These pulsars don't have RAJ, DECJ
					# I think they may all have ecliptic positions
					# which should be converted to ra,dec but I'm
					# going to ignore them for now. -- paulr
					#print "PSRCAT: No position for pulsar ",name
					bad = True
					
				# Add source to list if good.
				if not bad:                    
					sourcePos = astro.equ_posn(ra, dec) 
					entry = CatalogEntry(name, transform.CelestialPosition(sourcePos, name=name))
					self.source_map[name] = entry
					if debug: print('Added ',name)
					
					if alias is not None:
						alias = alias.rstrip()
						self.alias_map[alias] = entry
						entry.alias_list.append(alias)
						if debug: print('Alias : ',alias.rstrip())
						
				# Clear out vars for next source
				psrb = None
				psrj = None
				ra = None
				dec = None
				bad = False
				
		catFile.close()


class PKS_Catalog(Catalog):
	"""
	Specific definition for PKS source catalog.
	"""
	
	def __init__(self):
		"""
		Create a PKS catalog instance.
		"""
		
		Catalog.__init__(self, 'PKS')
		
	def parse_file(self):
		"""
		Read a source catalog data file.
		"""
		
		# open data file 
		fileName = os.path.join(self.get_directory(), 'pkscat.txt')
		catFile = open(fileName, 'r')
		
		# read source info
		lineNum = 1
		line = catFile.readline()
		while line:
			if line == '\n':
				line = catFile.readline()
				lineNum += 1
				continue
				
			try:
				name = line[0:8]
				alias = line[12:19]
				raHours = int(line[22:24])
				raMinutes = int(line[25:27])
				raSeconds = float(line[28:32])
				decSign = line[33]
				decDegrees = int(line[34:36])
				decMinutes = int(line[37:39])
				decSeconds = int(line[40:42])
			except ValueError:
				raise RuntimeError("file %s, line %d incorectly formated [%s]" % (fileName, lineNum, line))
				
			ra = astro.hms(raHours, raMinutes, raSeconds)
			if decSign == '-':
				sign = True
			else:
				sign = False
			dec = astro.dms(sign, decDegrees, decMinutes, decSeconds)
			sourcePos = astro.equ_posn(ra, dec)
			
			# precess coordinates from B1950
			entry = CatalogEntry(name, transform.CelestialPosition(sourcePos, epoch='B1950', name=name))
			self.source_map[name] = entry
			
			if len(alias.strip()):
				alias = alias.rstrip()
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			line = catFile.readline()
			lineNum += 1
			
		catFile.close()


class PKS90_Catalog(Catalog):
	"""
	Specific definition for PKS90 source catalogue data file.
	"""
	
	def __init__(self):
		"""
		Create a PKS90 catalog instance.
		"""
		
		Catalog.__init__(self, 'PKS90')
		
	def parse_file(self):
		"""
		Read a source catalog data file.
		"""
		
		# open data file 
		fileName = os.path.join(self.get_directory(), 'PKS90.txt')
		catFile = open(fileName, 'r')
		
		# read source info
		lineNum = 1
		line = catFile.readline()
		while line:
			if line == '\n':
				line = catFile.readline()
				lineNum += 1
				continue
				
			try:
				name = line[0:10]
				alias0 = line[139:148]
				alias1 = line[168:175]
				raHours = int(line[10:12])
				raMinutes = int(line[13:15])
				raSeconds = float(line[16:20])
				decSign = line[23]
				decDegrees = int(line[24:26])
				decMinutes = int(line[27:29])
				decSeconds = float(line[30:34])
			except ValueError:
				raise RuntimeError("file %s, line %d incorectly formated [%s]" % (fileName, lineNum, line))
				
			ra = astro.hms(raHours, raMinutes, raSeconds)
			if decSign == '-':
				sign = True
			else:
				sign = False
			dec = astro.dms(sign, decDegrees, decMinutes, decSeconds)
			sourcePos = astro.equ_posn(ra, dec)
			
			entry = CatalogEntry(name, transform.CelestialPosition(sourcePos, name=name))
			self.source_map[name] = entry
			
			if len(alias0.strip()): 
				alias = alias0.rstrip()
				entry.alias_list.append(alias)
				self.alias_map[alias] = entry
			if len(alias1.strip()):
				alias = alias1.rstrip()
				entry.alias_list.append(alias)
				self.alias_map[alias] = entry
				
			line = catFile.readline()
			lineNum += 1
			
		catFile.close()


class C3C_Catalog(Catalog):
	"""
	Specific definition for Cambridge 3C source catalogue data file.
	""" 
	
	def __init__(self):
		"""
		Create a 3C catalog instance.
		"""
		
		Catalog.__init__(self, '3C')
		
	def parse_file(self):
		"""
		Read a source catalog data file.
		"""
		
		# open data file 
		fileName = os.path.join(self.get_directory(), '3c.dat')
		catFile = open(fileName, 'r')
		
		# read source info
		lineNum = 1
		line = catFile.readline()
		while line:
			try:
				name = line[0:3]
				raHours = int(line[12:14])
				raMinutes = int(line[15:17])
				raSeconds = float(line[18:22])
				decSign = line[27]
				decDegrees = int(line[28:30])
				decMinutes = float(line[31:35])
			except ValueError:
				raise RuntimeError("file %s, line %d incorectly formated [%s]" % (fileName, lineNum, line))
				
			name = ('3C' + name.strip())         
			
			ra = astro.hms(raHours, raMinutes, raSeconds)
			(decSeconds, decMinutes) = math.modf(decMinutes)
			if decSign == '-':
				sign = True
			else:
				sign = False
			dec = astro.dms(sign, decDegrees, int(decMinutes), (60 * decSeconds))
			sourcePos = astro.equ_posn(ra, dec)
			
			# precess coordinates from B1950
			entry = CatalogEntry(name, transform.CelestialPosition(sourcePos, epoch='B1950', name=name))
			self.source_map[name] = entry        
			
			line = catFile.readline()
			lineNum += 1      
			
		catFile.close()


class C4C_Catalog(Catalog):
	"""
	Specific definition for Cambridge 4C source catalogue data file.
	""" 
	
	def __init__(self):
		"""
		Create a 4C catalog instance.
		"""
		
		Catalog.__init__(self, '4C')
		
	def parse_file(self):
		"""
		Read a source catalog data file.
		"""
		
		# open data file 
		fileName = os.path.join(self.get_directory(), '4c.dat')
		catFile = open(fileName, 'r')
		
		# read source info
		lineNum = 1
		line = catFile.readline()
		while line:
			try:
				name = line[0:8]
				raHours = int(line[10:12])
				raMinutes = int(line[13:15])
				raSeconds = float(line[16:20])
				decSign = line[22]
				decDegrees = int(line[23:25])
				decMinutes = float(line[26:30])
				alias = line[64:-1]
			except ValueError:
				raise RuntimeError("file %s, line %d incorectly formated [%s]" % (fileName, lineNum, line))
				
			name = name.strip()        
			name = ('4C' + name)
			
			alias = alias.strip()
			alias = alias.rstrip('*')         
			
			ra = astro.hms(raHours, raMinutes, raSeconds)
			(decSeconds, decMinutes) = math.modf(decMinutes)
			if decSign == '-':
				sign = True
			else:
				sign = False
			dec = astro.dms(sign, decDegrees, int(decMinutes), (60 * decSeconds))
			sourcePos = astro.equ_posn(ra, dec)
			
			# precess coordinates from B1950
			entry = CatalogEntry(name, transform.CelestialPosition(sourcePos, epoch='B1950', name=name))
			self.source_map[name] = entry
			
			if len(alias.strip()):
				alias = alias.rstrip()
				entry.alias_list.append(alias)
				self.alias_map[alias] = entry       
				
			line = catFile.readline()
			lineNum += 1      
			
		catFile.close()


class F1FGL_Catalog(Catalog):
	"""
	Specific definition for Fermi LAT first point source catalog.
	"""
	
	def __init__(self):
		"""
		Create a 1FGL catalog instance.
		"""
		
		Catalog.__init__(self, '1FGL')
		
	def parse_file(self):
		"""
		Read a source catalogue data file.
		"""
		
		import pyfits
		
		# open data file
		fileName = os.path.join(self.get_directory(), 'gll_psc_v02.fit')
		catFile = pyfits.open(fileName)
		
		# read source info
		sourceTable = catFile['LAT_POINT_SOURCE_CATALOG'].data
		
		for row in sourceTable:
			name = str(row.field('Source_Name')).replace(' ', '_')
			ra = float(row.field('RA'))
			dec = float(row.field('DEC'))
			entry = CatalogEntry(name, 
			transform.CelestialPosition((ra, dec), name=name))
			self.source_map[name] = entry
			
			alias = str(row.field('0FGL_Name'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC_GAM1'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC_GAM2'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC1'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC2'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
		catFile.close()


class F2FGL_Catalog(Catalog):
	"""
	Specific definition for Fermi LAT 2-year point source catalog.
	"""
	
	def __init__(self):
		"""
		Create a 2FGL catalog instance.
		"""
		
		Catalog.__init__(self, '2FGL')
		
	def parse_file(self):
		"""
		Read a source catalogue data file.
		"""
		
		import pyfits
		
		# open data file
		fileName = os.path.join(self.get_directory(), 'gll_psc_v08.fit')
		catFile = pyfits.open(fileName)
		
		# read source info
		sourceTable = catFile['LAT_POINT_SOURCE_CATALOG'].data
		
		for row in sourceTable:
			name = str(row.field('Source_Name')).replace(' ', '_')
			ra = float(row.field('RAJ2000'))
			dec = float(row.field('DEJ2000'))
			entry = CatalogEntry(name, 
			transform.CelestialPosition((ra, dec), name=name))
			self.source_map[name] = entry
			
			alias = str(row.field('0FGL_Name'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('1FGL_Name'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC_GAM1'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC_GAM2'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC1'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
			alias = str(row.field('ASSOC2'))
			if len(alias):
				alias = alias.replace(' ', '_')
				self.alias_map[alias] = entry
				entry.alias_list.append(alias)
				
		catFile.close()


class CatalogFactory(object):
	"""
	Get catalog objects by name.  Caches the catalog data so that
	the data file is parsed only once per session.
	"""
	
	# a mapping of catalog names to classes
	catalog_class_map = \
					{
						'LWA'   : LWA_Catalog,
						'PSR'   : PSR_Catalog,
						'PKS'   : PKS_Catalog,
						'PKS90' : PKS90_Catalog,
						'3C'    : C3C_Catalog,
						'4C'    : C4C_Catalog,
						'1FGL'  : F1FGL_Catalog,
						'2FGL'  : F2FGL_Catalog,
					}
					
	# a mapping of catalog names to instances
	catalog_instance_map = {}
	
	@classmethod
	def get_catalog(klass, name):
		"""
		Returns a Catalog object representing the catalog
		given by name.
		"""
		
		# check parameters
		if name not in list(klass.catalog_class_map.keys()):
			raise ValueError("unknown catalog \'%s\'" % name)
			
		# try to find an already created instance
		# if not found, create a new instance and cache
		
		try:
			catalog = klass.catalog_instance_map[name]
		except KeyError:
			catalogClass = klass.catalog_class_map[name]
			catalog = catalogClass()
			klass.catalog_instance_map[name] = catalog
			
		return catalog
		
	@classmethod
	def get_names(klass):
		"""
		Return a list of known catalog names.
		"""
		
		return list(klass.catalog_class_map.keys())
