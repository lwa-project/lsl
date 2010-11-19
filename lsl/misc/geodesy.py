#/usr/bin/env python

import os
import io
import sys
import gzip
import ephem
import urllib

import lsl.astro as astro
from lsl.common.paths import data as dataPath

class EOP(object):
	def __init__(self, mjd=0.0, x=0.0, y=0.0, utDiff=0.0, type='final'):
		self.mjd = mjd
		self.x = x
		self.y = y
		self.utDiff = utDiff
		self.type = type

	def fromMAIA(self, line):
		self.mjd = float(line[7:15])
		self.x = float(line[18:27])
		self.y = float(line[38:46])
		self.utDiff = float(line[58:68])
		if line[57] == 'I':
			self.type = 'final'
		else:
			self.type = 'prediction'

	def __str__(self):
		return "%.1f: x=%.6f y=%.6f UT1-UTC=%.6f (%s)" % (self.mjd, self.x, self.y, self.utDiff, self.type)
		
heopFH = gzip.GzipFile(os.path.join(dataPath, 'astro', 'eop-historical.dat.gz'), 'rb')
lines = heopFH.readlines()

heops = []
for line in lines:
	line = heopFH.read(185)
	newEOP = EOP()
        newEOP.fromMAIA(line)
	if newEOP.type == 'final':
	        heops.append(newEOP)
        	print "%s - h" % str(heops[-1])
heopFH.close()

sys.exit()

eopFH = urllib.urlopen('http://maia.usno.navy.mil/ser7/finals2000A.daily')
lines = eopFH.readlines()
eopFH.close()

eops = []
for line in lines:
	newEOP = EOP()
	newEOP.fromMAIA(line) 
	eops.append(newEOP)
	print eops[-1]


