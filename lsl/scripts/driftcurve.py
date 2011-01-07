#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Predict driftcurve for a given site using a given antenna model."""

import sys
import optparse
import logging
import math
import os

import numpy
import pylab

from lsl import skymap, astro
from lsl.common import stations
from lsl.sim import nec_util
        

__revision__ = "$Revision: 93 $"
__version__  = "0.1"
__author__    = "D.L.Wood"
__maintainer__ = "Jayce Dowell"


if __name__ == '__main__':
	# setup logger

	logging.basicConfig()
	log = logging.getLogger()
	log.setLevel(logging.INFO)
	
	# parse command line
	parser = optparse.OptionParser(usage="%prog [options] FREQ NECFILENAME SKYMAPFILENAME", \
		description="Generate a drift curve for a dipole at SITE observing at a given FREQ (MHz).  SITE must be one of the sites known by the station module in lwda_util.")
	parser.add_option("-v", "--verbose", action = "store_true", dest = "verbose",
		default = False, help = "enable debug messages")
	parser.add_option("-x", "--doplot", action = "store_true", dest = "doplot",
		default = False, help = "Make an X-windows plot")
	parser.add_option("-s", "--site", action = "store", dest = "site", type = "str",
		default = "LWA", help = "site name (default LWA)")
	parser.add_option("-p", "--polarization", action = "store", dest = "polarization",
		default = "NS", help = "antenna polarization orientation (NS or EW)")
	
	(opts, args) = parser.parse_args()    

	if len(args) != 3:
		parser.error("wrong number of arguments")
		
	if opts.verbose:
		log.setLevel(logging.DEBUG)
		
	if (opts.polarization != 'EW') and (opts.polarization != 'NS'):
		parser.error("only 'EW' and 'NS' allowed for polarization")
		
	# Extract the command line arguments
	freq = float(args[0])
	if freq < 1.0 or freq > 1000.0:
		log.error("Please specify a frequency in MHz!")
		sys.exit(1)
	necname = args[1]
	skymapname = args[2]
	
	# get site geo position
	try:
		sta = stations.lwa1()
	except:
		log.error("cannot get info for site %s", site)
		sys.exit(1)
	log.info("Site %s, Long %f", sta.name, sta.long*180.0/math.pi)
	
	# read skymap and scale to target frequency
	log.info("Reading skymap file")
	smap = skymap.SkyMap(skymapname)
	log.info("Read skymap of %d x %d pixels, min=%f, max=%f",smap.numPixelsX,smap.numPixelsY,smap._power.min(),smap._power.max())
	
	log.info("Reading antenna info")    
	# get user-supplied pattern
	antPatAnt = nec_util.NECPattern(necname, freq).antenna_pat_dB
	#pylab.figure(3)
	log.info("Read pattern.  Min %f Max %f", antPatAnt.min(), antPatAnt.max())
	#pylab.imshow(antPatAnt.transpose(), origin='lower', vmin=-10.0, cmap=pylab.cm.hot)
	antPatAnt = numpy.power(10.0, antPatAnt / 10.0)
	
	# calculate times in both site LST and UTC
	t0 = astro.get_julian_from_sys()
	lst = astro.get_local_sidereal_time(sta.long*180.0/math.pi, t0) / 24.0
	t0 -= lst*(23.933/24.0) # Compensate for shorter sidereal days
	times = numpy.arange(0, 1, 0.2/24) + t0
	
	lstList = []
	powListAnt = [] 
	
	for t in times:
		# project skymap to site location and observation time
		pmap = skymap.ProjectedSkyMap(smap, sta.lat*180.0/math.pi, sta.long*180.0/math.pi, t)
		lst = astro.get_local_sidereal_time(sta.long*180.0/math.pi, t)
		iaz = pmap.visibleAz.astype(numpy.int_)
		ialt = pmap.visibleAlt.astype(numpy.int_)
		lstList.append(lst)
		
		#log.info("iaz : %f - %f",iaz.min(),iaz.max())
		#log.info("ialt : %f - %f",ialt.min(),ialt.max())
		if opts.polarization == 'EW':
			iaz += 90.0
			iaz = (iaz >= 360).choose(iaz, iaz - 360)
		
		cdec = numpy.cos(pmap.visibleDec * smap.degToRad)
				
		# convolution of user antenna pattern with visible skymap
		gain = antPatAnt[iaz, ialt]
		powerAnt = (pmap.visiblePower * gain * cdec).sum() / (gain * cdec).sum()
		powListAnt.append(powerAnt)
		
		log.debug("LST=%f, power_ant=%f", lst, powerAnt) 
			
	# plot results
	if opts.doplot:
		pylab.figure(1)
		pylab.title("Driftcurve: Site %s Pol %s Ant %s Freq %0.2f MHz" % \
			(opts.site, opts.polarization, os.path.basename(necname), freq))
		pylab.plot(lstList, powListAnt, "ro",label="NEC Antenna Pattern")
		pylab.xlabel("LST")
		pylab.ylabel("Temp (K)")
		pylab.grid(1)
		pylab.show()
	
	mf = file("model_%s_%s_%s_%0.2f.txt" % (opts.site, os.path.basename(necname), opts.polarization, freq),"w")
	for lst,pow in zip(lstList,powListAnt):
		mf.write("%f  %f\n" % (lst,pow))
	mf.close()
	
	sys.exit(0)
