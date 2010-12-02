#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Application to calculate real-time ephemeris for a LWDA site."""

import curses
import math
import time
import traceback
import sys, os, string
import optparse

from lsl import astro
from lsl import transform
from lsl.common import stations as lwa_common

__revision__  = "$Revision: 94 $"
__version__   = "dev"
__author__    = "P.S.Ray"


def restorescreen():
    curses.nocbreak()
    curses.echo()
    curses.endwin()
    

display = string.Template("""==========================================================================
UTC : $utc          MJD(UTC) $mjdutc
Site: $loc                              LST: $lst
Long: $long                    Lat: $lat
==========================================================================
----------- Sun --------------          ------------ Jupiter -------------
RA(Cur)          $sun_ra        RA(Cur)               $jup_ra
Dec(Cur)        $sun_dec         Dec(Cur)             $jup_dec
Alt                   $sun_alt deg        Alt                  $jup_alt deg
Az                    $sun_az deg        Az                   $jup_az deg
---- Galactic Center ---------          ---------- Tau A (Crab) ----------
RA(Cur)          $gal_ra        RA(Cur)               $tau_ra
Dec(Cur)        $gal_dec         Dec(Cur)             $tau_dec
Alt                   $gal_alt deg        Alt                  $tau_alt deg
Az                    $gal_az deg        Az                   $tau_az deg
---------- Cas A -------------          ------------ Cyg A ---------------
RA(Cur)          $cas_ra        RA(Cur)               $cyg_ra
Dec(Cur)        $cas_dec         Dec(Cur)             $cyg_dec
Alt                   $cas_alt deg        Alt                   $cyg_alt deg
Az                    $cas_az deg        Az                    $cyg_az deg
==========================================================================""")


if __name__ == '__main__':
	usage = "astrostatus [options]"

	# check command line
	parser = optparse.OptionParser(usage = usage)
	parser.add_option("-s", "--site", action = "store", type = "string",
					dest = "site", default = 'LWDA', help = "site name (default LWDA)")
	(opts, args) = parser.parse_args()
			
	# setup transform objects
	lwa1 = lwa_common.lwa1()
	site = transform.GeographicalPosition((lwa1.long*180.0/math.pi, lwa1.lat*180.0/math.pi), name=lwa1.name)
	
	sun_pos = transform.PlanetaryPosition('Sun')
	sun_dir = transform.PointingDirection(sun_pos, site)
	
	jup_pos = transform.PlanetaryPosition('Jupiter')
	jup_dir = transform.PointingDirection(jup_pos, site)
	
	gal_pos = transform.CelestialPosition((astro.hms(17, 42, 0.1),
	astro.dms(True, 28, 55, 8)), name = 'SgrA')
	gal_dir = transform.PointingDirection(gal_pos, site)
	
	tau_pos = transform.CelestialPosition((astro.hms(5, 34, 0.5),
	astro.dms(False, 22, 0, 52)), name = 'TauA')
	tau_dir = transform.PointingDirection(tau_pos, site)
	
	cyg_pos = transform.CelestialPosition((astro.hms(19, 59, 0.8),
	astro.dms(False, 40, 44, 2)), name = 'CygA')
	cyg_dir = transform.PointingDirection(cyg_pos, site)
	
	cas_pos = transform.CelestialPosition((astro.hms(23, 23, 0.7),
	astro.dms(False, 58, 49, 16)), name = 'CasA')
	cas_dir = transform.PointingDirection(cas_pos, site)
    
	try:
		# setup display screen
		screen = curses.initscr()
		curses.noecho()
		curses.cbreak()
		screen.nodelay(1)

		strdict = {}
		while (True):
			screen.clear()

			# Build string to print
			# update current time and local sidereal time
			currentTime = transform.Time.from_system() 
						
			strdict['utc'] = "%-25s" % currentTime.utc_str
			strdict['mjdutc'] = "%.5f" % currentTime.utc_mjd
			strdict['loc'] = site.name
			(lng, lat) = site.geo.format()
			strdict['lat'] = lat
			strdict['long'] = "%-15s" % lng
			strdict['lst'] = astro.deg_to_dms(site.sidereal(currentTime))
			
			# calculate Sun position 
			sun_hrz = sun_dir.hrz(currentTime)
			(ra, dec) = sun_pos.apparent_equ(currentTime).format()
			strdict['sun_ra'] = "%-15s" % ra
			strdict['sun_dec'] = "%-15s" % dec
			strdict['sun_alt'] = "%6.2f" % sun_hrz.alt
			strdict['sun_az'] = "%6.2f" % sun_hrz.az

			# calculate Jupiter posn
			jup_hrz = jup_dir.hrz(currentTime)
			(ra, dec) = jup_pos.apparent_equ(currentTime).format()
			strdict['jup_ra'] = ra
			strdict['jup_dec'] = dec
			strdict['jup_alt'] = "%6.2f" % jup_hrz.alt
			strdict['jup_az'] = "%6.2f" % jup_hrz.az

			# Galactic Center
			gal_hrz = gal_dir.hrz(currentTime)
			(ra, dec) = gal_pos.apparent_equ(currentTime).format()
			strdict['gal_ra'] = "%-15s" % ra
			strdict['gal_dec'] = "%-15s" % dec
			strdict['gal_alt'] = "%6.2f" % gal_hrz.alt
			strdict['gal_az'] = "%6.2f" % gal_hrz.az

			# calculate Tau A (Crab) position
			tau_hrz = tau_dir.hrz(currentTime)
			(ra, dec) = tau_pos.apparent_equ(currentTime).format()
			strdict['tau_ra'] = ra
			strdict['tau_dec'] = dec
			strdict['tau_alt'] = "%6.2f" % tau_hrz.alt
			strdict['tau_az'] = "%6.2f" % tau_hrz.az

			# calculate Cygnus A position
			cyg_hrz = cyg_dir.hrz(currentTime)
			(ra, dec) = cyg_pos.apparent_equ(currentTime).format()
			strdict['cyg_ra'] = ra
			strdict['cyg_dec'] = dec
			strdict['cyg_alt'] = "%6.2f" % cyg_hrz.alt
			strdict['cyg_az'] = "%6.2f" % cyg_hrz.az

			# calculate Cassiopeia A position
			cas_hrz = cas_dir.hrz(currentTime)
			(ra, dec) = cas_pos.apparent_equ(currentTime).format()
			strdict['cas_ra'] = "%-15s" % ra
			strdict['cas_dec'] = "%-15s" % dec
			strdict['cas_alt'] = "%6.2f" % cas_hrz.alt
			strdict['cas_az'] = "%6.2f" % cas_hrz.az
		
			# Refresh screen and wait 1 second
			screen.addstr(0,0,display.safe_substitute(strdict))
			screen.refresh()
			time.sleep(1)
				
			# Check for keypress and exit if Q
			c = screen.getch()
			if (c > 0):
				if chr(c) == 'q': 
					break
				if chr(c) == 'Q': 
					break
			
		restorescreen()
		
	except KeyboardInterrupt:
		restorescreen()
