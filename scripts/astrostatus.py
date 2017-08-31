#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Application to calculate real-time ephemeris for a LWA site."""

import math
import time
import curses
import optparse

from lsl import astro
from lsl import transform
from lsl.common import stations

__revision__  = "$Revision: 94 $"
__version__   = "dev"
__author__    = "P.S.Ray"


def restorescreen():
    curses.nocbreak()
    curses.echo()
    curses.endwin()


if __name__ == '__main__':
	usage = "astrostatus.py [options]"

	# check command line
	parser = optparse.OptionParser(usage = usage)
	parser.add_option("-s", "--site", action = "store", type = "string",
					dest = "site", default = 'LWA-1', help = "site name (default LWA-1)")
	(opts, args) = parser.parse_args()
	
	# setup transform objects
	opts['site'] = opts['site'].lower().replace('-', '')
	if opts['site'] == 'lwa1':
		station = stations.lwa1
	elif opts['site'] == 'lwasv':
		station = stations.lwasv
	else:
		raise RuntimeError("Unknown site name: %s" % opts['site'])
	site = transform.GeographicalPosition((station.long*180.0/math.pi, station.lat*180.0/math.pi), name=station.name)
	
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
		curses.start_color()
		curses.noecho()
		curses.cbreak()
		screen.nodelay(1)
		screen.border()
		
		if curses.has_colors():
			curses. use_default_colors()
			
			curses.init_pair(1, -1, -1)
			curses.init_pair(2, curses.COLOR_RED, -1)
			curses.init_pair(3, curses.COLOR_GREEN, -1)
			curses.init_pair(4, curses.COLOR_BLUE, -1)
			curses.init_pair(5, curses.COLOR_YELLOW, -1)
			curses.init_pair(6, curses.COLOR_CYAN, -1)
			
			std = curses.color_pair(1)
			dwn = curses.color_pair(2)
			up  = curses.color_pair(3)
			st  = curses.color_pair(4) | curses.A_BOLD
			rs  = curses.color_pair(5) | curses.A_BOLD
			inf = curses.color_pair(6)
		else:
			std = curses.A_NORMAL
			dwn = curses.A_NORMAL
			up  = curses.A_NORMAL
			st  = curses.A_NORMAL
			rs  = curses.A_NORMAL
			inf = curses.A_NORMAL
			
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
			strdict['lat'] = "%-15s" % lat
			strdict['long'] = "%-15s" % lng
			strdict['lst'] = "%-15s" % astro.deg_to_dms(site.sidereal(currentTime))
			
			# calculate Sun position 
			sun_hrz = sun_dir.hrz(currentTime)
			(ra, dec) = sun_pos.apparent_equ(currentTime).format()
			strdict['sun_ra'] = "%-15s" % ra
			strdict['sun_dec'] = "%-15s" % dec
			strdict['sun_alt'] = "%6.2f" % sun_hrz.alt
			strdict['sun_az'] = "%6.2f" % sun_hrz.az
			try:
				if sun_hrz.alt > last_sun_hrz.alt:
					sun_rising = True
				else:
					sun_rising = False
			except NameError:
				pass
			last_sun_hrz = sun_hrz
			
			# calculate Jupiter posn
			jup_hrz = jup_dir.hrz(currentTime)
			(ra, dec) = jup_pos.apparent_equ(currentTime).format()
			strdict['jup_ra'] = "%-15s" % ra
			strdict['jup_dec'] = "%-15s" % dec
			strdict['jup_alt'] = "%6.2f" % jup_hrz.alt
			strdict['jup_az'] = "%6.2f" % jup_hrz.az
			try:
				if jup_hrz.alt > last_jup_hrz.alt:
					jup_rising = True
				else:
					jup_rising = False
			except NameError:
				pass
			last_jup_hrz = jup_hrz
			
			# Galactic Center
			gal_hrz = gal_dir.hrz(currentTime)
			(ra, dec) = gal_pos.apparent_equ(currentTime).format()
			strdict['gal_ra'] = "%-15s" % ra
			strdict['gal_dec'] = "%-15s" % dec
			strdict['gal_alt'] = "%6.2f" % gal_hrz.alt
			strdict['gal_az'] = "%6.2f" % gal_hrz.az
			try:
				if gal_hrz.alt > last_gal_hrz.alt:
					gal_rising = True
				else:
					gal_rising = False
			except NameError:
				pass
			last_gal_hrz = gal_hrz
			
			# calculate Tau A (Crab) position
			tau_hrz = tau_dir.hrz(currentTime)
			(ra, dec) = tau_pos.apparent_equ(currentTime).format()
			strdict['tau_ra'] = "%-15s" % ra
			strdict['tau_dec'] = "%-15s" % dec
			strdict['tau_alt'] = "%6.2f" % tau_hrz.alt
			strdict['tau_az'] = "%6.2f" % tau_hrz.az
			try:
				if tau_hrz.alt > last_tau_hrz.alt:
					tau_rising = True
				else:
					tau_rising = False
			except NameError:
				pass
			last_tau_hrz = tau_hrz
			
			# calculate Cassiopeia A position
			cas_hrz = cas_dir.hrz(currentTime)
			(ra, dec) = cas_pos.apparent_equ(currentTime).format()
			strdict['cas_ra'] = "%-15s" % ra
			strdict['cas_dec'] = "%-15s" % dec
			strdict['cas_alt'] = "%6.2f" % cas_hrz.alt
			strdict['cas_az'] = "%6.2f" % cas_hrz.az
			try:
				if cas_hrz.alt > last_cas_hrz.alt:
					cas_rising = True
				else:
					cas_rising = False
			except NameError:
				pass
			last_cas_hrz = cas_hrz
			
			# calculate Cygnus A position
			cyg_hrz = cyg_dir.hrz(currentTime)
			(ra, dec) = cyg_pos.apparent_equ(currentTime).format()
			strdict['cyg_ra'] = "%-15s" % ra
			strdict['cyg_dec'] = "%-15s" % dec
			strdict['cyg_alt'] = "%6.2f" % cyg_hrz.alt
			strdict['cyg_az'] = "%6.2f" % cyg_hrz.az
			try:
				if cyg_hrz.alt > last_cyg_hrz.alt:
					cyg_rising = True
				else:
					cyg_rising = False
			except NameError:
				pass
			last_cyg_hrz = cyg_hrz
			
			# Refresh screen and wait 1 second
			## Background
			screen.addstr( 0, 1, "==========================================================================")
			screen.addstr( 1, 1, "UTC :                                   MJD(UTC):                         ")
			screen.addstr( 2, 1, "Site:                                   LST:                              ")
			screen.addstr( 3, 1, "Long:                                   Lat:                              ")
			screen.addstr( 4, 1, "==========================================================================")
			screen.addstr( 5, 1, "----------- Sun --------------          ------------ Jupiter -------------")
			screen.addstr( 6, 1, "RA(Cur)                                 RA(Cur)                           ")
			screen.addstr( 7, 1, "Dec(Cur)                                Dec(Cur)                          ")
			screen.addstr( 8, 1, "Alt                                     Alt                               ")
			screen.addstr( 9, 1, "Az                                      Az                                ")
			screen.addstr(10, 1, "---- Galactic Center ---------          ---------- Tau A (Crab) ----------")
			screen.addstr(11, 1, "RA(Cur)                                 RA(Cur)                           ")
			screen.addstr(12, 1, "Dec(Cur)                                Dec(Cur)                          ")
			screen.addstr(13, 1, "Alt                                     Alt                               ")
			screen.addstr(14, 1, "Az                                      Az                                ")
			screen.addstr(15, 1, "---------- Cas A -------------          ------------ Cyg A ---------------")
			screen.addstr(16, 1, "RA(Cur)                                 RA(Cur)                           ")
			screen.addstr(17, 1, "Dec(Cur)                                Dec(Cur)                          ")
			screen.addstr(18, 1, "Alt                                     Alt                               ")
			screen.addstr(19, 1, "Az                                      Az                                ")
			screen.addstr(20, 1, "==========================================================================")
			## Date/time/site
			screen.addstr(1,  6, strdict['utc'], inf)
			screen.addstr(2,  6, strdict['loc'], inf)
			screen.addstr(3,  6, strdict['long'], inf)
			screen.addstr(1, 50, strdict['mjdutc'], inf)
			screen.addstr(2, 50, strdict['lst'], inf)
			screen.addstr(3, 50, strdict['lat'], inf)
			## Sun
			screen.addstr( 6, 10, strdict['sun_ra'],  std)
			screen.addstr( 7,  9, strdict['sun_dec'], std)
			screen.addstr( 8,  9, strdict['sun_alt'], dwn if sun_hrz.alt < 0 else up)
			try:
				screen.addstr( 8, 16, "^" if sun_rising else "v", rs if sun_rising else st)
			except NameError:
				pass
			screen.addstr( 9,  9, strdict['sun_az'],  dwn if sun_hrz.alt < 0 else up)
			## Jupiter
			screen.addstr( 6, 51, strdict['jup_ra'],  std)
			screen.addstr( 7, 50, strdict['jup_dec'], std)
			screen.addstr( 8, 50, strdict['jup_alt'], dwn if jup_hrz.alt < 0 else up)
			try:
				screen.addstr( 8, 57, "^" if jup_rising else "v", rs if jup_rising else st)
			except NameError:
				pass
			screen.addstr( 9, 50, strdict['jup_az'],  dwn if jup_hrz.alt < 0 else up)
			## GC
			screen.addstr(11, 10, strdict['gal_ra'],  std)
			screen.addstr(12,  9, strdict['gal_dec'], std)
			screen.addstr(13,  9, strdict['gal_alt'], dwn if gal_hrz.alt < 0 else up)
			try:
				screen.addstr(13, 16, "^" if gal_rising else "v", rs if gal_rising else st)
			except NameError:
				pass
			screen.addstr(14,  9, strdict['gal_az'],  dwn if gal_hrz.alt < 0 else up)
			## Tau A
			screen.addstr(11, 51, strdict['tau_ra'],  std)
			screen.addstr(12, 50, strdict['tau_dec'], std)
			screen.addstr(13, 50, strdict['tau_alt'], dwn if tau_hrz.alt < 0 else up)
			try:
				screen.addstr(13, 57, "^" if tau_rising else "v", rs if tau_rising else st)
			except NameError:
				pass
			screen.addstr(14, 50, strdict['tau_az'],  dwn if tau_hrz.alt < 0 else up)
			## Cas A
			screen.addstr(16, 10, strdict['cas_ra'],  std)
			screen.addstr(17,  9, strdict['cas_dec'], std)
			screen.addstr(18,  9, strdict['cas_alt'], dwn if cas_hrz.alt < 0 else up)
			try:
				screen.addstr(18, 16, "^" if cas_rising else "v", rs if cas_rising else st)
			except NameError:
				pass
			screen.addstr(19,  9, strdict['cas_az'],  dwn if cas_hrz.alt < 0 else up)
			## Cyg A
			screen.addstr(16, 51, strdict['cyg_ra'],  std)
			screen.addstr(17, 50, strdict['cyg_dec'], std)
			screen.addstr(18, 50, strdict['cyg_alt'], dwn if cyg_hrz.alt < 0 else up)
			try:
				screen.addstr(18, 57, "^" if cyg_rising else "v", rs if cyg_rising else st)
			except NameError:
				pass
			screen.addstr(19, 50, strdict['cyg_az'],  dwn if cyg_hrz.alt < 0 else up)
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
