#!/usr/bin/env python3

"""
LWA_Common.astro module informal unit test / demo: give LWA-1 observation data
at current system clock time

Moved out of the lsl.astro module into this script and updated for LWA-1
"""

import math    
import argparse

from lsl import astro
from lsl.common import stations

from lsl.misc import telemetry
telemetry.track_script()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-s', '--site', type=str, default='LWA-1',
                        help='site name')
    args = parser.parse_args()
    
    # setup transform objects
    site = args.site.lower().replace('-', '')
    if site == 'lwa1':
        station = stations.lwa1
    elif site == 'lwasv':
        station = stations.lwasv
    elif site == 'ovrolwa':
        station = stations.lwa1
        station.name = 'OVRO-LWA'
        station.lat, station.long, station.elev = ('37.23977727', '-118.2816667', 1183.48)
    else:
        raise RuntimeError(f"Unknown site name: {site}")
        
    nam = station.name
    lng = astro.deg_to_dms(station.long*180/math.pi)
    lat = astro.deg_to_dms(station.lat*180/math.pi)
    lwa_lnlat = astro.lnlat_posn(lng, lat)
    
    print('---------------------------------------------------------------')
    print('%s location' % nam)
    print('---------------------------------------------------------------')
    print('Longitude:         %s (%0.3f)' % (lng, lwa_lnlat.lng))
    print('Latitude:          %s (%0.3f)' % (lat, lwa_lnlat.lat)) 
        
    # calculate offset from GM
    
    lwa_gm_hms = astro.deg_to_hms(-lwa_lnlat.lng)
    lwa_gm_off = lwa_gm_hms.to_sec()
    print('GM Offset:         %s (%0.3f)' % (lwa_gm_hms, lwa_gm_off))
    
    # get current UTC time from system clock
    
    utc = astro.get_julian_from_sys()
    utcd = astro.get_date(utc)
    lcld = utcd.to_zone()
    unixt = astro.get_timet_from_julian(utc)
    
    # caculate sidereal times
    
    gm_sid = astro.get_apparent_sidereal_time(utc)
    lwa_sid = astro.get_local_sidereal_time(lwa_lnlat.lng, utc)
        
    print('---------------------------------------------------------------')
    print('Current time')
    print('---------------------------------------------------------------')
    print('UTC time:             %s (%0.3f)' % (utcd, utc))
    print('%s local time:     %s' % (nam, str(lcld)))
    print('GM sidereal time:     %0.3f' % gm_sid)
    print('%s sidereal time:  %0.3f' % (nam, lwa_sid))
    print('UNIX time:            %d' % unixt)   
    
    # calculate solar system phenomena
    sun_equ = astro.get_solar_equ_coords(utc)
    sun_ecl = sun_equ.to_ecl(utc)
    
    names = ('Sun', 'Moon', 'Venus', 'Mars', 'Jupiter', 'Saturn')
    rst_funcs = (astro.get_solar_rst, astro.get_lunar_rst, astro.get_venus_rst, astro.get_mars_rst, astro.get_jupiter_rst, astro.get_saturn_rst)
    equ_funcs = (astro.get_solar_equ_coords, astro.get_lunar_equ_coords, astro.get_venus_equ_coords, astro.get_mars_equ_coords, astro.get_jupiter_equ_coords, astro.get_saturn_equ_coords)
    for name,rst_func,equ_func in zip(names, rst_funcs, equ_funcs):
        bdy_rst = rst_func(utc, lwa_lnlat)
        (bdy_utc_rise, bdy_utc_set, bdy_utc_trans) = bdy_rst.format()
        bdy_lcl_rise = bdy_utc_rise.to_zone()
        bdy_lcl_trans = bdy_utc_trans.to_zone()
        bdy_lcl_set = bdy_utc_set.to_zone()
        
        bdy_equ = equ_func(utc) 
        (bdy_ra, bdy_dec) = bdy_equ.format()
        bdy_hrz = bdy_equ.to_hrz(lwa_lnlat, utc)
        bdy_ecl = bdy_equ.to_ecl(utc)
        (bdy_lng, bdy_lat) = bdy_ecl.format()
        
        if name != 'Sun':
            bdy_sun_ang = sun_equ.angular_separation(bdy_equ)
            
        if name == 'Venus':
            venus_elong = sun_ecl.lng - bdy_ecl.lng
            if venus_elong < 0:
                venus_dir = 'E'
            else:
                venus_dir = 'W'
            
        print('---------------------------------------------------------------')
        print(name)
        print('---------------------------------------------------------------')
        print('RA:                %s (%0.3f)' % (bdy_ra, bdy_equ.ra))
        print('DEC:               %s (%0.3f)' % (bdy_dec, bdy_equ.dec))
        print('Ecl longitude:     %s (%0.3f)' % (bdy_lng, bdy_ecl.lng))
        print('Ecl latitude:      %s (%0.3f)' % (bdy_lat, bdy_ecl.lat))      
        print('Rise:              %s (%0.3f) [%s]' % (bdy_utc_rise, bdy_rst.rise, bdy_lcl_rise))
        print('Transit:           %s (%0.3f) [%s]' % (bdy_utc_trans, bdy_rst.transit, bdy_lcl_trans))
        print('Set:               %s (%0.3f) [%s]' % (bdy_utc_set, bdy_rst.set, bdy_lcl_set))
        print('Azimuth:           %0.3f %s' % (bdy_hrz.az, astro.hrz_to_nswe(bdy_hrz)))
        print('Altitude:          %0.3f' % bdy_hrz.alt)
        print('Zenith:            %0.3f' % bdy_hrz.zen())
        if name != 'Sun':
            print('Sun angle:         %0.3f' % bdy_sun_ang)
        if name == 'Venus':
            print('Elongation:        %s %s (%0.3f)' % (astro.deg_to_dms(venus_elong), venus_dir, venus_elong))
            
    # calculate celestial phenomena
    names = ('SgrA', 'CasA', 'CygA', 'TauA', 'B1919+21')
    bdy_coords = ((astro.hms(17, 42, 48.1), astro.dms(True, 28, 55, 8)),
                  (astro.hms(23, 23, 22.7), astro.dms(False, 58, 49, 16)),
                  (astro.hms(19, 59, 27.8), astro.dms(False, 40, 44, 2)),
                  (astro.hms(5, 34, 32.0), astro.dms(False, 22, 00, 52.1)),
                  (astro.hms(19, 21, 44.80), astro.dms(False, 21, 53, 1.8)))
    for name,(bdy_ra,bdy_dec) in zip(names, bdy_coords):
        bdy_j2000_equ = astro.equ_posn(bdy_ra, bdy_dec)
        bdy_equ = astro.get_apparent_posn(bdy_j2000_equ, utc)
        bdy_rst = astro.get_object_rst(utc, lwa_lnlat, bdy_equ)
        try:
            (bdy_utc_rise, bdy_utc_set, bdy_utc_trans) = bdy_rst.format()
            bdy_lcl_rise = bdy_utc_rise.to_zone()
            bdy_lcl_trans = bdy_utc_trans.to_zone()
            bdy_lcl_set = bdy_utc_set.to_zone()
        except AttributeError:
            bdy_utc_rise = bdy_utc_set = bdy_utc_trans = None
            
        (bdy_ra, bdy_dec) = bdy_equ.format()
        bdy_hrz = bdy_equ.to_hrz(lwa_lnlat, utc)
        bdy_gal = bdy_equ.to_gal(utc)
        (bdy_l, bdy_b) = bdy_gal.format()
        
        bdy_sun_ang = sun_equ.angular_separation(bdy_equ)
        
        print('---------------------------------------------------------------')
        print(name)
        print('---------------------------------------------------------------')
        print('RA:                %s (%0.3f)' % (bdy_ra, bdy_equ.ra))
        print('DEC:               %s (%0.3f)' % (bdy_dec, bdy_equ.dec))
        print('Gal longitude:     %s (%0.3f)' % (bdy_l, bdy_gal.l))
        print('Gal latitude:      %s (%0.3f)' % (bdy_b, bdy_gal.b)) 
        if bdy_utc_rise is not None:
            print('Rise:              %s (%0.3f) [%s]' % (bdy_utc_rise, bdy_rst.rise, bdy_lcl_rise))
            print('Transit:           %s (%0.3f) [%s]' % (bdy_utc_trans, bdy_rst.transit, bdy_lcl_trans))
            print('Set:               %s (%0.3f) [%s]' % (bdy_utc_set, bdy_rst.set, bdy_lcl_set))
        else:
            print('Rise:              Circumpolar')
            print('Transit:           -----------')
            print('Set:               Circumpolar')
        print('Azimuth:           %0.3f %s' % (bdy_hrz.az, astro.hrz_to_nswe(bdy_hrz)))
        print('Altitude:          %0.3f' % bdy_hrz.alt)
        print('Zenith:            %0.3f' % bdy_hrz.zen())
        print('Sun angle:         %0.3f' % bdy_sun_ang)
