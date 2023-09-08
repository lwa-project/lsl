#!/usr/bin/env python

"""
Script for making and displaying images of correlated data files.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    input = raw_input
    
import os
import sys
import pytz
import numpy
import argparse
from astropy.io import fits as astrofits

from lsl import astro
from lsl.sim import vis as simVis
from lsl.imaging import utils, overlay
from lsl.writer.fitsidi import NUMERIC_STOKES
from lsl.misc import parser as aph

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from lsl.misc import telemetry
telemetry.track_script()


MST = pytz.timezone('US/Mountain')
UTC = pytz.UTC


NPIX_SIDE = 350


def main(args):
    filename = args.filename
    
    idi = utils.CorrelatedData(filename)
    aa = idi.get_antennaarray()
    lo = idi.get_observer()
    
    nStand = len(idi.stands)
    nchan = len(idi.freq)
    freq = idi.freq
    
    print("Raw Stand Count: %i" % nStand)
    print("Final Baseline Count: %i" % (nStand*(nStand-1)/2,))
    print("Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nchan, (freq[1] - freq[0])/1e3))
    try:
        print("Polarization Products: %s" % ' '.join([NUMERIC_STOKES[p] for p in idi.pols]))
    except KeyError:
        # Catch for CASA MS that use a different numbering scheme
        NUMERIC_STOKESMS = {1:'I', 2:'Q', 3:'U', 4:'V', 
                        9:'XX', 10:'XY', 11:'YX', 12:'YY'}
        print("Polarization Products: %s" % ' '.join([NUMERIC_STOKESMS[p] for p in idi.pols]))
        
    print("Reading in FITS IDI data")
    nSets = idi.integration_count
    for set in range(1, nSets+1):
        if args.dataset != -1 and args.dataset != set:
            continue
            
        print("Set #%i of %i" % (set, nSets))
        dataDict = idi.get_data_set(set, min_uv=args.uv_min)
        
        # Build a list of unique JDs for the data
        pols = dataDict.pols
        jdList = [dataDict.mjd + astro.MJD_OFFSET,]
        
        # Find the LST
        lo.date = jdList[0] - astro.DJD_OFFSET
        utc = str(lo.date)
        lst = str(lo.sidereal_time())   # pylint:disable=no-member
        
        # Pull out the right channels
        toWork = numpy.where( (freq >= args.freq_start) & (freq <= args.freq_stop) )[0]
        if len(toWork) == 0:
            raise RuntimeError("Cannot find data between %.2f and %.2f MHz" % (args.freq_start/1e6, args.freq_stop/1e6))
            
        # Integration report
        print("    Date Observed: %s" % utc)
        print("    LST: %s" % lst)
        print("    Selected Frequencies: %.3f to %.3f MHz" % (freq[toWork[0]]/1e6, freq[toWork[-1]]/1e6))
        
        # Prune out what needs to go
        if args.include is not None or args.exclude is not None:
            print("    Processing include/exclude lists")
            dataDict = dataDict.get_antenna_subset(include=args.include, 
                                                   exclude=args.exclude, 
                                                   indicies=False)
            
            ## Report
            for pol in dataDict.pols:
                print("        %s now has %i baselines" % (pol, len(dataDict.baselines)))
                
        # Build up the images for each polarization
        print("    Gridding")
        img1 = None
        lbl1 = 'XX'
        for p in ('XX', 'RR', 'I'):
            try:
                img1 = utils.build_gridded_image(dataDict, size=NPIX_SIDE//2, res=0.5, pol=p, chan=toWork)
                lbl1 = p.upper()
            except:
                pass
                
        img2 = None
        lbl2 = 'YY'
        for p in ('YY', 'LL', 'Q'):
            try:
                img2 = utils.build_gridded_image(dataDict, size=NPIX_SIDE//2, res=0.5, pol=p, chan=toWork)
                lbl2 = p.upper()
            except:
                pass
                
        img3 = None
        lbl3 = 'XY'
        for p in ('XY', 'RL', 'U'):
            try:
                img3 = utils.build_gridded_image(dataDict, size=NPIX_SIDE//2, res=0.5, pol=p, chan=toWork)
                lbl3 = p.upper()
            except:
                pass
                
        img4 = None
        lbl4 = 'YX'
        for p in ('YX', 'LR', 'V'):
            try:
                img4 = utils.build_gridded_image(dataDict, size=NPIX_SIDE//2, res=0.5, pol=p, chan=toWork)
                lbl4 = p.upper()
            except:
                pass
                
        # Plots
        print("    Plotting")
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        for ax, img, pol in zip([ax1, ax2, ax3, ax4], [img1, img2, img3, img4], [lbl1, lbl2, lbl3, lbl4]):
            # Skip missing images
            if img is None:
                ax.text(0.5, 0.5, 'Not found in file', color='black', size=12, horizontalalignment='center')
                
                ax.xaxis.set_major_formatter( NullFormatter() )
                ax.yaxis.set_major_formatter( NullFormatter() )
                
                if not args.utc:
                    ax.set_title("%s @ %s LST" % (pol, lst))
                else:
                    ax.set_title("%s @ %s UTC" % (pol, utc))
                continue
                
            # Display the image, save the limits, and label with the polarization/LST
            cb = utils.plot_gridded_image(ax, img)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            fig.colorbar(cb, ax=ax)
            if not args.utc:
                ax.set_title("%s @ %s LST" % (pol, lst))
            else:
                ax.set_title("%s @ %s UTC" % (pol, utc))
                
            junk = img.image(center=(img.shape[0]//2,img.shape[1]//2))
            print("%s: image is %.4f to %.4f with mean %.4f" % (pol, junk.min(), junk.max(), junk.mean()))
            
            # Turn off tick marks
            ax.xaxis.set_major_formatter( NullFormatter() )
            ax.yaxis.set_major_formatter( NullFormatter() )
            
            # Compute the positions of major sources and label the images
            overlay.sources(ax, aa, simVis.SOURCES, label=not args.no_labels)
            
            # Add in the horizon
            overlay.horizon(ax, aa)
            
            # Add lines of constant RA and dec.
            if not args.no_grid:
                if not args.topo:
                    overlay.graticule_radec(ax, aa)
                else:
                    overlay.graticule_azalt(ax, aa)

            # Reset the axes
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        plt.show()
        
        if args.fits is not None:
            ## Loop over the images to build up the FITS file
            hdulist = [astrofits.PrimaryHDU(),]
            for img,pol in zip((img1,img2,img3,img4), (lbl1,lbl2,lbl3,lbl4)):
                if img is None:
                    continue
                    
                ### Create the HDU
                try:
                    hdu = astrofits.ImageHDU(data=img.image(center=(img.shape[0]//2,img.shape[1]//2)), name=pol)
                except AttributeError:
                    hdu = astrofits.ImageHDU(data=img, name=pol)
                    
                ### Add in the coordinate information
                hdu.header['EPOCH'] = 2000.0 + (jdList[0] - 2451545.0) / 365.25
                hdu.header['CTYPE1'] = 'RA---SIN'
                hdu.header['CRPIX1'] = img.shape[0]//2+1
                hdu.header['CDELT1'] = -360.0/img.shape[0]/numpy.pi
                hdu.header['CRVAL1'] = lo.sidereal_time()*180/numpy.pi	# pylint:disable=no-member
                hdu.header['CTYPE2'] = 'DEC--SIN'
                hdu.header['CRPIX2'] = img.shape[1]//2+1
                hdu.header['CDELT2'] = 360.0/img.shape[1]/numpy.pi
                hdu.header['CRVAL2'] = lo.lat*180/numpy.pi
                hdu.header['LONPOLE'] = 180.0
                hdu.header['LATPOLE'] = 90.0
                
                ### Add the HDU to the list
                hdulist.append(hdu)
                
            ## Save the FITS file to disk
            hdulist = astrofits.HDUList(hdulist)
            overwrite = False
            if os.path.exists(args.fits):
                yn = input("WARNING: '%s' exists, overwrite? [Y/n]" % args.fits)
                if yn not in ('n', 'N'):
                    overwrite = True
            try:
                hdulist.writeto(args.fits, overwrite=overwrite)
            except IOError as e:
                print("WARNING: FITS image file not saved")
                
    print("...Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='create images from a FITS-IDI file', 
        epilog='NOTE:  If both -i/--include and -e/--exclude are specified the include list has priority.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to plot')
    parser.add_argument('-1', '--freq-start', type=aph.frequency, default='10.0', 
                        help='first frequency to analyze in MHz')
    parser.add_argument('-2', '--freq-stop', type=aph.frequency, default='88.0', 
                        help='last frequency to analyze in MHz')
    parser.add_argument('-s', '--dataset', type=int, default=-1, 
                        help='data set to image')
    parser.add_argument('-m', '--uv-min', type=float, default=0.0, 
                        help='minimun baseline uvw length to include in lambda at the midpoint frequency')
    parser.add_argument('-i', '--include', type=aph.csv_int_list, default='all', 
                        help='comma seperated list of dipoles to include')
    parser.add_argument('-e', '--exclude', type=aph.csv_int_list, default='none', 
                        help='comma seperated list of dipoles to exclude')
    parser.add_argument('-t', '--topo', action='store_true', 
                        help='display an az/el grid instead of a RA/Dec grid')
    parser.add_argument('-u', '--utc', action='store_true', 
                        help='display the time as UTC, instead of LST')
    parser.add_argument('-n', '--no-labels', action='store_true', 
                        help='disable source and grid labels')
    parser.add_argument('-g', '--no-grid', action='store_true', 
                        help='disable the coordinate grid')
    parser.add_argument('-f', '--fits', type=str, 
                        help='save the images to the specified FITS image file')
    args = parser.parse_args()
    if args.include == 'all':
        args.include = None
    if args.exclude == 'none':
        args.exclude = None
    main(args)
