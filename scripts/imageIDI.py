#!/usr/bin/env python3

"""
Script for making and displaying images of correlated data files.
"""

import os
import sys
import pytz
import numpy as np
import argparse

from astropy.io import fits as astrofits
from astropy import units as astrounits
from astropy.time import Time as AstroTime
from astropy.coordinates import FK5

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



def main(args):
    filename = args.filename
    
    idi = utils.CorrelatedData(filename)
    aa = idi.get_antennaarray()
    lo = idi.get_observer()
    
    nStand = len(idi.stands)
    nchan = len(idi.freq)
    freq = idi.freq
    
    chan_width = freq[1]-freq[0]
    chan = np.round(freq / chan_width)
    nif = len(np.where(np.diff(chan) > 1)[0]) + 1
    freq_if = freq*1.0
    freq_if = freq_if.reshape(nif, -1)  
    
    print(f"Raw Stand Count: {nStand}")
    print(f"Final Baseline Count: {nStand*(nStand-1)//2}")
    print(f"Spectra Coverage: {freq_if[0,0]/1e6:.3f} to {freq_if[0,-1]/1e6:.3f} MHz in {freq_if.shape[1]} channels ({(freq_if[0,1] - freq_if[0,0])/1e3:.2f} kHz/channel)")
    for i in range(1, nif):
        print(f"                  {freq_if[i,0]/1e6:.3f} to {freq_if[i,-1]/1e6:.3f} MHz in {freq_if.shape[1]} channels ({(freq_if[i,1] - freq_if[i,0])/1e3:.2f} kHz/channel)")
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
            
        print(f"Set #{set} of {nSets}")
        dataDict = idi.get_data_set(set, min_uv=args.uv_min)
        
        # Build a list of unique JDs for the data
        pols = dataDict.pols
        jdList = [dataDict.mjd + astro.MJD_OFFSET,]
        
        # Find the LST
        lo.date = jdList[0] - astro.DJD_OFFSET
        utc = str(lo.date)
        lst = str(lo.sidereal_time())   # pylint:disable=no-member
        
        # Pull out the right channels
        toWork = np.where( (freq >= args.freq_start) & (freq <= args.freq_stop) )[0]
        if len(toWork) == 0:
            raise RuntimeError(f"Cannot find data between {args.freq_start/1e6:.3f} and {args.freq_stop/1e6:.3f} MHz")
            
        # Integration report
        print(f"    Date Observed: {str(utc)}")
        print(f"    LST: {str(lst)}")
        print(f"    Selected Frequencies: {freq[toWork[0]]/1e6:.3f} to {freq[toWork[-1]]/1e6:.3f} MHz")
        
        # Prune out what needs to go
        if args.include != 'all' or args.exclude != 'none':
            print("    Processing include/exclude lists")
            dataDict = dataDict.get_antenna_subset(include=args.include, 
                                                   exclude=args.exclude, 
                                                   indicies=False)
            
            ## Report
            for pol in dataDict.pols:
                print(f"        {pol} now has {len(dataDict.baselines)} baselines")
                
        # Build up the images for each polarization
        print("    Gridding")
        img1 = None
        lbl1 = 'XX'
        for p in ('XX', 'RR', 'I'):
            try:
                img1 = utils.build_gridded_image(dataDict, im_size=args.npix_side, im_res=args.image_res, pol=p, chan=toWork)
                lbl1 = p.upper()
            except:
                pass
                
        img2 = None
        lbl2 = 'YY'
        for p in ('YY', 'LL', 'Q'):
            try:
                img2 = utils.build_gridded_image(dataDict, im_size=args.npix_side, im_res=args.image_res, pol=p, chan=toWork)
                lbl2 = p.upper()
            except:
                pass
                
        img3 = None
        lbl3 = 'XY'
        for p in ('XY', 'RL', 'U'):
            try:
                img3 = utils.build_gridded_image(dataDict, im_size=args.npix_side, im_res=args.image_res, pol=p, chan=toWork)
                lbl3 = p.upper()
            except:
                pass
                
        img4 = None
        lbl4 = 'YX'
        for p in ('YX', 'LR', 'V'):
            try:
                img4 = utils.build_gridded_image(dataDict, im_size=args.npix_side, im_res=args.image_res, pol=p, chan=toWork)
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
                    ax.set_title(f"{pol} @ {str(lst)} LST")
                else:
                    ax.set_title(f"{pol} @ {str(utc)} UTC")
                continue
                
            # Display the image, save the limits, and label with the polarization/LST
            cb = utils.plot_gridded_image(ax, img)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            fig.colorbar(cb, ax=ax)
            if not args.utc:
                ax.set_title(f"{pol} @ {str(lst)} LST")
            else:
                ax.set_title(f"{pol} @ {str(utc)} UTC")
                
            junk = img.image(center=(img.shape[0]//2,img.shape[1]//2))
            print(f"{pol}: image is {junk.min():.4f} to {junk.max():.4f} with mean {junk.mean():.4f}")
            
            # Turn off tick marks
            ax.xaxis.set_major_formatter( NullFormatter() )
            ax.yaxis.set_major_formatter( NullFormatter() )
            
            # Compute the positions of major sources and label the images
            overlay.sources(ax, img, simVis.SOURCES, label=not args.no_labels)
            
            # Add in the horizon
            overlay.horizon(ax, img)
            
            # Add lines of constant RA and dec.
            if not args.no_grid:
                if not args.topo:
                    overlay.graticule_radec(ax, img)
                else:
                    overlay.graticule_azalt(ax, img)

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
                wcs_hdr = img.wcs.to_header()
                for key in wcs_hdr:
                    hdu.header[key] = wcs_hdr[key]
                    
                ### Add the HDU to the list
                hdulist.append(hdu)
                
            ## Save the FITS file to disk
            hdulist = astrofits.HDUList(hdulist)
            fitsname = args.fits
            if args.dataset == -1 and nSets > 1:
                fitsparts = os.path.splitext(args.fits)
                fitsname = f"{fitsparts[0]}-t{set}{fitsparts[1]}"
                
            overwrite = False
            if os.path.exists(fitsname):
                yn = input(f"WARNING: '{fitsname}' exists, overwrite? [Y/n]")
                if yn not in ('n', 'N'):
                    overwrite = True
            try:
                hdulist.writeto(fitsname, overwrite=overwrite)
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
    parser.add_argument('-p', '--npix-side', type=int, default=350,
                        help='number of pixels per side in image plane')
    parser.add_argument('-r', '--image-res', type=float, default=0.37142,
                        help='resolution in the image plane as degrees/pix')
    parser.add_argument('-s', '--dataset', type=int, default=-1, 
                        help='data set to image')
    parser.add_argument('-m', '--uv-min', type=float, default=0.0, 
                        help='minimum baseline uvw length to include in lambda at the midpoint frequency')
    parser.add_argument('-i', '--include', type=aph.csv_int_list, default='all',
                        help='comma separated list of dipoles to include')
    parser.add_argument('-e', '--exclude', type=aph.csv_int_list, default='none',
                        help='comma separated list of dipoles to exclude')
    parser.add_argument('-t', '--topo', action='store_true', 
                        help='display an az/el grid instead of a RA/Dec grid')
    parser.add_argument('-u', '--utc', action='store_true', 
                        help='display the time as UTC, instead of LST')
    parser.add_argument('-n', '--no-labels', action='store_true', 
                        help='disable source and grid labels')
    parser.add_argument('-g', '--no-grid', action='store_true', 
                        help='disable the coordinate grid')
    parser.add_argument('-f', '--fits', type=str, 
                        help='save the images to the FITS image file; if multiple data sets are imaged the filename is updated to include a "-t<data_set_number>" tag')
    args = parser.parse_args()
    main(args)
