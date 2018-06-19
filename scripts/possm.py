#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy
import getopt

from lsl import astro
from lsl.imaging import utils
from lsl.common.progress import ProgressBar

import matplotlib.pyplot as plt


def usage(exitCode=None):
    print """possm.py - Create AIPS POSSM style plots of FITS IDI files

Usage: possm.py [OPTIONS] file

Options:
-h, --help             Display this help information
-1, --freq-start       First frequency to image in MHz (Default = 10 MHz)
-2, --freq-stop        Last frequency to image in MHz (Default = 88 MHz)
-s, --dataset          Data set to image (Default = All)
-m, --uv-min           Minimun baseline uvw length to include 
                    (Default = 0 lambda at midpoint frequency)
-i, --include          Comma seperated list of dipoles to include 
                    (Default = All)
-e, --exclude          Comma seperated list of dipoles to exclude
                    (Default = None)
-f, --frequency        Label the channels in frequency rather than channel
                    number (Default = no, use channel number)
-l, --log              Use a log scale for the visbility amplitudes
                    (Default = no, use linear)

NOTE:  If both -i/--include and -e/--exclude are specified the include list
    has priority.
"""

    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseConfig(args):
    config = {}
    # Command line flags - default values
    config['freq1'] = 10e6
    config['freq2'] = 88e6
    config['dataset'] = 0
    config['uvMin'] = 0.0
    config['include'] = None
    config['exclude'] = None
    config['labelFreq'] = False
    config['ampLog'] = False
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "h1:2:s:m:i:e:fl", ["help", "freq-start=", "freq-stop=", "dataset=", "uv-min=", "include=", "exclude=", "frequency", "log"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-1', '--freq-start'):
            config['freq1'] = float(value)*1e6
        elif opt in ('-2', '--freq-stop'):
            config['freq2'] = float(value)*1e6
        elif opt in ('-s', '--dataset'):
            config['dataset'] = int(value)
        elif opt in ('-m', '--uv-min'):
            config['uvMin'] = float(value)
        elif opt in ('-i', '--include'):
            config['include'] = [int(v) for v in value.split(',')]
        elif opt in ('-e', '--exclude'):
            config['exclude'] = [int(v) for v in value.split(',')]
        elif opt in ('-f', '--frequency'):
            config['labelFreq'] = True
        elif opt in ('-l', '--log'):
            config['ampLog'] = True
        else:
            assert False
    
    # Add in arguments
    config['args'] = arg
    
    # Return configuration
    return config


def main(args):
    # Parse the command line
    config = parseConfig(args)
    
    # Grab the filename
    filename = config['args'][0]
    
    idi = utils.CorrelatedData(filename)
    aa = idi.getAntennaArray()
    lo = idi.getObserver()
    lo.date = idi.dateObs.strftime("%Y/%m/%d %H:%M:%S")
    jd = lo.date + astro.DJD_OFFSET
    lst = str(lo.sidereal_time())

    nStand = len(idi.stands)
    nChan = len(idi.freq)
    freq = idi.freq
    
    print "Raw Stand Count: %i" % nStand
    print "Final Baseline Count: %i" % (nStand*(nStand-1)/2,)
    print "Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nChan, (freq[1] - freq[0])/1e3)
    print "Polarization Products: %i starting with %i" % (len(idi.pols), idi.pols[0])
    print "JD: %.3f" % jd

    print "Reading in FITS IDI data"
    nSets = idi.integrationCount
    for set in range(1, nSets+1):
        if config['dataset'] != 0 and config['dataset'] != set:
            continue
            
        print "Set #%i of %i" % (set, nSets)
        dataDict = idi.getDataSet(set, uvMin=config['uvMin'])
        
        # Prune out what needs to go
        if config['include'] is not None or config['exclude'] is not None:
            print "    Processing include/exclude lists"
            
            ## Create an empty output data dictionary
            newDict = {}
            for key in dataDict.keys():
                if key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
                    newDict[key] = {}
                    for pol in dataDict[key].keys():
                        newDict[key][pol] = []
                else:
                    newDict[key] = dataDict[key]
                    
            ## Fill it
            for pol in dataDict['bls'].keys():
                for b in xrange(len(dataDict['bls'][pol])):
                    a0,a1 = dataDict['bls'][pol][b]
                    
                    if config['include'] is not None:
                        if idi.stands[a0] in config['include'] and idi.stands[a1] in config['include']:
                            for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
                                newDict[key][pol].append( dataDict[key][pol][b] )
                        else:
                            continue
                            
                    if config['exclude'] is not None:
                        if idi.stands[a0] in config['exclude'] or idi.stands[a1] in config['exclude']:
                            continue
                        else:
                            for key in ['bls', 'uvw', 'vis', 'wgt', 'msk', 'jd']:
                                newDict[key][pol].append( dataDict[key][pol][b] )
                                
            ## Make the substitution so that we are working with the right data now
            dataDict = newDict
            
            ## Report
            for pol in dataDict['bls'].keys():
                print "        %s now has %i baselines" % (pol, len(dataDict['bls'][pol]))
                
        # Pull out the right channels
        toWork = numpy.where( (freq >= config['freq1']) & (freq <= config['freq2']) )[0]
        if len(toWork) == 0:
            raise RuntimeError("Cannot find data between %.2f and %.2f MHz" % (config['freq1']/1e6, config['freq2']/1e6))
        if config['labelFreq']:
            xValues = freq[toWork]/1e6
            xLabel = 'Frequency [MHz]'
        else:
            xValues = toWork
            xLabel = 'Channel'
            
        # Plot
        print "    Plotting the first 50 baselines"
        pols = dataDict['jd'].keys()
        nBL = len(dataDict['bls'][pols[0]])
        pb = ProgressBar(max=nBL)
        i = 0
        for k in range(2):
            fig = plt.figure()

            for j in range(25):
                try:
                    stnd1, stnd2 = dataDict['bls'][pols[0]][i]
                    stnd1 = idi.stands[stnd1]
                    stnd2 = idi.stands[stnd2]
                    vis = dataDict['vis'][pols[0]][i]
                    i += 1
                except IndexError:
                    plt.draw()
                    break

                amp = numpy.abs(vis)
                if config['ampLog']:
                    amp = numpy.log10(amp)
                phs = numpy.angle(vis)*180/numpy.pi

                ax = fig.add_subplot(10, 5, 2*(j/5)*5+j%5+1)
                if ((phs+360)%360).std() < phs.std():
                    ax.plot(xValues, (phs[toWork]+360)%360, linestyle=' ', marker='x')
                    ax.set_ylim([0, 360])
                else:
                    ax.plot(xValues, phs[toWork], linestyle=' ', marker='x')
                    ax.set_ylim([-180, 180])
                ax.set_title('%i - %i' % (stnd1, stnd2))
                if j % 5 == 0:
                    ax.set_ylabel('Phs')
                    
                ax = fig.add_subplot(10, 5, 2*(j/5)*5+j%5+1+5)
                ax.plot(xValues, amp[toWork], linestyle=' ', marker='x', color='green')
                ax.set_title('%i - %i' % (stnd1, stnd2))
                if j % 5 == 0:
                    ax.set_xlabel(xLabel)
                    ax.set_ylabel('%sAmp' % '' 'Log ' if config['ampLog'] else '')
                    
                pb.inc(amount=1)
                if pb.amount != 0 and pb.amount % 10 == 0:
                    sys.stdout.write(pb.show()+'\r')
                    sys.stdout.flush()
            plt.draw()

        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        plt.show()
    

if __name__ == "__main__":
    numpy.seterr(all='ignore')
    main(sys.argv[1:])
