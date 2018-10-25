#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given a TBW file, plot the time averaged spectra for each digitizer input."""

import os
import sys
import math
import numpy
import ephem
import getopt

from lsl.common import stations, metabundle
from lsl.reader.ldp import LWA1DataFile
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET

import matplotlib.pyplot as plt


def usage(exitCode=None):
    print """tbwSpectra.py - Read in TBW files and create a collection of 
time-averaged spectra.

Usage: tbwSpectra.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-m, --metadata              Name of SSMIF or metadata tarball file to use for 
                            mappings
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-q, --quiet                 Run tbwSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
-g, --gain-correct          Correct signals for the cable losses
-s, --stack                 Stack spectra in groups of 6 (if '-g' is enabled only)
-o, --output                Output file name for spectra image
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['metadata'] = ''
    config['LFFT'] = 4096
    config['maxFrames'] = 30000*260
    config['window'] = fxc.null_window
    config['applyGain'] = False
    config['stack'] = False
    config['output'] = None
    config['verbose'] = True
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hm:qtbnl:gso:", ["help", "metadata=", "quiet", "bartlett", "blackman", "hanning", "fft-length=", "gain-correct", "stack", "output="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-m', '--metadata'):
            config['metadata'] = value
        elif opt in ('-q', '--quiet'):
            config['verbose'] = False
        elif opt in ('-t', '--bartlett'):
            config['window'] = numpy.bartlett
        elif opt in ('-b', '--blackman'):
            config['window'] = numpy.blackman
        elif opt in ('-n', '--hanning'):
            config['window'] = numpy.hanning
        elif opt in ('-l', '--fft-length'):
            config['LFFT'] = int(value)
        elif opt in ('-g', '--gain-correct'):
            config['applyGain'] = True
        elif opt in ('-s', '--stack'):
            config['stack'] = True
        elif opt in ('-o', '--output'):
            config['output'] = value
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def main(args):
    # Parse command line options
    config = parseOptions(args)
    
    # Setup the LWA station information
    if config['metadata'] != '':
        try:
            station = stations.parse_ssmif(config['metadata'])
        except ValueError:
            station = metabundle.get_station(config['metadata'], apply_sdm=True)
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    # Length of the FFT
    LFFT = config['LFFT']
    
    # Make sure that the file chunk size contains is an integer multiple
    # of the FFT length so that no data gets dropped
    maxFrames = int(config['maxFrames']/float(LFFT))*LFFT
    # It seems like that would be a good idea, however...  TBW data comes one
    # capture at a time so doing something like this actually truncates data 
    # from the last set of stands for the first integration.  So, we really 
    # should stick with
    maxFrames = config['maxFrames']
    
    idf = LWA1DataFile(config['args'][0])
    
    nFrames = idf.get_info('nFrames')
    srate = idf.get_info('sample_rate')
    dataBits = idf.get_info('dataBits')
    # The number of ant/pols in the file is hard coded because I cannot figure out 
    # a way to get this number in a systematic fashion
    antpols = len(antennas)
    nChunks = int(math.ceil(1.0*nFrames/maxFrames))
    if dataBits == 12:
        nSamples = 400
    else:
        nSamples = 1200
        
    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    beginDate = ephem.Date(unix_to_utcjd(idf.get_info('tStart')) - DJD_OFFSET)
    
    # File summary
    print "Filename: %s" % config['args'][0]
    print "Date of First Frame: %s" % str(beginDate)
    print "Ant/Pols: %i" % antpols
    print "Sample Length: %i-bit" % dataBits
    print "Frames: %i" % nFrames
    print "Chunks: %i" % nChunks
    print "==="
    
    # Master loop over all of the file chunks
    nChunks = 1
    masterSpectra = numpy.zeros((nChunks, antpols, LFFT))
    masterWeight = numpy.zeros((nChunks, antpols, LFFT))
    
    readT, t, data = idf.read(0.061)
    
    # Calculate the spectra for this block of data and then weight the results by 
    # the total number of frames read.  This is needed to keep the averages correct.
    # NB:  The weighting is the same for the x and y polarizations because of how 
    # the data are packed in TBW
    freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'])
    for stand in xrange(masterSpectra.shape[1]):
        masterSpectra[0,stand,:] = tempSpec[stand,:]
        masterWeight[0,stand,:] = int(readT*srate/LFFT)
        
    # We don't really need the data array anymore, so delete it
    del(data)
    
    # Apply the cable loss corrections, if requested
    if config['applyGain']:
        for s in xrange(masterSpectra.shape[1]):
            currGain = antennas[s].cable.gain(freq)
            for c in xrange(masterSpectra.shape[0]):
                masterSpectra[c,s,:] /= currGain
                
    # Now that we have read through all of the chunks, perform the final averaging by
    # dividing by all of the chunks
    spec = masterSpectra.mean(axis=0)
    
    # The plots:  This is setup for the current configuration of 20 antpols
    if config['applyGain'] & config['stack']:
        # Stacked spectra - only if cable loss corrections are to be applied
        colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 
            'purple', 'salmon', 'olive', 'maroon', 'saddlebrown', 'yellowgreen', 
            'teal', 'steelblue', 'seagreen', 'slategray', 'mediumorchid', 'lime', 
            'dodgerblue', 'darkorange']
            
        for f in xrange(numpy.ceil(antpols/20)):
            fig = plt.figure()
            ax1 = fig.add_subplot(1, 1, 1)
            for i in xrange(f*20, f*20+20):
                currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
                ax1.plot(freq/1e6, currSpectra, label='%i,%i' % (antennas[i].stand.id, antennas[i].pol), color=colors[i % 20])
                
            ax1.set_xlabel('Frequency [MHz]')
            ax1.set_ylabel('P.S.D. [dB/RBW]')
            ax1.set_xlim([20,88])
            #ax1.set_ylim([10,90])
            leg = ax1.legend(loc=0, ncol=3)
            for l in leg.get_lines():
                l.set_linewidth(1.7)  # the legend line width
    else:
        for f in xrange(int(numpy.ceil(antpols/20))):
            # Normal plotting
            fig = plt.figure()
            figsY = 4
            figsX = 5
            fig.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.20, hspace=0.50)
            for i in xrange(f*20, f*20+20):
                ax = fig.add_subplot(figsX, figsY, (i%20)+1)
                try:
                    currSpectra = numpy.squeeze( numpy.log10(spec[i,:])*10.0 )
                except IndexError:
                    break
                ax.plot(freq/1e6, currSpectra, label='Stand: %i, Pol: %i (Dig: %i)' % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer))
                
                # If there is more than one chunk, plot the difference between the global 
                # average and each chunk
                if nChunks > 1:
                    for j in xrange(nChunks):
                        # Some files are padded by zeros at the end and, thus, carry no 
                        # weight in the average spectra.  Skip over those.
                        if masterWeight[j,i,:].sum() == 0:
                            continue
                            
                        # Calculate the difference between the spectra and plot
                        subspectra = numpy.squeeze( numpy.log10(masterSpectra[j,i,:])*10.0 )
                        diff = subspectra - currSpectra
                        ax.plot(freq/1e6, diff)
                        
                ax.set_title('Stand: %i (%i); Dig: %i [%i]' % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer, antennas[i].combined_status))
                ax.set_xlabel('Frequency [MHz]')
                ax.set_ylabel('P.S.D. [dB/RBW]')
                ax.set_xlim([10,90])
                ax.set_ylim([10,80])
                
            # Save spectra image if requested
            if config['output'] is not None:
                base, ext = os.path.splitext(config['output'])
                outFigure = "%s-%02i%s" % (base, f+1, ext)
                fig.savefig(outFigure)
                
        plt.draw()
        
    print "RBW: %.1f Hz" % (freq[1]-freq[0])
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])

