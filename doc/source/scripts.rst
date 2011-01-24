General Utilities
=================
astrostatus.py
  :Description: Application to calculate real-time ephemeris for a LWA site.

  :Usage: astrostatus.py [options]

  :Options: -h, --help            show this help message and exit
          
            -s SITE, --site=SITE  site name (default LWA-1)

driftcurve.py
  :Description: Generate a drift curve for a dipole at SITE observing at a given FREQ (MHz).
                SITE must be one of the sites known by the station module in lwda_util.

  :Usage: driftcurve.py [options] FREQ NECFILENAME SKYMAPFILENAME

  :Options: -h, --help            show this help message and exit

            -v, --verbose         enable debug messages

            -x, --doplot          Make an X-windows plot

            -s SITE, --site=SITE  site name (default LWDA)

            -p POLARIZATION, --polarization=POLARIZATION
                                 antenna polarization orientation (NS or EW)

lwa_cat_view.py
  :Description: Simple LWDA astronomical source catalogue display application.

  :Usage: lwa_cat_view.py [options]

  :Options: -h, --help            show this help message and exit

            -s SITE, --site=SITE  site name (default LWDA)

            -p PERIOD, --period=PERIOD
                        update period in seconds (default 5)

plotAntenna.py
  :Description: Example script to plot the relative response of an isolated LWA antenna
                as a function of azimuth and elevation using an NEC model at a particular
                frequency in MHz.

  :Usage: plotAntenna.py FREQ

  :Options: None

plotStands.py
  :Description: Example script to read in the positions of stands at LWA-1 and make a plot
                of the site.

  :Usage: plotStands.py

  :Options: None

plotUVCoverage.py
  :Description: Randomly select 20 antennae from LWA-1 and plot the uv-plane coverage for
                a zenith snapshot and the expected beam.  Alternatively, select some 
                FRACTION of the stands with installed FEEs to use or use the specified
                list of stands.

  :Usage: plotUVCoverage.py [FRACTION | STAND LIST]

  :Options: None

Data Reading and Writing
========================
splitTBW.py
  :Description: Split a TBW file containing multiple captures into several single capture files.

  :Usage: splitTBW.py [options] file

  :Options: -h, --help             	Display this help information

            -c, --count            	Number of capturs to split off
            -o, --offset           	Number of captures to skip before splitting
            -d, --date             	Label the split files with a date rather than a 
					sequence number

splitTBN.py
  :Description: Split a TBN file containing multiple seconds into several files

  :Usage: splitTBN.py [options] file

  :Options: -h, --help             	Display this help information

            -c, --count            	Number of seconds to keep
            -o, --offset           	Number of seconds to skip before splitting
            -d, --date             	Label the split files with a date rather than a 
					sequence number

readTBW.py
  :Description: Example script for reading in TBW data and writing it to a TSFITS file.

  :Usage: readTBW.py file

  :Options: None

readTBN.py
  :Description: Example script for reading in TBN data and writing it to a TSFITS file.

  :Usage: readTBN.py file

  :Options: None

readTBN_buffered.py
  :Description: Example script for reading in TBW data and writing it to a TSFITS file.
                This version differs from the regular readTBN script in that it uses a frame
                buffer to reorder out-of-order packets and dropped frames.

  :Usage: readTBN_buffered.py file

  :Options: None

readDRX.py
  :Description: Example script for reading in DRX data and writing it to a SD-FITS file.

  :Usage: readDRX.py file

  :Options: None

readS60.py
  :Description: Python script to read in a S60 file and average it in time.  The output is a
                npz file of the time-averaged spectra and a PNG of the bandpass/waterfall diagram.

  :Usage: readS60.py [OPTIONS] file

  :Options: -h, --help                 Display this help information

            -e, --enable-model          Use the CFTOOL bandpass model if
                                       it is present in the current
                                       directory

            -q, --quiet                 Run readS60 in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -t, --avg-time              Window to average spectra in time

plotMapper.py
  :Description: Read and plot the NOSTA_MAPPER table in a FITS IDI file writen by
                :mod:`lsl.writer.fitsidi` if it exists.

  :Usage: plotMapper.py file

  :Options: None

Data Analysis
=============
tbwSpectra.py
  :Description: Given a TBW file, plot the time averaged spectra for each digitizer input.

  :Usage: tbwSpectra.py [OPTIONS] file

  :Options: -h, --help                  Display this help information

            -b, --blackman              Apply a Blackman window to the data

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -o, --output                Output file name for spectra imag

tbnSpectra.py
  :Description: Given a TBN file, plot the time averaged spectra for each digitizer input.

  :Usage: tbnSpectra.py [OPTIONS] file

  :Options: -h, --help                  Display this help information

            -b, --blackman              Apply a Blackman window to the data

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -o, --output                Output file name for spectra image

correlateTBN.py
  :Description: Example script that reads in TBN data and runs a cross-correlation on it.
                The results are saved in the Miriad UV format.

  :Usage: correlateTBN.py [OPTIONS] file
  
  :Options: -h, --help             Display this help information

            -c, --central-freq     Central frequency of the observations in MHz

            -f, --fft-length       Set FFT length (default = 512)

            -t, --avg-time         Window to average visibilities in time (seconds;
                                   default = 6 s)

            -s, --samples          Number of average visibilities to generate
                                   (default = 10)

            -q, --quiet            Run correlateTBN in silent mode

possm.py
  :Description:  Script that takes a FITS IDI file and mimics the AIPS task POSSM by plotting
                 average cross-power spectra for all baselines in the FITS IDI file.

  :Usage: possm.py file

  :Options: None
