General Utilities
=================
astroevents.py
  :Description: Application to display rise, transit, and set times for various astronomical sources from LWA-1 for the current date.

  :Usage: astroevents.py

  :Options: None

astrostatus.py
  :Description: Application to calculate real-time ephemeris for a LWA site.

  :Usage: astrostatus.py [options]

  :Options: -h, --help            show this help message and exit
          
            -s SITE, --site=SITE  site name (default LWA-1)

driftcurve.py
  :Description: Generate a drift curve for a dipole at LWA-1 observing at a given frequency in MHz.

  :Usage: driftcurve.py [OPTIONS]

  :Options: -h, --help             Display this help information

            -f, --freq             Frequency of the observations in MHz (default = 74 MHz)

            -p, --polarization     Polarization of the observations (NS or EW; default = EW)

            -l, --lf-map           Use LF map instead of GSM

            -t, --time-step        Time step of simulations in minutes (default = 10)

            -x, --do-plot          Plot the driftcurve data

            -v, --verbose          Run driftcurve in vebose mode

lwa_cat_view.py
  :Description: Simple LWDA astronomical source catalogue display application.

  :Usage: lwa_cat_view.py [options]

  :Options: -h, --help            show this help message and exit

            -s SITE, --site=SITE  site name (default LWDA)

            -p PERIOD, --period=PERIOD
                        update period in seconds (default 5)

plotAntenna.py
  :Description: Plot the relative dipole response for both polarizations of an isolated LWA-1 antenna at a particular frequency.

  :Usage: plotAntenna.py [OPTIONS]

  :Options: -h, --help             Display this help information
            -f, --freq             Frequency of the observations in MHz (default = 74 MHz)
            -v, --verbose          Run plotAntenna in vebose mode

plotStands.py
  :Description: Plot the x, y, and z locations of stands at LWA-1.  Also, mark and label particular stands, if requested.

  :Usage: plotStands.py [OPTIONS] [stand1 [stand2 [...]]]

  :Options: -h, --help             Display this help information
            -l, --label            Label the stands with their ID numbers (default = No)
            -v, --verbose          Run plotStands in vebose mode

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

            -d, --date             	Label the split files with a date rather than a sequence number

  .. note::
	This script does not use a :mod:`lsl.reader.buffer` buffer to try to re-order or verify all
	packets and simply splits files based on size.

splitTBN.py
  :Description: Split a TBN file containing multiple seconds into several files

  :Usage: splitTBN.py [options] file

  :Options: -h, --help             	Display this help information

            -c, --count            	Number of seconds to keep

            -o, --offset           	Number of seconds to skip before splitting

            -d, --date             	Label the split files with a date rather than a squence number

  .. note::
	This script does not use a :mod:`lsl.reader.buffer` buffer to try to re-order or verify all
	packets and simply splits files based on size.

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

  :Options: -h, --help                  Display this help information

            -e, --enable-model          Use the CFTOOL bandpass model if it is present in the current directory

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

            -m, --metadata              Name of SSMIF file to use for mappings

            -t, --bartlett              Apply a Bartlett window to the data

            -b, --blackman              Apply a Blackman window to the data

            -n, --hanning               Apply a Hanning window to the data

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -g, --gain-correct          Correct signals for the cable losses

            -s, --stack                 Stack spectra in groups of 6 (if '-g' is enabled only)

            -d, --disable-chunks        Display plotting chunks in addition to the global average

            -o, --output                Output file name for spectra imag

  .. warning::
	tbwSpectra.py currently assumed that the system it is running on has enough memory to read in
	a full TBW capture.  Due to data representation and processing overheads this amounts to about
	16 GB.

tbnSpectra.py
  :Description: Given a TBN file, plot the time averaged spectra for each digitizer input.

  :Usage: tbnSpectra.py [OPTIONS] file

  :Options: -h, --help                  Display this help information

            -m, --metadata              Name of SSMIF file to use for mappings

            -t, --bartlett              Apply a Bartlett window to the data

            -b, --blackman              Apply a Blackman window to the data

            -n, --hanning               Apply a Hanning window to the data

            -s, --skip                  Skip the specified number of seconds at the beginning of the file (default = 0)

            -a, --average               Number of seconds of data to average for spectra (default = 10)

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -d, --disable-chunks        Display plotting chunks in addition to the global average

            -o, --output                Output file name for spectra image

drxTimeseries.py
  :Description: Read in DRX files and create a collection of timeseries (I/Q) plots.

  :Usage: drxTimeseries.py [OPTIONS] file

  :Options: -h, --help                  Display this help information

            -s, --skip                  Skip the specified number of seconds at the beginning of the file (default = 0)

            -p, --plot-range            Number of seconds of data to show in the I/Q plots (default = 0.0001)

            -q, --quiet                 Run drxSpectra in silent mode

            -o, --output                Output file name for time series image

drxSpectra.py
  :Description: Given a DRX file, plot the time averaged spectra for each beam output.

  :Usage: drxSpectra.py [OPTIONS] file

  :Options: -h, --help                  Display this help information

            -t, --bartlett              Apply a Bartlett window to the data

            -b, --blackman              Apply a Blackman window to the data

            -n, --hanning               Apply a Hanning window to the data

            -s, --skip                  Skip the specified number of seconds at the beginning of the file (default = 0)

            -a, --average               Number of seconds of data to average for spectra (default = 10)

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -d, --disable-chunks        Display plotting chunks in addition to the global average

            -o, --output                Output file name for spectra image

correlateTBW.py
  :Description: Cross-correlate data in a TBW file

  :Usage: correlateTBW.py [OPTIONS] file

  :Options: -h, --help             Display this help information

            -f, --fft-length       Set FFT length (default = 512)

            -s, --samples          Number of average visibilities to generate (default = 10)

            -q, --quiet            Run correlateTBW in silent mode

correlateTBN.py
  :Description: Example script that reads in TBN data and runs a cross-correlation on it.
                The results are saved in the Miriad UV format.

  :Usage: correlateTBN.py [OPTIONS] file
  
  :Options: -h, --help             Display this help information

            -c, --central-freq     Central frequency of the observations in MHz

            -f, --fft-length       Set FFT length (default = 512)

            -t, --avg-time         Window to average visibilities in time (seconds; default = 6 s)

            -s, --samples          Number of average visibilities to generate (default = 10)

            -q, --quiet            Run correlateTBN in silent mode

possm.py
  :Description:  Script that takes a FITS IDI file and mimics the AIPS task POSSM by plotting
                 average cross-power spectra for all baselines in the FITS IDI file.

  :Usage: possm.py file

  :Options: None


Post-IOC Compatible Scripts
============================
The scripts listed below are compatible with post-IOC observing conditions and use the meta-data 
tarballs provided with each observation to set all relevant parameters (frequency, stand mappings, 
etc.)  These scripts can be found in the scripts/IOC directory included with LSL 0.5.

.. seealso::
	`MCS0030 <http://www.ece.vt.edu/swe/lwavt/doc/MCS0030v5_Observing.pdf>`_ for a description 
	of the post-IOC observing procedure and the contents of the meta-data tarball.

plotStands.py
  :Description: Plot the x, y, and z locations of stands at LWA-1.  Also, mark and label particular stands, if requested.

  :Usage: plotStands.py [OPTIONS] metaData [stand1 [stand2 [...]]]

  :Options: -h, --help             Display this help information
            -l, --label            Label the stands with their ID numbers (default = No)
            -v, --verbose          Run plotStands in vebose mode

plotUVCoverage.py
  :Description: Randomly select 20 antennae from LWA-1 and plot the uv-plane coverage for
                a zenith snapshot and the expected beam.  Alternatively, select some 
                FRACTION of the stands with installed FEEs to use or use the specified
                list of stands.

  :Usage: plotUVCoverage.py [OPTIONS] metaData

  :Options: -h, --help             Display this help information

            -f, --frequency        Frequency in MHz to compute the uv coverage (default 50 MHz)


readTBN_buffered.py
  :Description: Example script for reading in TBW data and writing it to a TSFITS file.
                This version differs from the regular readTBN script in that it uses a frame
                buffer to reorder out-of-order packets and dropped frames.

  :Usage: readTBN_buffered.py metaData data

  :Options: None

tbwSpectra.py
  :Description: Given a TBW file, plot the time averaged spectra for each digitizer input.

  :Usage: tbwSpectra.py [OPTIONS] metaData data

  :Options: -h, --help                  Display this help information

            -t, --bartlett              Apply a Bartlett window to the data

            -b, --blackman              Apply a Blackman window to the data

            -n, --hanning               Apply a Hanning window to the data

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -g, --gain-correct          Correct signals for the cable losses

            -s, --stack                 Stack spectra in groups of 6 (if '-g' is enabled only)

            -d, --disable-chunks        Display plotting chunks in addition to the global average

            -o, --output                Output file name for spectra imag

  .. warning::
	tbwSpectra.py currently assumed that the system it is running on has enough memory to read in
	a full TBW capture.  Due to data representation and processing overheads this amounts to about
	16 GB.

tbnSpectra.py
  :Description: Given a TBN file, plot the time averaged spectra for each digitizer input.

  :Usage: tbnSpectra.py [OPTIONS] metaData data

  :Options: -h, --help                  Display this help information

            -t, --bartlett              Apply a Bartlett window to the data

            -b, --blackman              Apply a Blackman window to the data

            -n, --hanning               Apply a Hanning window to the data

            -s, --skip                  Skip the specified number of seconds at the beginning of the file (default = 0)

            -a, --average               Number of seconds of data to average for spectra (default = 10)

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -d, --disable-chunks        Display plotting chunks in addition to the global average

            -o, --output                Output file name for spectra image

drxSpectra.py
  :Description: Given a DRX file, plot the time averaged spectra for each beam output.

  :Usage: drxSpectra.py [OPTIONS] metaData data

  :Options: -h, --help                  Display this help information

            -t, --bartlett              Apply a Bartlett window to the data

            -b, --blackman              Apply a Blackman window to the data

            -n, --hanning               Apply a Hanning window to the data

            -s, --skip                  Skip the specified number of seconds at the beginning of the file (default = 0)

            -a, --average               Number of seconds of data to average for spectra (default = 10)

            -q, --quiet                 Run tbwSpectra in silent mode

            -l, --fft-length            Set FFT length (default = 4096)

            -d, --disable-chunks        Display plotting chunks in addition to the global average

            -o, --output                Output file name for spectra image

correlateTBW.py
  :Description: Cross-correlate data in a TBW file

  :Usage: correlateTBW.py [OPTIONS] metaData data

  :Options: -h, --help             Display this help information

            -f, --fft-length       Set FFT length (default = 512)

            -s, --samples          Number of average visibilities to generate (default = 10)

            -q, --quiet            Run correlateTBW in silent mode

correlateTBN.py
  :Description: Example script that reads in TBN data and runs a cross-correlation on it.
                The results are saved in the Miriad UV format.

  :Usage: correlateTBN.py [OPTIONS] metaData data
  
  :Options: -h, --help             Display this help information

            -f, --fft-length       Set FFT length (default = 512)

            -t, --avg-time         Window to average visibilities in time (seconds; default = 6 s)

            -s, --samples          Number of average visibilities to generate (default = 10)

            -q, --quiet            Run correlateTBN in silent mode