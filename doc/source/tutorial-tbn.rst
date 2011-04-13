TBN - How do I...
=================
.. highlight:: python
   :linenothreshold: 2

Use the Included Scripts
------------------------
The scripts included with LSL are primarily intended to be used as a guide to writing new scripts
to accomplished particular tasks.  However, some of the scripts are genuinely useful.

Paring Down TBN Files with splitTBN.py
++++++++++++++++++++++++++++++++++++++
At 200 kB per stand per polarization per second, it is easy to end up with large TBN files.  splitTBN.py
makes it easy to remove only a small section of a particular file for analysis.  The advantage of this is
that is makes it easier to perform a *quick look* using scripts like tbnSpectra.py.  All of the options to
splitTBN.py are safe and have no adverse side-effects.  To split off the first hour of a large TBN file::

	python splitTBN.py -c 3600 TBN_file.dat

For the second hour, use::
	
	python splitTBN.py -o 3600 -c 3600 TBN_file.dat

These commands create files with names of the form ``TBN_file_s<first frame number>.dat``.  A more convenient
way to name them is using the '-d' flag with appends the date and time of the first frame in the newly 
created file.

splitTBN.py is also useful for aligning TBN files with full frames.  In some situations, the first frame in a
TBN file does not correspond to the first stand.  splitTBN.py examines the start of every file and removes the 
first capture if it is missing some of the frames.

Generate Quick Spectra with tbnSpectra.py
+++++++++++++++++++++++++++++++++++++++++
tbnSpectra.py provides a way to plot integrated spectra for *all* data in a TBN file.  For short TBN captures
(up to a few minutes) it can be directly run on the file.  However, longer TBN files should be pared down with
splitTBN.py to ensure that tbnSpectra.py does not take hours to run.  To use tbnSpectra.py::

	python tbnSpectra.py -l 1024 -o spectra.png TBN_file.dat

The above command creates a 1,024 channel spectra of the given TBN file, displays the spectra to the screen, and
then saves the plots to the file ``spectra.png``.  

For files that contain more data than can fit in the machine's memory at once, tbnSpectra.py 'chunks' the data into
400,000 frame sections.  These chunks are averaged together to generate the global average of the file.  In addition, 
the relative difference between the average within each of the chunks and the global average is show in the spectra 
plot as variations around 0 dB.

Read in Data
------------
Here is a Python code snippet for reading in TBN data::

	>>> from lsl.reader import tbn, errors
	>>> fh = open('TBN_file.dat', 'rb')
	>>> frame = tbn.readFrame(fh)
	>>> print frame.parseID()
	(1, 0)
	>>> print frame.data.iq.mean()
	5.03+1.03j

In the above code, line 3 reads the raw TBN frame into a :mod:`lsl.reader.tbn.Frame` object.  Lines 4 and 6 access the Frame objects various attributes.  Line 4, for example, parses the TBN ID field and returns a two-element tuple containing the stand number and polarization.  Line 6 prints the mean value of the I/Q data associated with this frame.

.. warning::
	The TBN reader can throw various errors when reading in a TBN frame if the frame
	is truncated or corrupted.  These errors are defined in :mod:`lsl.reader.errors`.
	tbn.readFrame should be wrapped inside a try...except to check for these.


Plot Spectra
------------
After the TBN data have been read in, spectra can by computed and plotted using the function
:mod:`lsl.correlator.fx.calcSpectra`.  For example::

	>>> from lsl.correlator import fx as fxc
	>>> freq, spec = fxc.calcSpectra(data, LFFT=2048, SampleRate=1e5, CentralFreq=38e6, DisablePool=True)

Where data is a 2-D array of where the first dimension loops through stands  and the second samples.  Unlike TBW data,
the additional keywords 'SampleRate' and 'CentralFreq' are needed to create the correct frequencies associated with
the FFTs.  The sample rate can be obtained from the data using::

	>>> sampleRate = tbn.getSampleRate(fh)

which uses the time tags of sequetial frames to determine the sample rate.  Currently there is not a way to determine
the central frequency of the observations from the data frames.

Once the spectra have been computed, they can be plotted via *matplotlib* via::

	>>> import numpy
	>>> from matplotlib import pyplot as plt
	>>> fig = plt.figure()
	>>> ax = fig.gca()
	>>> ax.plot(freq/1e3, numpy.log10(spec[0,:])*10.0)
	>>> ax.set_xlabel('Frequency [kHz]')
	>>> ax.set_ylabel('PSD [Arb. dB]')

.. note::
	In the above example, the thread pool has been disabled for :mod:`lsl.correlator.fx.calcSpectra` which
	forces the function to run single-threaded.  By default, calcSpectra runs with 4 threads and this can
	cause problems if a Ctrl-C is issued.  Ctrl-C kills the main python thread but leaves the worker 
	threads running. 

Post-Acquisition Beam Form
--------------------------
For post-acquisition beam forming, you need need an azimuth (in degrees) and elevation 
(in degrees) to point the beam towards.  For planets, this can be accomplished using the
*pyephem* package that is required by lsl.  For example, compute the location of Jupiter
at LWA-1 on 12/17/2010 at 21:18 UTC (JD 2,455,548.38787)::

	>>> import math
	>>> import ephem
	>>> from lsl.common import stations
	>>> lwa1 = stations.lwa1
	>>> lwaObserver = lwa1.getObserver(2455548.38787, JD=True)
	>>> jove = ephem.Jupiter()
	>>> jove.compute(lwaObserver)
	>>> print "Jupiter:  az -> %.1f, el -> %.1f" % (jove.az*180/math.pi, 
	... jove.alt*180/math.pi)
	Jupiter:  az -> 112.4, el -> 24.4

Line 4 defines the location for LWA-1 as a :mod:`lsl.common.stations.LWAStation` object while line 5 create an ephem.Observer object that can be used to calculate the sky positions of various bodies.  The position of Jupiter is calculated using this Observer object on lines 6 and 7.

.. note::
	When working with positions from *pyephem* objects, all values are in radians.  For more
	information about pyehem, see http://rhodesmill.org/pyephem/

For fixed positions, use::

	>>> cyga = ephem.FixedBody()
	>>> cyga._ra = '19:59:28.30'
	>>> cyga._dec = '+40:44:02'
	>>> cyga.compute(lwaObserver)
	>>> print "Cygnus A:  az -> %.1f, el -> %.1f" % (cyga.az*180/math.pi, 
	... cyga.alt*180/math.pi)
	Cygnus A:  az -> 10.0, el -> 83.2

After TBN data have been read in and a pointing position has been found, a beam can be 
formed.  For example, forming a N-S beam via integer sample delay-and-sum on Cygnus A for 
data taken on JD 2,455,548.38787::

	>>> from lsl.misc import beamformer
	>>> antennas = []
	>>> for ant in lwa1.getAntennas():
	...     if ant.pol == 0:
	...         antennas.append(ant)
	...
	>>> beamdata = beamformer.intDelayAndSum(antennas, data, sampleRate=1e5, 
	... azimuth=10.0, elevation=83.2)

Lines 2 through 5 retrieves the list of antennas used for observations and selects only antennas with N-S polarization.  This information is needed in order to get the correct delays geometric and cable delays to use for the beam forming.

	




