TBW - How do I...
=================
.. highlight:: python
   :linenothreshold: 2

Use the Included Scripts
------------------------
The scripts included with LSL are primarily intended to be used as a guide to writing new scripts
to accomplished particular tasks.  However, some of the scripts are genuinely useful. 

Paring Down TBW Files with splitTBW.py
++++++++++++++++++++++++++++++++++++++
At 350 MB per 61 ms capture, it is easy to end up with large TBW files.  splitTBW.py makes it easy to split
a series of sequential TBW captures into appropriately named single captures.  To split a file into all of its
constituent captures, use::

	python splitTBW.py TBW_file.dat

For the first 10 captures, use::
	
	python splitTBW.py -c 10 TBW_file.dat

These commands create files with names of the form ``TBW_file_s<first frame number>.dat``.  A more convenient
way to name them is using the '-d' flag with appends the date and time of the first frame in the newly 
created file.  In addition, splitTBW.py accepts a '-o' flag which sets how many captures from the start of the
file should be skipped before splitting.

splitTBW.py is also useful for aligning TBW files with full frames.  In some situations, the first frame in a
TBW file does not correspond to the first stand.  splitTBW.py examines the start of every file and removes the 
first capture if it is missing some of the frames.

Quick Spectra with tbwSpectra.py
+++++++++++++++++++++++++++++++++
tbwSpectra.py provides a way to plot integrated spectra for *all* captures in a TBW file.  For files that contain
only a few TBW captures the script can be directly run on the file.  However, longer TBW files should be pared 
down with splitTBW.py to ensure that tbwSpectra.py does not take hours to run.  To use tbwSpectra.py::

	python tbwSpectra.py -l 1024 -o spectra.png TBW_file.dat

The above command creates a 1,024 channel spectra of the given TBW file, displays the spectra to the screen, and
then saves the plots to the file ``spectra.png``.  

For files that contain more data than can fit in the machine's memory at once, tbwSpectra.py 'chunks' the data into
300,000 frame sections (one complete capture for 10 dipole inputs).  These chunks are averaged together to generate 
the global average of the file.  In addition, the relative difference between the average within each of the chunks 
and the global average is show in the spectra plot as variations around 0 dB.

Read in Data
------------
Here is a Python code snippet for reading in TBW data::

	>>> from lsl.reader import tbw, errors
	>>> fh = open('TBW_file.dat', 'rb')
	>>> frame = tbw.readFrame(fh)
	>>> print frame.parseID()
	1
	>>> print frame.data.xy.mean()
	5.43

In the above code, line 3 reads the raw TBW frame into a :class:`lsl.reader.tbw.Frame` object.  Lines 4 and 6 access the Frame objects various attributes.  Line 4, for example, parses the TBW ID field and returns the stand number.  Line 6 prints the mean value of the X and Y polarizations.

.. warning::
	The TBW reader can throw various errors when reading in a TBW frame if the frame
	is truncated or corrupted.  These errors are defined in :mod:`lsl.reader.errors`.
	tbw.readFrame should be wrapped inside a try...except to check for these.


Plot Spectra
------------
After the TBW data have been read in, spectra can by computed and plotted using the function
:func:`lsl.correlator.fx.calcSpectra`.  For example::

	>>> from lsl.correlator import fx as fxc
	>>> freq, spec = fxc.calcSpectra(data, LFFT=2048, DisablePool=True)

LSL 0.4.0 introduces a new way to compute spectra with the :func:`lsl.correlator.fx.SpecMaster`
function.  This function uses a C extension and OpenMP to provide better overall performance.  SpecMaster
is called in the same way as the original calcSpectra function::

	>>> freq, spec = fxc.SpecMaster(data, LFFT=2048)

Where data is a 2-D array of where the first dimension loops through stands  and the second samples.
Once the spectra have been computed, they can be plotted via *matplotlib* via::

	>>> import numpy
	>>> from matplotlib import pyplot as plt
	>>> fig = plt.figure()
	>>> ax = fig.gca()
	>>> ax.plot(freq/1e6, numpy.log10(spec[0,:])*10.0)
	>>> ax.set_xlabel('Frequency [MHz]')
	>>> ax.set_ylabel('PSD [Arb. dB]')

.. note::
	In the above example, the thread pool has been disabled for :func:`lsl.correlator.fx.calcSpectra` which
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

Line 4 defines the location for LWA-1 as a :class:`lsl.common.stations.LWAStation` object while line 5 create an ephem.Observer object that can be used to calculate the sky positions of various bodies.  The position of Jupiter is calculated using this Observer object on lines 6 and 7.

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
formed.  For example, forming a beam via integer sample delay-and-sum on Cygnus A for 
data taken on JD 2,455,548.38787::

	>>> from lsl.misc import beamformer
	>>> antennas = []
	>>> for ant in lwa1.getAntennas():
	...     if ant.pol == 0:
	...         antennas.append(ant)
	...
	>>> beamdata = beamformer.intDelayAndSum(antennas, data, sampleRate=1e5, 
	... azimuth=10.0, elevation=83.2)

Line 2 retrieves the list of stands used for observations on the given date.  This information is needed in order to get the
correct delays geometric and cable delays to use for the beam forming.
