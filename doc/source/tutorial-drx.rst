DRX - How do I...
=================
.. highlight:: python
   :linenothreshold: 2

Use the Included Scripts
------------------------
The scripts included with LSL are primarily intended to be used as a guide to writing new scripts
to accomplished particular tasks.  However, some of the scripts are genuinely useful.

Generate Quick Spectra with drxSpectra.py
+++++++++++++++++++++++++++++++++++++++++
drxSpectra.py provides a way to plot integrated spectra for *all* data in a DRX file.  To use drxSpectra.py::

	python drxSpectra.py -l 1024 -o spectra.png DRX_file.dat

The above command creates a 1,024 channel spectra of the given DRX file, displays the spectra to the screen, and
then saves the plots to the file ``spectra.png``.  

For files that contain more data than can fit in the machine's memory at once, drxSpectra.py 'chunks' the data into
10,000 frame sections.  These chunks are averaged together to generate the global average of the file.  In addition, 
the relative difference between the average within each of the chunks and the global average is show in the spectra 
plot as variations around 0 dB.

Read in Data
------------
Here is a Python code snippet for reading in DRX data::

	>>> from lsl.reader import drx, errors
	>>> fh = open('DRX_file.dat', 'rb')
	>>> frame = drx.readFrame(fh)
	>>> print frame.parseID()
	(1, 1, 0)
	>>> print frame.data.iq.mean()
	0.03+1.03j

In the above code, line 3 reads the raw DRX frame into a :class:`lsl.reader.drx.Frame` object.  Lines 4 and 6 access the Frame objects various attributes.  Line 4, for example, parses the DRX ID field and returns a three-element tuple containing the beam number, tuning, and polarization.  Line 6 prints the mean value of the I/Q data associated with this frame.

.. warning::
	The DRX reader can throw various errors when reading in a DRX frame if the frame
	is truncated or corrupted.  These errors are defined in :mod:`lsl.reader.errors`.
	drx.readFrame should be wrapped inside a try...except to check for these.


Plot Spectra
------------
After the DRX data have been read in, spectra can by computed and plotted using the function
:func:`lsl.correlator.fx.calcSpectra`.  For example::

	>>> from lsl.correlator import fx as fxc
	>>> freq, spec = fxc.calcSpectra(data, LFFT=2048, SampleRate=19.6e6, CentralFreq=38e6, DisablePool=True)

Where data is a 2-D array of where the first dimension loops through stands  and the second samples.  Unlike TBW data,
the additional keywords 'SampleRate' and 'CentralFreq' are needed to create the correct frequencies associated with
the FFTs.  The sample rate can be obtained from the data using::

	>>> frame = drx.readFrame(fh)
	>>> sampleRate = frame.getSampleRate()

Currently there is not a way to determine the central frequency of the observations from the data frames.

LSL 0.4.0 introduces a new way to compute spectra with the :func:`lsl.correlator.fx.SpecMaster`
function.  This function uses a C extension and OpenMP to provide better overall performance.  SpecMaster
is called in the same way as the original calcSpectra function::

	>>> freq, spec = fxc.SpecMaster(data, LFFT=2048, SampleRate=1e5, CentralFreq=38e6)

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

