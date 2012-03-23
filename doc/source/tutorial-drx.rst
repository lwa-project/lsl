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
:func:`lsl.correlator.fx.SpecMaster`.  For example::

	>>> from lsl.correlator import fx as fxc
	>>> freq, spec = fxc.SpecMaster(data, LFFT=2048, SampleRate=19.6e6, CentralFreq=38e6)

Where data is a 2-D array of where the first dimension loops through stands  and the second samples.  Unlike TBW data,
the additional keywords 'SampleRate' and 'CentralFreq' are needed to create the correct frequencies associated with
the FFTs.  The sample rate can be obtained from the data using::

	>>> frame = drx.readFrame(fh)
	>>> sampleRate = frame.getSampleRate()

The central frequency of the observations from the data frames with::

	>>> frame = drx.readFrame(fh)
	>>> cFreq = frame.getCentralFreq()

.. note::
	Each of the DRX tunings may be at different frequecies so more than one frame will need to be inspected.  The following
	loop::

		>>> cFreq1, cFreq2 = 0, 0
		>>> for i in xrange(4):
		...      frame = drx.readFrame(fh)
		...      beam,tune,pol = frame.parseID()
		...      if tune == 1 and pol == 0:
		...           cFreq1 = frame.getCentralFreq()
		...      elif tune == 2 and pol == 1:
		...           cFreq1 = frame.getCentralFreq()
		...      else:
		...           pass
		...

	will determine both frequencies.

Once the spectra have been computed, they can be plotted via *matplotlib* via::

	>>> import numpy
	>>> from matplotlib import pyplot as plt
	>>> fig = plt.figure()
	>>> ax = fig.gca()
	>>> ax.plot(freq/1e6, numpy.log10(spec[0,:])*10.0)
	>>> ax.set_xlabel('Frequency [MHz]')
	>>> ax.set_ylabel('PSD [Arb. dB]')

Computing Stokes Parameters
---------------------------
The :func:`lsl.correlator.fx.SpecMaster` computes only linear polarization combination, e.g., XX and YY, for the data.  To compute
the Stokes parameters for the data use the :func:`lsl.correlator.fx.StokesMaster` function.  To generate the Stokes parameters for 
a dataset::

	>>> antennas = []
	>>> for i in xrange(4):
	...      if i / 2 == 0:
	...           newAnt = stations.Antenna(1)
	...      else:
	...           newAnt = stations.Antenna(2)
	...      if i % 2 == 0:
	...           newAnt.pol = 0
	...      else:
	...           newAnt.pol = 0
	...      antennas.append(newAnt)
	...
	>>> freq, spec = fxc.StokesMaster(data, antennas, LFFT=2048, SampleRate=19.6e6, CentralFreq=38e6)

This function differs from :func:`lsl.correlator.fx.SpecMaster` in that it requires a list of :class:`lsl.common.stations.Antenna` instances to match the polarization data for the two tunings.  This is accomplished on lines 1 through 11 where a list of psuedo-antennas is created.  The output of array
spec is three dimensional with the Stokes parameters on the first access.  The parameter order is I, Q, U, and V.  

Once the spectra have been computed, they can be plotted with *mathplotlib*::

	>>> import numpy
	>>> from matplotlib import pyplot as plt
	>>> fig = plt.figure()
	>>> ax = fig.gca()
	>>> for i,p in enumerate(('I', 'Q', 'U', 'V')):
	...      ax.plot(freq/1e6, numpy.log10(spec[i,0,:])*10.0, label='Stokes %i' % p)
	...
	>>> ax.set_xlabel('Frequency [MHz]')
	>>> ax.set_ylabel('PSD [Arb. dB]')

The loop in lines 5 and 6 plots all four parameters for tuning 1.
