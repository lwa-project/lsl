General - How do I...
=====================
.. highlight:: python
   :linenothreshold: 2

Get a List of Stands for a Particular Date
------------------------------------------
To obtain an ordered numpy array of stands connected to the DP system at a particular station::

	from lsl.common import stations
	lwa1 = stations.lwa1
	stands = lwa1.getStands()

This retrieves a :class:`lsl.common.stations.LWAStation` instance for LWA-1 built using the SSMIF 
file included with the LSL release.  To get the list of stands corresponding to a different SSMIF
file::

	from lsl.common import stations
	lwa = stations.parseSSMIF('your_ssmif.txt')
	stands = lwa1.getStands()

Although the stand list is useful, it is more useful to have a list of :class:`lsl.common.stations.Antenna`
instances.  This provide much more information about the station setup (stand, FEE, position, etc.) 
than the simple stands list.  To get this list::

	from lsl.common import stations
	lwa1 = stations.lwa1
	antennas = lwa1.getAntennas()

It should be noted that this function returned the antennas in order of DP digitizer number so 
that it is easier to convert digitizer numbers to something useful.

Get the Coordinates of a Particular Stand
-----------------------------------------
Once a list of :class:`lsl.common.stations.Antenna` instances has been retrieved, it is easy to get 
stand positions::
	
	xyz = numpy.zeros((3,))
	xyz[0] = antennas[0].stand.x
	xyz[1] = antennas[0].stand.y
	xyz[2] = antennas[0].stand.z

Where x refers to the east-west difference of the antenna from the station's center post, y refers to 
the north-south distance from the station's center post, and z to the height above the center post.  
All of the coordinates are in meters.

Get the Cable Delays and Gains for a Particular Stand
-----------------------------------------------------
The list of :class:`lsl.common.stations.Antenna` also makes it easy to find out the cable delay for a 
particular antennas::

	dly = antennas[0].cable.getDelay(49e6)

By default, delays are computed in seconds and the input frequencies are in Hz.  The delays can be
automatically converted to nanoseconds by setting the `ns` keyword to True.

Similarly, the cable loss for a particular antenna can be found::

	los = antennas[0].cable.getGain(49e6)

Compute the *uvw* Coordinates of a Baseline
-------------------------------------------
The :func:`lsl.correlator.uvUtils.computeUVW` function allows for *uvw* coordinates to be computed for a all baselines
formed from a collection of stands for a particular hour angle, declination, and frequency::

	from lsl.common import stations
	from lsl.correlator import uvUtils
	lwa1 = stations.lwa1
	antennas = lwa1.getAntennas()
	uvw = uvUtils.computeUVW(antennas, HA=1.0, dec=65.0, freq=50e6, IncludeAuto=True)

The above code block computes the *uvw* coordinates at 50 MHz for all stands used on JD 2,455,548.38 for an object at
an hour angle of 1 hour and a declination of +65 degrees.  The returned list of *uvw* coordinates will also include 
entries for autocorrelations.  The order of the baselines is given by the function :func:`lsl.correlator.uvUtils.getBaselines`.

Retrieve Earth Orientation Parameters
-------------------------------------
The :mod:`lsl.misc.geodesy` modules includes functions to retrieve earth orientation parameters (EOP; x, y, and UT1-UTC) for
a given modified Julian Date (MJD) or MJD range.  To retrieve the parameters as a :class:`lsl.misc.geodesy.EOP` object, use::

	from lsl import astro
	from lsl.misc import geodesy
	jd = 2455548.38787
	mjd = jd - astro.MJD_OFFSET
	eop = geodesy.getEOP(mjd)
	
For multiple MJDs, use :func:`lsl.misc.geodesy.getEOPRange` to return a list of EOPs, one for each day::
	
	from lsl import astro
	from lsl.misc import geodesy
	jd1 = 2455548.38787
	jd2 = 2455552.38787
	mjd1 = jd1 - astro.MJD_OFFSET
	mjd2 = jd2 - astro.MJD_OFFSET
	eops = geodesy.getEOPRange(mjd1, mdj2)
	for eop in eops:
		print str(eop)



