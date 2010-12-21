General - How do I...
=====================
.. highlight:: python
   :linenothreshold: 2

Get a List of Stands for a Particular Date
------------------------------------------
To obtain an ordered numpy array of stands connected to the DP system at a particular station
for a particular date::

	from lsl.common import stations
	lwa1 = stations.lwa1()
	stands = lwa1.getStands(2455548.38787, JD=True)

This retrieves a :mod:`lsl.common.stations.LWAStation` object for LWA-1 and queries the obejct for
a list of stands for JD 2,455,548.38.  Alternatively, the date used can be set as a string, i.e.,::

	stands = lwa1.getStands('2010/12/17 21:18:00')

This fetches the list of stands for UT December 17, 2010 at 21:18.

Get the Coordinates of a Particular Stand
-----------------------------------------
The :mod:`lsl.correlator.uvUtils` module contains functions to return the location of a particular stand.
Stands are numbered from 1 to 258 with stand #257 being the isolated stand in the south west corner of LWA-1
and stand #258 the outlier (RTA).  Two method exist to get stand position, :mod:`lsl.correlator.uvUtils.getXYZ`
and :mod:`lsl.correlator.uvUtils.PositionCache`.  The difference between these two methods is that the former 
reads in the list of stand coordinates each time it is called whereas the latter reads in the locations only 
once and caches the results.  To use either method::

	import numpy
	from lsl.correlator import uvUtils
	xyz = uvUtils.getXYZ( numpy.array([1, 2, 3]) )
	pCache = uvUtils.PositionCache()
	xyz = pCache.getXYZ( numpy.array([1, 2, 3]) )

Both lines 3 and 5 retrive the x, y, and z coordinates (in meters) of stands 1, 2, and 3.  

Get the Cable Delays for a Particular Stand
-------------------------------------------
Similar to the coordinates, the :mod:`lsl.correlator.uvUtils` module has functions to retrieve the the cable
delays as a function of frequency for a particular stand over a particular range of frequencies.  In addition,
there also exists both stand-alone and cached versions of the functions.  To calculate cable delays in seconds::

	import numpy
	from lsl.correlator import uvUtils
	dly = uvUtils.cableDelay(1, freq=numpy.array([50e6, 55e6, 60e6]))

which returns a numpy array of delays for stand #1 at teh frequencies 50, 55, and 60 MHz.  The cached version is
implemented differently so that the frequencies used are also stored in the object::

	import numpy
	from lsl.correlator import uvUtils
	dCache = uvUtils.CableCache(numpy.array([50e6, 55e6, 60e6]))
	dly = dCache.cableDelay(1)
	dCache.updateFreq(numpy.array([40e6, 45e6]))
	dly = dCache.cableDelay(1)

Line 3 creates the cache for frequencies of 50, 55, and 60 MHz and line 4 returns the delays at those frequencies
for stand #1.  Unless the frequencies stored in the cache are updated (i.e., line 5), all other calls to the cache
use the initial frequency set.

Compute the *uvw* Coordinates of a Baseline
-------------------------------------------
The :mod:`lsl.correlator.uvUtils.computeUVW` function allows for *uvw* coordinates to be computed for a all baselines
formed from a collection of stands for a particular hour angle, declination, and frequency::

	from lsl.common import stations
	from lsl.correlator import uvUtils
	lwa1 = stations.lwa1()
	stands = lwa1.getStands(2455548.38787, JD=True)
	uvw = uvUtils.computeUVW(stands, HA=1.0, dec=65.0, freq=50e6, IncludeAuto=True)

The above code block computes the *uvw* coordinates at 50 MHz for all stands used on JD 2,455,548.38 for an object at
an hour angle of 1 hour and a declination of +65 degrees.  The returned list of *uvw* coordinates will also include 
entries for autocorrelations.  The order of the baselines is given by the function :mod:`lsl.correlator.uvUtils.getBaselines`.

Retrieve Earth Orientation Parameters
-------------------------------------
The :mod:`lsl.misc.geodesy` modules includes functions to retrieve earth orientation parameters (EOP; x, y, and UT1-UTC) for
a given modified Julian Date (MJD) or MJD range.  To retrieve the parameters as a :mod:`lsl.misc.geodesy.EOP` object, use::

	from lsl import astro
	from lsl.misc import geodesy
	jd = 2455548.38787
	mjd = jd - astro.MJD_OFFSET
	eop = geodesy.getEOP(mjd)
	
For multiple MJDs, use :mod:`lsl.misc.geodesy.getEOPRange` to return a list of EOPs, one for each day::
	
	from lsl import astro
	from lsl.misc import geodesy
	jd1 = 2455548.38787
	jd2 = 2455552.38787
	mjd1 = jd1 - astro.MJD_OFFSET
	mjd2 = jd2 - astro.MJD_OFFSET
	eops = geodesy.getEOPRange(mjd1, mdj2)
	for eop in eops:
		print str(eop)



