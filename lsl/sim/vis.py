"""
Module for generating simulated arrays and visibility data.  The chief 
functions of this module are:

build_sim_array
  given a station object, a list of stands, and a list of frequencies, build 
  a AIPY AntennaArray-like object.  This module can also generate AntennaArray 
  objects with positional errors by setting the 'pos_error' keyword to a 
  positive value.

build_sim_data
  given a SimArray and a list of aipy.src sources, build up a collection of 
  visibilities for a given set of Julian dates

scale_data
  given a dictionary of simulated visibilities from build_sim_data, apply 
  antenna-based gains and delays to the visibilities

shift_data
  given a dictionary of simulated visibilities from build_sim_data, shift the uvw 
  coordinates of the visibilities.
  .. note::
    This only changes the uvw values and does not phase-shift the data.

The format of the is descrined in the :mod:`lsl.imaging.data` module.

In addition to simulation functions, this module includes buildGriddedImage
which takes a dictionary of visibilities and returns and aipy.im.ImgW object.

.. versionchanged:: 0.3.0
    This module was formerly called lsl.sim.sim
    
.. versionchanged:: 0.5.0
    Moved buildGriddedImage to the :mod:`lsl.imaging.utils` module.

.. versionchanged:: 1.0.1
    Switched over to a new C-based simulation package
    
.. versionchanged:: 1.0.2
    Added in a function, add_baseline_noise, to help bring noise into the 
    simulations
    
.. versionchanged:: 1.0.3
    Moved the calculateSEFD function into lsl.misc.rfutils
    Changed the meaning of the force_gaussian parameter of the build_sim_array()
    function to be the Gaussian full width at half maximum in degrees
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import aipy
import math
import ephem
import numpy
import warnings
from scipy.interpolate import interp1d
from astropy.constants import c as speedOfLight

from lsl import astro
from lsl.common.paths import DATA as DATA_PATH
from lsl.correlator import uvutils
from lsl.common.stations import lwa1
from lsl.imaging.data import PolarizationDataSet, VisibilityDataSet, VisibilityData
from lsl.sim._simfast import FastVis
from lsl.common.color import colorfy

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.6'
__all__ = ['SOURCES', 'RadioEarthSatellite', 'BeamAlm', 'Antenna', 'AntennaArray', 
           'build_sim_array', 'build_sim_data', 'scale_data', 'shift_data', 'add_baseline_noise']


speedOfLight = speedOfLight.to('m/s').value


#: A dictionary of bright sources in the sky to use for simulations
SOURCES = aipy.src.get_catalog(srcs=['Sun', 'Jupiter', 'cas', 'crab', 'cyg', 'her', 'sgr', 'vir'])


class RadioEarthSatellite(object):
    """
    Implement a aipy.amp.RadioBody-lime simulation object for an Earth-
    orbiting satellite using a two-line element set.
    """
    
    def __init__(self, tle, tfreq, tpower=0.0, tbw=1.0e6, ionref=(0.,0.)):
        """
        Initialize the class using:
          * a three-element list of strings that specify the two-line element
            set of the satellite,
          * a list of frequencies (in GHz) where the satellite transmits, 
          * the transmitter power in W, and
          * the transmission bandwidth in Hz.
        
        .. note::
            For calculating the power received on the ground this class 
            assumes that the signal is emitted isotropically.
        """
        
        # Location
        self.Body = ephem.readtle(*tle)
        self.src_name = self.Body.name
        
        # Transmitter - tuning, power, and bandwidth
        try:
            self.tfreq = list(tfreq)
        except TypeError:
            self.tfreq = [tfreq,]
        self.tpower = tpower
        self.tbw = tbw
        
        # Ionospheric distortion
        self.ionref = list(ionref)
        
        # Shape (unresolved)
        self.srcshape = list((0.,0.,0.))
        
    def __getattr__(self, nm):
        """
        First try to access attribute from this class, but if that fails, 
        try to get it from the underlying PyEphem object.
        """
        
        try:
            return object.__getattr__(self, nm)
        except AttributeError:
            return self.Body.__getattribute__(nm)
            
    def compute(self, observer):
        """
        Update coordinates relative to the provided observer.  Must be
        called at each time step before accessing information.
        """
        
        self.Body.compute(observer)
        self.update_jys(observer.get_afreqs())
        
    def get_crds(self, crdsys, ncrd=3):
        """
        Return the coordinates of this location in the desired coordinate
        system ('eq','top') in the current epoch.  If ncrd=2, angular
        coordinates (ra/dec or az/alt) are returned, and if ncrd=3,
        xyz coordinates are returned.
        """
        
        assert(crdsys in ('eq','top'))
        assert(ncrd in (2,3))
        if crdsys == 'eq':
            if ncrd == 2:
                return (self.ra, self.dec)
            return aipy.coord.radec2eq((self.ra, self.dec))
        else:
            if ncrd == 2:
                return (self.az, self.alt)
            return aipy.coord.azalt2top((self.az, self.alt))
            
    def update_jys(self, afreqs):
        """
        Update fluxes relative to the provided observer.  Must be
        called at each time step before accessing information.
        """
        
        # Setup
        r = self.Body.range				# m
        v = self.Body.range_velocity		# m/s
        self.jys = numpy.zeros_like(afreqs)
        
        # Compute the flux coming from the satellite assuming isotropic emission
        self._jys = self.tpower / (4*numpy.pi*r**2) / self.tbw		# W / m^2 / Hz
        self._jys /= 10**-26								# Jy
        
        try:
            dFreq = afreqs[1]-afreqs[0]
        except IndexError:
            dFreq = afreqs[0,1]-afreqs[0,0]
        for f in self.tfreq:
            ## Apply the Doppler shift
            fPrime = f / (1 + v/speedOfLight)
            
            ## Figure out which frequency bin is within one channel of the
            ## shifted signal.  If it's close, set it to a 
            diff = numpy.abs(fPrime - afreqs)*1e9
            good = numpy.where( diff <= self.tbw/2.0 )
            self.jys[good] = self._jys * min([1.0, dFreq*1e9/self.tbw])
                
    def get_jys(self):
        """
        Return the fluxes vs. freq that should be used for simulation.
        """
        
        return self.jys


class BeamAlm(aipy.amp.BeamAlm):
    """
    AIPY-based representation of a beam model where each pointing has a 
    response defined as a polynomial in frequency, and the spatial 
    distributions of these coefficients decomposed into spherical 
    harmonics.
    
    This differs from the AIPY version in that the response() method 
    accepts two and three-dimensions arrays of topocentric coordinates, 
    similar to what aipy.img.ImgW.get_top() produces, and computes the 
    beam response at all points.
    """
    
    def _response_primitive(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y').
        """
        
        return aipy.amp.BeamAlm.response(self, top)
        
    def response(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y').
        
        .. note::
            This function also accepts two and three-dimensions arrays of 
            topocentric coordinates, similar to what aipy.img.ImgW.get_top() 
            produces, and computes the beam response at all points
        """
        
        test = numpy.array(top)
        x,y,z = top
        
        if len(test.shape) == 1:
            temp = self._response_primitive((x,y,z))
            
        elif len(test.shape) == 2:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                temp[:,i] = numpy.squeeze(self._response_primitive((x[i],y[i],z[i])))
                
        elif len(test.shape) == 3:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                for j in range(temp.shape[2]):
                    temp[:,i,j] = numpy.squeeze(self._response_primitive((x[i,j],y[i,j],z[i,j])))
                    
        else:
            raise ValueError("Cannot compute response for %s" % str(test.shape))
            
        return temp


class Beam2DGaussian(aipy.amp.Beam2DGaussian):
    """
    AIPY-based representation of a 2-D Gaussian beam pattern, with default 
    setting for a flat beam.  The width parameters denote the 'sigma' of
    the Gaussian in radians.
    
    This differs from the AIPY version in that the response() method 
    accepts two and three-dimensions arrays of topocentric coordinates, 
    similar to what aipy.img.ImgW.get_top() produces, and computes the 
    beam response at all points.
    
    .. versionchanged:: 1.0.3
        Clarified what 'xwidth' and 'ywidth' are.
    """
    
    def _response_primitive(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates: (x=E,y=N,z=UP). x,y,z may be arrays of multiple 
        coordinates.  Returns 'x' linear polarization (rotate pi/2 for 'y').
        """
        
        return aipy.amp.Beam2DGaussian.response(self, top)
        
    def response(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y').
        
        .. note::
            This function also accepts two and three-dimensions arrays of 
            topocentric coordinates, similar to what aipy.img.ImgW.get_top() 
            produces, and computes the beam response at all points
        """
        
        test = numpy.array(top)
        x,y,z = top
        
        if len(test.shape) == 1:
            temp = self._response_primitive((x,y,z))
            
        elif len(test.shape) == 2:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                temp[:,i] = numpy.squeeze(self._response_primitive((x[i],y[i],z[i])))
                
        elif len(test.shape) == 3:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                for j in range(temp.shape[2]):
                    temp[:,i,j] = numpy.squeeze(self._response_primitive((x[i,j],y[i,j],z[i,j])))
                    
        else:
            raise ValueError("Cannot compute response for %s" % str(test.shape))
            
        return temp


class BeamPolynomial(aipy.amp.BeamPolynomial):
    """
    AIPY-based representation of a Gaussian beam model whose width varies 
    with azimuth angle and with frequency.
    
    This differs from the AIPY version in that the response() method 
    accepts two and three-dimensions arrays of topocentric coordinates, 
    similar to what aipy.img.ImgW.get_top() produces, and computes the 
    beam response at all points.
    """
    
    def _response_primitive(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y').
        """
        
        return aipy.amp.BeamPolynomial.response(self, top)
        
    def response(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y').
        
        .. note::
            This function also accepts two and three-dimensions arrays of 
            topocentric coordinates, similar to what aipy.img.ImgW.get_top() 
            produces, and computes the beam response at all points
        """
        
        test = numpy.array(top)
        x,y,z = top
        
        if len(test.shape) == 1:
            temp = self._response_primitive((x,y,z))
            
        elif len(test.shape) == 2:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                temp[:,i] = numpy.squeeze(self._response_primitive((x[i],y[i],z[i])))
                
        elif len(test.shape) == 3:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                for j in range(temp.shape[2]):
                    temp[:,i,j] = numpy.squeeze(self._response_primitive((x[i,j],y[i,j],z[i,j])))
                    
        else:
            raise ValueError("Cannot compute response for %s" % str(test.shape))
            
        return temp


class Beam(aipy.amp.Beam):
    """
    AIPY-based representation of a flat (gain=1) antenna beam pattern.
    
    This differs from the AIPY version in that the response() method 
    accepts two and three-dimensions arrays of topocentric coordinates, 
    similar to what aipy.img.ImgW.get_top() produces, and computes the 
    beam response at all points.
    """
    
    def _response_primitive(self, top):
        """
        Return the (unity) beam response as a function of position.
        """
        
        return aipy.amp.Beam.response(self, top)
        
    def response(self, top):
        """
        Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y').
        
        .. note::
            This function also accepts two and three-dimensions arrays of 
            topocentric coordinates, similar to what aipy.img.ImgW.get_top() 
            produces, and computes the beam response at all points
        """
        
        test = numpy.array(top)
        x,y,z = top
        
        if len(test.shape) == 1:
            temp = self._response_primitive((x,y,z))
            
        elif len(test.shape) == 2:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                temp[:,i] = numpy.squeeze(self._response_primitive((x[i],y[i],z[i])))
                
        elif len(test.shape) == 3:
            temp = numpy.zeros((self.afreqs.size,)+test.shape[1:])
            for i in range(temp.shape[1]):
                for j in range(temp.shape[2]):
                    temp[:,i,j] = numpy.squeeze(self._response_primitive((x[i,j],y[i,j],z[i,j])))
                    
        else:
            raise ValueError("Cannot compute response for %s" % str(test.shape))
            
        return temp


class Antenna(aipy.amp.Antenna):
    """
    Modification to the aipy.amp.Antenna class to also store the stand ID 
    number in the Antenna.stand attribute.  This also add a getBeamShape 
    attribute that pulls in the old vis.getBeamShape function.
    """

    def __init__(self, x, y, z, beam, phsoff=[0.,0.], bp_r=numpy.array([1]), bp_i=numpy.array([0]), amp=1, pointing=(0.,numpy.pi/2,0), stand=0, **kwargs):
        """
        New init function that include stand ID number support.  From aipy.amp.Antenna:
          * x,y z = antenna coordinates in equatorial (ns) coordinates
          * beam = Beam object (implements response() function)
          * phsoff = polynomial phase vs. frequency.  Phs term that is linear
                     with freq is often called 'delay'.
          * bp_r = polynomial (in freq) modeling real component of passband
          * bp_i = polynomial (in freq) modeling imaginary component of passband
          * amp = overall multiplicative scaling of gain
          * pointing = antenna pointing (az,alt).  Default is zenith.
        """
        
        aipy.amp.Antenna.__init__(self, x,y,z, beam=beam, phsoff=phsoff,
                                  bp_r=bp_r, bp_i=bp_i, amp=amp,
                                  pointing=pointing)
        self.stand = stand
        
    def bm_response(self, top, pol='x'):
        """
        Return response of beam for specified polarization.
        
        .. note::
            This differs from the AIPY implementation in that the LWA X-pol.
            is oriented N-S, not E-W.
            
        .. note::
            This function also accepts two and three-dimensions arrays of 
            topocentric coordinates, similar to what img.ImgW.get_top() 
            produces, and computes the beam response at all points.
        """
        top = numpy.array(top)
        pol = pol.lower()
        
        def _robust_dot(a, b):
            """
            Dot product that operations on multi-dimensional coordinate sets.
            """
            
            if len(b.shape) == 1:
                temp = numpy.dot(a, b)
                
            elif len(b.shape) == 2:
                temp = numpy.zeros_like(b)
                for i in range(b.shape[1]):
                    temp[:,i] = numpy.dot(a, b[:,i])
                    
            elif len(b.shape) == 3:
                temp = numpy.zeros_like(b)
                for i in range(b.shape[1]):
                    for j in range(top.shape[2]):
                        temp[:,i,j] = numpy.dot(a, b[:,i,j])
                        
            else:
                raise ValueError("Cannot dot a (%s) with b (%s)" % (str(a.shape), str(b.shape)))
            
            return temp
            
        top = {'y':_robust_dot(self.rot_pol_x, top), 
               'x':_robust_dot(self.rot_pol_y, top), 
               'l':_robust_dot(self.rot_pol_x, top), 
               'r':_robust_dot(self.rot_pol_y, top)}[pol]
        x,y,z = top
        
        return self.beam.response((x,y,z))
        
    def get_beam_shape(self, pol='x'):
        """
        Return a 360 by 90 by nFreqs numpy array showing the beam pattern of a
        particular antenna in the array.  The first two dimensions of the output 
        array contain the azimuth (from 0 to 359 degrees in 1 degree steps) and 
        altitude (from 0 to 89 degrees in 1 degree steps).
        """
        
        # Build azimuth and altitude arrays.  Be sure to convert to radians
        az = numpy.zeros((360,90))
        for i in range(360):
            az[i,:] = i*numpy.pi/180.0
        alt = numpy.zeros((360,90))
        for i in range(90):
            alt[:,i] = i*numpy.pi/180.0
            
        # The beam model is computed in terms of topocentric coordinates, so make that
        # conversion right quick using the aipy.coord module.
        xyz = aipy.coord.azalt2top(numpy.concatenate([[az],[alt]]))
        
        # I cannot figure out how to do this all at once, so loop through azimuth/
        # altitude pairs
        resp = numpy.zeros((360,90,len(self.beam.freqs)))
        for i in range(360):
            for j in range(90):
                resp[i,j,:] = numpy.squeeze( self.bm_response(numpy.squeeze(xyz[:,i,j]), pol=pol) )
                
        return resp


class AntennaArray(aipy.amp.AntennaArray):
    """
    Modification to the aipy.ant.AntennaArray class to add a function to 
    retrieve the stands stored in the AntennaArray.ants attribute.  Also add 
    a function to set the array time from a UNIX timestamp.
    
    .. versionchanged:: 1.0.1
        Added an option to set the ASP filter for simulation proposes.  
        This updates the bandpasses used by AIPY to include the antenna
        impedance mis-match and the mean ARX response.
    """
    
    def __str__(self):
        return "AntennaArray at lat: %.3f, lng: %.3f, elev: %.1f m with %i antennas" % (self.lat*180.0/numpy.pi, self.long*180.0/numpy.pi, self.elev, len(self.ants))
        
    def __repr__(self):
        return str(self)
        
    def __copy__(self):
        result = AntennaArray((self.lat, self.lon, self.elev), self.ants)
        result.__dict__.update(self.__dict__)
        return result
        
    def __reduce__(self):
        return (AntennaArray, ((self.lat, self.lon, self.elev), self.ants))
        
    def get_stands(self):
        """
        Return a numpy array listing the stands found in the AntennaArray 
        object.
        """

        stands = []
        for ant in self.ants:
            stands.append(ant.stand)
        
        return numpy.array(stands)

    def set_unixtime(self, timestamp):
        """
        Set the array time using a UNIX timestamp (epoch 1970).
        """
        
        self.set_jultime(astro.unix_to_utcjd(timestamp))
        
    def set_asp_filter(self, filter='split'):
        """
        Update the bandpasses for the antennas to include the effect of 
        the antenna impedance mis-match (IMM) and the mean LWA1 ARX 
        response.
        
        Valid filters are:
          * split
          * full
          * reduced
          * split@3MHz
          * full@3MHz
          * none
        
        None is a special case where both the IMM and ARX response are 
        removed, i.e., the bandpass is unity for all frequencies.
        
        .. versionchanged:: 1.2.1
            Added support for the 'split@3MHz' and 'full@3MHz' filters at
            LWA-SV.
        
        .. versionadded:: 1.0.1
        """
        
        # Get the frequencies we need to estimate the bandpass for
        freqs = self.get_afreqs()*1e9
        if len(freqs.shape) == 2:
            freqs = freqs[0,:]
            
        if filter == 'none' or filter is None:
            # Build up the bandpass - of ones
            resp = numpy.ones(freqs.size)
            
        else:
            # Load in the LWA antennas so we can grab some data.  If we don't know
            # what station we are at, assume LWA1.
            station = getattr(self, '_station', lwa1)
            ants = station.antennas
            
            # Antenna impedance mis-match
            immf, immr = ants[0].response(dB=False)
            immr /= immr.max()
            immIntp = interp1d(immf, immr, kind='cubic', bounds_error=False)
            
            # Mean ARX filter response
            arxf, arxr = ants[0].arx.response(filter, dB=False)
            for i in range(1, len(ants)):
                arxfi, arxri = ants[i].arx.response(filter, dB=False)
                arxr += arxri
            arxr /= arxr.max()
            arxIntp = interp1d(arxf, arxr, kind='cubic', bounds_error=False)
            
            # Build up the bandpass
            resp = numpy.ones(freqs.size)
            resp *= immIntp(freqs)
            resp *= arxIntp(freqs)
            
        # Update the AIPY passbands
        for i in range(len(self.ants)):
            self[i].amp = resp
            self[i].update()
            
    def get_baseline_fast(self, i, j, src='z', map=None):
        """Return the baseline corresponding to i,j in various coordinate 
        projections: src='e' for current equatorial, 'z' for zenith 
        topocentric, 'r' for unrotated equatorial, or a RadioBody for
        projection toward that source - fast."""
        bl = self[j] - self[i]
        
        if type(src) == str:
            if src == 'e':
                return numpy.dot(self._eq2now, bl)
            elif src == 'z':
                return numpy.dot(self._eq2zen, bl)
            elif src == 'r':
                return bl
            else:
                raise ValueError('Unrecognized source:' + src)
        try:
            if src.alt < 0:
                raise RuntimeError('%s below horizon' % src.src_name)
            m = src.map
        except AttributeError:
            if map is None:
                ra,dec = aipy.coord.eq2radec(src)
                m = aipy.coord.eq2top_m(self.sidereal_time() - ra, dec)
            else:
                m = map
        return numpy.dot(m, bl).transpose()
        
    def gen_uvw_fast(self, i, j, src='z', w_only=False, map=None):
        """Compute uvw coordinates of baseline relative to provided RadioBody, 
        or 'z' for zenith uvw coordinates.  If w_only is True, only w (instead
        of (u,v,w) will be returned) - fast."""
        
        x,y,z = self.get_baseline_fast(i,j, src=src, map=map)
        
        afreqs = self.get_afreqs()
        afreqs = numpy.reshape(afreqs, (1,afreqs.size))
        if len(x.shape) == 0:
            if w_only:
                return z*afreqs
            else:
                return numpy.array([x*afreqs, y*afreqs, z*afreqs])
                
        #afreqs = numpy.reshape(afreqs, (1,afreqs.size))
        x.shape += (1,); y.shape += (1,); z.shape += (1,)
        
        if w_only:
            out = numpy.dot(z,afreqs)
        else:
            out = numpy.array([numpy.dot(x,afreqs), numpy.dot(y,afreqs), numpy.dot(z,afreqs)])
            
        return out
        
    def gen_phs_fast(self, src, i, j, mfreq=.150, ionref=None, srcshape=None, resolve_src=False, u=None, v=None, w=None):
        """Return phasing that is multiplied to data to point to src - fast."""
        if ionref is None:
            try:
                ionref = src.ionref
            except AttributeError:
                pass
        if ionref is not None or resolve_src:
            if u is None or v is None or w is None:
                u,v,w = self.gen_uvw_fast(i,j,src=src)
        else:
            if w is None:
                w = self.gen_uvw_fast(i,j,src=src, w_only=True)
        if ionref is not None:
            w += self.refract(u, v, mfreq=mfreq, ionref=ionref)
        o = self.get_phs_offset(i,j)
        phs = numpy.exp(-1j*2*numpy.pi*(w + o))
        if resolve_src:
            if srcshape is None:
                try:
                    srcshape = src.srcshape
                except AttributeError:
                    pass
        if srcshape is not None:
            phs *= self.resolve_src(u, v, srcshape=srcshape)
        return phs.squeeze()
        
    def sim(self, i, j, pol='xx'):
        """
        Simulate visibilities for the specified (i,j) baseline and 
        polarization.  sim_cache() must be called at each time step before 
        this will return valid results.
        
        This function differs from aipy.amp.AntennaArray.sim in the fact that
        *ionref* and *srcshape* are both None in the call to gen_phs and that
        *resolve_src* is set to False.
        """
        
        assert(pol in ('xx','yy','xy','yx'))
        
        if self._cache is None:
            raise RuntimeError('sim_cache() must be called before the first sim() call at each time step.')
        elif self._cache == {}:
            return numpy.zeros_like(self.passband(i,j))
            
        s_eqs = self._cache['s_eqs']
        try:
            s_map = self._cache['s_map']
        except KeyError:
            s_map = None
        w = self.gen_uvw_fast(i, j, src=s_eqs, map=s_map, w_only=True)
        I_sf = self._cache['jys']
        Gij_sf = self.passband(i,j)
        try:
            self.set_active_pol(pol)
            self._cache['s_top'] = self._cache['s_top'].T
            Bij_sf = self.bm_response(i, j)
            self._cache['s_top'] = self._cache['s_top'].T
        except AttributeError:
            # Older versions of AIPY
            Bij_sf = self.bm_response(i, j, pol=pol)
        if len(Bij_sf.shape) == 2:
            Gij_sf = numpy.reshape(Gij_sf, (1, Gij_sf.size))
            
        # Get the phase of each src vs. freq, also does resolution effects
        E_sf = numpy.conjugate( self.gen_phs_fast(s_eqs, i, j, mfreq=self._cache['mfreq'], resolve_src=False, w=w) )
        try:
            E_sf.shape = I_sf.shape
        except(AttributeError):
            pass
            
        # Combine and sum over sources
        GBIE_sf = Gij_sf * Bij_sf * I_sf * E_sf
        Vij_f = GBIE_sf.sum(axis=0)
        
        return Vij_f


def build_sim_array(station, antennas, freq, jd=None, pos_error=0.0, force_flat=False, force_gaussian=False, verbose=False):
    """
    Build a AIPY AntennaArray for simulation purposes.  Inputs are a station 
    object defined from the lwa_common module, a numpy array of stand 
    numbers, and a numpy array of frequencies in either Hz of GHz.  Optional 
    inputs are a Julian Date to set the array to and a positional error terms 
    that perturbs each of the stands in x, y, and z.  The output of this 
    module is an AIPY AntennaArray object.
    
    The shape of the antenna response is either flat (gain of 1 in all 
    directions), modeled by a 2-D Gaussian with the specified full width at
    half maximum in degrees, or modeled by a collection of spherical 
    harmonics that are polynomials in frequency.  The spherical harmonics 
    are used if the file 'beam_shape.npz' is found in the current directory.
    
    .. versionchanged:: 1.0.3
        Changed the meaning of the force_gaussian parameters so that the
        Gaussian full width at half maximum in degrees is passed in.
        
    .. versionchanged:: 1.0.1
        Moved the simulation code over from AIPY to the new _simFast module.  
        This should be much faster but under the caveats that the bandpass
        and antenna gain patterns are the same for all antennas.  This 
        should be a reasonable assumption for large-N arrays.
        
        Added an option to use a 2-D Gaussian beam pattern via the force_gaussian
        keyword.
        
    .. versionchanged:: 0.4.0
        Switched over to passing in Antenna instances generated by the
        :mod:`lsl.common.station` module instead of a list of stand ID numbers.
    """
    
    # If the frequencies are in Hz, we need to convert to GHz
    try:
        freqs = freq.copy()
    except AttributeError:
        freqs = numpy.array(freq)
        if freqs.shape == ():
            freqs.shape = (1,)
    if freqs.min() > 1e6:
        freqs /= 1.0e9
        
    # If the beam Alm coefficient file is present, build a more realistic beam 
    # response.  Otherwise, assume a flat beam
    if force_gaussian:
        try:
            xw, yw = force_gaussian
            xw, yw = float(xw), float(yw)
        except (TypeError, ValueError) as e:
            xw = float(force_gaussian)
            yw = 1.0*xw
            
        # FWHM to sigma
        xw /= 2.0*numpy.sqrt(2.0*numpy.log(2.0))
        yw /= 2.0*numpy.sqrt(2.0*numpy.log(2.0))
        
        # Degrees to radians
        xw *= numpy.pi/180
        yw *= numpy.pi/180
        
        if verbose:
            print("Using a 2-D Gaussian beam with sigmas %.1f by %.1f degrees" % (xw*180/numpy.pi, yw*180/numpy.pi))
        beam = Beam2DGaussian(freqs, xw, yw)
        
    elif force_flat:
        if verbose:
            print("Using flat beam model")
        beam = Beam(freqs)
        
    else:
        if os.path.exists(os.path.join(DATA_PATH, 'beam-shape.npz')):
            dd = numpy.load(os.path.join(DATA_PATH, 'beam-shape.npz'))
            coeffs = dd['coeffs']
            
            deg = coeffs.shape[0]-1
            lmax = int((math.sqrt(1+8*coeffs.shape[1])-3)/2)
            beamShapeDict = {}
            for i in range(deg+1):
                beamShapeDict[i] = numpy.squeeze(coeffs[-1-i,:])
            try:
                dd.close()
            except AttributeError:
                pass
                
            if verbose:
                print("Using Alm beam model with %i-order freq. polynomial and %i-order sph. harmonics" % (deg, lmax))
            beam = BeamAlm(freqs, lmax=lmax, mmax=lmax, deg=deg, nside=128, coeffs=beamShapeDict)
        else:
            if verbose:
                print("Using flat beam model")
            beam = Beam(freqs)
            
    if pos_error != 0:
        warnings.warn(colorfy("{{%%yellow Creating array with positional errors between %.3f and %.3f m}}" % (-pos_error, pos_error)), RuntimeWarning)

    # Build an array of AIPY Antenna objects
    ants = []
    for antenna in antennas:
        top = numpy.array([antenna.stand.x, antenna.stand.y, antenna.stand.z])
        top += (2*pos_error*numpy.random.rand(3)-pos_error)	# apply a random positional error if needed
        top.shape = (3,)
        eq = numpy.dot( aipy.coord.top2eq_m(0.0, station.lat), top )
        eq /= speedOfLight	# m -> s
        eq *= 1e9	# s -> ns
        
        delayCoeff = numpy.zeros(2)
        
        amp = 0*antenna.cable.gain(freqs*1e9) + 1
        
        ants.append( Antenna(eq[0], eq[1], eq[2], beam, phsoff=delayCoeff, amp=amp, stand=antenna.stand.id) )

    # Combine the array of antennas with the array's location to generate an
    # AIPY AntennaArray object
    simAA = AntennaArray(station.aipy_location, ants)
    simAA._station = station
    
    # Set the Julian Data for the AntennaArray object if it is provided.  The try...except
    # clause is used to deal with people who may want to pass an array of JDs in rather than
    # just one.  If one isn't provided, use the date set for the input 'station'.
    if jd is None:
        simAA.date = station.date
    else:
        try:
            simAA.set_jultime(jd[0])
        except TypeError:
            simAA.set_jultime(jd)
            
    return simAA


def __build_sim_data(aa, srcs, pols=['xx', 'yy', 'xy', 'yx'], jd=None, chan=None, phase_center='z', baselines=None, mask=None, verbose=False, count=None, max=None, flat_response=False, resolve_src=False):
    """
    Helper function for build_sim_data so that build_sim_data can be called with 
    a list of Julian Dates and reconstruct the data appropriately.
    
    .. versionchanged:: 1.0.1
        * Moved the simulation code over from AIPY to the new _simFast module.  
          This should be much faster but under the caveats that the bandpass
          and antenna gain patterns are the same for all antennas.  This 
          should be a reasonable assumption for large-N arrays.
        * Added a 'flat_response' keyword to make it easy to toggle on and off
          the spectral and spatial response of the array for the simulation
        * Added a 'resolve_src' keyword to turn on source resolution effects
    """
    
    rawFreqs = aa.get_afreqs()
    
    if len(rawFreqs.shape) == 2:
        nFreq = rawFreqs.shape[1]
    else:
        nFreq = rawFreqs.size
    if chan is None:
        chanMin = 0
        chanMax = -1
    else:
        chanMin = chan[0]
        chanMax = chan[-1]
        
    # Update the JD if necessary
    if jd is None:
        jd = aa.get_jultime()
    else:
        if verbose:
            if count is not None and max is not None:
                print("Setting Julian Date to %.5f (%i of %i)" % (jd, count, max))
            else:
                print("Setting Julian Date to %.5f" % jd)
        aa.set_jultime(jd)
    Gij_sf = aa.passband(0,1)
    def Bij_sf(xyz, pol):
        Bi = aa[0].bm_response(xyz, pol=pol[0]).transpose()
        Bj = aa[1].bm_response(xyz, pol=pol[1]).transpose()
        Bij = numpy.sqrt( Bi*Bj.conj() )
        return Bij.squeeze()
    if flat_response:
        Gij_sf *= 0.0
        Gij_sf += 1.0
        
    # Compute the source parameters
    srcs_tp = []
    srcs_eq = []
    srcs_ha = []
    srcs_dc = []
    srcs_jy = []
    srcs_fq = []
    srcs_sh = []
    if verbose:
        print("Sources Used for Simulation:")
    for name in srcs:
        ## Update the source's coordinates
        src = srcs[name]
        src.compute(aa)
        
        ## Remove sources below the horizon
        srcTop = src.get_crds(crdsys='top', ncrd=3)
        srcAzAlt = aipy.coord.top2azalt(srcTop)
        if srcAzAlt[1] < 0:
            if verbose:
                print("  %s: below horizon" % name)
            continue
            
        ## Topocentric coordinates for the gain pattern calculations
        srcTop.shape = (1,3)
        srcs_tp.append( srcTop )
        
        ## RA/dec -> equatorial 
        srcEQ = src.get_crds(crdsys='eq', ncrd=3)
        srcEQ.shape = (3,1)
        srcs_eq.append( srcEQ )
        
        ## HA/dec
        srcRA,srcDec = aipy.coord.eq2radec(srcEQ)
        srcHA = aa.sidereal_time() - srcRA
        srcs_ha.append( srcHA )
        srcs_dc.append( srcDec )
        
        ## Source flux over the bandpass - corrected for the bandpass
        jys = src.get_jys()
        jys.shape = (1,nFreq)
        jys = Gij_sf * jys
        srcs_jy.append( jys )
        
        ## Frequencies that the fluxes correspond to
        frq = aa.get_afreqs()
        frq.shape = (1,nFreq)
        srcs_fq.append( frq )
        
        ## Source shape parameters
        shp = numpy.array(src.srcshape)
        shp.shape = (3,1)
        srcs_sh.append( shp )
        
    # Build the simulation cache
    aa.sim_cache( numpy.concatenate(srcs_eq, axis=1), 
                jys=numpy.concatenate(srcs_jy, axis=0),
                mfreqs=numpy.concatenate(srcs_fq, axis=0),
                srcshapes=numpy.concatenate(srcs_sh, axis=1) )
    aa._cache['s_top'] = numpy.concatenate(srcs_tp, axis=0)
    aa._cache['s_ha'] = numpy.concatenate(srcs_ha, axis=0)
    aa._cache['s_dec'] = numpy.concatenate(srcs_dc, axis=0)
    
    # Build the simulated data.  If no baseline list is provided, build all 
    # baselines available
    if baselines is None:
        baselines = uvutils.get_baselines(numpy.zeros(len(aa.ants)), indicies=True)
        
    # Define output data structure
    freq = aa.get_afreqs()*1e9
    if len(freq.shape) == 2:
        freq = freq[0,:]
    UVData = VisibilityDataSet(jd, freq, baselines, [], antennaarray=aa, phase_center=phase_center)
    
    # Go!
    if phase_center != 'z':
        phase_center.compute(aa)
        pcAz = phase_center.az*180/numpy.pi
        pcEl = phase_center.alt*180/numpy.pi
    else:
        pcAz = 0.0
        pcEl = 90.0
        
    for p,pol in enumerate(pols):
        ## Apply the antenna gain pattern for each source
        if not flat_response:
            if p == 0:
                for i in range(aa._cache['jys'].shape[0]):
                    aa._cache['jys'][i,:] *= Bij_sf(aa._cache['s_top'][i,:], pol)
            else:
                for i in range(aa._cache['jys'].shape[0]):
                    aa._cache['jys'][i,:] /= Bij_sf(aa._cache['s_top'][i,:], pols[p-1])	# Remove the old pol
                    aa._cache['jys'][i,:] *=  Bij_sf(aa._cache['s_top'][i,:], pol)
                    
        ## Simulate
        if not flat_response:
            uvw1, vis1 = FastVis(aa, baselines, chanMin, chanMax, pcAz, pcEl, resolve_src=resolve_src)
        else:
            currentVars = locals().keys()
            if 'uvw1' not in currentVars or 'vis1' not in currentVars:
                uvw1, vis1 = FastVis(aa, baselines, chanMin, chanMax, pcAz, pcEl, resolve_src=resolve_src)
                
        ## Unpack the data and add it to the data set
        if p == 0:
            UVData.uvw = uvw1
        pds = PolarizationDataSet(pol.upper(), data=vis1)
        UVData.append( pds )
        
    # Cleanup
    if not flat_response:
        for i in range(aa._cache['jys'].shape[0]):
            aa._cache['jys'][i,:] /= Bij_sf(aa._cache['s_top'][i,:], pols[-1])	# Remove the old pol
            aa._cache['jys'][i,:] /= Gij_sf								# Remove the bandpass
            
    # Return
    return UVData


def build_sim_data(aa, srcs, pols=['xx', 'yy', 'xy', 'yx'], jd=None, chan=None, phase_center='z', baselines=None, mask=None,  flat_response=False, resolve_src=False, verbose=False):
    """
    Given an AIPY AntennaArray object and a dictionary of sources from 
    aipy.src.get_catalog, returned a :class:`lsl.imaging.data.VisibilityDataSet` 
    object of simulated data taken at zenith.  Optinally, the data can be 
    masked using some referenced (observed) data set or only a specific sub-set 
    of baselines.
    
    .. versionchanged:: 1.0.1
        * Added a 'flat_response' keyword to make it easy to toggle on and off
          the spectral and spatial response of the array for the simulation
        * Added a 'resolve_src' keyword to turn on source resolution effects
        
    ..versionchanged:: 0.4.0
        Added the 'pols' keyword to only compute certain polarization components
    """
    
    # Update the JD if necessary
    if jd is None:
        jd = [aa.get_jultime()]
    else:
        try:
            len(jd)
        except TypeError:
            jd = [jd]
            
    # Build up output dictionary
    freq = aa.get_afreqs()*1e9
    if len(freq.shape) == 2:
        freq = freq[0,:]
    UVData = VisibilityData()
    
    # Loop over Julian days to fill in the simulated data set
    jdCounter = 1
    for juldate in jd:
        oBlk = __build_sim_data(aa, srcs, pols=pols, jd=juldate, chan=chan, phase_center=phase_center, baselines=baselines, mask=mask, verbose=verbose, count=jdCounter, max=len(jd), flat_response=flat_response, resolve_src=resolve_src)
        jdCounter = jdCounter + 1
        UVData.append( oBlk )
        
    # Trim
    if len(UVData) == 1:
        UVData = UVData.pop()
        
    return UVData


def scale_data(dataSet, amps, delays, phase_offsets=None):
    """
    Apply a set of antenna-based real gain values and phase delays in ns to a 
    :class:`lsl.imaging.data.VisibilityDataSet` object.  Returned the new
    scaled and delayed VisibilityDataSet.
    
    ..versionchanged:: 0.6.3
        Added a keyword so that phase offsets (in radians) can also be specified
    
    ..versionchanged:: 0.4.0
        The delays are now expected to be in nanoseconds rather than radians.
    """
    
    # Make sure we have the right kind of object
    if not isinstance(dataSet, (VisibilityDataSet, VisibilityData)):
        raise TypeError("Expected data to be stored in a VisibilityData or VisibilityDataSet object")
        
    # Build the VisibilityDataSet to hold the scaled and delayed data
    sclData = dataSet.copy(include_pols=True)
    fq = dataSet.freq / 1e9
    
    cGains = []
    for i in range(len(amps)):
        gain = 2j*numpy.pi*fq*delays[i]
        if phase_offsets is not None:
            gain += 1j*phase_offsets[i]
        cGains.append( amps[i]*numpy.exp(gain) )

    # Apply the scales and delays for all polarization pairs found in the original data
    for pds in sclData:
        if isinstance(pds, VisibilityDataSet):
            for ppds in pds:
                for b,(i,j) in enumerate(sclData.baselines):
                    ppds.data[b,:] *= cGains[j].conj()*cGains[i] 
        else:
            for b,(i,j) in enumerate(sclData.baselines):
                pds.data[b,:] *= cGains[j].conj()*cGains[i]
                
    return sclData
    

def shift_data(dataSet, aa):
    """
    Shift the uvw coordinates in one :class:`lsl.imaging.data.VisibilityDataSet` 
    object to a new set of uvw coordinates that correspond to a new 
    AntennaArray object.  This is useful for looking at how positional errors 
    in the array affect the data.
    """
    
    # Make sure we have the right kind of object
    if not isinstance(dataSet, (VisibilityData, VisibilityDataSet)):
        raise TypeError("Expected data to be stored in a VisibilityData or VisibilityDataSet object")
        
    # Build the VisibilityDataSet to hold the scaled and delayed data
    shftData = dataSet.copy(include_pols=True)
    
    # Apply the coordinate shift
    if isinstance(shftData, VisibilityData):
        for data_set in shftData:
            for b,(i,j) in enumerate(data_set.baselines):
                crds = aa.gen_uvw(j, i, src=data_set.phase_center)[:,0,:]
                data_set.uvw[b,:] = crds
    else:
        for b,(i,j) in enumerate(shftData.baselines):
            crds = aa.gen_uvw(j, i, src=shftData.phase_center)[:,0,:]
            shftData.uvw[b,:] = crds
            
    return shftData


def add_baseline_noise(dataSet, SEFD, tInt, bandwidth=None, efficiency=1.0):
    """
    Given a :class:`lsl.imaging.data.VisibilityDataSet` object of visibilities, 
    an SEFD or array SEFDs in Jy, and an integration time in seconds, add noise 
    to the visibilities assuming that the "weak source" limit.
    
    This function implements Equation 9-15 from Chapter 9 of "Synthesis 
    Imaging in Radio Astronomy II".
    
    .. versionadded:: 1.0.2
    """
    
    # Make sure we have the right kind of object
    if not isinstance(dataSet, (VisibilityData, VisibilityDataSet)):
        raise TypeError("Expected data to be stored in a VisibilityData or VisibilityDataSet object")
        
    # Figure out the bandwidth from the frequency list in the 
    if bandwidth is None:
        try:
            bandwidth = dataSet.freq[1] - dataSet.freq[0]
        except IndexError:
            raise RuntimeError("Too few frequencies to determine the bandwidth, use the 'bandwidth' keyword")
            
    # Make sure that we have an SEFD for each antenna
    ants = []
    for i,j in dataSet.baselines:
        if i not in ants:
            ants.append( i )
        if j not in ants:
            ants.append( j )
    nAnts = len(ants)
    try:
        if len(SEFD) != nAnts:
            raise RuntimeError("Mis-match between the number of SEFDs supplied and the number of antennas in the data")
            
    except TypeError:
        SEFD = numpy.ones(nAnts)*SEFD
        
    # Calculate the standard deviation of the real/imaginary noise
    visNoiseSigma = numpy.zeros((len(dataSet.baselines), len(dataSet.freq)))
    for k,(i,j) in enumerate(dataSet.baselines):
        visNoiseSigma[k,:] = numpy.sqrt(SEFD[i]*SEFD[j])
    visNoiseSigma *= 1 / efficiency / numpy.sqrt(2.0*bandwidth*tInt)
    
    # Build the VisibilityDataSet to hold the data with noise added
    bnData = dataSet.copy(include_pols=True)
    
    # Apply the scales and delays for all polarization pairs found in the original data
    for pds in bnData:
        if isinstance(pds, VisibilityDataSet):
            for ppds in pds:
                ## Calculate the expected noise
                visNoise = visNoiseSigma * (numpy.random.randn(*ppds.data.shape) \
                           + 1j*numpy.random.randn(*ppds.data.shape))
                
                ## Apply
                ppds.data += visNoise
        else:
            ## Calculate the expected noise
            visNoise = visNoiseSigma * (numpy.random.randn(*pds.data.shape) \
                       + 1j*numpy.random.randn(*pds.data.shape))
            
            ## Apply
            pds.data += visNoise
            
    return bnData
