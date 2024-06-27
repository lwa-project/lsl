"""
Module for virtually converting/storing correlator output so that it can be
directly worked with in the :mod:`lsl.imaging.utils` module.

.. versionadded:: 2.1.3
"""

import aipy
import ephem
import numpy as np
from datetime import datetime

from astropy import units as astrounits
from astropy.time import Time as AstroTime
from astropy.constants import c as speedOfLight
from astropy.coordinates import AltAz, ITRS, FK5

from lsl import astro
from lsl.reader.base import FrameTimestamp
from lsl.common.mcs import datetime_to_mjdmpm
from lsl.imaging.data import VisibilityData, VisibilityDataSet, PolarizationDataSet
from lsl.writer.fitsidi import WriterBase
from lsl.sim.vis import build_sim_array
from lsl.correlator.uvutils import compute_uvw


__version__ = '0.1'
__all__ = ['VirtualWriter',]


speedOfLight = speedOfLight.to('m/s').value


class VirtualWriter(WriterBase):
    """
    Class for storing/converting visibility data in memory and provided access
    to it in a way that is compatible with the :mod:`lsl.imaging.utils` module.
    """
    
    def __init__(self, ref_time=0.0, verbose=False):
        """
        Initialize a new virtual writer object using a reference time given in
        seconds since the UNIX 1970 epoch, a python datetime object, or a
        string in the format of 'YYYY-MM-DDTHH:MM:SS'.
        """
        
        WriterBase. __init__(self, '', ref_time=ref_time, verbose=verbose)
        self.visibility_data = VisibilityData()
        
    def set_frequency(self, freq):
        """
        Given a numpy array of frequencies, set the relevant common observation
        parameters.
        """
        
        self.freq = freq
        self._build_antenna_array()
        
    def set_geometry(self, site, antennas, bits=8):
        """
        Given a station and an array of stands, set the relevant common observation
        parameters.
        """
        
        self.site = site
        self.antennas = antennas
        self._build_antenna_array()
        
    def _build_antenna_array(self):
        """
        Update the internal observer and antenna array information if and only
        if both the set_frequency and set_geometry methods have been called.
        """
        
        try:
            assert(len(self.freq) > 0)
            self.site
            self.antennas
        except (AssertionError, AttributeError):
            return False
            
        # Convert the reference time to a JD
        ref_time = datetime.strptime(self.ref_time, "%Y-%m-%dT%H:%M:%S")
        mjd, mpm = datetime_to_mjdmpm(ref_time)
        ref_jd = mjd + mpm/1000.0 / 86400 + astro.MJD_OFFSET
        
        # Create the observer and antenna array
        self.el= self.site.earth_location
        self.antenna_array = build_sim_array(self.site, self.antennas,
                                             self.freq/1e9, jd=ref_jd)
        
        return True
        
    def _build_uvw(self, jd, baselines, source):
        """
        Generate and return a 3-D array of (u,v,w) coordiantes for the specified
        baselines and source at the given UTC JD.
        """
        
        date = AstroTime(jd, format='jd', scale='utc')
        
        if source == 'z':
            tc = AltAz('0deg', '90deg', location=self.el, obstime=date)
            equ = tc.transform_to(FK5(equinox=date))
        else:
            equ = FK5(source.a_ra*astrounits.rad, source.a_dec*astrounits.rad,
                      equinox=date)
            
        # Phase center coordinates
        it = equ.transform_to(ITRS(location=self.el, obstime=date))
        HA = ((self.el.lon - it.spherical.lon).wrap_at('180deg')).deg
        dec = it.spherical.lat.deg
        
        # (u,v,w) coordinates
        uvw = compute_uvw(baselines, HA=HA, dec=dec, freq=self.freq, site=self.el)
        
        return uvw
        
    def add_comment(self, comment):
        """
        Dummy method for compatibilty with the :class:`lsl.writer.fitsidi.WriterBase`
        class.
        """
        
        raise NotImplementedError
            
    def add_history(self, history):
        """
        Dummy method for compatibilty with the :class:`lsl.writer.fitsidi.WriterBase`
        class.
        """
        
        raise NotImplementedError
        
    def add_header_keyword(self, name, value, comment=None):
        """
        Dummy method for compatibilty with the :class:`lsl.writer.fitsidi.WriterBase`
        class.
        """
        
        raise NotImplementedError
        
    def convert_to_data_set(self,  obsTime, intTime, baselines, visibilities, weights=None, pol='XX', source='z'):
        """
        Create a :class:`lsl.imaging.data.VisibilityDataSet` object to store a
        collection of visibilities for the specified TAI MJD time.  This is
        similar to the add_data_set method but the VisbilityDataSet is returned
        instead of appended to the internal data structure.
        """
        
        if isinstance(obsTime, FrameTimestamp):
            obsTime = obsTime.tai_mjd
        elif isinstance(obsTime, AstroTime):
            obsTime = obsTime.tai.mjd
            
        numericBaselines = []
        for (a1,a2) in baselines:
            numericBaselines.append((self.antennas.index(a1), self.antennas.index(a2)))
            
        if source == 'z':
            ot = AstroTime(obsTime, format='mjd', scale='tai').utc
            tc = AltAz('0deg', '90deg', location=self.el, obstime=ot)
            pc = tc.transform_to(FK5(equinox=ot))
        
            phase_center = aipy.amp.RadioFixedBody(pc.ra.rad, pc.dec.rad,
                                                   name=f"ZA{pc.ra.to_string(sep='')}",
                                                   epoch=ot.jd - astro.DJD_OFFSET)
            old_obs_date = self.antenna_array.date
            self.antenna_array.date = ot.iso
            phase_center.compute(self.antenna_array)
            self.antenna_array.date = old_obs_date
        else:
            phase_center = source
                                                      
        obsJD = astro.taimjd_to_utcjd(obsTime)
        uvw = self._build_uvw(obsJD, baselines, source)
        pds = PolarizationDataSet(pol, visibilities, weight=weights)
        vds = VisibilityDataSet(obsJD, self.freq, numericBaselines, uvw,
                                antennaarray=self.antenna_array,
                                phase_center=phase_center)
        vds.append(pds)
        return vds
        
    def add_data_set(self, obsTime, intTime, baselines, visibilities, weights=None, pol='XX', source='z'):
        """
        Create a class:`lsl.imaging.data.VisibilityDataSet` object to store a
        collection of visibilities for the specified TAI MJD time.
        """
        
        vds = self.convert_to_data_set(obsTime, intTime, baselines, visibilities,
                                       weights=weights, pol=pol, source=source)
        self.visibility_data.append(vds)
        
    def write(self):
        """
        Dummy method for compatibilty with the :class:`lsl.writer.fitsidi.WriterBase`
        class.
        """
        
        raise NotImplementedError
            
    def close(self):
        """
        Dummy method for compatibilty with the :class:`lsl.writer.fitsidi.WriterBase`
        class.
        """
        
        raise NotImplementedError
        
    def get_data_set(self, sets, include_auto=False, sort=True, min_uv=0, max_uv=np.inf):
        """
        Return a :class:`lsl.imaging.data.VisibilityDataSet` or 
        :class:`lsl.imaging.data.VisibilityData` object for all 
        baselines for a given set of observations for the specified data set.  
        By default this excludes the autocorrelations.  To include 
        autocorrelations set the value of 'include_auto' to True.  Setting the
        'sort' keyword to False will disable the baseline sorting.  Optionally,
        baselines with lengths between min_uv and max_uv can only be returned.

        .. note::
            min_uv and max_uv should be specified in lambda
        """
        
        dataSets = VisibilityData()
        try:
            len(sets)
        except TypeError:
            sets = range(sets, sets+1)
        for set in sets:
            dataSets.append(self.visibility_data[set-1])
            
        # Sort
        if sort:
            dataSets.sort()
            
        # Prune
        if not include_auto and min_uv == 0:
            min_uv = 1e-3
        if min_uv != 0 or max_uv != np.inf:
            dataSets = dataSets.get_uv_range(min_uv=min_uv, max_uv=max_uv)
            
        # Prune a different way
        if len(dataSets) == 1:
            dataSets = dataSets.pop()
            
        # Return
        return dataSets
