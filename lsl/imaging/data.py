"""
Classes to hold visibility data that are used internally with LSL.
"""

import aipy
import copy
import numpy

from lsl import astro

__version__ = "0.1"
__revision__ = "$Rev$"
__all__ = ['PolarizationDataSet', 'VisibilityDataSet', 'VisibilityData', '__version__', '__revision__', '__all__']


class PolarizationDataSet(object):
    """
    Class to hold the visibility data/weight/mask for a single polarization.
    """
    
    def __init__(self, polarization, data=None, weight=None, mask=None):
        self.polarization = polarization.upper()
        if data is not None:
            self.data = data
        else:
            self.data = numpy.array([], dtype=numpy.complex64)
            
        if weight is None:
            weight = numpy.ones(self.data.shape, dtype=numpy.float32)
        self.weight = weight
        
        if mask is None:
            mask = numpy.zeros(self.data.shape, dtype=numpy.bool)
        self.mask = mask
            
    @property
    def nbaseline(self):
        return len(self.data)
        
    def copy(self):
        """
        Return a deep copy of the object.
        """
        
        pol_copy = PolarizationDataSet(self.polarization, 
                                       data=self.data.copy(), 
                                       weight=self.weight.copy(), 
                                       mask=self.mask.copy())
        return pol_copy
        
    def append(self, data, weight=None, mask=None):
        """
        Append a new visibility to the object.
        """
        
        if self.data.size == 0:
            ## Deal with the first batch of data being added
            self.data = data
            
            if weight is None:
                weight = numpy.ones(data.shape, dtype=numpy.float32)
            if self.weight.shape != data.shape:
                raise ValueError("weight shape does not match the data shape")
            self.weight = weight
            
            if mask is None:
                mask = numpy.zeros(data.shape, dtype=numpy.bool)
            if self.mask.shape != data.shape:
                raise ValueError("mask shape does not match the data shape")
            self.mask = mask
            
        else:
            self.data = numpy.vstack([self.data, data])
            
            if weight is None:
                weight = numpy.ones(data.shape, dtype=numpy.float32)
            if self.weight.shape != data.shape:
                raise ValueError("weight shape does not match the data shape")
            self.weight = numpy.vstack([self.weight, weight])
            
            if mask is None:
                mask = numpy.zeros(data.shape, dtype=numpy.bool)
            if self.mask.shape != data.shape:
                raise ValueError("mask shape does not match the data shape")
            self.mask = numpy.vstack([self.mask, mask])
            
    def extend(self, data, weight=None, mask=None):
        """
        Extend the object with a collection of visibilites.
        """
        
        if weight is None:
            weight = [None for i in data]
        if mask is None:
            mask = [None for i in data]
        for d,w,m in zip(data, weight, mask):
            self.append(d, weight=w, mask=m)
            
    def sort(self, order):
        """
        Apply the specified ordering on the object.
        """
        
        if len(order) != self.data.shape[0]:
            raise ValueError("order length does not match data length")
        self.data = numpy.take(self.data, order, axis=0)
        self.weight = numpy.take(self.weight, order, axis=0)
        self.mask = numpy.take(self.mask, order, axis=0)
        
    def subset(self, selection):
        """
        Return a subset of the data.
        """
        
        new_pol = PolarizationDataSet(self.polarization, 
                                      data=numpy.take(self.data, selection, axis=0),
                                      weight=numpy.take(self.weight, selection, axis=0), 
                                      mask=numpy.take(self.mask, selection, axis=0))
        return new_pol


class VisibilityDataSet(object):
    """
    Class to wrap a single integration of visibility data.
    """
    
    def __init__(self, jd, freq, baselines, uvw, antennaarray=None, phase_center=None, **kwds):
        self.jd = jd
        self.freq = freq
        self.baselines = baselines
        self.uvw = uvw
        if antennaarray is not None and not isinstance(antennaarray, aipy.amp.AntennaArray):
            raise TypeError("Expected antennaarray to be either None or AntennaArray")
        self.antennaarray = antennaarray.copy()
        if self.antennaarray is not None:
            self.antennaarray.set_jultime(self.jd)
        self.phase_center = phase_center
        self.pols = [key for key in kwds]
        for key in kwds:
            setattr(self, key.upper(), kwds[key])
            
    def __contains__(self, value):
        return True if value in self.pols else False
        
    def __iter__(self):
        i = 0
        while i < self.npol:
            pol = self.pols[i]
            yield getattr(self, pol)
            i += 1
           
    @property
    def nbaseline(self):
        return len(self.baselines)
        
    @property
    def nchan(self):
        return len(self.freq)
        
    @property
    def npol(self):
        return len(self.pols)
        
    @property
    def mjd(self):
        return self.jd - astro.MJD_OFFSET
        
    def copy(self, include_pols=False):
        set_copy = VisibilityDataSet(self.jd*1.0, 
                                     self.freq.copy(), 
                                     copy.deepcopy(self.baselines), 
                                     self.uvw.copy(), 
                                     antennaarray=self.antennaarray.copy(), 
                                     phase_center=self.phase_center)
        if include_pols:
            for pol in self:
                pol_copy = pol.copy()  
                set_copy.append(pol_copy)
        return set_copy               
        
    def append(self, value):
        if not isinstance(value, PolarizationDataSet):
            raise TypeError("Must be a VisibilityDataSet instance")
        if self.nbaseline != value.nbaseline:
            raise RuntimeError("Baseline count mis-match")
            
        pol = value.polarization.upper()
        if pol not in self.pols:
            self.pols.append(pol)
        setattr(self, pol, value)
        
    @property
    def _baseline_order(self):
        def __cmpBaseline(bl):
            return 4096*bl[0] + bl[1]
            
        return [i for (v, i) in sorted((v, i) for (i, v) in enumerate([__cmpBaseline(bl) for bl in self.baselines]))]
        
    def sort(self, order=None):
        if order is None:
            order = self._baseline_order
            
        self.baselines = [self.baselines[i] for i in order]
        self.uvw = numpy.take(self.uvw, order, axis=0)
        for pds in self:
            pds.sort(order)
            
    def rephase(self, new_phase_center):
        if self.phase_center is None:
            raise AttributeError("No phase center defined for this data set")
        if self.antennaarray is None:
            raise AttributeError("No anntennaarray defined for this data set")
            
        # Update the time for the AntennaArray
        self.antennaarray.set_jultime(self.jd)
        
        # Recompute
        if self.phase_center is not 'z':
            self.phase_center.compute(self.antennaarray)
        if new_phase_center is not 'z':
            new_phase_center.compute(self.antennaarray)
            
        # Go!
        for pds in self:
            pds.data = pds.data.astype(numpy.complex128)
            ## Loop over baselines
            for k in xrange(self.nbaseline):
                ### Load in the data
                i,j = self.baselines[k]
                vis = pds.data[k,:]
                
                ### Compute the uvw coordinates and the new phasing
                try:
                    crd = self.antennaarray.gen_uvw(j, i, src=new_phase_center)[:,0,:]
                    vis = self.antennaarray.unphs2src(vis, self.phase_center, j, i)
                    vis = self.antennaarray.phs2src(vis, new_phase_center, j, i)
                except aipy.phs.PointingError:
                    raise RuntimeError("Rephasing center is below the horizon")
                    
                ### Save
                self.uvw[k,:,:] = crd
                pds.data[k,:] = vis
                
        # Update the phase center
        self.phase_center = new_phase_center
        
    def get_uv_range(self, min_uv=0.0, max_uv=numpy.inf):
        # Force min to be less than max
        if min_uv > max_uv:
            temp = min_uv
            min_uv = max_uv
            max_uv = temp
            
        # Find out the baseline lengths and create a list of good ones
        center_uvw = self.uvw[:,:2,self.nchan/2]
        center_uvw = numpy.sqrt( (center_uvw**2).sum(axis=1) ).ravel()
        selection = numpy.where( (center_uvw >= min_uv) & (center_uvw < max_uv) )[0]
        
        # Create the new data set
        new_baselines = [self.baselines[b] for b in selection]
        new_data = self.copy(include_pols=False)
        new_data.baselines = new_baselines
        
        for pds in self:
            new_data.append( pds.subset(selection) )
            
        # Done
        return new_data


class VisibilityData(object):
    """
    Class to hold multiple integrations of visibility data.
    """
    
    def __init__(self, data=None):
        if data is not None:
            self._data = data
        else:
            self._data = []
            
    def __iter__(self):
        return self._data.__iter__()
        
    def __len__(self):
        return len(self._data)
        
    def __getitem__(self, value):
        return self._data[value]
        
    def __setitem__(self, index, value):
        if not isinstance(value, VisibilityDataSet):
            raise TypeError("Expected type to be VisibilityDataSet")
        self._data[index] = value
        
    def append(self, value):
        if not isinstance(value, VisibilityDataSet):
            raise TypeError("Expected type to be VisibilityDataSet")
        self._data.append(value)
        
    def extend(self, values):
        for value in values:
            if not isinstance(value, VisibilityDataSet):
                raise TypeError("Expected type to be VisibilityDataSet")
            self._data.append(value)
            
    def pop(self, index=-1):
        return self._data.pop(index)
        
    def sort(self):
        self._data.sort(key=lambda x: x.jd)
        for data_set in self._data:
            data_set.sort()
            
    def rephase(self, new_phase_center, ignore_errors=False):
        for data_set in self._data:
            try:
                data_set.rephase(new_phase_center)
            except Exception as e:
                if not ignore_errors:
                    raise e
                    
    def get_uv_range(self, min_uv=0.0, max_uv=numpy.inf):
        new_data = VisibilityData()
        for data_set in self:
            new_data.append( data_set.get_uv_range(min_uv=min_uv, max_uv=max_uv) )
        return new_data
