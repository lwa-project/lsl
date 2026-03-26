"""
Module to simulate observations made with the DP system.
"""

import time
import numpy as np
from aipy import coord as aipycoord
from astropy.constants import c as speedOfLight

from lsl import astro
from lsl.common import ndp as ndp_common
from lsl.common import stations as lwa_common
from lsl.sim import drx
from lsl.reader.drx import FILTER_CODES as DRXFilters

from lsl.misc import telemetry
telemetry.track_module()
speedOfLight = speedOfLight.to('m/s').value


__version__ = '0.6'
__all__ = ['basic_signal',]


def _basic_drx(fh, stands, nframes, **kwargs):
    """
    Private function for generating a basic TBN signal.
    """

    start_time = kwargs['start_time']
    filter = kwargs['filter']
    ntuning = kwargs['ntuning']
    verbose = kwargs['verbose']
    noise_strength = kwargs['noise_strength']
    sample_rate = DRXFilters[filter]
    decimation = int(round(ndp_common.fS / sample_rate))
    
    maxValue = 7
    samplesPerFrame = 4096
    upperSpike1 = sample_rate / 4.0
    lowerSpike1 = -sample_rate / 4.0
    upperSpike2 = sample_rate / 3.0
    lowerSpike2 = -sample_rate / 3.0

    if verbose:
        print("Simulating %i frames of DRX Data @ %.2f MHz for %i beams, %i tunings each:" % \
            (nframes, sample_rate/1e6, len(stands), ntuning))

    beams = stands
    for i in range(nframes):
        if i % 1000 == 0 and verbose:
            print(" frame %i" % i)
        t = int(start_time*ndp_common.fS) + int(i*ndp_common.fS*samplesPerFrame/sample_rate)
        tFrame = t/ndp_common.fS - start_time + np.arange(samplesPerFrame, dtype=np.float32) / sample_rate
        for beam in beams:
            for tune in range(1, ntuning+1):
                for pol in (0, 1):
                    if tune == 1:
                        if pol == 0:
                            spike = upperSpike1
                        else:
                            spike = lowerSpike1
                    else:
                        if pol == 0:
                            spike = lowerSpike2
                        else:
                            spike = upperSpike2
                            
                    cFrame = drx.SimFrame(beam=beam, tune=tune, pol=pol, frame_count=i+1, decimation=decimation, time_offset=0, obs_time=t, flags=0)
                    cFrame.data = np.zeros(samplesPerFrame, dtype=np.complex64)
                    cFrame.data += np.random.randn(samplesPerFrame) + 1j*np.random.randn(samplesPerFrame)
                    cFrame.data *= maxValue*noise_strength
                    cFrame.data += maxValue*np.exp(2j*np.pi*spike*tFrame)
                    cFrame.write_raw_frame(fh)


def basic_signal(fh, stands, nframes, station=lwa_common.lwa1, mode='DRX', filter=6, ntuning=2, start_time=0, noise_strength=0.1, verbose=False):
    """
    Generate a collection of frames with a basic test signal for DRX.  The 
    signals for the three modes are:

    DRX
     * noise + (sample_rate/4) kHz signal for x-pol. and noise + 
        (-sample_rate/4) for y-pol. -> tuning 1
     * noise + (-sample_rate/3) kHz signal for x-pol. and noise + 
        (sample_rate/3) for y-pol. -> tuning 2
        
    All modes need to have stands (a list of integer beams numbers for DRX)
    and number of frames to generate.  The DRX frames need the 'filter'
    keyword set to specify the filter width.  In addition, the 'stands' 
    argument is interpreted as beam numbers for DRX.
    
    .. versionchanged:: 0.6
        Dropped support for TBN.
        
    .. versionchanged:: 0.4.4
        Added the `noise_strength` keyword to control how much noise is added to 
        the data.
        
    .. versionchanged:: 2.0.0
        Removed support for generating TBW data.
        
    .. versionchanged:: 2.1.8:
        Add the `station` keyword and documentation cleanup
        `stands` is now a list of :class:`lsl.common.stations.Antenna`
        instances for TBN
    """

    if start_time == 0:
        start_time = time.time()

    if mode == 'DRX':
        _basic_drx(fh, stands, nframes, filter=filter, ntuning=ntuning, start_time=start_time, noise_strength=noise_strength, verbose=verbose)
    else:
        raise RuntimeError(f"Unknown observations mode: {mode}")
