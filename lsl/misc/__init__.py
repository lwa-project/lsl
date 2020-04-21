"""
lsl.misc - Miscellanous modules including:
 * mathutil       - math utilities that were part of the lwa_user package, 
 * beamformer     - post data aquisition beam former for TBW and TBN data,
 * ionosphere     - access to ionospheric and geomagnetic models, 
 * dedispersion   - dedispersion module,
 * scattering     - multi-path scattering module, 
 * rfutil         - convert RF engineering lingo into radio astrnomy lingo, and
 * wisdom         - build LSL-specific FFTW wisdom.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    