import os
import sys

from skyfield import api

from lsl.config import LSL_CONFIG
EPHEM_CONFIG = LSL_CONFIG.view('skyephem')


__all__ = ['load_planetary_ephemeris',]


# Create the cache directory
if not os.path.exists(os.path.join(os.path.expanduser('~'), '.lsl')):
    os.mkdir(os.path.join(os.path.expanduser('~'), '.lsl'))
_CACHE_DIR = os.path.join(os.path.expanduser('~'), '.lsl', 'skyfield_cache')
if not os.path.exists(_CACHE_DIR):
    os.mkdir(_CACHE_DIR)


def load_planetary_ephemeris():
    load = api.Loader(_CACHE_DIR)
    solar_system = load(EPHEM_CONFIG.get('ephemeris'))
    return solar_system
