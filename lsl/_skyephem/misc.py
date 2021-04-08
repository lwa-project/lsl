"""
Python module that provides miscellanous functions from ephem related to the
equinoxes, solstices, and the phases of the Moon.
"""

from __future__ import print_function, division

import numpy
from scipy.optimize import minimize_scalar

from skyfield import api

from lsl._skyephem.dates import Date
from lsl._skyephem.cache import load_planetary_ephemeris


__all__ = ['previous_equinox', 'next_equinox', 'previous_solstice', 'next_solstice',
           'previous_new_moon', 'next_new_moon', 'previous_first_quarter_moon',
           'next_first_quarter_moon', 'previous_full_moon', 'next_full_moon',
           'previous_last_quarter_moon', 'next_last_quarter_moon']


_solar_system = load_planetary_ephemeris()
_sol = _solar_system['sun']
_ter = _solar_system['earth']
_lun = _solar_system['moon']


def _seasons(djd, value, specific):
    """
    Private function to help determine when the Sun is near a particular
    ecliptic longitude.  specific controlls whether or not you are looking
    for a specific longitude (1) or if you looking for any equinox or any
    solstice (0).
    """
    
    d = Date(djd)
    e = _ter.at(d)
    s = e.observe(_sol)
    lat, lon, _ = s.ecliptic_latlon(epoch=d)
    
    return abs(lon.radians - value) % (numpy.pi+specific*numpy.pi)


def previous_equinox(date):
    """
    Return the date of the previous equinox relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi, 0), method='bounded',
                          bounds=(date-365.25, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_equinox(date):
    """
    Return the date of the next equinox relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi, 0), method='bounded',
                          bounds=(date*1.0, date+365.25),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_solstice(date):
    """
    Retun the date of the previous solstice relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi/2, 0), method='bounded',
                          bounds=(date-365.25, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_solstice(date):
    """
    Retun the date of the next solstice relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi/2, 0), method='bounded',
                          bounds=(date*1.0, date+365.25),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def _phases(djd, value):
    """
    Private function to help determine when the Moon is near a particular phase.
    """
    
    d = Date(djd)
    e = _ter.at(d)
    s = e.observe(_sol)
    m = e.observe(_lun)
    s_elon = s.ecliptic_latlon(epoch=d)[0]
    m_elon = m.ecliptic_latlon(epoch=d)[0]
    return abs(m_elon.radians - s_elon.radians - value) % (2*numpy.pi)


def previous_new_moon(date):
    """
    Return the date of the previous new moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(0.0,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_new_moon(date):
    """
    Return the date of the next new moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(0.0,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_first_quarter_moon(date):
    """
    Return the date of the previous first quarter moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi/2,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_first_quarter_moon(date):
    """
    Return the date of the next first quarter moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi/2,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_full_moon(date):
    """
    Return the date of the previous full moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_full_moon(date):
    """
    Return the date of the next full moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_last_quarter_moon(date):
    """
    Return the date of the previous last quarter moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(3*numpy.pi/2,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_last_quarter_moon(date):
    """
    Return the date of the next last quarter moon relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_phases, args=(3*numpy.pi/2,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)
