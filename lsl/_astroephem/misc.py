"""
Python module that provides miscellanous functions from ephem related to the
equinoxes, solstices, and the phases of the Moon.
"""

from __future__ import print_function, division

import numpy
from scipy.optimize import minimize_scalar

from astropy import units as u
from astropy.coordinates import get_body, GeocentricTrueEcliptic

from dates import Date

__all__ = ['previous_equinox', 'next_equinox', 'previous_solstice', 'next_solstice',
           'previous_new_moon', 'next_new_moon', 'previous_first_quarter_moon',
           'next_first_quarter_moon', 'previous_full_moon', 'next_full_moon',
           'previous_last_quarter_moon', 'next_last_quarter_moon']


def _seasons(djd, value):
    """
    Private function to help determine when the Sun is near a particular
    ecliptic longitude.
    """
    
    d = Date(djd)
    s = get_body('sun', d)
    ecl = GeocentricTrueEcliptic(equinox=d)
    s_ecl = s.transform_to(ecl)
    return abs(s_ecl.lon.rad - value) % numpy.pi


def previous_equinox(date):
    """
    Return the date of the previous equinox relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi,), method='bounded',
                          bounds=(date-365.25, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_equinox(date):
    """
    Return the date of the next equinox relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi,), method='bounded',
                          bounds=(date*1.0, date+365.25),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_solstice(date):
    """
    Retun the date of the previous solstice relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi/2,), method='bounded',
                          bounds=(date-365.25, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_solstice(date):
    """
    Retun the date of the next solstice relative to the date provided.
    """
    
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi/2,), method='bounded',
                          bounds=(date*1.0, date+365.25),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def _phases(djd, value):
    """
    Private function to help determine when the Moon is near a particular phase.
    """
    
    d = Date(djd)
    s = get_body('sun', d)
    m = get_body('moon', d)
    ecl = GeocentricTrueEcliptic(equinox=d)
    s_ecl = s.transform_to(ecl)
    m_ecl = m.transform_to(ecl)
    return abs(m_ecl.lon.rad - s_ecl.lon.rad - value) % (2*numpy.pi)


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
