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
    d = Date(djd)
    s = get_body('sun', d)
    ecl = GeocentricTrueEcliptic(equinox=d)
    s_ecl = s.transform_to(ecl)
    return abs(s_ecl.lon.rad - value) % numpy.pi


def previous_equinox(date):
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi,), method='bounded',
                          bounds=(date-365.25, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_equinox(date):
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi,), method='bounded',
                          bounds=(date*1.0, date+365.25),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_solstice(date):
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi/2,), method='bounded',
                          bounds=(date-365.25, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_solstice(date):
    date = Date(date)
    sol = minimize_scalar(_seasons, args=(numpy.pi/2,), method='bounded',
                          bounds=(date*1.0, date+365.25),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def _phases(djd, value):
    d = Date(djd)
    s = get_body('sun', d)
    m = get_body('moon', d)
    ecl = GeocentricTrueEcliptic(equinox=d)
    s_ecl = s.transform_to(ecl)
    m_ecl = m.transform_to(ecl)
    return abs(m_ecl.lon.rad - s_ecl.lon.rad - value) % (2*numpy.pi)


def previous_new_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(0.0,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_new_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(0.0,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_first_quarter_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi/2,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_first_quarter_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi/2,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_full_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_full_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(numpy.pi,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def previous_last_quarter_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(3*numpy.pi/2,), method='bounded',
                          bounds=(date-29.53, date*1.0),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)


def next_last_quarter_moon(date):
    date = Date(date)
    sol = minimize_scalar(_phases, args=(3*numpy.pi/2,), method='bounded',
                          bounds=(date*1.0, date+29.53),
                          options={'xatol': 1/86400.0})
    return Date(sol.x)
