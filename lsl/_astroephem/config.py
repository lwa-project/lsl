"""
Simple configuration module that helps toggle how the __repr__ and __str__
methods in this class of this module behave.  Setting 'PYEPHEM_REPR' to True
causes them work similar to their counterparts in PyEphem.  False cause them to
work more like their counterparts in AstroPy.
"""

__all__ = ['PYEPHEM_REPR',]

PYEPHEM_REPR = True
