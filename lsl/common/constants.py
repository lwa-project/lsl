# -*- coding: utf-8 -*

"""
Module that stores various useful constants in one convenient location.
The constants defined in this file are:

c
  the speed of light in m/s

deg2rad
  the conversion factor for degrees to radians

tpi
  :math: 2 \\pi \\sqrt{-1}
"""

import math

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['c', 'deg2rad', 'tpi', '__version__', '__revision__', '__all__']

c = 2.9979245800e8			# speed of light in m/s
deg2rad = math.pi / 180.0	# degrees to radians conversion
tpi = 2j*math.pi			# two pi i