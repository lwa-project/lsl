"""
Module to set up path information for the LWA Software Library.  Two 
variables are defined by this module:
  * MODULE
  * DATA

MODULE
    specifies the absolute path to the module

DATA
    the absolute path to the data directory where data files are stored
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
try:
    from imp import find_module as find_spec
except ImportError:
    from importlib.util import find_spec
    
from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.2'
__all__ = ['MODULE', 'DATA', 'WISDOM', 'MODULE_BUILD', 'DATA_BUILD', 'WISDOM_BUILD']


modInfo = find_spec('lsl')

#: Absolute path to the LSL intall location
MODULE = os.path.abspath(modInfo[1])

#: Absolute path to the data directory where data files for LSL are stored
DATA = os.path.join(MODULE, 'data')

#: Absolute path to where the LSL-specific FFTW wisdom file is stored
WISDOM = os.path.join(os.path.expanduser('~'), '.lsl')


# If we seem to be in the building directory, make the module and 
# data build paths point to the right place.  This is done so that 
# the testing scripts run via:
#	python setup.py test
# access the right files and tests.
# 
# If we don't seem to be in the 
# source directory, MODULE_BUILD points to module and DATA_BUILD 
# points to data.
currentDir = os.path.abspath(os.getcwd())
if os.path.exists(os.path.join(currentDir, 'setup.py')) and os.path.exists(os.path.join(currentDir, 'lsl')):
    modInfoBuild = find_spec('lsl', [currentDir])
    MODULE_BUILD =  os.path.abspath(modInfoBuild[1])
    DATA_BUILD = os.path.join(MODULE_BUILD, 'data')
else:
    MODULE_BUILD = MODULE
    DATA_BUILD = DATA
WISDOM_BUILD = WISDOM
