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

import os
import importlib.util

from lsl.config import LSL_CONFIG

from lsl.misc import telemetry
telemetry.track_module()


__version__ = '0.3'
__all__ = ['MODULE', 'DATA', 'WISDOM', 'MODULE_BUILD', 'DATA_BUILD', 'WISDOM_BUILD']


modInfo = importlib.util.find_spec('lsl')

#: Absolute path to the LSL intall location
MODULE = os.path.abspath(modInfo.origin)
MODULE = os.path.dirname(MODULE)

#: Absolute path to the data directory where data files for LSL are stored
DATA = os.path.join(MODULE, 'data')

#: Absolute path to where the LSL-specific FFTW wisdom file is stored
WISDOM = LSL_CONFIG.dirname


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
    modInfoBuild = importlib.util.find_spec('lsl', [currentDir])
    MODULE_BUILD =  os.path.abspath(modInfoBuild.origin)
    MODULE_BUILD =  os.path.dirname(MODULE_BUILD)
    DATA_BUILD = os.path.join(MODULE_BUILD, 'data')
else:
    MODULE_BUILD = MODULE
    DATA_BUILD = DATA
WISDOM_BUILD = WISDOM
