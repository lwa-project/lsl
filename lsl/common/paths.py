# -*- coding: utf-8 -*

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
import imp

__version__ = '0.1'
__revision__ = '$Rev$'
__all__ = ['MODULE', 'DATA', 'MODULE_BUILD', 'DATA_BUILD']

modInfo = imp.find_module('lsl')
MODULE = os.path.abspath(modInfo[1])
DATA = os.path.join(MODULE, 'data')

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
    modInfoBuild = imp.find_module('lsl', [currentDir])
    MODULE_BUILD =  os.path.abspath(modInfoBuild[1])
    DATA_BUILD = os.path.join(MODULE_BUILD, 'data')
else:
    MODULE_BUILD = MODULE
    DATA_BUILD = DATA
