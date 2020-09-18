"""
Utilities for helping to build over extensions that want complex integer support.
"""

import os

__version__ = "0.1"
__all__ = ['get_include',]

def get_include():
    """
    Return the directory that contains the LSL complex integer types \*.h 
    header files.
    
    Extension modules that need to compile against the LSL complex integer
    types should use this function to locate the appropriate include directory.
    """
    
    return os.path.dirname(os.path.abspath(__file__))
