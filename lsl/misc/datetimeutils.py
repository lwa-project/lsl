"""
Module that provides backwards compatible naive UTC datetime generator
functions to deal with Python's move to default timezone-aware objects
in Python 3.12+.

.. versionadded:: 4.0.0
"""

from datetime import datetime, timezone

__all__ = ['utcfromtimestamp', 'utcnow']


def utcfromtimestamp(value):
    """
    Replacement for `datetime.datetime.utcfromtimestamp()` that avoids the
    deprecation warning introduced in Python 3.12.
    """
    
    dt = datetime.fromtimestamp(value, tz=timezone.utc)
    return dt.replace(tzinfo=None)


def utcnow():
    """
    Replacement for `datetime.datetime.utcnow()` that avoids the deprecation
    warning introduced in Python 3.12.
    """
    
    dt = datetime.now(tz=timezone.utc)
    return dt.replace(tzinfo=None)
