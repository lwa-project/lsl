RF Antenna Parameters
=====================
.. automodule:: lsl.misc.rfutils
   :members:

Simulated Stand Response
========================
.. automodule:: lsl.sim.necutils
   :members:

Ionospheric and Geomagnetic Field Models
========================================
.. versionadded:: 1.1.0

.. automodule:: lsl.misc.ionosphere
   :members:

Pulsar Scattering Removal
=========================
.. automodule:: lsl.misc.scattering
   :members:

File Cache Management
=====================
.. versionadded:: 2.1.3

.. automodule:: lsl.misc.file_cache
   :members:

File Locking
============
.. versionadded:: 2.1.3

.. automodule:: lsl.misc.file_lock
   :members:

FFTW Wisdom Management
======================
.. versionadded:: 1.0.1

.. automodule:: lsl.misc.wisdom
   :members:

UTC Datetime Utilities
======================
.. versionadded:: 4.0.0

Python 3.12 deprecated ``datetime.utcfromtimestamp()`` and ``datetime.utcnow()``
as part of the transition toward timezone-aware datetime objects.  The
``lsl.misc.datetimeutils`` module provides drop-in replacements that avoid the
deprecation warnings while maintaining compatibility with existing code that
uses naive UTC datetimes.

.. automodule:: lsl.misc.datetimeutils
   :members:

Conversion Helper Functions for argparse
========================================
.. versionadded:: 1.2.4

.. automodule:: lsl.misc.parser
   :members:
