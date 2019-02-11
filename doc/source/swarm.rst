LWA Swarm Observations
======================

Interferometer Definition Files
-------------------------------
.. automodule:: lsl.common.idf

Run Structure
+++++++++++++
Similar to the :mod:`lsl.common.sdf` module, :mod:`lsl.common.idf` provides a means to represent a set of scans using the LWA Swarm interferometer as Python objects.  For each
:class:`lsl.common.idf.Project`, there is:

  1) An observer (:class:`lsl.common.idf.Observer`)
  2) The project office comments (:class:`lsl.common.idf.ProjectOffice`)
  3) A single run that defines the IDF (:class:`lsl.common.idf.Run`)
  
The run contains one or more scans (:class:`lsl.common.idf.Scan`).  Each observing mode supported
by the LWA Swarm is sub-classed (see below).

.. autoclass:: lsl.common.idf.Project
   :members:
.. autoclass:: lsl.common.idf.Observer
   :members:
.. autoclass:: lsl.common.idf.ProjectOffice
   :members:
.. autoclass:: lsl.common.idf.Run
   :members:
.. autoclass:: lsl.common.idf.Scan
   :members:

Observing Modes
+++++++++++++++
.. autoclass:: lsl.common.idf.DRX
   :members:
.. autoclass:: lsl.common.idf.Solar
   :members:
.. autoclass:: lsl.common.idf.Jovian
   :members:
