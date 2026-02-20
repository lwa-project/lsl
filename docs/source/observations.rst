Defining Single Station Observations and Observation Metadata
=============================================================

Session Definition Files
------------------------
.. automodule:: lsl.common.sdf

Session Structure
+++++++++++++++++
:mod:`lsl.common.sdf` provides means to represent a set of observations as Python objects.  For each
:class:`lsl.common.sdf.Project`, there is:

  1) An observer (:class:`lsl.common.sdf.Observer`)
  2) The project office comments (:class:`lsl.common.sdf.ProjectOffice`)
  3) A single session that defines the SDF (:class:`lsl.common.sdf.Session`)

The session contains one or more observerions (:class:`lsl.common.sdf.Observation`).  Each observing mode supported
by the LWA is sub-classed (see below).

.. autoclass:: lsl.common.sdf.Project
   :members:
.. autoclass:: lsl.common.sdf.Observer
   :members:
.. autoclass:: lsl.common.sdf.ProjectOffice
   :members:
.. autoclass:: lsl.common.sdf.Session
   :members:
.. autoclass:: lsl.common.sdf.Observation
   :members:

Observing Modes
+++++++++++++++
.. autoclass:: lsl.common.sdf.TBT
   :members:
.. autoclass:: lsl.common.sdf.TBS
   :members:
.. autoclass:: lsl.common.sdf.DRX
   :members:
.. autoclass:: lsl.common.sdf.Solar
   :members:
.. autoclass:: lsl.common.sdf.Jovian
   :members:
.. autoclass:: lsl.common.sdf.Stepped
   :members:
.. autoclass:: lsl.common.sdf.BeamStep
   :members:

MCS Metadata Tarball Utilities
------------------------------
.. automodule:: lsl.common.metabundle
   :members:

MCS Supporting Functions
------------------------
.. automodule:: lsl.common.mcs
   :members:

