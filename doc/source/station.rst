Station Meta-Data
=================

Station-Level Data
------------------
.. automodule:: lsl.common.stations
   :members:

Stand-Level Data
----------------
.. automodule:: lsl.correlator.uvUtils
   :members: validateStand

Positions and Cable Model
+++++++++++++++++++++++++
.. autofunction:: lsl.correlator.uvUtils.getXYZ
.. autofunction:: lsl.correlator.uvUtils.getRelativeXYZ
.. autoclass:: lsl.correlator.uvUtils.PositionCache
   :members:

Cable Model and Delays
++++++++++++++++++++++
.. autofunction:: lsl.correlator.uvUtils.cableDelay
.. autofunction:: lsl.correlator.uvUtils.cableAttenuation
.. autoclass:: lsl.correlator.uvUtils.CableCache
   :members:
.. autofunction:: lsl.correlator.uvUtils.signalDelay
.. autoclass:: lsl.correlator.uvUtils.SignalCache
   :show-inheritance:
   :members:

Exceptions
++++++++++
.. autoexception:: lsl.correlator.uvUtils.uvUtilsError
   :members:

