C Extensions
============
.. warning::
	The modules documented before are compiled C extensions and are subject 
	to change.  It is recommended that users use the standard modules detailed
	in the previous sections over directly using these modules.

Go-Fast! Readers
----------------
.. automodule:: lsl.reader._gofast
   :members:
.. autoexception:: lsl.reader._gofast.EOFError
.. autoexception:: lsl.reader._gofast.SyncError

Power Spectral Density Calculation
----------------------------------

Linear Polarization
+++++++++++++++++++
.. automodule:: lsl.correlator._spec
   :members:

Stokes Parameters
+++++++++++++++++
.. automodule:: lsl.correlator._stokes

.. autofunction:: lsl.correlator._stokes.FPSD
.. autofunction:: lsl.correlator._stokes.PFBPSD
.. autofunction:: lsl.correlator._stokes.XEngine3

FX Correlator Core
------------------
.. automodule:: lsl.correlator._core
   :members:
.. autofunction:: lsl.correlator._core.FEngine
.. autofunction:: lsl.correlator._core.PFBEngine
.. autofunction:: lsl.correlator._core.XEngine2
.. autofunction:: lsl.correlator._core.XEngine3

DP-Style Signal Processing
--------------------------
.. automodule:: lsl.common._fir
   :members:
