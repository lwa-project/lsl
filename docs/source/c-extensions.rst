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
   :members:

FX Correlator Core
------------------
.. automodule:: lsl.correlator._core
   :members:

DP-Style Signal Processing
--------------------------
.. automodule:: lsl.common._fir
   :members:
