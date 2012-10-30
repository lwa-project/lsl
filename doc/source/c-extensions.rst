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
.. autoexception:: lsl.reader._gofast.eofError
.. autoexception:: lsl.reader._gofast.syncError

Power Spectral Density Calculation
----------------------------------

Linear Polarization
+++++++++++++++++++
.. automodule:: lsl.correlator._spec
   :members:

Stokes Parameters
+++++++++++++++++
.. automodule:: lsl.correlator._stokes

.. autofunction:: lsl.correlator._stokes.FPSDR2
.. autofunction:: lsl.correlator._stokes.FPSDR3
.. autofunction:: lsl.correlator._stokes.FPSDC2
.. autofunction:: lsl.correlator._stokes.FPSDC3

FX Correlator Core
------------------
.. automodule:: lsl.correlator._core
   :members:
.. autofunction:: lsl.correlator._stokes.XEngine2

