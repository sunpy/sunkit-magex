*********************
Improving performance
*********************

numba
=====

`sunkit_magex.pfss` automatically detects an installation of `numba`_, which compiles some of the numerical code to speed up the pfss calculations.
To enable this simply `install numba`_  and use `sunkit_magex.pfss` as normal.

Streamline tracing
==================

`sunkit_magex.pfss` has two streamline tracers: a pure python `sunkit_magex.pfss.tracing.PythonTracer` and a complied tracer `sunkit_magex.pfss.tracing.FortranTracer`.
The FORTRAN version is significantly faster, using the `streamtracer`_ package.

.. _numba: https://numba.pydata.org
.. _install numba: http://numba.pydata.org/numba-doc/latest/user/installing.html
.. _streamtracer: https://github.com/sunpy/streamtracer
