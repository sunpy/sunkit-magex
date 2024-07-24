*********************
Improving performance
*********************

numba
=====

`sunkit_magex.pfss` automatically detects an installation of `numba`_, which compiles some of the numerical code to speed up the pfss calculations.
To enable this simply `install numba`_  and use `sunkit_magex.pfss` as normal.

Streamline tracing
==================

`sunkit_magex.pfss` uses a complied tracer `sunkit_magex.pfss.tracing.PerformanceTracer` using the `streamtracer`_ package.

.. _numba: https://numba.pydata.org
.. _install numba: http://numba.pydata.org/numba-doc/latest/user/installing.html
.. _streamtracer: https://docs.sunpy.org/projects/streamtracer/en/stable/
