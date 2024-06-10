**********
Installing
**********

If you need to install Python, you should start by following the :external+sunpy:ref:`sunpy core installation guide <sunpy-tutorial-installing>`.

``sunkit-magex`` can then be installed with conda:

.. code:: shell

    conda install -c conda-forge sunkit-magex

Alternatively ``sunkit-magex`` can be installed from PyPI using:

.. code:: shell

    pip install sunkit-magex

This will install sunkit_magex and all of its dependencies.
In addition to the core dependencies, there are optional dependencies (numba, streamtracer) that can improve code performance.
These can be installed with

.. code:: shell

    pip install sunkit-magex[performance]
