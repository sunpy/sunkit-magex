*************
Code Heritage
*************

``pfsspy``
==========

As of initial release, the main component of ``sunkit-magex`` (`sunkit_magex.pfss`) is forked from `pfsspy` package (archived as of August 24, 2023).
`pfsspy` implemented the Potential Field Source Surface (PFSS) model, a widely used method to extrapolate the magnetic fields of the solar corona.
`pfsspy` was developed, integrated with `sunpy` and maintained by `David Stansby <https://www.davidstansby.com/>`__, based on an `original PFSS implementation <https://doi.org/10.5281/zenodo.1472183>`__  by `Anthony Yeates <https://www.maths.dur.ac.uk/users/anthony.yeates/>`__.

Details of the numerical methods underlying the solver can be found in :any:`numerical_methods_pfss/index`.

Citing
------

If you use `sunkit_magex.pfss` in work that results in publication, please cite the original ``pfsspy`` Journal of Open Source Software paper at https://doi.org/10.21105/joss.02732.
A ready made bibtex entry is

.. code:: bibtex

  @article{Stansby2020,
    doi = {10.21105/joss.02732},
    url = {https://doi.org/10.21105/joss.02732},
    year = {2020},
    publisher = {The Open Journal},
    volume = {5},
    number = {54},
    pages = {2732},
    author = {David Stansby and Anthony Yeates and Samuel T. Badman},
    title = {pfsspy: A Python package for potential field source surface modelling},
    journal = {Journal of Open Source Software}
  }
