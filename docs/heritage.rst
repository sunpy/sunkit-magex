.. _code-heritage:

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

``outflowpy``
=============

The outflow field extrapolation (`sunkit_magex.pfss.outflow`) is adapted from `outflowpy <https://github.com/oekrice/outflowpy>`__, developed by `Oliver Rice <https://www.durham.ac.uk/staff/oliver-e-k-rice/>`__ at Durham University.
`outflowpy` implements the outflow field model, which generalises PFSS by including a solar wind outflow speed profile, and is itself based on `pfsspy`.
Both `outflowpy` and `sunkit-magex` are distributed under the GNU General Public License v3.

Citing
------

If you use `sunkit_magex.pfss.outflow` in work that results in publication, please cite Rice & Yeates (2021) at https://doi.org/10.3847/1538-4357/ac2c71.
A ready made bibtex entry is

.. code:: bibtex

  @article{Rice2021,
    doi = {10.3847/1538-4357/ac2c71},
    url = {https://doi.org/10.3847/1538-4357/ac2c71},
    year = {2021},
    publisher = {American Astronomical Society},
    volume = {923},
    number = {1},
    pages = {57},
    author = {Oliver E. K. Rice and Anthony R. Yeates},
    title = {Global Coronal Equilibria with Solar Wind Outflow},
    journal = {The Astrophysical Journal}
  }
