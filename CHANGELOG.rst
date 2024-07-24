<<<<<<<< HEAD:CHANGELOG.rst
1.0.0 (2024-05-31)
==================

This is the first release and aims to keep the API the same from `pfsspy` to
`sunkit_magex`.  The main difference is that you will need to replace the
import like so:

  .. code-block:: python

    # Before
    import pfsspy
    from pfsspy import tracing

    # After
    from sunkit_magex import pfss as pfsspy
    from sunkit_magex.pfss import tracing

The main changes from the previous release of `pfsspy` are as follows:

* The ``GongSynopticMap`` class has moved into `sunpy`, note that the new
  class in ``sunpy`` is slightly different in that it does not modify the
  header.
|||||||| parent of cc08f79 (Handle old functions and update examples):docs/changelog.rst
.. _sunkit-magex-changelog:

**************
Full Changelog
**************

1.0.0 (2024-05-31)
==================

This is the first release and aims to keep the API the same from `pfsspy` to
`sunkit_magex`.  The main difference is that you will need to replace the
import like so:

  .. code-block:: python

    # Before
    import pfsspy
    from pfsspy import tracing

    # After
    from sunkit_magex import pfss as pfsspy
    from sunkit_magex.pfss import tracing

The main changes from the previous release of `pfsspy` are as follows:

* The ``GongSynopticMap`` class has moved into `sunpy`, note that the new
  class in ``sunpy`` is slightly different in that it does not modify the
  header.
========
.. _sunkit-magex-changelog:

**************
Full Changelog
**************

.. changelog::
   :towncrier: ../
   :towncrier-skip-if-empty:
   :changelog_file: ../CHANGELOG.rst
>>>>>>>> cc08f79 (Handle old functions and update examples):docs/changelog.rst
