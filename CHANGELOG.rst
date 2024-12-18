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
