.. _sunkit-magex-changelog:

**************
Full Changelog
**************

1.0.0 (2024-05-31)
==================

* First release that keeps API the same from `pfsspy` to `sunkit-magex`.
  The main difference is that you will need to replace the import like so:

  .. code-block:: python

    # Before
    import pfsspy
    from pfsspy import tracing

    # After
    from sunkit_magex import pfss as pfsspy
    from sunkit_magex.pfss import tracing

The main breaking changes is that the old sunpy Map classes were removed and added to sunpy.
You will need to make sure you have sunpy 6.0.0 installed.
