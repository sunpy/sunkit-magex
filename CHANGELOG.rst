1.1.0 (2025-05-07)
==================

Breaking Changes
----------------

- Increased the minimum required version of ``sunpy``  to v6.0.0. (`#63 <https://github.com/sunpy/sunkit-magex/pull/63>`__)
- Made streamtracer a hard dependency of sunkit-magex. (`#63 <https://github.com/sunpy/sunkit-magex/pull/63>`__)
- Deprecated the following pfss utils:

  - pfss.utils.fix_hmi_meta
  - pfss.utils.load_adapt

  `If you want to fix metadata, read this how-to guide. <https://docs.sunpy.org/en/latest/how_to/fix_map_metadata.html>`__
  `If you want to load ADAPT maps, read this example. <https://docs.sunpy.org/en/latest/generated/gallery/saving_and_loading_data/load_adapt_fits_into_map.html>`__ (`#63 <https://github.com/sunpy/sunkit-magex/pull/63>`__)


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
