import astropy.io

import sunpy.map

import sunkit_magex.pfss.map


def test_adapt_map(adapt_test_file):
    with astropy.io.fits.open(adapt_test_file) as adapt_fits:
        for map_slice in adapt_fits[0].data:
            m = sunpy.map.Map((map_slice, adapt_fits[0].header))
            assert isinstance(m, sunkit_magex.pfss.map.ADAPTMap)
