import numpy as np
import pytest

from astropy.time import Time

from sunpy.map import Map

import sunkit_magex.pfss
from sunkit_magex.tests.helpers import get_dummy_map_from_header, get_fitsfile_from_header


@pytest.fixture
def zero_map():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((nphi, ns))
    header = sunkit_magex.pfss.utils.carr_cea_wcs_header(Time('1992-12-21'), br.shape)
    input_map = Map((br.T, header))

    input = sunkit_magex.pfss.Input(input_map, nr, rss)
    output = sunkit_magex.pfss.pfss(input)
    return input, output


@pytest.fixture
def dipole_map():
    ntheta = 30
    nphi = 20

    phi = np.linspace(0, 2 * np.pi, nphi)
    theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta)
    theta, phi = np.meshgrid(theta, phi)

    def dipole_Br(r, theta):
        return 2 * np.sin(theta) / r**3

    br = dipole_Br(1, theta)
    header = sunkit_magex.pfss.utils.carr_cea_wcs_header(Time('1992-12-21'), br.shape)
    header['bunit'] = 'nT'
    return Map((br.T, header))


@pytest.fixture
def dipole_result(dipole_map):
    nr = 10
    rss = 2.5

    input = sunkit_magex.pfss.Input(dipole_map, nr, rss)
    output = sunkit_magex.pfss.pfss(input)
    return input, output


@pytest.fixture
def dipole_result_closed(dipole_map):
    nr = 10
    rss = 2.5

    br_zeros = np.zeros(dipole_map.data.shape)
    header_zeros = sunkit_magex.pfss.utils.carr_cea_wcs_header(Time('2020-1-1'), br_zeros.shape[::-1])
    map_zeros = Map((br_zeros, header_zeros))

    input = sunkit_magex.pfss.Input(dipole_map, nr, rss, map_zeros)
    output = sunkit_magex.pfss.pfss(input)
    return input, output


@pytest.fixture
def adapt_test_file(tmp_path):
    """
    Return a fake
    """
    return get_fitsfile_from_header(
        "adapt_header.header",
        tmp_path / "adapt.fits",
        package="sunkit_magex.pfss.tests.data"
    )

@pytest.fixture
def adapt_map(adapt_test_file):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Map(adapt_test_file, hdus=0)


@pytest.fixture
def gong_map():
    return get_dummy_map_from_header(
        "gong_header.header",
        package="sunkit_magex.pfss.tests.data"
    )
