import pathlib
from datetime import timedelta

import numpy as np
import pytest

import astropy.constants as const
import astropy.coordinates as coord
import astropy.units as u
from astropy.tests.helper import quantity_allclose

import sunpy.map
import sunpy.util.exceptions
from sunpy.coordinates import frames

import sunkit_magex.pfss
import sunkit_magex.pfss.coords
from sunkit_magex.pfss import tracing

R_sun = const.R_sun
test_data = pathlib.Path(__file__).parent / 'data'


@pytest.mark.array_compare(reference_dir="arraydiff_reference")
def test_pfss(gong_map):
    # Regression test to check that the output of pfss doesn't change
    m = sunpy.map.Map(gong_map)
    # Resample to lower res for easier comparisons
    m = m.resample([30, 15] * u.pix)

    expected = np.loadtxt(test_data / 'br_in.txt')
    m = sunpy.map.Map(expected, m.meta)

    pfss_in = sunkit_magex.pfss.Input(m, 50, 2)
    pfss_out = sunkit_magex.pfss.pfss(pfss_in)

    br = pfss_out.source_surface_br.data
    assert br.shape == m.data.shape
    return br


@pytest.mark.array_compare(reference_dir="arraydiff_reference")
def test_pfss_br_outer(gong_map):
    # Regression test to check that the output of pfss doesn't change
    # when an equivalent custom outer boundary condition is provided
    m = sunpy.map.Map(gong_map)
    # Resample to lower res for easier comparisons
    m = m.resample([30, 15] * u.pix)

    expected = np.loadtxt(test_data / 'br_in.txt')
    m = sunpy.map.Map(expected, m.meta)

    # Get a source surface
    br_outer_in = sunkit_magex.pfss.Input(m, 50, 2)
    br_outer_out = sunkit_magex.pfss.pfss(br_outer_in)
    br_outer = br_outer_out.source_surface_br

    # Use calculated source surface as outer boundary condition
    pfss_in = sunkit_magex.pfss.Input(m, 50, 2, br_outer)
    pfss_out = sunkit_magex.pfss.pfss(pfss_in)

    br = pfss_out.source_surface_br.data
    assert br.shape == m.data.shape
    return br


def test_bunit(gong_map):
    # Regression test to check that the output of pfss doesn't change
    m = sunpy.map.Map(gong_map)
    pfss_in = sunkit_magex.pfss.Input(m, 2, 2)
    pfss_out = sunkit_magex.pfss.pfss(pfss_in)
    assert pfss_out.bunit == m.unit

    pfss_out.input_map.meta['bunit'] = 'notaunit'
    with pytest.warns(UserWarning, match='Could not parse unit string "notaunit"'):
        assert pfss_out.bunit == u.dimensionless_unscaled

    pfss_out.input_map.meta.pop('bunit')
    assert pfss_out.bunit == u.dimensionless_unscaled


def test_expansion_factor(dipole_result):
    _, out = dipole_result
    out_frame = out.coordinate_frame

    tracer = tracing.PythonTracer()
    seed = coord.SkyCoord(0 * u.deg, 80 * u.deg, 1.1 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    assert field_line.expansion_factor > 1

    seed = coord.SkyCoord(0 * u.deg, -80 * u.deg, 1.1 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    assert field_line.expansion_factor > 1

    # This is a closed field line
    seed = coord.SkyCoord(0 * u.deg, 0 * u.deg, 1.1 * R_sun, frame=out_frame)
    eq_field_line = tracer.trace(seed, out)[0]
    assert np.isnan(eq_field_line.expansion_factor)

    # Check that a field line near the equator has a bigger expansion
    # factor than one near the pole
    seed = coord.SkyCoord(0 * u.deg, 10 * u.deg, 2.49 * R_sun, frame=out_frame)
    pil_field_line = tracer.trace(seed, out)[0]
    assert pil_field_line.expansion_factor > field_line.expansion_factor


def test_field_line_polarity(dipole_result):
    _, out = dipole_result
    out_frame = out.coordinate_frame

    tracer = tracing.PythonTracer()
    seed = coord.SkyCoord(0 * u.deg, 90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)
    assert field_line[0].polarity == 1

    seed = coord.SkyCoord(0 * u.deg, -90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)
    assert field_line[0].polarity == -1

    # This is a closed field line
    seed = coord.SkyCoord(0 * u.deg, 0 * u.deg, 1.01 * R_sun, frame=out_frame)
    eq_field_line = tracer.trace(seed, out)[0]
    assert eq_field_line.polarity == 0


def test_footpoints(dipole_result):
    _, out = dipole_result
    out_frame = out.coordinate_frame

    tracer = tracing.PythonTracer(atol=1e-8, rtol=1e-8)

    def check_radius(coord, r):
        coord.representation_type = 'spherical'
        assert quantity_allclose(coord.radius, r)

    def check_open_fline(fline):
        check_radius(fline.solar_footpoint, const.R_sun)
        check_radius(fline.source_surface_footpoint, 2.5 * const.R_sun)

    seed = coord.SkyCoord(0 * u.deg, 90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    check_open_fline(field_line)

    seed = coord.SkyCoord(0 * u.deg, -90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    check_open_fline(field_line)

    seed = coord.SkyCoord(0 * u.deg, 0 * u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    check_radius(field_line.solar_footpoint, const.R_sun)
    check_radius(field_line.source_surface_footpoint, const.R_sun)


def test_shape(zero_map):
    # Test output map shapes
    _input, out = zero_map
    nr = _input.grid.nr
    nphi = _input.grid.nphi
    ns = _input.grid.ns

    out = sunkit_magex.pfss.pfss(_input)
    alr, als, alp = out._al
    for comp in (alr, als, alp):
        assert np.all(comp == 0)

    assert alr.shape == (nphi + 1, ns + 1, nr)
    assert als.shape == (nphi + 1, ns, nr + 1)
    assert alp.shape == (nphi, ns + 1, nr + 1)

    br, bs, bp = out.bc
    for comp in (br, bs, bp):
        assert np.all(comp == 0)

    assert br.shape == (nphi, ns, nr + 1)
    assert bs.shape == (nphi, ns + 1, nr)
    assert bp.shape == (nphi + 1, ns, nr)

    bg = out.bg
    assert np.all(bg == 0)
    assert bg.shape == (nphi + 1, ns + 1, nr + 1, 3)
    assert (bg[0, ...] == bg[-1, ...]).all()


def test_wrong_projection_error(dipole_map):
    dipole_map.meta['ctype1'] = 'HGLN-CAR'
    with pytest.raises(ValueError, match='must be CEA'):
        sunkit_magex.pfss.Input(dipole_map, 5, 2.5)


def test_nan_value(dipole_map):
    dipole_map.data[0, 0] = np.nan
    with pytest.raises(
            ValueError, match='At least one value in the input is NaN'):
        sunkit_magex.pfss.Input(dipole_map, 5, 2.5)


def test_non_map_input():
    with pytest.raises(ValueError, match='br must be a sunpy Map'):
        sunkit_magex.pfss.Input(np.random.rand(2, 2), 1, 1)


def test_bvec_interpolator(dipole_result):
    _, out = dipole_result
    test_coord = coord.SkyCoord(
        x=np.array([1 / out.grid.rss, 0.5, 0.75, 1]) * out.grid.rss * u.R_sun,
        y=np.zeros(4) * u.R_sun,
        z=np.zeros(4) * u.R_sun,
        frame=out.coordinate_frame,
        representation_type="cartesian"
    )
    b_cart = out.get_bvec(test_coord, out_type="cartesian")
    b_sph = out.get_bvec(test_coord, out_type="spherical")

    # Check the output shape matches is [N,3] where
    # N is the length of test_coord, the unit is nT
    assert b_cart.shape == (len(test_coord), 3)
    assert b_sph.shape == (len(test_coord), 3)
    assert b_cart.unit == u.nT

    # Test interpolation matches the analytic expectation
    # of the boundary condition at the solar surface.
    # The first test coordinate is located at [1,0,0] Rs
    # i.e. on the solar surface at the equator of the dipole,
    # where we expect Br to vanish.
    assert np.isclose(b_cart[0, 0], 0.0)

    # The test coordinates are all located at
    # rhat = [1,0,0], thetahat = [0,0,1], phihat = [0,1,0]
    # Use this to sanity check the rotation to spherical coordinates
    assert np.all(b_sph[:, 0] == b_cart[:, 0])
    assert np.all(b_sph[:, 1] == -b_cart[:, 2])
    assert np.all(b_sph[:, 2] == b_cart[:, 1])

    # Test warnings and errors are raised correctly
    with pytest.raises(ValueError, match="coords must be of type astropy.coordinates.SkyCoord"):
        out.get_bvec([])
    with pytest.raises(ValueError, match="out_type must be 'cartesian' or 'spherical'"):
        out.get_bvec(test_coord, out_type='')
    # Create test coord with different datetime
    wrong_datetime = coord.SkyCoord(
        x=np.array([1]) * out.grid.rss * u.R_sun,
        y=np.zeros(1) * u.R_sun,
        z=np.zeros(1) * u.R_sun,
        frame=frames.HeliographicCarrington,
        obstime=out.dtime.datetime + timedelta(days=1),
        observer='Earth',
        representation_type="cartesian"
    )
    with pytest.warns(UserWarning) as record:
        out.get_bvec(wrong_datetime)
    assert record[0].message.args[0] == "The obstime of one of more input coordinates do not match the pfss model obstime."
