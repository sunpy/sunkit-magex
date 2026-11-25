import numpy as np
import pytest

import astropy.constants as const
import astropy.coordinates as acoord
import astropy.units as u

import sunkit_magex.pfss
from sunkit_magex.pfss import tracing
from sunkit_magex.pfss.fieldline import FieldLine, FieldLines
from sunkit_magex.pfss.grid import Grid
from sunkit_magex.pfss.outflow import (
    OutflowGrid,
    OutflowInput,
    _azimuthal_eigenmodes,
    _boundary_coefficient,
    _latitudinal_eigenmodes,
    _radial_g_function,
    _radial_h_function,
    outflow,
)
from sunkit_magex.pfss.output import OutflowOutput


@pytest.fixture
def outflow_grid():
    return OutflowGrid(ns=10, nphi=8, nr=5, rss=2.5)


@pytest.fixture
def outflow_result(dipole_map):
    # Same resolution as the `dipole_result` PFSS fixture in conftest.py,
    # so the two can be compared directly.
    nr = 10
    rss = 2.5
    input = OutflowInput(dipole_map, nr, rss, mf_constant=0.0)
    output = outflow(input)
    return input, output


def test_azimuthal_eigenmodes_shapes():
    nphi = 8
    dp = 2 * np.pi / nphi
    pc = np.linspace(0.5 * dp, 2 * np.pi - 0.5 * dp, nphi)

    ms, trigs = _azimuthal_eigenmodes(pc, dp)

    assert ms.shape == (nphi,)
    assert trigs.shape == (nphi, nphi)
    assert np.all(np.isfinite(ms))
    assert np.all(np.isfinite(trigs))
    # The m=0 (constant) mode should have eigenvalue zero.
    assert np.isclose(np.min(ms), 0.0, atol=1e-7)


def test_latitudinal_eigenmodes_shapes():
    ns = 10
    grid = Grid(ns=ns, nphi=8, nr=5, rss=2.5)

    ls, legs = _latitudinal_eigenmodes(0.0, grid.sc, grid.sg, grid.ds, ns)

    assert ls.shape == (ns,)
    assert legs.shape == (ns, ns)
    assert np.all(np.isfinite(ls))
    assert np.all(np.isfinite(legs))
    # Eigenvalues of this problem are non-negative.
    assert np.all(ls > -1e-8)


def test_boundary_coefficient_constant_mode():
    # A constant lower boundary projected onto the constant (all-ones)
    # basis function should just recover that constant.
    ns, nphi = 6, 4
    br0 = np.full((ns, nphi), 3.0)
    q = np.ones(ns)
    p = np.ones(nphi)
    assert np.isclose(_boundary_coefficient(br0, q, p), 3.0)


def test_radial_h_and_g_functions_zero_flow():
    nr = 8
    rss = 2.5
    dr = np.log(rss) / nr
    rcx = np.linspace(-0.5 * dr, np.log(rss) + 0.5 * dr, nr + 2)
    rg = np.linspace(0, np.log(rss), nr + 1)
    vcx = np.zeros(nr + 2)
    vdcx = np.zeros(nr + 2)

    hcx = _radial_h_function(2.0, rcx, vcx, vdcx, dr)
    gg = _radial_g_function(hcx, 2.0, rg)

    assert hcx.shape == (nr + 2,)
    assert gg.shape == (nr + 1,)
    assert np.all(np.isfinite(hcx))
    assert np.all(np.isfinite(gg))
    # G is normalised to 1 at the lower boundary by construction.
    assert np.isclose(gg[0], 1.0)


class TestOutflowGrid:
    def test_inherits_base_grid_properties(self, outflow_grid):
        base = Grid(ns=10, nphi=8, nr=5, rss=2.5)
        np.testing.assert_allclose(outflow_grid.rc, base.rc)
        np.testing.assert_allclose(outflow_grid.sc, base.sc)
        np.testing.assert_allclose(outflow_grid.pc, base.pc)
        np.testing.assert_allclose(outflow_grid.rg, base.rg)

    def test_eigenmode_shapes(self, outflow_grid):
        ns, nphi = outflow_grid.ns, outflow_grid.nphi
        assert outflow_grid.ms.shape == (nphi,)
        assert outflow_grid.trigs.shape == (nphi, nphi)
        assert outflow_grid.ls.shape == (nphi, ns)
        assert outflow_grid.legs.shape == (nphi, ns, ns)
        assert np.all(np.isfinite(outflow_grid.ms))
        assert np.all(np.isfinite(outflow_grid.ls))
        assert np.all(np.isfinite(outflow_grid.legs))

    def test_rcx_shape_and_bounds(self, outflow_grid):
        rcx = outflow_grid.rcx
        assert rcx.shape == (outflow_grid.nr + 2,)
        # rcx should be rc extended by one ghost cell either side.
        np.testing.assert_allclose(rcx[1:-1], outflow_grid.rc)


class TestOutflowInput:
    def test_zero_flow(self, dipole_map):
        input = OutflowInput(dipole_map, nr=5, rss=2.5, mf_constant=0.0)
        assert isinstance(input.grid, OutflowGrid)
        np.testing.assert_array_equal(input.vg, 0.0)
        np.testing.assert_array_equal(input.vcx, 0.0)
        np.testing.assert_array_equal(input.vdcx, 0.0)
        assert input.vg.shape == (input.grid.nr + 1,)
        assert input.vcx.shape == (input.grid.nr + 2,)
        assert input.vdcx.shape == (input.grid.nr + 2,)

    def test_default_profile(self, dipole_map):
        input = OutflowInput(dipole_map, nr=5, rss=2.5)
        assert np.all(np.isfinite(input.vg))
        assert np.all(np.isfinite(input.vcx))
        assert np.all(np.isfinite(input.vdcx))
        assert np.all(input.vg >= 0)
        assert np.all(input.vcx >= 0)
        assert input.vg.shape == (input.grid.nr + 1,)

    @pytest.mark.parametrize('polynomial_type', ['abs', 'clip', 'raw', 'smooth', 'smooth_monotonic'])
    def test_polynomial_profile(self, dipole_map, polynomial_type):
        input = OutflowInput(dipole_map, nr=5, rss=2.5,
                             polynomial_coeffs=[0.1, 0.2],
                             polynomial_type=polynomial_type)
        assert np.all(np.isfinite(input.vg))
        assert np.all(np.isfinite(input.vcx))
        assert np.all(np.isfinite(input.vdcx))

    def test_polynomial_profile_bad_type(self, dipole_map):
        with pytest.raises(ValueError, match='polynomial_type not recognised'):
            OutflowInput(dipole_map, nr=5, rss=2.5,
                        polynomial_coeffs=[0.1], polynomial_type='nonsense')

    def test_parker_profile_sound_speed(self, dipole_map):
        input = OutflowInput(dipole_map, nr=5, rss=2.5, mf_constant=1e-3,
                             sound_speed=150 * u.km / u.s)
        assert np.all(np.isfinite(input.vg))
        assert np.all(np.isfinite(input.vcx))
        assert np.all(np.isfinite(input.vdcx))

    def test_parker_profile_corona_temp(self, dipole_map):
        input = OutflowInput(dipole_map, nr=5, rss=2.5, mf_constant=1e-3,
                             corona_temp=1.5e6 * u.K)
        assert np.all(np.isfinite(input.vg))
        assert np.all(np.isfinite(input.vcx))
        assert np.all(np.isfinite(input.vdcx))

    def test_sound_speed_and_corona_temp_both_given_raises(self, dipole_map):
        with pytest.raises(ValueError, match='not both'):
            OutflowInput(dipole_map, nr=5, rss=2.5,
                        sound_speed=150 * u.km / u.s, corona_temp=1.5e6 * u.K)

    def test_sound_speed_requires_quantity(self, dipole_map):
        with pytest.raises(TypeError):
            OutflowInput(dipole_map, nr=5, rss=2.5, mf_constant=1e-3, sound_speed=150)

    def test_corona_temp_requires_quantity(self, dipole_map):
        with pytest.raises(TypeError):
            OutflowInput(dipole_map, nr=5, rss=2.5, mf_constant=1e-3, corona_temp=1.5e6)

    def test_sound_speed_wrong_unit(self, dipole_map):
        with pytest.raises(u.UnitsError):
            OutflowInput(dipole_map, nr=5, rss=2.5, mf_constant=1e-3, sound_speed=150 * u.K)


class TestOutflow:
    def test_returns_outflow_output(self, outflow_result):
        _, output = outflow_result
        assert isinstance(output, OutflowOutput)

    def test_default_profile_runs(self, dipole_map):
        input = OutflowInput(dipole_map, nr=6, rss=2.5)
        output = outflow(input)
        br, bth, bph = output.bc
        assert np.all(np.isfinite(br.value))
        assert np.all(np.isfinite(bth.value))
        assert np.all(np.isfinite(bph.value))

    def test_zero_flow_matches_pfss(self, dipole_map, outflow_result):
        # outflow() with mf_constant=0 solves the same physical problem as
        # pfss(), via a completely different numerical method (eigenfunction
        # expansion + radial finite-difference recursion, vs FFT + a
        # closed-form radial quadratic). They are independent
        # discretisation's of the same continuous problem, so they agree up
        # to the (shrinking, as resolution increases) discretisation error
        # of each method rather than to machine precision.
        # Checked convergence at several grid resolutions during
        # development. At this fixture's resolution (ns=30, nphi=20, nr=10)
        # the two methods agree to within ~1%.
        nr = 10
        rss = 2.5
        pfss_input = sunkit_magex.pfss.Input(dipole_map, nr, rss)
        pfss_output = sunkit_magex.pfss.pfss(pfss_input)

        _, outflow_output = outflow_result

        bg_pfss = pfss_output.bg.value
        bg_outflow = outflow_output.bg.value
        max_val = np.max(np.abs(bg_pfss))
        np.testing.assert_allclose(bg_outflow, bg_pfss, atol=0.02 * max_val)


class TestOutflowOutputTracing:
    @pytest.fixture(params=[tracing.PythonTracer(), tracing.PerformanceTracer()],
                    ids=['python', 'compiled'])
    def flines(self, outflow_result, request):
        tracer = request.param
        _, output = outflow_result
        out_frame = output.coordinate_frame
        seed = acoord.SkyCoord(2 * u.deg, -45 * u.deg, 1.01 * const.R_sun, frame=out_frame)
        return tracer.trace(seed, output)

    def test_field_lines(self, flines):
        # Proves the existing tracers (including the rust-backed
        # PerformanceTracer) work unmodified on OutflowOutput, since they
        # only ever touch output.bg / output.grid / output.coordinate_frame.
        assert isinstance(flines, FieldLines)
        assert isinstance(flines[0], FieldLine)

    def test_fline_in_bounds(self, flines):
        # The tracers terminate a field line by root-finding for the
        # (r - 1) * (r - rss) = 0 crossing, so the last point can land a
        # tiny floating-point amount outside [1, rss] depending on
        # platform/BLAS-dependent rounding in the field reconstruction.
        tol = 1e-9  # ~1 meter
        radius = flines[0].coords.radius
        np.testing.assert_array_less(const.R_sun * (1 - tol), radius)
        np.testing.assert_array_less(radius, 2.5 * const.R_sun * (1 + tol))
