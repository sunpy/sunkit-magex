import numpy as np
import pytest

import astropy.constants as const
import astropy.coordinates as coord
import astropy.units as u

import sunkit_magex.pfss
from sunkit_magex.pfss import tracing


@pytest.fixture(params=[tracing.PythonTracer(),
                        tracing.PerformanceTracer()],
                ids=['python', 'compiled'])
def flines(dipole_result, request):
    tracer = request.param
    _, out = dipole_result
    out_frame = out.coordinate_frame

    seed = coord.SkyCoord(2*u.deg, -45*u.deg, 1.01*const.R_sun, frame=out_frame)
    return tracer.trace(seed, out)


def test_field_lines(flines):
    assert isinstance(flines[0],
                      sunkit_magex.pfss.fieldline.FieldLine)
    assert isinstance(flines.open_field_lines.solar_feet,
                      coord.SkyCoord)
    assert isinstance(flines.open_field_lines.source_surface_feet,
                      coord.SkyCoord)
    assert isinstance(flines.polarities,
                      np.ndarray)


def test_fline_in_bounds(flines):
    assert np.all(flines[0].coords.radius >= const.R_sun)
    assert np.all(flines[0].coords.radius <= 2.5 * const.R_sun)


def test_fline_step_size(dipole_result):
    input, out = dipole_result
    out_frame = out.coordinate_frame
    seed = coord.SkyCoord(2*u.deg, -45*u.deg, 1.01*const.R_sun,
                          frame=out_frame)

    tracer = tracing.PerformanceTracer(step_size=0.5)
    flines = tracer.trace(seed, out)
    assert out.grid.nr == 10
    # With a step size of 0.5, this should be ~20
    assert len(flines[0]) == 21

    tracer = tracing.PerformanceTracer(step_size=0.2)
    flines = tracer.trace(seed, out)
    assert out.grid.nr == 10
    # With a step size of 0.2, this should be ~50
    assert len(flines[0]) == 52


def test_rot_warning(dipole_result):
    tracer = tracing.PerformanceTracer(max_steps=2)
    input, out = dipole_result
    out_frame = out.coordinate_frame
    seed = coord.SkyCoord(0*u.deg, -45*u.deg, 1.01*const.R_sun,
                          frame=out_frame)

    with pytest.warns(UserWarning, match='ran out of steps'):
        tracer.trace(seed, out)
