"""
Code for calculating an outflow field extrapolation.

An outflow field generalises the PFSS model by including a prescribed,
purely radial solar wind outflow speed profile :math:`v(r)` in the
force-free-like solve. Setting :math:`v(r) \\equiv 0` everywhere recovers
the PFSS solution exactly.

The numerical method implemented here (the eigenfunction expansion in
:func:`outflow`, and the ``OutflowGrid``/``OutflowInput`` velocity profile
machinery) is adapted from `outflowpy <https://github.com/oekrice/outflowpy>`__
(Copyright (C) Oliver E. Rice, Durham University), which is itself based on
`pfsspy <https://github.com/dstansby/pfsspy>`__, the package
`sunkit_magex.pfss` is itself forked from (see :ref:`code-heritage`). Both
``outflowpy`` and ``sunkit_magex`` are distributed under the GNU General
Public License v3.

See Rice & Yeates, 2021, *Global Coronal Equilibria with Solar Wind
Outflow*, ApJ 923, 57, https://doi.org/10.3847/1538-4357/ac2c71 for the
underlying physics and numerical method.
"""
from pathlib import Path

import numpy as np
from scipy import interpolate
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize_scalar, root_scalar

import astropy.units as u

from sunkit_magex.pfss import OutflowOutput
from sunkit_magex.pfss.grid import Grid
from sunkit_magex.pfss.input import Input

HAS_NUMBA = False
try:
    import numba
    HAS_NUMBA = True
except Exception:
    pass

__all__ = ['OutflowGrid', 'OutflowInput', 'outflow']

_DATA_DIR = Path(__file__).resolve().parent / 'data'


def _azimuthal_eigenmodes(pc, dp):
    """
    Eigenvalues and eigenvectors of the azimuthal (phi) direction.

    Unlike PFSS (where the radial ODE has a closed-form solution and a
    plain FFT can be used), the outflow radial ODE has no closed form. So
    the azimuthal basis has to be built explicitly by combining the
    eigenvectors of the equivalent sine and cosine tridiagonal eigenvalue
    problems.
    """
    num = len(pc)
    dvals = 2 * np.ones(num)
    evals = -np.ones(num - 1)

    dvals[0] = 1
    dvals[-1] = 1
    w1, v1 = eigh_tridiagonal(dvals, evals)
    ms1 = np.sqrt(np.abs(w1)) / dp

    dvals[0] = 3
    dvals[-1] = 3
    w2, v2 = eigh_tridiagonal(dvals, evals)
    ms2 = np.sqrt(np.abs(w2)) / dp

    v = np.zeros_like(v1)
    ms = np.zeros(num)
    for i in range(num):
        if i % 2 == 0:
            ms[i] = ms1[i]
            v[:, i] = v1[:, i]
        else:
            ms[i] = ms2[i]
            v[:, i] = v2[:, i]

    return ms, v


def _latitudinal_eigenmodes(m, sc, sg, ds, ns):
    """
    Eigenvalues and eigenvectors (Legendre-like functions) of the
    latitudinal (:math:`s = \\cos\\theta`) direction, for a given azimuthal
    mode ``m``.
    """
    sigc = np.sqrt(1 - sc**2)
    sigs = np.sqrt(1 - sg**2)
    evals = -sigs[1:-1]**2
    dvals = (sigs[1:]**2 + sigs[:-1]**2
             - (ds * m**2 / sigc) * (np.arcsin(sg[:-1]) - np.arcsin(sg[1:])))
    w, v = eigh_tridiagonal(dvals, evals)

    for j in range(ns):
        v[:, j] = v[:, j] / np.sign(v[0, j] + 1e-10)

    return w / ds**2, v


def _boundary_coefficient(br0, q, p):
    """
    Boundary-fitting coefficient for a single (l, m) mode: the orthogonal
    projection of the lower boundary condition ``br0`` onto the basis
    function formed from the latitudinal eigenvector ``q`` and azimuthal
    eigenvector ``p``.
    """
    lhs = np.sum(br0 * q[:, np.newaxis] * p[np.newaxis, :])
    rhs = np.sum((q * q)[:, np.newaxis] * (p * p)[np.newaxis, :])
    return lhs / rhs


def _radial_h_function(l_mode, rcx, vcx, vdcx, dr):
    """
    Radial function H(rho), found via a finite-difference recursion
    (shooting from the source surface inward) that includes the outflow
    speed terms ``vcx`` (v(r)) and ``vdcx`` (dv/drho). Normalised to
    satisfy the lower boundary condition.
    """
    hcx = np.zeros(len(rcx))
    hcx[-1] = 1.0
    hcx[-2] = 0.0

    for i in range(len(rcx) - 3, -1, -1):
        A = 1.0
        B = 3 - vcx[i + 1] * np.exp(rcx[i + 1])
        C = 2 - l_mode - 3 * vcx[i + 1] * np.exp(rcx[i + 1]) - vdcx[i + 1] * np.exp(rcx[i + 1])
        top = hcx[i + 1] * (2 * A / dr**2 - C) + hcx[i + 2] * (-A / dr**2 - B / (2 * dr))
        bottom = A / dr**2 - B / (2 * dr)
        hcx[i] = top / bottom

    grad = ((hcx[1] * np.exp(rcx[1]) - hcx[0] * np.exp(rcx[0])) / dr
            - 0.5 * (vcx[0] * np.exp(rcx[0]) * hcx[0] + vcx[1] * np.exp(rcx[1]) * hcx[1]))

    return hcx / grad


def _radial_g_function(hcx, l_mode, rg):
    """
    Radial function G(rho), derived from H via the divergence-free
    coupling between the two.
    """
    gs = np.zeros(len(rg))
    gs[0] = 1.0
    for i in range(1, len(gs)):
        gs[i] = np.exp(-2 * rg[i]) * (
            0.5 * l_mode * hcx[i] * (np.exp(2 * rg[i]) - np.exp(2 * rg[i - 1]))
            + gs[i - 1] * np.exp(2 * rg[i - 1]))
    return gs


def _radial_h_g_batch(l_modes, rcx, vcx, vdcx, dr, rg):
    """
    Compute the H and G radial functions for every l-eigenvalue in
    ``l_modes`` in a single call, instead of one Python-level call per
    mode. When jitted, nested calls to `_radial_h_function`/
    `_radial_g_function` from within this function compile to direct
    calls with no Python<->numba boundary crossing, avoiding the
    per-call dispatch overhead that adds up when there are many
    non-negligible modes (e.g. for realistically-scaled magnetogram
    data, see the performance note on :func:`outflow`).
    """
    n = len(l_modes)
    hcx_all = np.empty((n, len(rcx)))
    gg_all = np.empty((n, len(rg)))
    for idx in range(n):
        l_mode = l_modes[idx]
        hcx = _radial_h_function(l_mode, rcx, vcx, vdcx, dr)
        hcx_all[idx] = hcx
        gg_all[idx] = _radial_g_function(hcx, l_mode, rg)
    return hcx_all, gg_all


if HAS_NUMBA:
    _radial_h_function = numba.jit(nopython=True)(_radial_h_function)
    _radial_g_function = numba.jit(nopython=True)(_radial_g_function)
    _boundary_coefficient = numba.jit(nopython=True)(_boundary_coefficient)
    _radial_h_g_batch = numba.jit(nopython=True)(_radial_h_g_batch)


class OutflowGrid(Grid):
    r"""
    Grid on which an outflow field solution is calculated.

    In addition to everything provided by `~sunkit_magex.pfss.grid.Grid`,
    this pre-computes the azimuthal and latitudinal eigenfunctions needed
    by the outflow eigenfunction expansion (see :func:`outflow`), and the
    ghost-extended radial coordinate ``rcx`` needed for the radial
    finite-difference recursion.
    """
    def __init__(self, ns, nphi, nr, rss):
        super().__init__(ns, nphi, nr, rss)
        self.ms, self.trigs = _azimuthal_eigenmodes(self.pc, self.dp)
        n_modes = len(self.ms)
        self.ls = np.zeros((n_modes, self.ns))
        self.legs = np.zeros((n_modes, self.ns, self.ns))
        for i in range(n_modes):
            self.ls[i], self.legs[i] = _latitudinal_eigenmodes(self.ms[i], self.sc, self.sg, self.ds, self.ns)

    @property
    def rcx(self):
        """
        Location of the centre of cells in log(r), plus ghost cells either
        side.
        """
        return np.linspace(-0.5 * self.dr, np.log(self.rss) + 0.5 * self.dr, self.nr + 2)


class OutflowInput(Input):
    r"""
    Input to outflow field modelling.

    Parameters
    ----------
    br : `sunpy.map.GenericMap`
        Boundary condition of radial magnetic field at the inner surface.
        Note that the data *must* have a cylindrical equal area projection.
    nr : `int`
        Number of cells in the radial direction on which to calculate the
        3D solution.
    rss : `float`
        Radius of the source surface, in units of solar radii.
    corona_temp : `astropy.units.Quantity`, optional
        Temperature of the corona for the implicit Parker wind solution, in
        units convertible to Kelvin (e.g. ``1.5e6 * u.K``).
    mf_constant : `float`, optional
        Magnetofrictional constant. ``mf_constant=0`` gives the zero-outflow
        (PFSS-equivalent) solution.
    sound_speed : `astropy.units.Quantity`, optional
        Sound speed of the Parker wind solution, in units convertible to
        m/s (e.g. ``150 * u.km / u.s``).

        Note this is a case where using a `~astropy.units.Quantity` removes
        an ambiguity present in the upstream ``outflowpy`` reference
        implementation: its docstring describes the equivalent (bare
        float) parameter as being in km/s, but its code does not perform
        any unit conversion internally, so the value it actually uses is
        in m/s.
    polynomial_coeffs : array-like, optional
        Coefficients of a polynomial (in :math:`r`) outflow speed profile.
    polynomial_type : `str`, optional
        One of ``'abs'``, ``'clip'``, ``'raw'``, ``'smooth'`` or
        ``'smooth_monotonic'`` (the default), specifying how the polynomial
        given by ``polynomial_coeffs`` is transformed to ensure it is
        non-negative.

    Notes
    -----
    If none of ``corona_temp``, ``mf_constant``, ``sound_speed`` or
    ``polynomial_coeffs`` are given, the default outflow speed profile
    (optimised to match eclipse observations, see Rice & Yeates 2021) is
    used. If ``mf_constant=0``, the zero-outflow (PFSS-equivalent) solution
    is used. If ``mf_constant`` and either of ``corona_temp``/
    ``sound_speed`` are given, a Parker wind profile is used. If
    ``polynomial_coeffs`` is given, a polynomial profile is used.
    """
    @u.quantity_input
    def __init__(self, br, nr, rss, corona_temp: u.K = None, mf_constant=None,
                 sound_speed: u.m / u.s = None, polynomial_coeffs=None,
                 polynomial_type='smooth_monotonic'):
        super().__init__(br, nr, rss, br_outer='radial')

        base_grid = self._grid
        self._grid = OutflowGrid(base_grid.ns, base_grid.nphi, nr, rss)

        if polynomial_coeffs is not None:
            self._set_polynomial_profile(polynomial_coeffs, polynomial_type)
        elif mf_constant is not None and sound_speed is not None:
            self._set_parker_profile(mf_constant, sound_speed=sound_speed.to_value(u.m / u.s))
        elif sound_speed is not None and corona_temp is not None:
            raise ValueError(
                'Please specify either a sound speed or a corona '
                'temperature, but not both.')
        elif mf_constant is not None and corona_temp is not None:
            self._set_parker_profile(mf_constant, corona_temp=corona_temp.to_value(u.K))
        elif mf_constant == 0.0:
            self.vg = np.zeros(len(self._grid.rg))
            self.vcx = np.zeros(len(self._grid.rcx))
            self.vdcx = np.zeros(len(self._grid.rcx))
        else:
            self._set_default_profile()

    def _extended_rg(self):
        rgx = np.zeros(len(self._grid.rg) + 2)
        rgx[1:-1] = self._grid.rg
        rgx[0] = 2 * rgx[1] - rgx[2]
        rgx[-1] = 2 * rgx[-2] - rgx[-3]
        return rgx

    def _set_polynomial_profile(self, polynomial_coeffs, polynomial_type):
        def poly_at_pt(r):
            res = 0
            for i, c in enumerate(polynomial_coeffs):
                res = res + c * np.exp(r)**i
            return res

        rgx = self._extended_rg()
        vgx = poly_at_pt(rgx)
        vcx = poly_at_pt(self._grid.rcx)

        if polynomial_type == 'abs':
            vgx = np.abs(vgx)
            vcx = np.abs(vcx)
        elif polynomial_type == 'clip':
            vgx = np.clip(vgx, a_min=0.0, a_max=None)
            vcx = np.clip(vcx, a_min=0.0, a_max=None)
        elif polynomial_type == 'raw':
            pass
        elif polynomial_type in ('smooth', 'smooth_monotonic'):
            vgx = vgx + 1e-6
            vcx = vcx + 1e-6
            vgx = (vgx * np.exp(vgx)) / (np.exp(vgx) - 1)
            vcx = (vcx * np.exp(vcx)) / (np.exp(vcx) - 1)
            if polynomial_type == 'smooth_monotonic':
                vgmax_ind = np.argmax(vgx)
                if vgmax_ind != len(vgx) - 1:
                    vgx[vgmax_ind:] = vgx[vgmax_ind]
                vcmax_ind = np.argmax(vcx)
                if vcmax_ind != len(vcx) - 1:
                    vcx[vcmax_ind:] = vcx[vcmax_ind]
        else:
            raise ValueError(
                "polynomial_type not recognised. Currently allowed types "
                "are 'abs', 'clip', 'raw', 'smooth' and 'smooth_monotonic'.")

        self.vg = vgx[1:-1]
        self.vcx = vcx
        self.vdcx = (vgx[1:] - vgx[:-1]) / (rgx[1:] - rgx[:-1])

    def _parker_implicit_fn(self, r, v, r_c):
        """
        Implicit Parker solar wind function; its root in ``v`` gives the
        (dimensionless) wind speed at radius ``r``.
        """
        if np.abs(v) < 1e-12:
            return 1e12
        res = v**2
        res -= 2 * np.log(np.abs(v))
        res -= 4 * (np.log(np.abs(r / r_c)) + r_c / r)
        res += 3
        return res

    def _get_parker_wind_speed(self, r_c):
        min_r = -1.0
        max_r = self._grid.rg[-1] * 2.0
        vtest_min = 1e-6
        dr = (max_r - min_r) / 2000

        vfinals = []
        r0s = []
        r0 = min_r

        def implicit_fn(r, v):
            return self._parker_implicit_fn(r, v, r_c)

        while r0 <= max_r:
            minimum = minimize_scalar(lambda v: implicit_fn(np.exp(r0), v))
            p0, p1, p2 = vtest_min, minimum.x, 10.0 * minimum.x
            if (implicit_fn(np.exp(r0), p0) * implicit_fn(np.exp(r0), p1) < 0.0
                    and implicit_fn(np.exp(r0), p1) * implicit_fn(np.exp(r0), p2) < 0.0):
                vslow = root_scalar(lambda v: implicit_fn(np.exp(r0), v), bracket=[p0, p1]).root
                vfast = root_scalar(lambda v: implicit_fn(np.exp(r0), v), bracket=[p1, p2]).root
                if len(vfinals) < 2:
                    vfinals.append(vslow)
                    r0s.append(r0)
                else:
                    prediction = 2 * vfinals[-1] - vfinals[-2]
                    if np.abs(vslow - prediction) < np.abs(vfast - prediction):
                        vfinals.append(vslow)
                        r0s.append(r0)
                    else:
                        vfinals.append(vfast)
                        r0s.append(r0)
            else:
                if r0 < np.log(2.5):
                    vfinals.append(0.0)
                    r0s.append(r0)
                else:
                    raise RuntimeError(
                        'A sensible solution to the implicit Parker wind '
                        'speed equation could not be found.')
            r0 = r0 + dr

        vfinals = np.array(vfinals)
        r0s = np.array(r0s)
        vf = interpolate.interp1d(r0s, vfinals, bounds_error=False, fill_value='extrapolate')

        rgx = self._extended_rg()
        vgx = vf(rgx)
        vcx = vf(self._grid.rcx)
        vdcx = (vgx[1:] - vgx[:-1]) / (rgx[1:] - rgx[:-1])
        return vgx[1:-1], vcx, vdcx

    def _set_parker_profile(self, mf_constant, sound_speed=None, corona_temp=None):
        """
        ``sound_speed``/``corona_temp`` are plain floats in SI units here
        (m/s, K) -- unit handling and validation happens in ``__init__``,
        where the public API accepts `~astropy.units.Quantity`.
        """
        # Physical constants in SI units (see Rice & Yeates 2021 for the
        # isothermal Parker wind solution used here).
        gravitational_constant = 6.67408e-11
        solar_mass = 1.98847542e30
        boltzmann_constant = 1.38064852e-23
        proton_mass = 1.67262192e-27
        r_sun_m = 6.957e8
        r_sun_cm = 6.957e10

        if sound_speed is not None:
            sound_speed_m_s = sound_speed
        else:
            sound_speed_m_s = np.sqrt(boltzmann_constant * corona_temp / proton_mass)

        mf_in_sensible_units = mf_constant * r_sun_cm**2
        r_c = (gravitational_constant * solar_mass / (2 * sound_speed_m_s**2)) / r_sun_m
        c_s = mf_in_sensible_units * sound_speed_m_s / r_sun_m

        vg, vcx, vdcx = self._get_parker_wind_speed(r_c)
        self.vg = vg * c_s
        self.vcx = vcx * c_s
        self.vdcx = vdcx * c_s

    def _set_default_profile(self):
        flow_data = np.loadtxt(_DATA_DIR / 'opt_flow.txt', delimiter=',')
        rs_interp = flow_data[:, 0]
        vs_interp = flow_data[:, 1]

        rgx = self._extended_rg()
        f = interpolate.interp1d(rs_interp, vs_interp, fill_value='extrapolate')
        vgx = f(np.exp(rgx))
        vcx = f(np.exp(self._grid.rcx))

        vgx = np.clip(vgx, a_min=0.0, a_max=np.max(vgx))
        vcx = np.clip(vcx, a_min=0.0, a_max=np.max(vcx))

        self.vg = vgx[1:-1]
        self.vcx = vcx
        self.vdcx = (vgx[1:] - vgx[:-1]) / (rgx[1:] - rgx[:-1])


def outflow(input, mode_tol=1e-10):
    r"""
    Compute an outflow field extrapolation.

    Extrapolates a 3D magnetic field including a radial solar wind outflow
    speed profile, using an eigenfunction expansion in :math:`r, s, p`
    coordinates on the same dumfric grid used by `sunkit_magex.pfss.pfss`
    (equally spaced in :math:`\rho = \ln(r/r_{sun})`,
    :math:`s = \cos(\theta)`, and :math:`p = \phi`). Setting the outflow
    speed to zero everywhere (``OutflowInput(..., mf_constant=0)``)
    recovers the PFSS solution.

    Parameters
    ----------
    input : `OutflowInput`
        Input parameters, including the outflow speed profile.
    mode_tol : float, optional
        Modes whose contribution to the lower boundary condition is
        smaller than ``mode_tol`` times the peak ``|Br|`` of the input map
        are treated as negligible and skipped. This is a *relative*
        tolerance (unlike outflowpy's reference implementation, which
        uses a hardcoded absolute threshold of ``1e-10`` which never
        really skips anything so all modes use expensive radial recursion).

    Returns
    -------
    ~sunkit_magex.pfss.OutflowOutput

    See Also
    --------
    sunkit_magex.pfss.pfss : Compute a PFSS model.

    Notes
    -----
    See Rice & Yeates, 2021, *Global Coronal Equilibria with Solar Wind
    Outflow*, ApJ 923, 57, https://doi.org/10.3847/1538-4357/ac2c71, and
    https://github.com/oekrice/outflowpy for the original implementation
    this is adapted from.

    **Performance note.** Direct implementation of this method is
    ``O((ns*nphi)^2)`` for the boundary-fitting sweep and
    ``O(n_active * nr * ns * nphi)`` for the field reconstruction, in the
    total number of modes. For a normalised analytic test field (e.g. a
    dipole), only a handful of modes are ever non-negligible, so this is
    ok a real magnetogram data, essentially every mode is non-negligible,
    and the cost becomes prohibitive. Three changes address this:

    1. **Boundary-fitting coefficients (new)** ``C_{l,m}`` (`_boundary_coefficient`
       is the reference per-mode formula) are computed for every ``(i, j)``
       pair at once: the azimuthal projection ``T = Br @ Phi`` is one
       matmul, and the remaining latitudinal sum is a batched
       matrix-vector product (one small matmul per ``i``, executed as a
       single batched ``@`` call), instead of ``ns*nphi`` separate calls
       each redoing an ``O(ns*nphi)`` reduction from scratch.
    2. **Field reconstruction (in orig Fortran but not python)** :
       the radial (H, G) and latitudinal (Q) profiles of a mode depend
       on *both* azimuthal index ``i`` and latitudinal index ``j``, but
       the azimuthal profile (P) depends only on ``i``. That asymmetry
       is exploited to reconstruct br/bs/bp with two batched matrix
       multiplications instead of one broadcast-and-add per mode:
       for each ``i``, the contributions of its non-negligible ``j`` modes
       are combined with a single small matmul (at most ``ns`` terms), and
       the per-``i`` results are then combined with one final matmul that
       sums over ``i`` (at most``nphi`` terms).
    3. **Radial recursion (new numba optim)** (`_radial_h_g_batch`):
       all of a given ``i``'s non-negligible ``j`` modes are processed
       in a single jitted call rather than one Python-level call per mode,
       avoiding Python<->numba dispatch overhead when there are many
       non-negligible modes.

    Together these keep the solve cost close to
    ``O(n_active * nr * ns)`` even when nearly every mode is
    non-negligible, rather than scaling with the square of the total
    mode count. Measured on this machine:

    - Full-resolution HMI synoptic map (360x180 pixels, nr=35, ~65,000
      total modes, effectively all non-negligible): several minutes
      (naive per-mode python loop) down to ~0.3s.
    """
    grid = input.grid
    br0 = input.br

    nr = grid.nr
    ns = grid.ns
    nphi = grid.nphi

    sigs = np.sqrt(1 - grid.sg**2)
    sigc = np.sqrt(1 - grid.sc**2)

    br_scale = np.max(np.abs(br0))
    threshold = mode_tol * br_scale if br_scale > 0 else mode_tol

    # Vectorised boundary-fitting coefficients C_{l,m} for every (m, l)
    # pair at once (see `_boundary_coefficient` for the per-mode formula
    # this is derived from). T_{s,m} = sum_p Br(s,p) Phi_m(p) is one
    # matmul, and the l-sum is then a batched matrix-vector product (one
    # small matmul per m, executed as a single batched ``@`` call), instead
    # of ns*nphi separate calls each redoing an O(ns*nphi) reduction from
    # scratch (which made the boundary-fitting sweep O((ns*nphi)^2)
    # overall).
    T = br0 @ grid.trigs  # (ns, n_m)
    numerator = (T.T[:, np.newaxis, :] @ grid.legs).squeeze(1)  # (n_m, ns)
    q_sq_sum = np.sum(grid.legs**2, axis=1)  # (n_m, ns)
    p_sq_sum = np.sum(grid.trigs**2, axis=0)  # (n_m,)
    all_cmls = numerator / (q_sq_sum * p_sq_sum[:, np.newaxis])  # (nphi, ns)

    # Azimuthal (P) profiles only depend on i, so build them once for
    # every mode up front.
    p_pad = np.zeros((nphi, nphi + 2))
    p_pad[:, 1:-1] = grid.trigs.T
    p_pad[:, 0] = p_pad[:, -2]
    p_pad[:, -1] = p_pad[:, 1]
    p_central = p_pad[:, 1:-1]
    dp_ = (p_pad[:, 1:] - p_pad[:, :-1]) / grid.dp

    # Per-i partial reconstructions (summed over each i's active j modes).
    partial_br = np.zeros((nphi, nr + 1, ns))
    partial_bs = np.zeros((nphi, nr, ns + 1))
    partial_bp = np.zeros((nphi, nr, ns))

    for i in range(nphi):
        cmls_all_j = all_cmls[i]
        active_j = np.nonzero(np.abs(cmls_all_j) >= threshold)[0]
        if len(active_j) == 0:
            continue

        n_i = len(active_j)
        cmls_i = cmls_all_j[active_j]

        # One jitted call for all of this i's active modes, instead of
        # one Python-level call per mode (see `_radial_h_g_batch`).
        l_modes_i = grid.ls[i, active_j]
        hcx_all_i, gg_i = _radial_h_g_batch(l_modes_i, grid.rcx, input.vcx, input.vdcx, grid.dr, grid.rg)
        hcx_central_i = hcx_all_i[:, 1:-1]

        q_central_i = np.empty((n_i, ns))
        b_bs_i = np.empty((n_i, ns + 1))
        b_bp_i = np.empty((n_i, ns))

        for idx, j in enumerate(active_j):
            q_pad = np.zeros(ns + 2)
            q_pad[1:-1] = grid.legs[i, :, j]
            q_pad[0] = q_pad[1]
            q_pad[-1] = q_pad[-2]

            q_central_i[idx] = q_pad[1:-1]
            b_bs_i[idx] = sigs * (q_pad[1:] - q_pad[:-1]) / grid.ds
            b_bp_i[idx] = q_pad[1:-1] / sigc

        a_br_i = cmls_i[:, np.newaxis] * gg_i
        a_bsbp_i = cmls_i[:, np.newaxis] * hcx_central_i

        partial_br[i] = a_br_i.T @ q_central_i
        partial_bs[i] = a_bsbp_i.T @ b_bs_i
        partial_bp[i] = a_bsbp_i.T @ b_bp_i

    # Final batched contraction over i (the azimuthal mode index).
    br = (partial_br.reshape(nphi, -1).T @ p_central).reshape(nr + 1, ns, nphi)
    bs = (partial_bs.reshape(nphi, -1).T @ p_central).reshape(nr, ns + 1, nphi)
    bp = (partial_bp.reshape(nphi, -1).T @ dp_).reshape(nr, ns, nphi + 1)

    br = np.swapaxes(br, 0, 2)
    bs = np.swapaxes(bs, 0, 2)
    bp = np.swapaxes(bp, 0, 2)

    return OutflowOutput(br, bs, bp, grid, input.map)
