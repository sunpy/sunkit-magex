"""
Code for calculating a PFSS extrapolation.
"""
import numpy as np

from sunkit_magex.pfss.input import Input
from sunkit_magex.pfss.output import Output
from ._pfss_internal import _A_diag, _als_alp, _compute_r_term


def pfss(input: Input):
    r"""
    Compute PFSS model.

    Extrapolates a 3D PFSS using an eigenfunction method in :math:`r,s,p`
    coordinates, on the dumfric grid
    (equally spaced in :math:`\rho = \ln(r/r_{sun})`,
    :math:`s= \cos(\theta)`, and :math:`p=\phi`).

    Parameters
    ----------
    input : ~sunkit_magex.pfss.Input
        Input parameters.

    Returns
    -------
    ~sunkit_magex.pfss.Output

    Notes
    -----
    In order to avoid numerical issues, the monopole term (which should be zero
    for a physical magnetic field anyway) is explicitly excluded from the
    solution.

    The output should have zero current to machine precision,
    when computed with the DuMFriC staggered discretization.
    """
    br0 = input.br
    nr = input.grid.nr
    ns = input.grid.ns
    nphi = input.grid.nphi
    rss = input.grid.rss

    # Coordinates:
    ds = input.grid.ds
    dp = input.grid.dp
    dr = input.grid.dr

    input.grid.rg
    input.grid.rc

    sg = input.grid.sg
    sc = input.grid.sc

    br_outer = input.br_outer

    k = np.linspace(0, nr, nr + 1)

    Fp = sg * 0  # Lp/Ls on p-ribs
    Fp[1:-1] = np.sqrt(1 - sg[1:-1]**2) / (np.arcsin(sc[1:]) - np.arcsin(sc[:-1])) * dp
    Vg = Fp / ds / dp
    Fs = (np.arcsin(sg[1:]) - np.arcsin(sg[:-1])) / np.sqrt(1 - sc**2) / dp  # Ls/Lp on s-ribs
    Uc = Fs / ds / dp

    # FFT in phi of photospheric distribution at each latitude:
    brt = np.fft.rfft(br0, axis=1)
    brt = brt.astype(np.complex128)

    if isinstance(br_outer, np.ndarray):
        brt_outer = np.fft.rfft(br_outer, axis=1)
        brt_outer = brt_outer.astype(np.complex128)
    else:
        brt_outer = None

    # Prepare tridiagonal matrix:
    # - create off-diagonal part of the matrix:
    A = np.zeros((ns, ns))
    for j in range(ns - 1):
        A[j, j + 1] = -Vg[j + 1]
        A[j + 1, j] = A[j, j + 1]
    # - term required for m-dependent part of matrix:
    mu = np.fft.fftfreq(nphi)
    mu = 4 * np.sin(np.pi * mu)**2
    # - initialise:
    psir = np.zeros((nr + 1, ns), dtype='complex')
    psi = np.zeros((nr + 1, ns, nphi), dtype='complex')
    e1 = np.exp(dr)
    fact = np.sinh(dr) * (e1 - 1)

    # Loop over azimuthal modes (positive m):
    for m in range(nphi // 2 + 1):
        # - set diagonal terms of matrix:
        A = _A_diag(A, ns, Vg, Uc, mu, m)

        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that A is symmetric so use special solver)
        lam, Q = np.linalg.eigh(A)
        Q = Q.astype(np.complex128)
        # - solve quadratic:
        Flm = 0.5 * (1 + e1 + lam * fact)
        ffp = Flm + np.sqrt(Flm**2 - e1)
        ffm = e1 / ffp

        # - compute radial term for each l (for this m):
        psi, psir = _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir, rss, brt_outer)

        if (m > 0):
            psi[:, :, nphi - m] = np.conj(psi[:, :, m])

    # Past this point only psi, Fs, Fp are needed
    # Compute psi by inverse fft:
    psi = np.real(np.fft.ifft(psi, axis=2))

    # Hence compute vector potential [note index order, for netcdf]:
    # (note that alr is zero by definition)
    alr = np.zeros((nphi + 1, ns + 1, nr))
    als = np.zeros((nphi + 1, ns, nr + 1))
    alp = np.zeros((nphi, ns + 1, nr + 1))

    als, alp = _als_alp(nr, nphi, Fs, psi, Fp, als, alp)

    return Output(alr, als, alp, input.grid, input.map)
