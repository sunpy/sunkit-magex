r"""
Analytic inputs and solutions to the PFSS equations.

All angular quantities must be passed as astropy quantities. All radial
quantities are passed normalised to the source surface radius, and therefore
can be passed as normal scalar values.

Angular definitions
-------------------
- ``theta`` is the polar angle, in the range :math:`0, \pi`
  (ie. the co-latitude).
- ``phi`` is the azimuthal angle, in the range :math:`0, 2\pi`.

Using this module requires `sympy` to be installed.
"""
from functools import partial

import astropy.units as u
import numpy as np
import scipy.special
import sympy


@u.quantity_input
def _normalise_angles(theta: u.deg, phi: u.deg):
    """
    - Strip units and return radians
    """
    theta = theta.to_value(u.rad)
    phi = phi.to_value(u.rad)
    return theta, phi


def _Ynm(l, m, theta, phi):
    """
    Return values of spherical harmonic with numbers l, m at coordiantes
    theta, phi.
    """
    # Note swapped arguments phi, theta, as scipy has a different
    # definition of these
    return scipy.special.sph_harm(m, l, phi, theta)


def _cot(theta):
    return 1 / np.tan(theta)


_extras = {'Ynm': _Ynm, 'cot': _cot}


def _spherical_harmonic_sympy(l, m):
    """
    Return a complex spherical harmonic with numbers ``l, m``.
    """
    L, M = sympy.symbols('l, m')
    theta, phi = sympy.symbols('theta, phi')
    harm = sympy.Ynm(L, M, theta, phi)
    harm = harm.subs([(L, l), (M, m)])
    return harm, theta, phi


def _real_spherical_harmonic_sympy(l, m):
    """
    Return a real spherical harmonic.
    """
    sph, theta, phi = _spherical_harmonic_sympy(l, m)
    if m == 0:
        return sph, theta, phi
    elif m < 0:
        # Multiply by i to get imaginary part later
        return sympy.sqrt(2) * (-1)**m * 1j * sph, theta, phi
    elif m > 0:
        return sympy.sqrt(2) * (-1)**m * sph, theta, phi


def _c(l, zss):
    """
    """
    def cl(z):
        return (z**(-l - 2) *
                (l + 1 + l * (z / zss)**(2 * l + 1)) /
                (l + 1 + l * zss**(-2 * l - 1)))

    return cl


def _d(l, zss):
    """
    """
    def dl(z):
        return (z**(-l - 2) *
                (1 - (z / zss)**(2 * l + 1)) /
                (l + 1 + l * zss**(-2 * l - 1)))

    return dl


def Br(l, m, zss):
    """
    Returns
    -------
    function :
        ``Br(r, theta, phi)``, which takes coordiantes and returns the radial
        magnetic field component.
    """
    sph, t, p = _real_spherical_harmonic_sympy(l, m)
    sph = sympy.lambdify((t, p), sph, _extras)

    def f(r, theta, phi):
        theta, phi = _normalise_angles(theta, phi)
        return _c(l, zss)(r) * np.real(sph(theta, phi))

    return f


def Btheta(l, m, zss):
    """
    Returns
    -------
    function :
        ``Btheta(z, theta, phi)``, which takes coordiantes and returns the
        radial magnetic field component.
    """
    sph, t, p = _real_spherical_harmonic_sympy(l, m)
    sph = sympy.diff(sph, t)
    sph = sympy.lambdify((t, p), sph, [_extras, 'numpy'])

    @u.quantity_input
    def f(z, theta: u.deg, phi: u.deg):
        theta, phi = _normalise_angles(theta, phi)
        return _d(l, zss)(z) * np.real(sph(theta, phi))

    return f


def Bphi(l, m, zss):
    """
    Returns
    -------
    function :
        ``Btheta(z, theta, phi)``, which takes coordiantes and returns the
        radial magnetic field component.
    """
    sph, t, p = _real_spherical_harmonic_sympy(l, m)
    sphi = sympy.diff(sph, p)
    sph = sympy.lambdify((t, p), sph, [_extras, 'numpy'])

    @u.quantity_input
    def f(z, theta: u.deg, phi: u.deg):
        theta, phi = _normalise_angles(theta, phi)
        return _d(l, zss)(z) * np.real(sph(theta, phi)) / np.sin(theta)

    return f
