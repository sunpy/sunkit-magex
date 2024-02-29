==========================
Numerical methods for PFSS
==========================

Basic Equations
===============

The ``sunkit_magex.pfss`` code solves the for a magnetic field satisfying

.. math::

    \nabla\times\boldsymbol{B}=0,\qquad
    \nabla\cdot\boldsymbol{B}= 0

in a spherical shell :math:`1 \leq r \leq R_{ss}`, given boundary conditions

.. math::

    &B_r(\theta,\phi) = g(\theta,\phi) \quad \textrm{on} \quad r=1

    &B_\theta=B_\phi=0 \quad \textrm{on} \quad r=R_{ss}.

The function :math:`g(\theta,\phi)` is specified by the user.

Numerical Grid
==============

.. figure:: grid.jpg
    :name: fig:grid

    The numerical grid used in the solver.

The solver uses a rectilinear grid that is equally spaced in :math:`\rho`, :math:`s`, :math:`\phi`, where

.. math::

    \rho = \ln(r), \quad s=\cos\theta

in terms of spherical coordinates :math:`(r,\theta\,\phi)`.
The coordinate scale factors :math:`|\mathrm{d}\boldsymbol r/\mathrm{d}\rho|`, :math:`|\mathrm{d}\boldsymbol r/\mathrm{d}s|`, :math:`|\mathrm{d}\boldsymbol r/\mathrm{d}\phi|` are

.. math::

    h_\rho = r = \mathrm{e}^\rho,\quad h_s = \frac{r}{\sin\theta} = \frac{\mathrm{e}^\rho}{\sqrt{1-s^2}}, \quad h_\phi = r\sin\theta = \mathrm{e}^\rho\sqrt{1-s^2}

The grid is illustrated in Figure `1 <#fig:grid>`__.
Note that the longitudinal cell size goes to zero at the poles; these points are treated specially in the calculation of :math:`\boldsymbol{B}`.
Note also that, since :math:`s` is a decreasing function of :math:`\theta`, the components of a vector :math:`\boldsymbol{v}` in :math:`(\rho,s,\phi)` are :math:`v_\rho = v_r` but :math:`v_s = -v_\theta`.

We define the number of grid cells :math:`n_\rho`, :math:`n_s`, :math:`n_\phi`, with corresponding uniform spacings

.. math::

    \Delta_\rho= \frac{\ln(R_{ss})}{n_\rho}, \quad \Delta_s= \frac{2}{n_s}, \quad \Delta_\phi= \frac{2\pi}{n_\phi}

Note that the boundaries in :math:`s` are at the poles :math:`s=\pm1`, at which points :math:`h_\phi` is not defined. The solution is periodic in the longitudinal (:math:`\phi`) direction.

Numerical method
================

Overall strategy
----------------

Rather than writing :math:`\boldsymbol{B}= \nabla\chi` and solving :math:`\nabla^2\chi=0`, we write instead :math:`\boldsymbol{A}= \nabla\times\big(\psi \,\mathrm{e}_\rho\big)`.
Then accounting for the unusual coordinates we get

.. math::

    \boldsymbol{B}&=  \frac{1}{h_\rho h_sh_\phi}

    h_\rho\,\mathrm{e}_\rho& h_\phi\,\mathrm{e}_\phi& h_s\,\mathrm{e}_s
    \partial_\rho & \partial_\phi & \partial_s
    0 & \frac{h_\phi}{h_s}\partial_s\psi & -\frac{h_s}{h_\phi}\partial_\phi\psi

    &= -\frac{1}{h_sh_\phi}\left[\partial_s\left(\frac{h_\phi}{h_s}\partial_s\psi\right) + \partial_\phi\left(\frac{h_s}{h_\phi}\partial_\phi\psi\right) \right]\,\mathrm{e}_\rho+ \frac{1}{h_\rho h_\phi}\partial_\phi\partial_\rho\psi\,\mathrm{e}_\phi+ \frac{1}{h_\rho h_s}\partial_s\partial_\rho\psi\,\mathrm{e}_s

    &= -\Delta_\perp\psi\,\mathrm{e}_\rho+ \frac{1}{h_\phi}\partial_\phi\left(\frac{1}{h_\rho}\partial_\rho\psi\right)\,\mathrm{e}_\phi+ \frac{1}{h_s}\partial_s\left(\frac{1}{h_\rho}\partial_\rho\psi\right)\,\mathrm{e}_s

This will take the curl-free form :math:`\boldsymbol{B}= \nabla\big(\tfrac1{h_\rho}\partial_\rho\psi\big)`  provided that

.. math::

    \nabla^2_\perp\psi = -\frac{1}{h_\rho}\partial_\rho\left(\frac{1}{h_\rho}\partial_\rho\psi\right)

so our strategy is to solve Equation `[eqn:psi] <#eqn:psi>`__ for :math:`\psi`, then reconstruct :math:`\boldsymbol{A}` and :math:`\boldsymbol{B}`.
The reason for doing it this way is that it allows us to compute :math:`\boldsymbol{A}` as well as :math:`\boldsymbol{B}` (again, for legacy reasons).

Numerical solution
------------------

We follow the method described in [vanballe2000]_, except that we modify the finite-difference discretization to suit our particular coordinates.

The discretization is chosen so that we will have :math:`\nabla\times\boldsymbol{B}=0` to machine precision on a staggered grid, when the curl is taken using central differences.
This property of essentially zero current density is required when using the PFSS solution to, e.g., initialize non-potential simulations.
It would not typically be achieved by interpolating a spherical harmonic solution onto the numerical grid. However, we will see that the discrete solution effectively computes discrete approximations of the spherical harmonics, tailored to our particular difference formula.

In the following subsections, we describe the numerical solution in more detail.

Variables
~~~~~~~~~

Let the coordinate grid points be defined by

.. math::

    &\rho^k = k\Delta_\rho, \qquad k=0,\ldots, n_\rho

    &s^j = j\Delta_s- 1, \qquad j=0,\ldots, n_s

    &\phi^i = i\Delta_\phi, \qquad i=0,\ldots, n_\phi

In the code the first two arrays are called ``rg`` and ``sg`` (that for ``pg`` is not required). There are also arrays ``rc``, ``sc`` and ``pc`` corresponding to the cell centres, i.e. :math:`\rho^{k+1/2}`, :math:`s^{j+1/2}` and :math:`\phi^{i+1/2}`.

To deal with the curvilinear coordinates, we define the edge lengths

.. math::

    &L_\rho^{k+\nicefrac{1}{2},j,i} = \int_{\rho^k}^{\rho^{k+1}} h_\rho\,\mathrm{d}\rho = \,\mathrm{e}^{\rho^{k+1}} - \,\mathrm{e}^{\rho^k}

    &L_s^{k,j+\nicefrac{1}{2},i} = \int_{s^j}^{s^{j+1}} h_s\,\mathrm{d}s = \,\mathrm{e}^{\rho^k}\big(\arcsin(s^{j+1}) - \arcsin(s^j)\big)

    &L_\phi^{k,j,i+\nicefrac{1}{2}} = \int_{\phi^i}^{\phi^{i+1}} h_\phi\,\mathrm{d}\phi = \,\mathrm{e}^{\rho^k}\sigma^j\Delta_\phi.

Here we used the fact that :math:`\Delta_\rho`, :math:`\Delta_s` and :math:`\Delta_\phi` are constant, and used the shorthand

.. math::

    \sigma^j := \sqrt{1 - (s^j)^2}.

Similarly we define the areas of the cell faces

.. math::

    &S_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} =  \int_{\phi^i}^{\phi^{i+1}}\int_{s^j}^{s^{j+1}} h_s h_\phi\,\mathrm{d}s\mathrm{d}\phi = \,\mathrm{e}^{2\rho^k}\Delta_s\Delta_\phi

    &S_s^{k+\nicefrac{1}{2},j,i+\nicefrac{1}{2}} = \int_{\rho^k}^{\rho^{k+1}}\int_{\phi^i}^{\phi^{i+1}} h_\rho h_\phi\,\mathrm{d}\phi\mathrm{d}\rho = \tfrac12\big(\,\mathrm{e}^{2\rho^{k+1}} - \,\mathrm{e}^{2\rho^{k}}\big)\sigma^j\,\Delta_\phi

    &S_\phi^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i} = \int_{\rho^k}^{\rho^{k+1}}\int_{s^j}^{s^{j+1}}h_\rho h_s\,\mathrm{d}s\mathrm{d}\rho = \tfrac12\big(\,\mathrm{e}^{2\rho^{k+1}}- \,\mathrm{e}^{2\rho^k}\big)\big(\arcsin(s^{j+1}) - \arcsin(s^j)\big)

In the code the face areas are stored in arrays ``Sbr``, ``Sbs`` and ``Sbp`` (with only the dimensions required).

In the code the magnetic field :math:`\boldsymbol{B}` is defined staggered on the face centres, so :math:`B_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}`, :math:`B_s^{k+\nicefrac{1}{2},j,i+\nicefrac{1}{2}}`, :math:`B_\phi^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i}`.
These variables are called ``br``, ``bs`` and ``bp``.

The vector potential is located on the corresponding cell edges, so :math:`A_\rho^{k+\nicefrac{1}{2},j,i}`,
:math:`A_s^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}`, :math:`A_\phi^{k,j,i+\nicefrac{1}{2}}`.
In fact, these values are never required on their own, only multiplied by the corresponding edge lengths. So the variables ``alr``, ``als`` and ``alp`` correspond to the products :math:`L_\rho A_\rho`, :math:`L_sA_s` and :math:`L_\phi A_\phi`, respectively.

Finally, the potential :math:`\psi` is located on the :math:`\rho`-faces (like :math:`B_\rho`), so :math:`\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}`. It is stored in the variable ``psi``.

Derivatives
~~~~~~~~~~~

Firstly, we have :math:`\boldsymbol{A}= \nabla\times\big(\psi\,\mathrm{e}_\rho\big)`.
Numerically, this is approximated by

.. math::

    A_s^{k,j+\nicefrac{1}{2},i} = -\frac{\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k,j+\nicefrac{1}{2},i-\nicefrac{1}{2}}}{L_\phi^{k,j+\nicefrac{1}{2},i}}, \qquad A_\phi^{k,j,i+\nicefrac{1}{2}} = \frac{\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k,j-\nicefrac{1}{2},i+\nicefrac{1}{2}}}{L_s^{k,j,i+\nicefrac{1}{2}}}

The magnetic field :math:`\boldsymbol{B}= \nabla\times\boldsymbol{A}` is
then approximated by

.. math::

    &(S_\rho B_\rho)^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = (L_s A_s)^{k,j+\nicefrac{1}{2},i+1} - (L_s A_s)^{k,j+\nicefrac{1}{2},i} - (L_\phi A_\phi)^{k,j+1,i+\nicefrac{1}{2}} + (L_\phi A_\phi)^{k,j,i+\nicefrac{1}{2}},

    &(S_s B_s)^{k+\nicefrac{1}{2},j,i+\nicefrac{1}{2}} = (L_\phi A_\phi)^{k+1,j,i+\nicefrac{1}{2}} - (L_\phi A_\phi)^{k,j,i+\nicefrac{1}{2}},

    &(S_\phi B_\phi)^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i} =  - (L_s A_s)^{k+1,j+\nicefrac{1}{2},i} + (L_s A_s)^{k,j+\nicefrac{1}{2},i}.

These formulae correspond to Stokes' Theorem applied to the cell face.
The condition :math:`\nabla\times\boldsymbol{B}=0` may be expressed similarly as

.. math::

    &0 = (L_s B_s)^{k+\nicefrac{1}{2},j,i-\nicefrac{1}{2}} - (L_s B_s)^{k+\nicefrac{1}{2},j,i+\nicefrac{1}{2}} - (L_\phi B_\phi)^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i} + (L_\phi B_\phi)^{k+\nicefrac{1}{2},j-\nicefrac{1}{2},i}
    &0 = (L_\phi B_\phi)^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i} - (L_\phi B_\phi)^{k-\nicefrac{1}{2},j+\nicefrac{1}{2},i} - (L_\rho B_\rho)^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} + (L_\rho B_\rho)^{k,j+\nicefrac{1}{2},i-\nicefrac{1}{2}}
    &0 = (L_\rho B_\rho)^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - (L_\rho B_\rho)^{k,j-\nicefrac{1}{2},i+\nicefrac{1}{2}} - (L_s B_s)^{k+\nicefrac{1}{2},j,i+\nicefrac{1}{2}} + (L_s B_s)^{k-\nicefrac{1}{2},j,i+\nicefrac{1}{2}}

Note that the factors :math:`L_\rho`, :math:`L_s`, :math:`L_\phi` here are defined normal to the cell faces, not on the edges.
But they have the same formulae.

In fact, condition `[eqn:j1] <#eqn:j1>`__ is automatically satisfied.
This may be shown using equations `[eqn:as] <#eqn:as>`__ to `[eqn:bp] <#eqn:bp>`__, together with our formulae for :math:`L_s`, :math:`L_\phi`, :math:`S_s` and :math:`S_\phi`.

Below, we will discretize `[eqn:psi] <#eqn:psi>`__ in such a way that conditions `[eqn:j2] <#eqn:j2>`__ and `[eqn:j3] <#eqn:j3>`__ are also satisfied exactly (up to rounding error).

Boundary conditions for :math:`\boldsymbol{B}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boundary conditions are needed when ``br``, ``bs``, ``bp`` are averaged to the grid points for output.
We use a layer of "ghost cells", whose values are set by the following boundary conditions:

#. In :math:`\phi`, ``br`` and ``bs`` are simply periodic.

#. At the outer boundary :math:`\rho=\log(R_{ss})`, ghost values of ``bs`` and ``bp`` are set assuming constant gradient in :math:`\rho`.

#. At the inner boundary, :math:`\rho=0`, ghost values of ``bs`` and ``bp`` are set using equations `[eqn:j2] <#eqn:j2>`__ and `[eqn:j3] <#eqn:j3>`__ (effectively assuming zero horizontal current density).

#. At the polar boundaries, the ghost value of ``br`` is set to the polemost interior value from the opposite side of the grid.
   Similarly, ``bp`` is set to minus the polemost interior value from the opposite side of the grid. The values of ``bs`` at the poles are not meaningful as the cell faces have zero area.
   However, they are set to the average of the two neighboring interior values at that longitude (with the opposite one being reversed in sign).

Some of these conditions are chosen for compatibility with other codes, and are not necessarily the most straightforward option for a pure PFSS solver.

Discretization of Equation `[eqn:psi] <#eqn:psi>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, we approximate the two-dimensional Laplacian :math:`\nabla^2_\perp\psi` by

.. math::

    &\nabla^2_\perp\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = \frac{1}{S_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}\left[
    \frac{L_s^{k,j+\nicefrac{1}{2},i+1}}{L_\phi^{k,j+\nicefrac{1}{2},i+1}}\big(\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{3}{2}} - \psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}\big) -
    \frac{L_s^{k,j+\nicefrac{1}{2},i}}{L_\phi^{k,j+\nicefrac{1}{2},i}}\big(\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k,j+\nicefrac{1}{2},i-\nicefrac{1}{2}}\big) \right.

    &\left. +
    \frac{L_\phi^{k,j+1,i+\nicefrac{1}{2}}}{L_s^{k,j+1,i+\nicefrac{1}{2}}}\big(\psi^{k,j+\nicefrac{3}{2},i+\nicefrac{1}{2}} - \psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}\big) - \frac{L_\phi^{k,j,i+\nicefrac{1}{2}}}{L_s^{k,j,i+\nicefrac{1}{2}}}\big(\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k,j-\nicefrac{1}{2},i+\nicefrac{1}{2}}\big)
    \right]

As shorthand we define the quantities

.. math::

    U^{j+\nicefrac{1}{2}} = \left(\frac{L_s}{\Delta_s\Delta_\phi L_\phi}\right)^{j+\nicefrac{1}{2}}, \qquad V^j = \left(\frac{L_\phi}{\Delta_s\Delta_\phi L_s}\right)^j,

noting that these both depend on :math:`j` only.
In the code these are called ``Uc`` and ``Vg``. Then we can write our discretization as

.. math::

    \nabla^2_\perp\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}  = \frac{1}{\,\mathrm{e}^{2\rho^k}}\Big[U^{j+\nicefrac{1}{2}}\big(\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{3}{2}} - \psi^{k,j+\nicefrac{1}{2},i-\nicefrac{1}{2}} \big) + V^{j+1}\psi^{k,j+\nicefrac{3}{2},i+\nicefrac{1}{2}} + V^{j}\psi^{k,j-\nicefrac{1}{2},i+\nicefrac{1}{2}}

    - \Big(2U^{j+\nicefrac{1}{2}} + V^{j+1} + V^{j}\Big)\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}\Big].


This is the left-hand side of `[eqn:psi] <#eqn:psi>`__.

To discretize the right-hand side of `[eqn:psi] <#eqn:psi>`__, we use the approximation

.. math::

    -\frac{1}{h_\rho}\partial_\rho\left(\frac{1}{h_\rho}\partial_\rho\psi\right)^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = -\frac{c(\Delta_\rho)}{L_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}\left(\frac{\psi^{k+1,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}{L_\rho^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i+\nicefrac{1}{2}}} - \frac{\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k-1,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}{L_\rho^{k-\nicefrac{1}{2},j+\nicefrac{1}{2},i+\nicefrac{1}{2}}} \right),

where

.. math::

    c(\Delta_\rho) = \frac{2\,\mathrm{e}^{\Delta_\rho/2}}{\,\mathrm{e}^{\Delta_\rho} + 1} = \mathrm{sech}\left(\frac{\Delta_\rho}{2}\right).

Combining this with `[eqn:lapl] <#eqn:lapl>`__, we arrive at

.. math::

    U^{j+\nicefrac{1}{2}}\big(\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{3}{2}} - \psi^{k,j+\nicefrac{1}{2},i-\nicefrac{1}{2}} \big) + V^{j+1}\psi^{k,j+\nicefrac{3}{2},i+\nicefrac{1}{2}} + V^{j}\psi^{k,j-\nicefrac{1}{2},i+\nicefrac{1}{2}}
    - \Big(2U^{j+\nicefrac{1}{2}} + V^{j+1} + V^{j}\Big)\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}

    = -\frac{c(\Delta_\rho)\,\mathrm{e}^{2\rho^k}}{L_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}\left(\frac{\psi^{k+1,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}{L_\rho^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i+\nicefrac{1}{2}}} - \frac{\psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} - \psi^{k-1,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}{L_\rho^{k-\nicefrac{1}{2},j+\nicefrac{1}{2},i+\nicefrac{1}{2}}} \right).

The reader may verify algebraically that conditions `[eqn:j2] <#eqn:j2>`__ and `[eqn:j3] <#eqn:j3>`__ follow if this finite-difference equation is satisfied.

Method of solution
~~~~~~~~~~~~~~~~~~

Equation `[eqn:main] <#eqn:main>`__, together with the appropriate boundary conditions, yields a large (but sparse) system of :math:`n_\rho n_sn_\phi\times n_\rho n_sn_\phi` linear equations to solve.
Fortunately, following [vanballe2000]_, we can reduce this to a series of symmetric tridiagonal eigenvalue problems, if we look for eigenfunctions of the form

.. math::

    \psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = f^kQ_{lm}^{j+\nicefrac{1}{2}}\,\,\mathrm{e}^{2\pi I mi/n_\phi}.
    \label{eqn:eig}

Here the :math:`k` in :math:`f^k` is a power, not an index, and :math:`I` is the square root of :math:`-1` (since we already used :math:`i` and :math:`j` for indices).
This reduction will enable very efficient solution of the linear system.

Substituting `[eqn:eig] <#eqn:eig>`__ in Equation `[eqn:main] <#eqn:main>`__ gives

.. math::

    -V^{j}Q^{j-\nicefrac{1}{2}}_{lm} + \left(V^{j} + V^{j+1}+ 4U^{j+\nicefrac{1}{2}}\sin^2\left(\tfrac{\pi m}{n_\phi}\right) \right)Q^{j+\nicefrac{1}{2}}_{lm} - V^{j+1}Q^{j+\nicefrac{3}{2}}_{lm}

    = \frac{c(\Delta_\rho)\mathrm{e}^{2\rho^k}}{L_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}}\left(\frac{f - 1}{L_\rho^{k+\nicefrac{1}{2},j+\nicefrac{1}{2},i+\nicefrac{1}{2}}} - \frac{1 - f^{-1}}{L_\rho^{k-\nicefrac{1}{2},j+\nicefrac{1}{2},i+\nicefrac{1}{2}}} \right)Q_{lm}^{j+\nicefrac{1}{2}}.

The right-hand side can be simplified since the dependence on :math:`k`
cancels out. This leaves the tridiagonal eigenvalue problem

.. math:: -V^{j}Q^{j-\nicefrac{1}{2}}_{lm} + \left(V^{j} + V^{j+1}+ 4U^{j+\nicefrac{1}{2}}\sin^2\left(\tfrac{\pi m}{n_\phi}\right) \right)Q^{j+\nicefrac{1}{2}}_{lm} - V^{j+1}Q^{j+\nicefrac{3}{2}}_{lm} = \lambda_{lm}Q_{lm}^{j+\nicefrac{1}{2}},

where :math:`f` is determined from the eigenvalues :math:`\lambda_{lm}`
by solving the quadratic relation

.. math:: \lambda_{lm} = \frac{c(\Delta_\rho)}{\mathrm{e}^{\Delta_\rho/2} - \mathrm{e}^{-\Delta_\rho/2}} \left(\frac{1-f^{-1}}{1-\mathrm{e}^{-\Delta_\rho}} - \frac{f-1}{\mathrm{e}^{\Delta_\rho} - 1}\right).

This may be rearranged into the form

.. math:: f^2 - \left[1 + \mathrm{e}^{\Delta_\rho} + \mathrm{sech}\left(\frac{\Delta_\rho}{2}\right)\lambda_{lm}(\mathrm{e}^{\Delta_\rho}-1)(\mathrm{e}^{\Delta_\rho/2} - \mathrm{e}^{-\Delta_\rho/2}) \right]f + \mathrm{e}^{\Delta_\rho} = 0,

with two solutions for each :math:`l`, :math:`m` given by

.. math:: f_{lm}^+, f_{lm}^- = F_{lm} \pm \sqrt{F_{lm}^2 - \mathrm{e}^{\Delta_\rho}}, \quad \textrm{where} \quad F_{lm} = \tfrac12\Big[1 + \mathrm{e}^{\Delta_\rho} + \lambda_{lm}(\mathrm{e}^{\Delta_\rho}-1)\sinh(\Delta_\rho) \Big].

In the code, the eigenvalues are called ``lam`` and the corresponding
matrix of eigenvectors is ``Q``. The solutions :math:`f_{lm}^+` and
:math:`f_{lm}^-` are called ``ffp`` and ``ffm`` respectively.

The solution may then be written as a sum of these two sets of radial
eigenfunctions:

.. math::

    \psi^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = \sum_{l=0}^{n_s-1}\sum_{m=0}^{n_\phi-1}\Big[c_{lm}(f_{lm}^+)^k + d_{lm}(f_{lm}^-)^k) \Big] Q_{lm}^{j+\nicefrac{1}{2}}\,\mathrm{e}^{2\pi I mi/n_\phi}.
    \label{eqn:psisum}

The coefficients :math:`c_{lm}` and :math:`d_{lm}` are then determined
by the radial boundary conditions:

#. At the inner boundary :math:`\rho = 0`, where :math:`k=0`, we want
   :math:`B_\rho = -\nabla^2_\perp\psi` to match our given distribution
   :math:`g^{j+\nicefrac{1}{2},i+\nicefrac{1}{2}}`. We have

   .. math::

        B_\rho^{k,j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = \sum_{l=0}^{n_s-1}\sum_{m=0}^{n_\phi-1}\frac{\lambda_{lm}}{\mathrm{e}^{2\rho^k}}\Big[c_{lm}(f_{lm}^+)^k + d_{lm}(f_{lm}^-)^k) \Big] Q_{lm}^{j+\nicefrac{1}{2}}\mathrm{e}^{2\pi I mi/n_\phi},

   so

   .. math::

        g^{j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = \sum_{l=0}^{n_s-1}\sum_{m=0}^{n_\phi-1}\frac{\lambda_{lm}}{\mathrm{e}^{2\rho^0}}\Big[c_{lm} + d_{lm}\Big] Q_{lm}^{j+\nicefrac{1}{2}}\mathrm{e}^{2\pi I mi/n_\phi}.

   We take the discrete Fourier transform of :math:`g^{j+\nicefrac{1}{2},i+\nicefrac{1}{2}}` in :math:`i`, so that (noting :math:`\,\mathrm{e}^{\rho^0}=1`),

   .. math::

        \sum_{l=0}^{n_s-1}\sum_{m=0}^{n_\phi-1}\lambda_{lm}
        \Big[c_{lm} + d_{lm}\Big] Q_{lm}^{j+\nicefrac{1}{2}}\,\mathrm{e}^{2\pi I mi/n_\phi} = \sum_{m=0}^{n_\phi-1}b_m^{j+\nicefrac{1}{2}}\,\mathrm{e}^{2\pi I mi/n_\phi}.

   Then the orthonormality of :math:`Q_{lm}^{j+\nicefrac{1}{2}}` for different :math:`l` allows us to determine

   .. math::

        c_{lm} + d_{lm} = \frac{1}{\lambda_{lm}}\sum_{j=0}^{n_s-1}b_m^{j+\nicefrac{1}{2}}Q_{lm}^{j+\nicefrac{1}{2}}.
        \label{eqn:dc2}

#. At the source (outer) surface :math:`\rho=\ln(R_{ss})`, where :math:`k=n_\rho`, there are two options.

   #. *Radial field.*
      We impose :math:`\partial_\rho\psi = 0`, in the form
      :math:`\psi^{n_\rho,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}=\psi^{n_\rho-1,j+\nicefrac{1}{2},i+\nicefrac{1}{2}}`, which gives

      .. math::

            \frac{d_{lm}}{c_{lm}} = \frac{(f_{lm}^+)^{n_\rho} - (f_{lm}^+)^{n_\rho-1}}{(f_{lm}^-)^{n_\rho-1} - (f_{lm}^-)^{n_\rho}}.
            \label{eqn:dc1}

      (Numerically it is better to compute this ratio the other way up, to prevent overflow.)

   #. *Imposed :math:`B_r`.*
      In this case the boundary condition is treated similarly to the inner boundary.
      We require

      .. math::

            \hat{g}^{j+\nicefrac{1}{2},i+\nicefrac{1}{2}} = \sum_{l=0}^{n_s-1}\sum_{m=0}^{n_\phi-1}\frac{\lambda_{lm}}{\mathrm{e}^{2\rho^{n_\rho}}}\Big[c_{lm}(f_{lm}^+)^{n_\rho} + d_{lm}(f_{lm}^-)^{n_\rho}\Big] Q_{lm}^{j+\nicefrac{1}{2}}\mathrm{e}^{2\pi I mi/n_\phi},

      so we may again take the discrete Fourier transform to end up with

      .. math::

            c_{lm}(f_{lm}^+)^{n_\rho} + d_{lm}(f_{lm}^-)^{n_\rho} = \frac{\mathrm{e}^{2\rho^{n_\rho}}}{\lambda_{lm}}\sum_{j=0}^{n_s-1}\hat{b}_m^{j+\nicefrac{1}{2}}Q_{lm}^{j+\nicefrac{1}{2}},
            \label{eqn:dc3}

Solving `[eqn:dc1] <#eqn:dc1>`__ simultaneously with either `[eqn:dc2] <#eqn:dc2>`__ or `[eqn:dc3] <#eqn:dc3>`__ gives :math:`c_{lm}` and :math:`d_{lm}`.
These are called ``clm`` and ``dlm`` in the code.

Remark: as we increase the grid resolution, the eigenfunctions :math:`Q_{lm}^{j+\nicefrac{1}{2}}`, which are functions of :math:`\theta`, should converge to the corresponding associated Legendre polynomials :math:`P_l^m(\cos\theta)`, up to normalization.
This is illustrated in Figure `2 <#fig:Q>`__.

.. figure:: Q.png
    :name: fig:Q
    :width: 70.0%

    Comparison of :math:`P_l^m(\cos\theta)` (coloured lines) with the discrete eigenfunctions :math:`Q_{lm}^{j+\nicefrac{1}{2}}` (black dots), for :math:`m=6` and :math:`l=0,\ldots,4`, at resolution :math:`n_s=60` and :math:`n_\phi=120`.

.. [vanballe2000] `Mean Field Model for the Formation of Filament Channels on the Sun, ApJ, 539, 983-994, 2000. <http://adsabs.harvard.edu/abs/2000ApJ...539..983V>`__
