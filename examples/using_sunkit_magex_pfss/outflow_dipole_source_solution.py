"""
Outflow field dipole source solution
=====================================

A simple example showing how to use an outflow field extrapolation to
compute the solution to a dipole source field, and comparing it to the
equivalent PFSS solution.

Outflow fields (Rice & Yeates 2021) generalise PFSS by including a solar
wind outflow speed profile in the solve. Setting the outflow speed to zero
everywhere recovers the PFSS solution exactly; here we use the default
outflow speed profile (optimised to match eclipse observations) instead,
to see how it changes the solution compared to a plain PFSS extrapolation.
See :ref:`code-heritage` for more on how this method relates to PFSS.
"""
# sphinx_gallery_thumbnail_number = 3

import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

import sunpy.map

from sunkit_magex import pfss

###############################################################################
# As in the PFSS dipole example, set up a regular 2D grid in (phi, s), where
# s = cos(theta), and calculate a dipole boundary condition on it.

nphi = 180
ns = 90
phi = np.linspace(0, 2 * np.pi, nphi)
s = np.linspace(-1, 1, ns)
s, phi = np.meshgrid(s, phi)


def dipole_Br(r, s):
    return 2 * s / r**3


br = dipole_Br(1, s)

###############################################################################
# Choose the number of radial grid points and the source surface radius.

nrho = 30
rss = 2.5

###############################################################################
# Construct the input map, and from it a `~sunkit_magex.pfss.OutflowInput`
# (using the default outflow speed profile) and a plain
# `~sunkit_magex.pfss.Input` for the equivalent PFSS solution.

header = pfss.utils.carr_cea_wcs_header(Time('2020-1-1'), br.shape)
input_map = sunpy.map.Map((br.T, header))

outflow_in = pfss.OutflowInput(input_map, nrho, rss)
pfss_in = pfss.Input(input_map, nrho, rss)

###############################################################################
# The outflow speed profile used above can be inspected directly: it is
# stored (in code units) on ``outflow_in.vg``, as a function of altitude
# ``outflow_in.grid.rg`` (in units of solar radii).

fig, ax = plt.subplots()
ax.plot(np.exp(outflow_in.grid.rg), outflow_in.vg)
ax.set_xlabel(r'Altitude $r$ ($R_\odot$)')
ax.set_ylabel('Outflow speed (code units)')
ax.set_title('Default outflow speed profile')

###############################################################################
# Calculate both solutions.

outflow_out = pfss.outflow(outflow_in)
pfss_out = pfss.pfss(pfss_in)

###############################################################################
# Using the Output objects we can plot the source surface field, and the
# polarity inversion line, for both solutions side by side.

fig, axs = plt.subplots(1, 2, figsize=(10, 5),
                        subplot_kw={'projection': pfss_out.source_surface_br})

for ax, out, title in zip(axs, [pfss_out, outflow_out], ['PFSS', 'Outflow']):
    ss_br = out.source_surface_br
    ss_br.plot(axes=ax)
    ax.plot_coord(out.source_surface_pils[0])
    ax.set_title(f'{title} source surface field')

fig.tight_layout()

###############################################################################
# Finally, trace some field lines from both solutions and plot them in the
# meridional plane, as in the PFSS dipole example. 32 points equally spaced
# in theta are chosen and traced from just above the solar surface.

r = 1.01 * const.R_sun
lon = np.pi / 2 * u.rad
lat = np.linspace(-np.pi / 2, np.pi / 2, 33) * u.rad

tracer = pfss.tracing.PerformanceTracer()

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

for ax, out, title in zip(axs, [pfss_out, outflow_out], ['PFSS', 'Outflow']):
    ax.set_aspect('equal')
    seeds = SkyCoord(lon, lat, r, frame=out.coordinate_frame)
    field_lines = tracer.trace(seeds, out)

    for field_line in field_lines:
        coords = field_line.coords
        coords.representation_type = 'cartesian'
        color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
        ax.plot(coords.y / const.R_sun, coords.z / const.R_sun, color=color)

    ax.add_patch(mpatch.Circle((0, 0), 1, color='k', fill=False))
    ax.add_patch(mpatch.Circle((0, 0), rss, color='k', linestyle='--', fill=False))
    ax.set_title(f'{title} solution for a dipole source field')

fig.tight_layout()

plt.show()
