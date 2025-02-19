"""
GONG PFSS extrapolation
=======================

Calculating PFSS solution for a GONG synoptic magnetic field map.
"""
# sphinx_gallery_thumbnail_number = 4

import matplotlib.pyplot as plt
import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map

from sunkit_magex import pfss

###############################################################################
# Load a GONG magnetic field map.

gong_fname = pfss.sample_data.get_gong_map()
gong_map = sunpy.map.Map(gong_fname)

###############################################################################
# The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of rho grid points, and the source surface radius.

nrho = 35
rss = 2.5

###############################################################################
# From the boundary condition, number of radial grid points, and source
# surface, we now construct an Input object that stores this information.

pfss_in = pfss.Input(gong_map, nrho, rss)

###############################################################################
# Using the Input object, plot the input field.

input_map = pfss_in.map
fig = plt.figure()
ax = fig.add_subplot(projection=input_map)
input_map.plot(axes=ax)
plt.colorbar()
ax.set_title('Input field')

###############################################################################
# Now calculate the PFSS solution.

pfss_out = pfss.pfss(pfss_in)

###############################################################################
# Using the Output object we can plot the source surface field, and the
# polarity inversion line.

ss_br = pfss_out.source_surface_br

fig = plt.figure()
ax = fig.add_subplot(projection=ss_br)

ss_br.plot(axes=ax)
# Plot the polarity inversion line
ax.plot_coord(pfss_out.source_surface_pils[0])
plt.colorbar()
ax.set_title('Source surface magnetic field')

###############################################################################
# It is also easy to plot the magnetic field at an arbitrary height within
# the PFSS solution.

# Get the radial magnetic field at a given height
ridx = 15
br = pfss_out.bc[0][:, :, ridx]
# Create a sunpy Map object using output WCS
br = sunpy.map.Map(br.T, pfss_out.source_surface_br.wcs)
# Get the radial coordinate
r = np.exp(pfss_out.grid.rc[ridx])

fig = plt.figure()
ax = fig.add_subplot(projection=br)

br.plot(cmap='RdBu')
plt.colorbar()
ax.set_title('$B_{r}$ ' + f'at r={r:.2f}' + '$r_{\\odot}$')

###############################################################################
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case 64 points equally gridded in theta and phi are chosen and
# traced from the source surface outwards.

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

tracer = pfss.tracing.PerformanceTracer()
r = 1.2 * const.R_sun
lat = np.linspace(-np.pi / 2, np.pi / 2, 8, endpoint=False)
lon = np.linspace(0, 2 * np.pi, 8, endpoint=False)
lat, lon = np.meshgrid(lat, lon, indexing='ij')
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad

seeds = SkyCoord(lon, lat, r, frame=pfss_out.coordinate_frame)
field_lines = tracer.trace(seeds, pfss_out)
for field_line in field_lines:
    color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
    coords = field_line.coords
    coords.representation_type = 'cartesian'
    ax.plot(coords.x / const.R_sun,
            coords.y / const.R_sun,
            coords.z / const.R_sun,
            color=color, linewidth=1)
ax.set_title('PFSS solution')

plt.show()
