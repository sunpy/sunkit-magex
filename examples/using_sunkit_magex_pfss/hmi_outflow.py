"""
HMI outflow field solutions
============================

Calculating an outflow field solution from a HMI synoptic map, and
comparing it to the equivalent PFSS solution.

Outflow fields (Rice & Yeates 2021) generalise the PFSS model by
including a solar wind outflow speed profile in the solve. Setting the
outflow speed to zero everywhere recovers the PFSS solution exactly, which
this example uses as a point of comparison. See :ref:`code-heritage` for
more on how this method relates to PFSS.
"""
import os

import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
import numpy as np

import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

from sunkit_magex import pfss

###############################################################################
# Set up the search.
#
# The synoptic maps are labelled by Carrington rotation number instead of time.

series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
crot = a.jsoc.PrimeKey('CAR_ROT', 2210)

###############################################################################
# Do the search.
#
# If you use this code, please replace this email address
# with your own one, registered here:
# http://jsoc.stanford.edu/ajax/register_email.html

result = Fido.search(series, crot, a.jsoc.Notify(os.environ["JSOC_EMAIL"]))
files = Fido.fetch(result)

###############################################################################
# Read in a file, and resample it down to a smaller size so that the
# solution can be calculated in a reasonable time (see the PFSS HMI
# example for more on why this is needed).

hmi_map = sunpy.map.Map(files[0])
hmi_map = hmi_map.resample([360, 180] * u.pix)

###############################################################################
# Calculate the outflow field solution, using the default outflow speed
# profile (optimised to match eclipse observations, see Rice & Yeates 2021).
# For comparison, also calculate the equivalent zero-outflow (PFSS) solution
# on the same grid.

nrho = 35
rss = 2.5

outflow_in = pfss.OutflowInput(hmi_map, nrho, rss)
outflow_out = pfss.outflow(outflow_in)

pfss_in = pfss.Input(hmi_map, nrho, rss)
pfss_out = pfss.pfss(pfss_in)

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
# Using the Output objects we can plot the source surface field for both
# solutions, and the polarity inversion line, side by side.

fig, axs = plt.subplots(1, 2, figsize=(10, 5),
                        subplot_kw={'projection': pfss_out.source_surface_br})

for ax, out, title in zip(axs, [pfss_out, outflow_out], ['PFSS', 'Outflow']):
    ss_br = out.source_surface_br
    ss_br.plot(axes=ax)
    ax.plot_coord(out.source_surface_pils[0])
    ax.set_title(f'{title} source surface field')

fig.tight_layout()

###############################################################################
# Finally, trace a grid of field lines from both solutions using the
# rust-backed `~sunkit_magex.pfss.tracing.PerformanceTracer`. This tracer
# only ever touches ``output.bg`` and ``output.grid``, so it works
# unmodified on outflow field results, exactly as it does on PFSS ones. As
# in the open/closed field map example, plot the footpoint polarities:
# +/- 1 for open field regions, and 0 for closed field regions.

tracer = pfss.tracing.PerformanceTracer()

nsteps = 45
lon_1d = np.linspace(0, 2 * np.pi, nsteps * 2 + 1)
lat_1d = np.arcsin(np.linspace(-1, 1, nsteps + 1))
lon, lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')
lon, lat = lon * u.rad, lat * u.rad
seeds = SkyCoord(lon.ravel(), lat.ravel(), const.R_sun, frame=outflow_out.coordinate_frame)

outflow_field_lines = tracer.trace(seeds, outflow_out)
pfss_field_lines = tracer.trace(seeds, pfss_out)

print(f'Outflow: {len(outflow_field_lines.open_field_lines)} open field lines '
     f'out of {len(outflow_field_lines)}')
print(f'PFSS: {len(pfss_field_lines.open_field_lines)} open field lines '
     f'out of {len(pfss_field_lines)}')

fig, axs = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
cmap = mcolor.ListedColormap(['tab:red', 'black', 'tab:blue'])
norm = mcolor.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ncolors=3)

for ax, field_lines, title in zip(axs, [pfss_field_lines, outflow_field_lines], ['PFSS', 'Outflow']):
    pols = field_lines.polarities.reshape(2 * nsteps + 1, nsteps + 1).T
    ax.contourf(np.rad2deg(lon_1d), np.sin(lat_1d), pols, norm=norm, cmap=cmap)
    ax.set_ylabel('sin(latitude)')
    ax.set_title(f'{title}: open (blue/red) and closed (black) field')
    ax.set_aspect(0.5 * 360 / 2)

axs[-1].set_xlabel('Longitude (degrees)')
fig.tight_layout()

plt.show()
