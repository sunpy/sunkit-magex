"""
Open/closed field map
=====================

Creating an open/closed field map on the solar surface.
"""
import matplotlib.colors as mcolor
import matplotlib.pyplot as plt

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
# Set the model parameters.

nrho = 40
rss = 2.5

###############################################################################
# Construct the input, and calculate the output solution.

pfss_in = pfss.Input(gong_map, nrho, rss)
pfss_out = pfss.pfss(pfss_in)

###############################################################################
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case a grid of 90 x 180 points equally gridded in theta and phi are
# chosen and traced from the source surface outwards.
#
# First, set up the tracing seeds.

# Here we will use the gong map coordinates to set the seed points
coords = sunpy.map.all_coordinates_from_map(gong_map)
ny, nx = gong_map.data.shape

lon = coords.lon.wrap_at(360*u.deg)
lat = coords.lat

seeds = SkyCoord(lon.ravel(), lat.ravel(), const.R_sun, frame=pfss_out.coordinate_frame)

###############################################################################
# Trace the field lines.

tracer = pfss.tracing.PerformanceTracer()
field_lines = tracer.trace(seeds, pfss_out)

###############################################################################
# Create a Map object for the open/closed fields.

pols = field_lines.polarities.reshape(ny,nx)

wcs_header = gong_map.wcs.to_header(relax=True)
wcs_header['NAXIS']  = 2
wcs_header['NAXIS1'] = pols.shape[1]
wcs_header['NAXIS2'] = pols.shape[0]
wcs_header['BTYPE']  = 'polarity'
wcs_header['BUNIT']  = ''           # dimensionless

pols_map = sunpy.map.Map(pols.astype(float), wcs_header)
# Adjust plot settings
pols_map.plot_settings['cmap'] =  mcolor.ListedColormap(['tab:red', 'black', 'tab:blue'])
pols_map.plot_settings['norm'] = mcolor.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ncolors=3)

###############################################################################
# Plot the result. The top plot is the input magnetogram, and the bottom plot
# shows a contour map of the the footpoint polarities, which are +/- 1 for open
# field regions and 0 for closed field regions.

fig = plt.figure()
input_map = pfss_in.map
ax = fig.add_subplot(2, 1, 1, projection=input_map)
input_map.plot(axes=ax)
ax.set_title('Input GONG magnetogram')

ax = fig.add_subplot(2, 1, 2)

ax = fig.add_subplot(2, 1, 2, projection=pols_map)
pols_map.plot()
ax.set_title('Open (blue/red) and closed (black) field')


fig.tight_layout()

plt.show()
