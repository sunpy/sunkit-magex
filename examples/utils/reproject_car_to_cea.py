"""
Re-projecting from CAR to CEA
=============================

The `sunkit_magex.pfss` solver takes a cylindrical-equal-area (CEA) projected magnetic field
map as input, which is equally spaced in sin(latitude). Some synoptic field
maps are equally spaced in latitude, a plate carée (CAR) projection, and need
reprojecting.

This example shows how to use the `sunkit_magex.pfss.utils.car_to_cea` function to
reproject a CAR projection to a CEA projection that `sunkit_magex.pfss` can take as input.
"""
import matplotlib.pyplot as plt

import sunpy.map

from sunkit_magex import pfss

###############################################################################
# Load a sample ADAPT map, which has a CAR projection.

adapt_map_car = sunpy.map.Map(pfss.sample_data.get_adapt_map(), hdus=0)

###############################################################################
# Re-project into a CEA projection.

adapt_map_cea = pfss.utils.car_to_cea(adapt_map_car)

###############################################################################
# Plot the original map and the reprojected map.

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 2, 1, projection=adapt_map_car)
ax1 = fig.add_subplot(1, 2, 2, projection=adapt_map_cea)

adapt_map_car.plot(axes=ax)
adapt_map_cea.plot(axes=ax1)

fig.tight_layout()

plt.show()
