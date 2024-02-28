"""
Re-projecting from CAR to CEA
-----------------------------

The sunkit_magex.pfss solver takes a cylindrical-equal-area (CEA) projected magnetic field
map as input, which is equally spaced in sin(latitude). Some synoptic field
maps are equally spaced in latitude, a plate carée (CAR) projection, and need
reprojecting.

This example shows how to use the `sunkit_magex.pfss.utils.car_to_cea` function to
reproject a CAR projection to a CEA projection that sunkit_magex.pfss can take as input.
"""
import matplotlib.pyplot as plt

from sunkit_magex.pfss import sample_data, utils

###############################################################################
# Load a sample ADAPT map, which has a CAR projection
adapt_maps = utils.load_adapt(sample_data.get_adapt_map())
adapt_map_car = adapt_maps[0]

###############################################################################
# Re-project into a CEA projection
adapt_map_cea = utils.car_to_cea(adapt_map_car)

###############################################################################
# Plot the original map and the reprojected map
plt.figure()
adapt_map_car.plot()
plt.figure()
adapt_map_cea.plot()

plt.show()
