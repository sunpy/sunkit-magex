import copy

import numpy as np

import sunpy.map

import sunkit_magex.pfss.utils
from sunkit_magex.pfss.grid import Grid


class Input:
    r"""
    Input to PFSS modelling.


    Parameters
    ----------
    br : sunpy.map.GenericMap
        Boundary condition of radial magnetic field at the inner surface.
        Note that the data *must* have a cylindrical equal area projection.
    nr : int
        Number of cells in the radial direction to calculate the PFSS solution
        on.
    rss : float
        Radius of the source surface, as a fraction of the solar radius.
    br_outer : sunpy.map.GenericMap, String
        Boundary condition of radial magnetic field at the outer surface.

    Notes
    -----
    The input must be on a regularly spaced grid in :math:`\phi` and
    :math:`s = \cos (\theta)`. See `sunkit_magex.pfss.grid` for more
    information on the coordinate system.
    """
    def __init__(self, br, nr, rss, br_outer="radial"):
        if not isinstance(br, sunpy.map.GenericMap):
            raise ValueError('br must be a sunpy Map')
        if np.any(~np.isfinite(br.data)):
            raise ValueError('At least one value in the input is NaN or '
                             'non-finite. The input must consist solely of '
                             'finite values.')

        sunkit_magex.pfss.utils.is_cea_map(br, error=True)
        sunkit_magex.pfss.utils.is_full_sun_synoptic_map(br, error=True)


        if isinstance(br_outer, sunpy.map.GenericMap):
            if np.any(~np.isfinite(br_outer.data)):
                raise ValueError('At least one value in the input is NaN or '
                                 'non-finite. The input must consist solely of '
                                 'finite values.')
            if br.dimensions != br_outer.dimensions:
                raise ValueError('br and br_outer must have the same dimensions')

        self._map_in = copy.deepcopy(br)
        self.br = self.map.data

        if isinstance(br_outer, sunpy.map.GenericMap):
            self.br_outer = br_outer.data
        else:
            self.br_outer = br_outer

        # Force some nice defaults
        self._map_in.plot_settings['cmap'] = 'RdBu'
        lim = np.nanmax(np.abs(self._map_in.data))
        self._map_in.plot_settings['vmin'] = -lim
        self._map_in.plot_settings['vmax'] = lim

        ns = self.br.shape[0]
        nphi = self.br.shape[1]
        self._grid = Grid(ns, nphi, nr, rss)

    @property
    def map(self):
        """
        :class:`sunpy.map.GenericMap` representation of the input.
        """
        return self._map_in

    @property
    def grid(self):
        """
        `~sunkit_magex.pfss.grid.Grid` that the PFSS solution for this input is
        calculated on.
        """
        return self._grid
