"""
Custom `sunpy.map.GenericMap` sub-classes for different magnetogram sources.
"""
import astropy.units as u
import sunpy.map

__all__ = ['ADAPTMap']


class ADAPTMap(sunpy.map.GenericMap):
    def __init__(self, data, header, **kwargs):
        if 'date-obs' not in header:
            header['date-obs'] = header['maptime']
        # Fix CTYPE
        if header['ctype1'] == 'Long':
            header['ctype1'] = 'CRLN-CAR'
        if header['ctype2'] == 'Lat':
            header['ctype2'] = 'CRLT-CAR'

        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an ADAPT map."""
        return header.get('model') == 'ADAPT'


def _observer_coord_meta(observer_coord):
    """
    Convert an observer coordinate into FITS metadata.
    """
    new_obs_frame = sunpy.coordinates.HeliographicStonyhurst(
        obstime=observer_coord.obstime)
    observer_coord = observer_coord.transform_to(new_obs_frame)

    new_meta = {'hglt_obs': observer_coord.lat.to_value(u.deg)}
    new_meta['hgln_obs'] = observer_coord.lon.to_value(u.deg)
    new_meta['dsun_obs'] = observer_coord.radius.to_value(u.m)
    return new_meta


def _earth_obs_coord_meta(obstime):
    """
    Return metadata for an Earth obeserver coordinate.
    """
    return _observer_coord_meta(sunpy.coordinates.get_earth(obstime))
