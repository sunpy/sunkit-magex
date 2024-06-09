"""
Custom `sunpy.map.GenericMap` sub-classes for different magnetogram sources.
"""
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
