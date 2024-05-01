import importlib.resources
from os import PathLike

import numpy as np

from astropy.io import fits

import sunpy.map


def get_test_data(filename, package="sunkit_magex.tests.data") -> PathLike:
   with importlib.resources.as_file(importlib.resources.files(package) / filename) as file_path:
       return file_path


def get_header_data_from_header(filename, package="sunkit_magex.tests.data"):
    filepath = get_test_data(filename, package=package)
    header = fits.Header.fromtextfile(filepath)
    data = np.random.rand(*(header[f'NAXIS{n}'] for n in range(header["NAXIS"], 0, -1)))
    if 'BITPIX' in header:
        data = data.astype(fits.BITPIX2DTYPE[header['BITPIX']])
    return header, data


def get_fitsfile_from_header(filename, outputfile, *, package="sunkit_magex.tests.data", hdutype=fits.PrimaryHDU):
    """
    Generate a FITS file from a test header.
    """
    header, data = get_header_data_from_header(filename, package=package)
    hdu = hdutype(header=header, data=data)
    hdu.writeto(outputfile)
    return outputfile


def get_dummy_map_from_header(filename, package="sunkit_magex.tests.data"):
    """
    Generate a dummy `~sunpy.map.Map` from header-only test data.

    The "image" will be random numbers with the correct shape
    as specified by the header.
    """
    header, data = get_header_data_from_header(filename, package=package)
    # NOTE: by reading straight from the data header pair, we are skipping
    # the fixes that are applied in sunpy.io._fits, e.g. the waveunit fix
    # Would it be better to write this whole map back to a temporary file and then
    # read it back in by passing in the filename instead?
    return sunpy.map.Map(data, header)
