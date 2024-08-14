"""
This package contains all of SunPy's test data.
"""
import os
import re
import glob
import fnmatch
from pathlib import Path

import numpy as np

import astropy.io.fits
from astropy.utils.data import get_pkg_data_filename

import sunpy
import sunpy.io._fits as _fits

__all__ = [
    'rootdir',
    'file_list',
    'get_test_filepath',
    'get_test_data_filenames',
    'get_dummy_map_from_header',
    'write_header_file_from_image_file',
]

rootdir = Path(os.path.dirname(sunpy.__file__)) / "data" / "test"
file_list = glob.glob(os.path.join(rootdir, '*.[!p]*'))


def get_test_filepath(filename, package="sunpy.data.test", **kwargs):
    """
    Return the full path to a test file in the ``data/test`` directory.

    Parameters
    ----------
    filename : `str`
        The name of the file inside the ``data/test`` directory.
    package : `str`, optional
        The package in which to look for the file. Defaults to "sunpy.data.test".

    Returns
    -------
    filepath : `str`
        The full path to the file.

    Notes
    -----
    This is a wrapper around `astropy.utils.data.get_pkg_data_filename` which
    sets the ``package`` kwarg to be 'sunpy.data.test`.
    """
    if isinstance(filename, Path):
        # NOTE: get_pkg_data_filename does not accept Path objects
        filename = filename.as_posix()
    return get_pkg_data_filename(filename, package=package, **kwargs)


def get_test_data_filenames():
    """
    Return a list of all test files in ``data/test`` directory.

    This ignores any ``py``, ``pyc`` and ``__*__`` files in these directories.

    Returns
    -------
    `list`
        The name of all test files in ``data/test`` directory.
    """
    get_test_data_filenames_list = []
    excludes = ['*.pyc', '*'+os.path.sep+'__*__', '*.py']
    excludes = r'|'.join([fnmatch.translate(x) for x in excludes]) or r'$.'

    for root, _, files in os.walk(rootdir):
        files = [Path(root) / f for f in files]
        files = [f for f in files if not re.match(excludes, str(f))]
        get_test_data_filenames_list.extend(files)

    return get_test_data_filenames_list


def write_image_file_from_header_file(header_file, fits_directory):
    """
    Given a header-only file ``header_file``, write a dummy FITS file
    with the same name to ``fits_directory``

    Parameters
    ----------
    header_file: `~pathlib.Path`
    fits_directory: `~pathlib.Path`

    Returns
    -------
    fits_file: `~pathlib.Path`
        Path to dummy FITS file
    """
    header = astropy.io.fits.Header.fromtextfile(header_file)
    shape = [header[f'naxis{i+1}'] for i in range(header['naxis'])]
    data = np.random.rand(*shape[::-1])
    if 'BITPIX' in header:
        data = data.astype(astropy.io.fits.BITPIX2DTYPE[header['BITPIX']])
    hdu = astropy.io.fits.PrimaryHDU(data=data, header=header, do_not_scale_image_data=True, scale_back=True)
    fits_file = os.fspath(fits_directory.joinpath(header_file.with_suffix('.fits').name))
    hdu.writeto(fits_file)
    return fits_file


def get_dummy_map_from_header(filename, package="sunpy.data.test"):
    """
    Generate a dummy `~sunpy.map.Map` from header-only test data.

    The "image" will be random numbers with the correct shape
    as specified by the header.
    """
    import sunpy.map

    filepath = get_test_filepath(filename, package=package)
    header = _fits.format_comments_and_history(astropy.io.fits.Header.fromtextfile(filepath))
    data = np.random.rand(header['NAXIS2'], header['NAXIS1'])
    if 'BITPIX' in header:
        data = data.astype(astropy.io.fits.BITPIX2DTYPE[header['BITPIX']])
    # NOTE: by reading straight from the data header pair, we are skipping
    # the fixes that are applied in sunpy.io._fits, e.g. the waveunit fix
    # Would it be better to write this whole map back to a temporary file and then
    # read it back in by passing in the filename instead?
    return sunpy.map.Map(data, header)


def write_header_file_from_image_file(image_path, header_path, hdu=0, **kwargs):
    """
    Convert a FITS file containing an image and a header
    to a plaintext file containing only the string representation
    of the header.

    This is to be used in the context of generating header-only test data
    files.

    Parameters
    -----------
    image_path : path-like
        Path to original image data FITS file
    header_path : path-like
        Path to header-only test data plaintext file
    hdu : `int`, optional
        HDU index of the header to write the header for

    Additional Parameters
    ---------------------
    Keyword arguments accepted by `~astropy.io.fits.Header.totextfile`

    See Also
    --------
    get_dummy_map_from_header
    """
    with astropy.io.fits.open(image_path) as hdul:
        hdul[hdu].header.totextfile(header_path, **kwargs)
