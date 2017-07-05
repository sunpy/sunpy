"""SunPy test data files"""
from __future__ import absolute_import

import os
import glob

from astropy.utils.data import get_pkg_data_filename

import sunpy

__all__ = ['rootdir', 'file_list', 'get_test_filepath',
           'print_test_filenames']

rootdir = os.path.join(os.path.dirname(sunpy.__file__), "data", "test")


def get_test_filepath(filename, **kwargs):
    """
    Return the full path to a test file in the ``data/test`` directory.

    Parameters
    ----------
    filename : `str`
        The name of the file inside the ``data/test`` directory.

    Return
    ------
    filepath : `str`
        The full path to the file.

    See Also
    --------

    astropy.utils.data.get_pkg_data_filename : Get package data filename

    Notes
    -----

    This is a wrapper around `astropy.utils.data.get_pkg_data_filename` which
    sets the ``package`` kwarg to be 'sunpy.data.test`.

    """
    return get_pkg_data_filename(filename, package="sunpy.data.test", **kwargs)


def print_test_filenames(directory=rootdir):
    """
    Prints the file names of all of the test data files under the ``data/test``
    directory.

    Parameters
    ----------
    directory : 'str'
        The root directory to use to find and print all files within.

    """

    for (dirpath, _, filenames) in os.walk(directory):
        for fname in filenames:
            if not fname.startswith("__"):
                if dirpath == rootdir:
                    print(fname)
                else:
                    print(os.path.join(os.path.basename(dirpath), fname))


file_list = glob.glob(os.path.join(rootdir, '*.[!p]*'))
