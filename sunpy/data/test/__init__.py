"""SunPy test data files"""
from __future__ import absolute_import

import os
import glob

from astropy.utils.data import get_pkg_data_filename

import sunpy

import fnmatch
import re

__all__ = ['rootdir', 'file_list', 'get_test_filepath', 'test_data_filenames']

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


file_list = glob.glob(os.path.join(rootdir, '*.[!p]*'))


def test_data_filenames():
    """
    Return a list of all test files in ``data/test`` directory.

    Return
    ------
    get_all_test_filepath : `list`
        The name of all test files in ``data/test`` directory.

    """
    test_data_filenames_list = []
    excludes = ['*.pyc', '*/__*__', '*.py']
    excludes = r'|'.join([fnmatch.translate(x) for x in excludes]) or r'$.'

    for root, dirs, files in os.walk(rootdir):

        dirs[:] = [os.path.join(root, d) for d in dirs]
        dirs[:] = [d for d in dirs if not re.match(excludes, d)]

        files = [os.path.join(root, f) for f in files]
        files = [f for f in files if not re.match(excludes, f)]

        files = [file.replace(rootdir + '/', '') for file in files]

        test_data_filenames_list.extend(files)

    return test_data_filenames_list
