"""
This package contains all of SunPy's test data.
"""
import os
import re
import glob
import fnmatch
from pathlib import Path

from astropy.utils.data import get_pkg_data_filename

import sunpy

__all__ = ['rootdir', 'file_list', 'get_test_filepath', 'test_data_filenames']

rootdir = Path(os.path.dirname(sunpy.__file__)) / "data" / "test"
file_list = glob.glob(os.path.join(rootdir, '*.[!p]*'))


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

    Notes
    -----
    This is a wrapper around `astropy.utils.data.get_pkg_data_filename` which
    sets the ``package`` kwarg to be 'sunpy.data.test`.
    """
    return get_pkg_data_filename(filename, package="sunpy.data.test", **kwargs)


def test_data_filenames():
    """
    Return a list of all test files in ``data/test`` directory.

    This ignores any ``py``, ``pyc`` and ``__*__`` files in these directories.

    Return
    ------
    `list`
        The name of all test files in ``data/test`` directory.
    """
    test_data_filenames_list = []
    excludes = ['*.pyc', '*'+os.path.sep+'__*__', '*.py']
    excludes = r'|'.join([fnmatch.translate(x) for x in excludes]) or r'$.'

    for root, dirs, files in os.walk(rootdir):
        files = [Path(root) / f for f in files]
        files = [f for f in files if not re.match(excludes, str(f))]
        test_data_filenames_list.extend(files)

    return test_data_filenames_list
