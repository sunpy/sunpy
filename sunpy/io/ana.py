"""
This module provides an ANA file Reader.

This is a modified version of `pyana <https://github.com/tvwerkhoven/pyana>`__.

.. warning::

    The reading and writing of ana files is not supported under Windows.

    By default, this module is not installed on platforms other than Linux (x86-64) and macOS (x86-64 and ARM64).
    See the installation guide for more info.
"""
import os

from sunpy.io._header import FileHeader
from sunpy.util.decorators import deprecated
from sunpy.util.io import HDPair

try:
    from sunpy.io import _pyana
except ImportError:
    _pyana = None

ANA_DEPRECATION_MESSAGE = (
    "The ANA reader may be removed in a future version of sunpy, "
    "please comment here if you are using this code: "
    "https://community.openastronomy.org/t/possible-deprecation-of-ana-file-readers-and-writers-in-sunpy"
)

__all__ = ['read', 'get_header', 'write']


def check_ana_installed(func):
    ana_not_installed = (
        "C extension for ANA is missing. For more details see: "
        "https://docs.sunpy.org/en/stable/installation.html#installing-without-conda"
    )
    if _pyana is None and os.environ.get("SUNPY_NO_BUILD_ANA_EXTENSION") is None:
        raise ImportError(ana_not_installed)
    return func


@deprecated(since="6.0", message=ANA_DEPRECATION_MESSAGE)
@check_ana_installed
def read(filename, debug=False, **kwargs):
    """
    Loads an ANA file and returns the data and a header in a list of (data,
    header) tuples.

    Parameters
    ----------
    filename : `str`
        Name of file to be read.
    debug : `bool`, optional
        Prints verbose debug information.
    **kwargs : `dict`
        Unused, provided for compatibility with the unified IO layer.

    Returns
    -------
    `list`
        A list of (data, header) tuples
    """
    if not os.path.isfile(filename):
        raise OSError(f"File {filename} does not exist!")
    data = _pyana.fzread(filename, debug)
    return [HDPair(data['data'], FileHeader(data['header']))]


@deprecated(since="6.0", message=ANA_DEPRECATION_MESSAGE)
@check_ana_installed
def get_header(filename, debug=False):
    """
    Loads an ANA file and only return the header consisting of the dimensions,
    size (defined as the product of all dimensions times the size of the
    datatype, this not relying on actual filesize) and comments.

    Parameters
    ----------
    filename : `str`
        Name of file to be read.
    debug : `bool`, optional
        Prints verbose debug information.

    Returns
    -------
    `list`
        A list of `~sunpy.io._header.FileHeader` headers.
    """
    data = _pyana.fzread(filename, debug)
    return [FileHeader(data['header'])]


@deprecated(since="6.0", message=ANA_DEPRECATION_MESSAGE)
@check_ana_installed
def write(filename, data, comments=None, compress=True, debug=False):
    """
    Saves a 2D `numpy.array` as an ANA file and returns the bytes written or
    ``NULL``.

    Parameters
    ----------
    filename : `str`
        Name of file to be created.
    data : `numpy.ndarray`
        The data to be stored.
    comments : `~sunpy.io._header.FileHeader`, optional
        The comments to be stored as a header.
    compress : `bool`, optional
        Compress the data with `True` (the default).
    debug : `bool`, optional
        Prints verbose debug information, defaults to `False`.

    Returns
    -------
    `str`
        A new ANA compressed archive containing the data and header.
    """
    comments = comments or ""
    return _pyana.fzwrite(filename, data, int(compress), comments, debug)
