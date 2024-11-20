"""
This module provides a generic file reader for internal use.
"""
import re
import gzip
import pathlib

import fsspec

from sunpy.util.io import is_uri

try:
    from . import _fits as fits
except ImportError:
    fits = None

try:
    from . import _jp2
except ImportError:
    _jp2 = None

try:
    from . import ana
except ImportError:
    ana = None


__all__ = ['read_file', 'read_file_header', 'write_file', 'detect_filetype']

# File formats supported by sunpy.io
_KNOWN_EXTENSIONS = {
    ('fts', 'fits'): 'fits',
    ('jp2', 'j2k', 'jpc', 'jpt'): 'jp2',
    ('fz', 'f0'): 'ana'
}


class Readers(dict):
    # Define a dict which raises a custom error message if the value is None
    def __init__(self, *args):
        dict.__init__(self, *args)

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        if val is None:
            raise ReaderError(
                f"The Reader sunpy.io.{key} is not available, "
                "please check that you have the required dependencies "
                "installed."
            )
        return val

# File readers supported by sunpy.io
_READERS = Readers({
    'fits': fits,
    'jp2': _jp2,
    'ana': ana
})


def _read(filepath, function_name, filetype=None, **kwargs):
    """
    This functions provides the logic paths for reading a file.

    It is trying to handle, reading files with known extensions and reading files with known readers.

    We also want the detection logic to first try to detect the filetype based on the content and then
    fallback to using extension.

    Parameters
    ----------
    filepath : pathlib.Path, str
        The file to be read.
    function_name : {'read' | 'get_header'}
        The name of the function to call on the reader.
    filetype : {'jp2' | 'fits' | 'ana'}, optional
        Supported reader or extension to manually specify the filetype.
        Supported readers are ('jp2', 'fits', 'ana')
    **kwargs : `dict`
        Additional keyword arguments are handed to file specific reader.

    Returns
    -------
    pairs : `list`
        A list of (data, header) tuples.
    """
    filepath = str(filepath)
    if filetype is not None:
        return getattr(_READERS[filetype], function_name)(filepath, **kwargs)
    try:
        readername = detect_filetype(filepath)
        if readername in _READERS.keys():
            return getattr(_READERS[readername], function_name)(filepath, **kwargs)
        readername = None
    except UnrecognizedFileTypeError:
        readername = None
    for extension, name in _KNOWN_EXTENSIONS.items():
        if filepath.endswith(extension) or filetype in extension:
            readername = name
    if readername is not None:
        return getattr(_READERS[readername], function_name)(filepath, **kwargs)
    raise UnrecognizedFileTypeError("The requested filetype is not currently supported by sunpy.")


def read_file(filepath, filetype=None, **kwargs):
    """
    Automatically determine the filetype and read the file.

    Parameters
    ----------
    filepath : `str`, path-like
        The file to be read.
    filetype : `str`, optional
        Supported reader or extension to manually specify the filetype.
        Supported readers are ('jp2', 'fits', 'ana')
    memmap : `bool`, optional
        Should memory mapping be used, i.e. keep data on disk rather than in RAM.
        This is currently only supported by the FITS reader.
    **kwargs : `dict`
        All extra keyword arguments are passed to ``.read`` for the file specific reader.

    Returns
    -------
    pairs : `list`
        A list of (data, header) tuples.
    """
    return _read(filepath, 'read', filetype, **kwargs)


def read_file_header(filepath, filetype=None, **kwargs):
    """
    Reads the header from a given file.

    This should always return a instance of `~sunpy.io._header.FileHeader`.

    Parameters
    ----------
    filepath : `str`
        The file from which the header is to be read.
    filetype : `str`
        Supported reader or extension to manually specify the filetype.
        Supported readers are ('jp2', 'fits').
    **kwargs : `dict`
        All extra keyword arguments are passed to ``.get_header`` for the file specific reader.

    Returns
    -------
    headers : `list`
        A list of headers.
    """
    return _read(filepath, 'get_header', filetype, **kwargs)


def write_file(fname, data, header, filetype='auto', **kwargs):
    """
    Write a file from a data & header pair using one of the defined file types.

    Parameters
    ----------
    fname : `str`
        Filename of file to save.
    data : `numpy.ndarray`
        Data to save to a fits file.
    header : `collections.OrderedDict`
        Meta data to save with the data.
    filetype : `str`, {'auto', 'fits', 'jp2'}, optional
        Filetype to save if ``auto`` the  filename extension will
        be detected, else specify a supported file extension.
    **kwargs : `dict`
        All extra keyword arguments are passed to ``.write`` for the file specific reader.

    Notes
    -----
    * This routine currently only supports saving a single HDU.
    """
    if filetype == 'auto':
        # Get the extension without the leading dot
        filetype = pathlib.Path(fname).suffix[1:]
    for extension, readername in _KNOWN_EXTENSIONS.items():
        if filetype in extension:
            return _READERS[readername].write(fname, data, header, **kwargs)
    raise ValueError(f"The filetype provided ({filetype}) is not supported")


def detect_filetype(filepath, **kwargs):
    """
    Attempts to determine the type of file a given filepath is.

    Parameters
    ----------
    filepath : `str`, `pathlib.Path`
        Where the file is.
    **kwargs : `dict`
        Additional keyword arguments.

    Returns
    -------
    filetype : `str`
        The type of file.
    """
    if str(filepath).startswith('http') or  str(filepath).startswith('ftp'):
        return None
    if is_uri(filepath):
        fsspec_kw = kwargs.get("fsspec_kwargs", {})
        try:
            fileobj = fsspec.open(filepath, 'rb', **fsspec_kw).open()
        except Exception:
            return None
    else:
        fileobj = open(filepath, 'rb')
    with fileobj as fp:
        line1 = fp.readline()
        line2 = fp.readline()
        # Some FITS files do not have line breaks at the end of header cards.
        fp.seek(0)
        first80 = fp.read(80)
        # First 8 bytes of netcdf4/hdf5 to determine filetype as have same sequence
        fp.seek(0)
        first_8bytes = fp.read(8)
        # First 4 bytes of CDF
        fp.seek(0)
        cdf_magic_number = fp.read(4).hex()

    # For ASDF files
    if first80.startswith(b"#ASDF"):
        return "asdf"

    # FITS
    # Checks for gzip signature.
    # If found, decompresses first few bytes and checks for FITS
    if first80[:3] == b"\x1f\x8b\x08":
        with gzip.open(filepath, 'rb') as fp:
            first80 = fp.read(80)

    # Check for "KEY_WORD  =" at beginning of file
    match = re.match(br"[A-Z0-9_]{0,8} *=", first80)
    if match is not None:
        return 'fits'

    # JPEG 2000
    # Checks for one of two signatures found at beginning of all JP2 files.
    # Adapted from ExifTool
    # [1] https://www.sno.phy.queensu.ca/~phil/exiftool/
    # [2] http://www.hlevkin.com/Standards/fcd15444-2.pdf
    # [3] http://www.hlevkin.com/Standards/fcd15444-1.pdf
    jp2_signatures = [b"\x00\x00\x00\x0cjP  \x0d\x0a\x87\x0a",
                      b"\x00\x00\x00\x0cjP\x1a\x1a\x0d\x0a\x87\x0a"]
    for sig in jp2_signatures:
        if line1 + line2 == sig:
            return 'jp2'

    # netcdf4 and hdf5 files
    if first_8bytes == b'\x89HDF\r\n\x1a\n':
        return 'hdf5'

    if cdf_magic_number in ['cdf30001', 'cdf26002', '0000ffff']:
        return 'cdf'

    raise UnrecognizedFileTypeError("The requested filetype is not currently supported by sunpy.")


class UnrecognizedFileTypeError(OSError):
    """
    Exception to raise when an unknown file type is encountered.
    """


class ReaderError(ImportError):
    """
    Exception to raise when a reader errors.
    """
