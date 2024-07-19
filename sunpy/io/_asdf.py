from collections import OrderedDict

import asdf
from astropy.utils.introspection import minversion

from sunpy.io._header import FileHeader

__all__ = ["write", "read", "get_header"]

# use memmap for and after asdf 3.1.0 to avoid a warning
# for the deprecated copy_arrays
if minversion(asdf, "3.1.0"):
    _NO_MEMMAP_KWARGS = {"memmap": False, "lazy_load": False}
else:
    _NO_MEMMAP_KWARGS = {"copy_arrays": True, "lazy_load": False}


def write(filepath, data, header, **kwargs):
    """
    Take ``(data, header)`` pairs and save it to an ASDF file.

    Parameters
    ----------
    filepath : pathlib.Path, str
        File name, with extension.
    data : `numpy.ndarray`
        N-dimensional data array.
    header : `dict`
        A header dictionary.
    """
    asdf.AsdfFile({"object": {"meta": OrderedDict(header), "data": data}}).write_to(str(filepath), **kwargs)


def read(filepath, **kwargs):
    """
    Read an ASDF file.

    Parameters
    ----------
    filepath : pathlib.Path, `str`
        The ASDF file to be read.

    Returns
    -------
    `list`
        A list of ``(data, header)`` tuples.
    """
    # Provide the required options to prevent asdf from memory mapping
    # or lazy loading the data as the file will be closed at the end of
    # _read_obj.
    obj = _read_obj(str(filepath), **_NO_MEMMAP_KWARGS)
    return [(obj["data"][:], FileHeader(obj["meta"]))]


def get_header(filepath):
    """
    Read an ASDF file and return just the header(s).

    Parameters
    ----------
    filepath : pathlib.Path, str
        The ASDF file to be read.

    Returns
    -------
    `list`
        A list of `sunpy.io._header.FileHeader` headers.
    """
    return [FileHeader(_read_obj(str(filepath))["meta"])]


def _read_obj(fname, **kwargs):
    with asdf.open(fname, **kwargs) as af:
        # TODO as asdf files can be structured in many ways some tests
        # of the structure and appropriate errors are needed
        # - does the file contain an "object" key?
        # - does "object" contain "meta" and "data"?
        # - is meta a dict?
        # - is data a asdf.tags.core.NDArrayType or ndarray?
        #   (NDArrayType is an asdf-specific type used for arrays).
        return af["object"]
