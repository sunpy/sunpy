import os
from collections import OrderedDict

import numpy as np

import asdf
from asdf.tags.core import NDArrayType
from astropy.utils.introspection import minversion

from sunpy.io._header import FileHeader

__all__ = ["write", "read", "get_header"]

# use memmap for and after asdf 3.1.0 to avoid a warning
# for the deprecated copy_arrays
if minversion(asdf, "3.1.0"):
    _NO_MEMMAP_KWARGS = {"memmap": False, "lazy_load": False}
else:
    _NO_MEMMAP_KWARGS = {"copy_arrays": True, "lazy_load": False}


def write(filepath, data, header, overwrite=False, **kwargs):
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
    overwrite : bool, optional
        If True, overwrite the file if it already exists. Default is False.
    """
    if os.path.exists(filepath) and not overwrite:
        raise FileExistsError(f"The file '{filepath}' already exists. Set 'overwrite=True' to overwrite it.")

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
        if "object" not in af.tree:
            raise KeyError("The ASDF does not contain 'object' key")
        obj = af.tree["object"]
        # - does "object" contain "meta" and "data"?
        if 'meta' not in obj or 'data' not in obj:
            raise TypeError("The object does not have any meta and data")
        meta = obj["meta"]
        data = obj["data"]
        # - is meta a dict?
        if not isinstance(meta, dict):
            raise TypeError(f"meta must be a dict not {type(meta)}")
        # - is data a asdf.tags.core.ANDArrayType or ndarray?
        if not isinstance(data, (NDArrayType| np.ndarray)):
            raise TypeError(f"data must be a NDArrayType or numpy ndarray not {type(data)}")
        return obj
