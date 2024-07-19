from pathlib import Path
from collections import OrderedDict

import numpy as np

import asdf

from sunpy.io._header import FileHeader

__all__ = ["write", "read", "get_header", "get_keys_name"]

def write(fname, data, header, **kwargs):
    """
    Take ``(data, header)`` pairs and save it to an ASDF file.

    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        N-dimensional data array.
    header : `dict`
        A header dictionary.
    """
    map_name = Path(fname)
    map_name = map_name.name
    ty = data.dtype
    data = data.astype(ty)
    meta = dict(header)
    asdf.AsdfFile({map_name : {"meta" : meta, "data" : data}},lazy_load=False,**kwargs).write_to(fname)


def read(fname, **kwargs):
    """
    Read an ASDF file.

    Parameters
    ----------
    filepath : `str`
        The ASDF file to be read.

    Returns
    -------
    `list`
        A list of ``(data, header)`` tuples.
    """

    with asdf.open(fname) as af:
        map_name = get_keys_name(fname)
        if isinstance(af[map_name],dict):
            data = af[map_name]["data"].data
            data_array = np.asarray(data)
            meta_data= af[map_name]["meta"]
            meta_data = OrderedDict(meta_data)
            meta_data = FileHeader(meta_data)
            return [(data_array,meta_data)]
        else:
            data = af[map_name].data
            meta_data = af[map_name].meta
            meta_data = FileHeader(meta_data)
            return [(data,meta_data)]

def get_header(fname):
    """
    Read an ASDF file and return just the header(s).

    Parameters
    ----------
    fname : `str`
        The ASDF file to be read.

    Returns
    -------
    `list`
        A list of `sunpy.io._header.FileHeader` headers.
    """
    with asdf.open(fname) as af:
        map_name = get_keys_name(fname)
        if isinstance(af[map_name], dict):
            meta_data= af[map_name]["meta"]
            meta_data = OrderedDict(meta_data)
            meta_data = FileHeader(meta_data)
            return [meta_data]
        else:
            meta_data = af[map_name].meta
            meta_data = FileHeader(meta_data)
            return [meta_data]


def get_keys_name(fname):
    """
    Returns the keys of primary tree (excluding the "asdf" and "history" trees).

    Parameters
    ----------
    fname : `str`
        The ASDF file to be read.

    Returns
    -------
    `str`
        Name of primary tree keys (excluding "asdf" and "history").
    """
    with asdf.open(fname) as af:
        root_keys = af.tree.keys()
        main_data_keys = [key for key in root_keys if key not in ['asdf_library', 'history']]
        return main_data_keys[0]
