import asdf
import asdf_astropy
from pathlib import Path
import numpy as np 
from sunpy.io.header import FileHeader
from collections import OrderedDict
__all__ = ["write","read","get_header"]

def write(fname, data, header, **kwargs):
    """
    Take data and header pairs and save it to asdf file 
    inspired from https://docs.sunpy.org/en/stable/generated/gallery/saving_and_loading_data/genericmap_in_asdf.html

    format:
    map_name.asdf (main tree)
        meta: (header info)

        data : data blocks
    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        n-dimensional data array.
    header : `dict`
        A header dictionary.

    """
    map_name = Path(fname)
    map_name = map_name.name
    meta = dict(header)
    with asdf.AsdfFile() as af:
        af.tree={map_name:{"meta":meta,"data":data}}

        af.write_to(fname)


def read(fname,**kwargs):
    """
    A function to read asdf_file
    parameters
    ----------
    fname : Str

    returns : list of (data,meta) 
    
    """
    with asdf.open(fname) as af:

        map_name = Path(fname)
        map_name = map_name .name

        data = af[map_name]["data"].data
        data_array = np.asarray(data)

        meta_data= af[map_name]["meta"]
        meta_data = OrderedDict(meta_data)
        meta_data = FileHeader(meta_data)

        return [(data_array,meta_data)]
def get_header(fname):
    """
    read an asdf file and return the meta_data (meta) from the file
    
    """
    with asdf.open(fname) as af:

        map_name = Path(fname)
        map_name = map_name .name



        meta_data= af[map_name]["meta"]
        meta_data = OrderedDict(meta_data)
        meta_data = FileHeader(meta_data)

        return [meta_data]
    
