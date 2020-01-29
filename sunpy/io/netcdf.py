"""
This module provides a netCDF file reader.
"""
import collections
import numpy.ma as ma

from sunpy.io.header import FileHeader

__all__ = ['read', 'get_header', 'write']

HDPair = collections.namedtuple('HDPair', ['data', 'header'])


def read(filepath, **kwargs):
    """
    Reads a netCDF file.

    Parameters
    ----------
    filepath : `str`
        The file to be read.

    Returns
    -------
    pairs : `list`
        A list of (data, header) tuples.
    """
    # Put import here to speed up sunpy.io import time
    from netCDF4 import Dataset
    header = get_header(filepath)

    nc = Dataset(filepath, 'r')
    data = ma.masked_array(nc.variables['RAD'][:], mask=nc.variables['DQF'][:])

    return [HDPair(data, header[0])]


def get_header(filepath):
    """
    Reads the header from the file.

    Parameters
    ----------
    filepath : `str`
        The file to be read.

    Returns
    -------
    headers : list
        A list of headers read from the file.
    """
    # Put import here to speed up sunpy.io import time
    from netCDF4 import Dataset
    import datetime
    nc = Dataset(filepath, 'r')
    pydict = {}
    for k in nc.variables.keys():
        if not k == 'RAD' and not k == 'DQF':
            if nc.variables[k].dtype == '|S1':
                pydict[k.lower()] = nc.variables[k][:].tostring().decode()
            else:
                if not nc.variables[k][:].shape:
                    pydict[k.lower()] = nc.variables[k][:].item()
                else:
                    pydict[k.lower()] = nc.variables[k][:].filled()
    dateobs = datetime.datetime(2000, 1, 1, 12, 0) \
                + datetime.timedelta(seconds=float(nc.variables['DATE-OBS'][:]))
    pydict['date-obs'] = dateobs
    pydict['instrume'] = 'GOES-R Series Solar Ultraviolet Imager'

    return [FileHeader(pydict)]


def write(fname, data, header):
    """
    Place holder for required file writer.
    """
    raise NotImplementedError("No netCDF writer is implemented.")
