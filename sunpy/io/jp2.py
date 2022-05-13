"""
This module provides a JPEG 2000 file reader.
"""
from xml.etree import cElementTree as ET

from sunpy.io.header import FileHeader
from sunpy.util.io import HDPair, is_float
from sunpy.util.xml import xml_to_dict

__all__ = ['read', 'get_header', 'write']


def read(filepath, **kwargs):
    """
    Reads a JPEG2000 file.

    Parameters
    ----------
    filepath : `str`
        The file to be read.
    **kwargs : `dict`
        Unused.

    Returns
    -------
    `list`
        A list of (data, header) tuples.
    """
    # Put import here to speed up sunpy.io import time
    from glymur import Jp2k

    header = get_header(filepath)
    data = Jp2k(filepath)[...][::-1]
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
    `list`
        A list of one header read from the file.
    """
    # Put import here to speed up sunpy.io import time
    from glymur import Jp2k
    jp2 = Jp2k(filepath)
    xml_box = [box for box in jp2.box if box.box_id == 'xml ']
    xmlstring = ET.tostring(xml_box[0].xml.find('fits'))
    pydict = xml_to_dict(xmlstring)["fits"]

    # Fix types
    for k, v in pydict.items():
        if v.isdigit():
            pydict[k] = int(v)
        elif is_float(v):
            pydict[k] = float(v)

    # Remove newlines from comment
    if 'comment' in pydict:
        pydict['comment'] = pydict['comment'].replace("\n", "")

    # Is this file a Helioviewer Project JPEG2000 file?
    pydict['helioviewer'] = xml_box[0].xml.find('helioviewer') is not None

    return [FileHeader(pydict)]


def write(fname, data, header):
    raise NotImplementedError("No jp2 writer is implemented.")
