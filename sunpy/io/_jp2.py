"""
This module provides a JPEG 2000 file reader.
"""
import os
from xml.etree import cElementTree as ET
import warnings
import numpy as np

from sunpy.io.header import FileHeader
from sunpy.util.io import HDPair, string_is_float
from sunpy.util.xml import xml_to_dict

__all__ = ['_read', '_get_header', '_write']


def _read(filepath, **kwargs):
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
    warnings.warn("_read is depricated and it is meant to be used for internal use only",DeprecationWarning)

    # Put import here to speed up sunpy.io import time
    from glymur import Jp2k

    header = _get_header(filepath)
    data = Jp2k(filepath)[...][::-1]
    return [HDPair(data, header[0])]


def _get_header(filepath):
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
    warnings.warn("_get_header is depricated and it is meant to be used for internal use only",DeprecationWarning)
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
        elif string_is_float(v):
            pydict[k] = float(v)

    # Remove newlines from comment
    if 'comment' in pydict:
        pydict['comment'] = pydict['comment'].replace("\n", "")

    # Is this file a Helioviewer Project JPEG2000 file?
    pydict['helioviewer'] = xml_box[0].xml.find('helioviewer') is not None

    return [FileHeader(pydict)]


def _header_to_xml(header):
    """
    Converts image header metadata into an XML Tree that can be inserted into
    a JP2 file header.

    Parameters
    ----------
    header : `MetaDict`
        A header dictionary to convert to xml.

    Returns
    ----------
    `lxml.etree._Element`
        A fits element where each child is an xml element
        in the form <key>value</key> derived from the key/value
        pairs in the given header dictionary
    """
    warnings.warn("_header_to_xml is depricated and it is meant to be used for internal use only",DeprecationWarning)
    # glymur uses lxml and will crash if trying to use
    # python's builtin xml.etree
    import lxml.etree as ET

    fits = ET.Element("fits")

    already_added = set()
    for key in header:
        # Some headers span multiple lines and get duplicated as keys
        # header.get will appropriately return all data, so if we see
        # a key again, we can assume it was already added to the xml tree.
        if (key in already_added):
            continue

        # Add to the set so we don't duplicate entries
        already_added.add(key)

        el = ET.SubElement(fits, key)
        data = header.get(key)
        if type(data) == bool:
            data = "1" if data else "0"
        else:
            data = str(data)

        el.text = data

    return fits


def _generate_jp2_xmlbox(header):
    """
    Generates the JP2 XML box to be inserted into the jp2 file.

    Parameters
    ----------
    header : `MetaDict`
        A header dictionary.

    Returns
    ----------
    `XMLBox`
        XML box containing FITS metadata to be used in jp2 headers
    """
    warnings.warn("_get_header is depricated and it is meant to be used for internal use only",DeprecationWarning)
    # glymur uses lxml and will crash if trying to use
    # python's builtin xml.etree
    import lxml.etree as ET
    from glymur import jp2box

    header_xml = _header_to_xml(header)
    meta = ET.Element("meta")
    meta.append(header_xml)
    tree = ET.ElementTree(meta)
    return jp2box.XMLBox(xml=tree)


def _write(fname, data, header, **kwargs):
    """
    Take a data header pair and write a JP2 file.

    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        n-dimensional data array.
    header : `dict`
        A header dictionary.
    kwargs :
        Additional keyword args are passed to the glymur.Jp2k constructor

    Notes
    -----
    Saving as a JPEG2000 will cast the data array to
    uint8 values to support the JPEG2000 format.
    """
    from glymur import Jp2k
    warnings.warn("_write is depricated and it is meant to be used for internal use only",DeprecationWarning)

    tmpname = fname + "tmp.jp2"
    jp2_data = np.uint8(data)
    jp2 = Jp2k(tmpname, jp2_data, **kwargs)

    # Append the XML data to the header information stored in jp2.box
    meta_boxes = jp2.box
    target_index = len(meta_boxes) - 1
    fits_box = _generate_jp2_xmlbox(header)
    meta_boxes.insert(target_index, fits_box)

    # Rewrites the jp2 file on disk with the xml data in the header
    jp2.wrap(fname, boxes=meta_boxes)

    os.remove(tmpname)
