"""
This module provides a JPEG 2000 file reader for internal use.
"""
import os

# We have to use lxml as lxml can not  serialize xml from the standard library
import lxml.etree as ET
import numpy as np

from sunpy.io._header import FileHeader
from sunpy.util.io import HDPair, string_is_float

__all__ = ['read', 'get_header', 'write']


def _sanative_value(value):
    if value is None:
        return
    if value.isdigit() or value.isnumeric():
        value = int(value)
    elif string_is_float(value):
        value = float(value)
    # We replace newlines with spaces so when we join the history
    # it won't be so ugly
    elif "\n" in value:
        value = value.replace("\n", " ")
    return value

def _parse_xml_metadata(xml_data):
    keycomments = {}
    final_dict = {}
    history = []
    for node in xml_data.xml.getroot().iter():
        if node.tag in ['HISTORY']:
            history.append(node.attrib.get('comment', ''))
            continue
        if node.text is not None:
            final_dict[node.tag] = _sanative_value(node.text)
            keycomments[node.tag] = node.attrib.get('comment', '')
    return {**final_dict, "HISTORY": "".join(history), 'KEYCOMMENTS': keycomments}

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
    # We assume that the header is the first XMLBox in the file
    xml_box = [box for box in jp2.box if box.box_id == 'xml '][0]
    pydict = _parse_xml_metadata(xml_box)
    return [FileHeader(pydict)]


def header_to_xml(header):
    """
    Converts image header metadata into an XML Tree that can be inserted into
    a JPEG2000 file header.

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
        if isinstance(data, bool):
            data = "1" if data else "0"
        else:
            data = str(data)
        el.text = data
    return fits


def generate_jp2_xmlbox(header):
    """
    Generates the JPEG2000 XML box to be inserted into the JPEG2000 file.

    Parameters
    ----------
    header : `MetaDict`
        A header dictionary.

    Returns
    ----------
    `XMLBox`
        XML box containing FITS metadata to be used in JPEG2000 headers
    """
    from glymur import jp2box

    header_xml = header_to_xml(header)
    meta = ET.Element("meta")
    meta.append(header_xml)
    tree = ET.ElementTree(meta)
    return jp2box.XMLBox(xml=tree)


def write(fname, data, header, **kwargs):
    """
    Take a data header pair and write a JPEG2000 file.

    Parameters
    ----------
    fname : `str`
        File name, with extension.
    data : `numpy.ndarray`
        N-Dimensional data array.
    header : `dict`
        A header dictionary.
    kwargs :
        Additional keyword args are passed to glymur's Jp2k constructor.

    Notes
    -----
    Saving as a JPEG2000 will cast the data array to uint8 values to support the JPEG2000 format.
    """
    from glymur import Jp2k

    tmp_filename = f"{fname}tmp.jp2"
    jp2_data = np.uint8(data)
    # The jp2 data is flipped when read in, so we have to flip it back before
    # saving. See https://github.com/sunpy/sunpy/pull/768 for context.
    flipped = np.flip(jp2_data, 0)
    jp2 = Jp2k(tmp_filename, flipped, **kwargs)
    # Append the XML data to the header information stored in jp2.box
    meta_boxes = jp2.box
    target_index = len(meta_boxes) - 1
    fits_box = generate_jp2_xmlbox(header)
    meta_boxes.insert(target_index, fits_box)
    # Rewrites the jp2 file on disk with the xml data in the header
    jp2.wrap(fname, boxes=meta_boxes)
    os.remove(tmp_filename)
